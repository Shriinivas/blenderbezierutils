# bezier_utils/drawing/primitives.py

import bpy
from mathutils import Matrix
from math import sqrt, radians, degrees, atan, cos, sin, tan, pi
from ..utils.bezier_math import (
    get3DVector,
    getSegsForArc,
    getWSDataForSegs,
    getInterpBezierPts,
    getInterpolatedVertsCo,
)
from ..utils.view_utils import getClosestPlaneToView, getCoordFromLoc
from ..core.snap import SnapParams
from ..core.hotkeys import FTHotKeys
from ..utils.math_utils import vectCmpWithMargin


class BaseDraw:
    def __init__(self, parent):
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.curvePts = []

    def newPoint(self, loc, ltype, rtype, faceIdx=None, is_intersection=False):
        self.curvePts.append([loc, loc, loc, ltype, rtype, faceIdx, is_intersection])

    def setCurvePt(self, idx, pt):
        self.curvePts[idx] = pt

    def popCurvePt(self):
        self.curvePts.pop()

    def setCurvePts(self, curvePts):
        self.curvePts = curvePts


class Primitive2DDraw(BaseDraw):
    def __init__(self, parent):
        super(Primitive2DDraw, self).__init__(parent)
        self.shapeSegCnt = 4
        self.curveObjOrigin = None

    def initialize(self):
        super(Primitive2DDraw, self).initialize()
        self.bbStart = None
        self.bbEnd = None
        self.XYstart = None
        self.freeAxesN = None

    def set2DAxes(self):
        params = bpy.context.window_manager.bezierToolkitParams
        self.freeAxesN = self.parent.snapper.getFreeAxesNormalized()

        if params.snapOrient == 'SURFACE':
            self.freeAxesN = [0, 1]
        elif len(self.freeAxesN) > 2 and params.snapOrient in {
            "GLOBAL",
            "REFERENCE",
            "CURR_POS",
        }:
            self.freeAxesN = getClosestPlaneToView(self.parent.rmInfo.rv3d)

    def getNumSegsLimits(self):
        return 2, 100

    def getShapePts(
        self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, theta, axisIdxs, z
    ):
        raise NotImplementedError("Call to abstract method.")

    dynamicParams = {
        0: (["UP_ARROW", "DOWN_ARROW"], "Up (incr) / Down (decr)"),
        1: (["RIGHT_ARROW", "LEFT_ARROW"], "Left (incr) / Right (decr)"),
        2: (["PAGE_UP", "PAGE_DOWN"], "Page Up (incr) / Page Down (decr)"),
        3: (["W", "S"], "W (incr) / S (decr)"),
        4: (["A", "D"], "A (incr) / D (decr)"),
        5: (["LEFT_BRACKET", "RIGHT_BRACKET"], "[ (incr) / ] (decr)"),
    }

    def getParamCnt():
        return len(Primitive2DDraw.dynamicParams)

    def getParamHotKeyDescriptions():
        return [
            Primitive2DDraw.dynamicParams[p][1]
            for p in range(len(Primitive2DDraw.dynamicParams))
        ]

    def getParamPropDefs():
        """
        Returns a list of tuples defining dynamic properties.
        Each tuple: (prop_name, prop_type_class, keywords_dict)
        """
        from bpy.props import FloatProperty # Ensure we have FloatProperty available or import it here
        from .math_fn import MathFnDraw # Import locally if needed or assume available
        
        # We need access to MathFnDraw constants but they are imported at top level.
        # But MathFnDraw is in math_fn.py which is imported in params.py... 
        # primitives.py imports nothing from math_fn. 
        # Wait, params.py imported MathFnDraw. 
        # primitives.py does NOT import MathFnDraw. 
        # The original code in params.py used MathFnDraw.startPrefix etc.
        # So we should put this logic in params.py or import MathFnDraw here.
        # Let's import it here.

        defs = []
        hks = Primitive2DDraw.getParamHotKeyDescriptions()

        for i in range(Primitive2DDraw.getParamCnt()):
            char = chr(ord('A') + i)
            
            # Start Value Property
            propName = MathFnDraw.startPrefix + str(i)
            defs.append((
                propName,
                FloatProperty,
                {
                    "name": 'Constant ' + char + ' Value', 
                    "description": 'Value of ' + char + ' used in equation', 
                    "default": MathFnDraw.defConstStart
                }
            ))

            # Incr Value Property
            propName = MathFnDraw.incrPrefix + str(i)
            defs.append((
                propName,
                FloatProperty,
                {
                    "name": 'Constant ' + char + ' Step', 
                    "description": 'Constant ' + char + ' increment / decrement step ' + \
                        ' (hot keys: ' + hks[i]+ ')', 
                    "default": MathFnDraw.defConstIncr
                }
            ))
        return defs

    for i in range(len(dynamicParams)):
        exec("def updateParam" + str(i) + "(self, event, rmInfo, isIncr):\n\tpass")

    def afterShapeSegCnt(self):
        pass

    def getCurvePts(self, numSegs, axisIdxs, z=None):
        params = bpy.context.window_manager.bezierToolkitParams
        if params.snapOrient == 'SURFACE':
            tm = self.parent.rmInfo.rv3d.view_matrix
        else:
            tm = self.parent.snapper.tm if self.parent.snapper.tm is not None else Matrix()

        idx0, idx1, idx2 = axisIdxs
        startAngle = params.drawStartAngle
        sweep = params.drawAngleSweep

        bbEnd = tm @ self.bbEnd

        mode = params.drawObjMode

        if mode == "CENTER":
            diffV = self.bbEnd - self.bbStart
            bbStart = tm @ (self.bbStart - diffV)
        else:
            bbStart = tm @ self.bbStart

        cX = (bbEnd[idx0] - bbStart[idx0]) / 2
        cY = (bbEnd[idx1] - bbStart[idx1]) / 2
        if cX == 0 and cY == 0:
            return None

        if z is None:
            z = bbStart[idx2]
        center2d = complex(cX, cY)

        orig = complex(bbStart[idx0], bbStart[idx1])
        self.curveObjOrigin = tm.inverted_safe() @ get3DVector(
            orig + center2d, axisIdxs, z
        )

        curvePts = self.getShapePts(
            mode, numSegs, bbStart, bbEnd, center2d, startAngle, sweep, axisIdxs, z
        )

        if curvePts is None:
            return None

        curvePts = [
            [tm.inverted_safe() @ p if not isinstance(p, str) else p for p in pts]
            for pts in curvePts
        ]

        return curvePts

    def updateCurvePts(self):
        axisIdxs = self.freeAxesN + sorted(list({0, 1, 2} - set(self.freeAxesN)))

        curvePts = self.getCurvePts(axisIdxs=axisIdxs, numSegs=self.shapeSegCnt)

        if curvePts is not None:
            self.setCurvePts(curvePts)
        else:
            self.setCurvePts(
                [
                    [self.bbStart, self.bbStart, self.bbStart],
                    [self.bbEnd, self.bbEnd, self.bbEnd],
                ]
            )

    def getPtLoc(self):
        parent = self.parent
        rmInfo = parent.rmInfo
        if self.freeAxesN is None:
            self.set2DAxes()

        params = bpy.context.window_manager.bezierToolkitParams
        if params.snapOrient == 'SURFACE' and self.bbStart is not None:
            return self.parent.snapper.get3dLocSnap(
                rmInfo,
                SnapParams(
                    parent.snapper,
                    transType='VIEW',
                    origType='REFERENCE',
                    refLineOrig=self.bbStart,
                    lastCo1Axis=True,
                    freeAxesN=[0, 1],
                    snapToPlane=True,
                ),
            )

        return self.parent.snapper.get3dLocSnap(
            rmInfo,
            SnapParams(
                parent.snapper,
                lastCo1Axis=True,
                freeAxesN=self.freeAxesN,
                snapToPlane=(len(self.freeAxesN) == 2),
            ),
        )

    def procDrawEvent(self, context, event, snapProc):
        parent = self.parent
        rmInfo = parent.rmInfo
        parent.snapper.getMetakeys()

        # Local import to avoid circular dependency
        from ..operators.modal_ops import ModalDrawBezierOp

        if len(self.curvePts) > 0 and not snapProc:
            if event.type in {
                "WHEELDOWNMOUSE",
                "WHEELUPMOUSE",
                "NUMPAD_PLUS",
                "NUMPAD_MINUS",
                "PLUS",
                "MINUS",
            }:
                if (
                    event.type in {"NUMPAD_PLUS", "NUMPAD_MINUS", "PLUS", "MINUS"}
                    and event.value == "PRESS"
                ):
                    return {"RUNNING_MODAL"}
                self.updateSegCount(
                    event,
                    rmInfo,
                    (event.type == "WHEELUPMOUSE" or event.type.endswith("PLUS")),
                )
                return {"RUNNING_MODAL"}

            for i in range(len(self.dynamicParams)):
                hotkeys = self.dynamicParams[i][0]
                if event.type in hotkeys:
                    if event.value == "RELEASE":
                        exec(
                            "self.updateParam"
                            + str(i)
                            + "(event, rmInfo, (event.type == hotkeys[0]))"
                        )
                    return {"RUNNING_MODAL"}

            if event.type == "H" or event.type == "h":
                if event.value == "RELEASE":
                    ModalDrawBezierOp.h = not ModalDrawBezierOp.h
                    self.parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
                return {"RUNNING_MODAL"}

        if self.bbStart is not None and (event.type == "RET" or event.type == "SPACE"):
            if event.value == "RELEASE":
                self.updateCurvePts()
                self.parent.confirm(context, event, self.curveObjOrigin)
                self.parent.snapper.resetSnap()
                self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if event.type == "ESC":
            if event.value == "RELEASE":
                self.parent.initialize()
                self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and event.type == "LEFTMOUSE" and event.value == "RELEASE":
            if self.bbStart is None:
                self.set2DAxes()
                loc = self.getPtLoc()
                self.XYstart = self.parent.rmInfo.xy
                self.bbStart = loc
                self.bbEnd = loc
                self.setCurvePts(
                    [[loc, loc, loc, "FREE", "FREE"], [loc, loc, loc, "FREE", "FREE"]]
                )
            else:
                loc = self.getPtLoc()
                self.bbEnd = loc
                self.updateCurvePts()
                parent.confirm(context, event, self.curveObjOrigin)
            return {"RUNNING_MODAL"}

        if snapProc or event.type == "MOUSEMOVE":
            if self.bbStart is not None:
                self.bbEnd = self.getPtLoc()
                self.updateCurvePts()
            parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
            return {"RUNNING_MODAL"}

        return {"PASS_THROUGH"}

    def getRefLine(self):
        if self.bbStart is not None:
            if self.bbEnd is not None:
                return [self.bbStart, self.bbEnd]
            else:
                return [self.bbStart]

        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[0] if len(refLine) > 0 else None


class ClosedShapeDraw(Primitive2DDraw):
    def __init__(self, parent, star=False):
        super(ClosedShapeDraw, self).__init__(parent)

    def updateSegCount(self, event, rmInfo, isIncr):
        minSegs, maxSegs = self.getNumSegsLimits()
        if isIncr and self.shapeSegCnt < maxSegs:
            self.shapeSegCnt += 1
        if not isIncr and self.shapeSegCnt > minSegs:
            self.shapeSegCnt -= 1
        self.afterShapeSegCnt()
        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
        return True

    def updateParam1(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams
        theta = params.drawAngleSweep

        if isIncr:
            if theta >= -10 and theta <= 0:
                params.drawAngleSweep = 10
            elif theta > 350:
                params.drawAngleSweep = -350
            else:
                params.drawAngleSweep = theta + 10
        else:
            if theta <= 10 and theta >= 0:
                params.drawAngleSweep = -10
            elif theta < -350:
                params.drawAngleSweep = 350
            else:
                params.drawAngleSweep = theta - 10
        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
        return True


class RectangleDraw(ClosedShapeDraw):
    def __init__(self, parent, star=False):
        super(RectangleDraw, self).__init__(parent)

    def getShapePts(
        self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, theta, axisIdxs, z
    ):
        idx0, idx1, idx2 = axisIdxs
        pt0 = complex(bbStart[idx0], bbStart[idx1])
        pt1 = complex(bbEnd[idx0], bbStart[idx1])
        pt2 = complex(bbEnd[idx0], bbEnd[idx1])
        pt3 = complex(bbStart[idx0], bbEnd[idx1])

        curvePts = []
        for pt2d in [pt0, pt1, pt2, pt3, pt0]:
            pt = get3DVector(pt2d, axisIdxs, z)
            curvePts.append([pt, pt, pt, "VECTOR", "VECTOR"])

        return curvePts


class PolygonDraw(ClosedShapeDraw):
    def updateParam0(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams

        params.drawStarOffset += 0.1 if (isIncr) else -0.1

        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
        return True

    def getNumSegsLimits(self):
        return 3, 100

    def __init__(self, parent, star=False):
        super(PolygonDraw, self).__init__(parent)
        self.star = star

    def getShapePts(
        self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, theta, axisIdxs, z
    ):
        params = bpy.context.window_manager.bezierToolkitParams
        params.drawSides = numSegs
        idx0, idx1, idx2 = axisIdxs
        cX, cY = center2d.real, center2d.imag
        radius = sqrt(cX * cX + cY * cY)
        offset = params.drawStarOffset
        orig = complex(bbStart[idx0], bbStart[idx1])

        if cX == 0:
            startAngle = 90 * ((cY / abs(cY)) if cY != 0 else 1)
        else:
            startAngle = degrees(atan(cY / cX))
        if cX < 0:
            startAngle += 180
        thetaIncr = theta / numSegs
        curvePts = []
        for segCnt in range(numSegs + 1):
            if self.star:
                angle = startAngle + thetaIncr * segCnt - thetaIncr / 2
                shift = complex(
                    offset * radius * cos(radians(angle)),
                    offset * radius * sin(radians(angle)),
                )
                pt = get3DVector(orig + center2d + shift, axisIdxs, z)
                curvePts.append([pt, pt, pt, "VECTOR", "VECTOR"])
            if not self.star or segCnt < numSegs:
                angle = startAngle + thetaIncr * segCnt
                shift = complex(
                    radius * cos(radians(angle)), radius * sin(radians(angle))
                )
                pt = get3DVector(orig + center2d + shift, axisIdxs, z)
                curvePts.append([pt, pt, pt, "VECTOR", "VECTOR"])

        return curvePts


class EllipseDraw(ClosedShapeDraw):
    # https://math.stackexchange.com/questions/22064/calculating-a-point-that-lies-on-an-ellipse-given-an-angle
    def getPtAtAngle(a, b, theta):
        if theta < 0:
            theta += 2 * pi
        denom = sqrt(b * b + a * a * tan(theta) * tan(theta))
        num = a * b
        x = num / denom
        if pi / 2 < theta <= 3 * pi / 2:
            x = -x
        y = x * tan(theta)
        return complex(x, y)

    def __init__(self, parent):
        super(EllipseDraw, self).__init__(parent)

    def updateParam0(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams
        theta = params.drawStartAngle

        if isIncr:
            if theta > 350:
                params.drawStartAngle = 0
            else:
                params.drawStartAngle = theta + 10
        else:
            if theta < -350:
                params.drawStartAngle = 0
            else:
                params.drawStartAngle = theta - 10

        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs={}, hltEndSeg=False)
        return True

    def getShapePts(
        self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, theta, axisIdxs, z
    ):
        idx0, idx1, idx2 = axisIdxs
        cX, cY = center2d.real, center2d.imag

        if cX == 0 or cY == 0:
            return None

        radius = complex(cX, cY)  # Actually same as center2d
        orig = complex(bbStart[idx0], bbStart[idx1])


        sweep = 1
        startIdx = 0
        endIdx = 1
        rvs = False

        if theta < 0:
            sweep = 0
            startIdx = 1
            endIdx = 0
            rvs = True

        pt1 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(startAngle))
        a1 = startAngle + theta / 2
        if a1 > 360:
            a1 = a1 - 360
        if a1 < -360:
            a1 = a1 + 360
        a2 = startAngle + theta
        if a2 > 360:
            a2 = a2 - 360
        if a2 < -360:
            a2 = a2 + 360
        pt2 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(a1))
        pt3 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(a2))

        pt1 = orig + pt1 + center2d
        pt2 = orig + pt2 + center2d
        pt3 = orig + pt3 + center2d

        endPts = [pt1, pt2]
        segs1 = getSegsForArc(
            endPts[startIdx], radius, 1, endPts[endIdx], 10, axisIdxs, z
        )

        endPts = [pt2, pt3]
        segs2 = getSegsForArc(
            endPts[startIdx], radius, 1, endPts[endIdx], 10, axisIdxs, z
        )

        segElems = [segs1, segs2]
        segs = segElems[startIdx] + segElems[endIdx]

        curvePts = getWSDataForSegs(segs)
        if len(curvePts) < 2:
            return None

        pts = getInterpBezierPts(curvePts, subdivPerUnit=100, segLens=None)
        if len(pts) < 2:
            return None

        vertCos = getInterpolatedVertsCo(pts, numSegs)

        # TODO: A more efficient approach for dividing ellipse uniformly
        newSegs = []
        for i in range(1, len(vertCos)):
            segStart = complex(vertCos[i - 1][idx0], vertCos[i - 1][idx1])
            segEnd = complex(vertCos[i][idx0], vertCos[i][idx1])
            endPts = [segStart, segEnd]
            segs1 = getSegsForArc(
                endPts[startIdx], radius, sweep, endPts[endIdx], 1, axisIdxs, z
            )

            newSegs += segs1

        if rvs:
            newSegs = reversed(newSegs)
        curvePts = getWSDataForSegs(newSegs)

        if len(curvePts) < 2:
            return None

        if vectCmpWithMargin(curvePts[0][1], curvePts[-1][1]):
            ldiffV = curvePts[0][1] - curvePts[0][0]
            rdiffV = curvePts[-1][2] - curvePts[-1][1]
            if vectCmpWithMargin(ldiffV, rdiffV):
                curvePts[0][3] = "ALIGNED"
                curvePts[0][4] = "ALIGNED"
                curvePts[-1][3] = "ALIGNED"
                curvePts[-1][4] = "ALIGNED"

        return curvePts


class BezierDraw(BaseDraw):
    def __init__(self, parent):
        super(BezierDraw, self).__init__(parent)
        self.reset()

    def reset(self):
        self.capture = False
        self.grabRepos = False
        self.dissociateHdl = False

    def moveBezierPt(self, loc, faceIdx=None):
        if len(self.curvePts) > 0:
            pt = self.curvePts[-1][:]
            # Make sure list is long enough for faceIdx
            if len(pt) < 6:
                pt.append(faceIdx)
            else:
                pt[5] = faceIdx
            pt[0] = loc
            pt[1] = loc
            pt[2] = loc
            self.setCurvePt(-1, pt)

    def movePointByDelta(self, delta):
        if len(self.curvePts) > 0:
            pt = self.curvePts[-1][:]
            pt[1] += delta
            self.setCurvePt(-1, pt)

    def moveBptElem(self, handle, loc):
        idx = {"left": 0, "pt": 1, "right": 2}[handle]
        if len(self.curvePts) > 0:
            pt = self.curvePts[-1][:]
            pt[idx] = loc
            self.setCurvePt(-1, pt)

    def moveBptElemByDelta(self, handle, delta):
        idx = {"left": 0, "pt": 1, "right": 2}[handle]
        if len(self.curvePts) > 0:
            pt = self.curvePts[-1][:]
            pt[idx] += delta
            self.setCurvePt(-1, pt)

    def resetHandle(self, handle):
        if len(self.curvePts) > 0:
            self.moveBptElem(handle, self.curvePts[-1][1])

    def isHandleSet(self):
        if len(self.curvePts) == 0:
            return False
        co = self.curvePts[-1][1]
        lh = self.curvePts[-1][0]
        rh = self.curvePts[-1][2]
        if not vectCmpWithMargin(co, lh) or not vectCmpWithMargin(co, rh):
            return True
        return False

    def procDrawEvent(self, context, event, snapProc):
        rmInfo = self.parent.rmInfo
        snapper = self.parent.snapper
        metakeys = snapper.getMetakeys()

        if self.capture and FTHotKeys.isHotKey(
            FTHotKeys.hkGrabRepos, event.type, metakeys
        ):
            if event.value == "RELEASE":
                self.dissociateHdl = False
                self.grabRepos = not self.grabRepos
            return {"RUNNING_MODAL"}

        if self.capture and FTHotKeys.isHotKey(
            FTHotKeys.hkDissociateHdl, event.type, metakeys
        ):
            if event.value == "RELEASE":
                self.grabRepos = False
                self.dissociateHdl = not self.dissociateHdl
            return {"RUNNING_MODAL"}

        # This can happen only when space was entered and something was there
        # for Snapper to process
        if snapProc and snapper.digitsConfirmed:
            snapper.resetSnap()

            # Because resetSnap sets this to False (TODO: Refactor resetSnap)
            snapper.digitsConfirmed = True

            # First space / enter is equivalent to mouse press without release
            if not self.capture:
                self.capture = True
                snapper.setStatus(rmInfo.area, None)
                return {"RUNNING_MODAL"}
            else:
                # Second space / enter means it should be processed here,
                # set snapProc to False so this modal will process it
                snapProc = False

        if not snapProc and event.type == "ESC":
            if event.value == "RELEASE":
                if self.grabRepos:
                    self.grabRepos = False
                elif self.dissociateHdl:
                    self.dissociateHdl = False
                elif self.capture and self.isHandleSet():
                    self.resetHandle("left")
                    self.resetHandle("right")

                    # Needed to indicate next space / entered to be processed here
                    snapper.digitsConfirmed = True
                    snapper.setStatus(rmInfo.area, None)
                else:
                    self.parent.initialize()  # TODO: should not access parent.initialize
                self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and FTHotKeys.isHotKey(
            FTHotKeys.hkResetLastHdl, event.type, metakeys
        ):
            if event.value == "RELEASE":
                if len(self.curvePts) > 1:
                    if len(self.curvePts) > 2:
                        idx = len(self.curvePts) - 2
                        pt = self.curvePts[idx][:]
                        pt[2] = pt[1].copy()
                        self.setCurvePt(idx, pt)
                    self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and FTHotKeys.isHotKey(
            FTHotKeys.hkUndoLastSeg, event.type, metakeys
        ):
            if event.value == "RELEASE":
                snapper.resetSnap()
                if not self.capture:
                    if len(self.curvePts) > 0:
                        self.popLastSegment()

                # Because there is an extra point (the current one)
                if len(self.curvePts) <= 1:
                    self.initialize()
                    self.reset()
                else:
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc, snapper.lastFaceIdx)
                self.capture = False
                self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and (event.type == "RET" or event.type == "SPACE"):
            if event.value == "RELEASE":
                if snapper.digitsConfirmed:
                    self.reset()
                    snapper.digitsConfirmed = False
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.newPoint(loc, "ALIGNED", "ALIGNED", snapper.lastFaceIdx)
                    self.resolve_segment_on_click(context)
                    self.parent.redrawBezier(rmInfo)
                else:
                    if len(self.curvePts) > 0:
                        self.popCurvePt()
                    self.parent.confirm(context, event)
                    snapper.resetSnap()
                    self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and event.type == "LEFTMOUSE" and event.value == "PRESS":
            if len(self.curvePts) == 0:
                loc = snapper.get3dLocSnap(rmInfo)
                self.newPoint(loc, "ALIGNED", "ALIGNED", snapper.lastFaceIdx)

            # Special condition for hot-key single axis lock (useful)
            if len(snapper.freeAxes) == 1 and len(self.curvePts) > 1:
                snapper.resetSnap()

            if self.capture:
                self.parent.pressT = None  # Lock capture
            else:
                self.capture = True
            return {"RUNNING_MODAL"}

        if not snapProc and event.type == "LEFTMOUSE" and event.value == "RELEASE":
            # ~ if(snapper.isLocked()):
            # ~ if(len(self.curvePts) == 1):
            # ~ self.moveBptElem('right', \
            # ~ snapper.get3dLocSnap(rmInfo))# changes only rt handle
            # ~ return {'RUNNING_MODAL'}

            # See the special condition above in event.value == 'PRESS'
            # ~ if(len(snapper.freeAxes) > 1):
            snapper.resetSnap()
            self.reset()

            # Rare condition: This happens e. g. when user clicks on header menu
            # like Object->Transform->Move. These ops consume press event but not release
            # So update the snap locations anyways if there was some transformation
            if len(self.curvePts) == 0:
                self.parent.updateSnapLocs()  # Subclass (TODO: have a relook)

            elif self.parent.doubleClick:
                if len(self.curvePts) > 0:
                    self.popCurvePt()
                self.parent.confirm(context, event)
                self.parent.redrawBezier(rmInfo)

            else:
                if self.parent.click:
                    loc = self.curvePts[-1][1]
                    self.moveBptElem("left", loc)
                    self.moveBptElem("right", loc)
                else:
                    loc = snapper.get3dLocSnap(rmInfo)

                # ~ if(len(self.curvePts) == 1):
                # ~ self.moveBptElem('right', loc)# changes only rt handle

                self.newPoint(loc, "ALIGNED", "ALIGNED", snapper.lastFaceIdx)
                self.resolve_segment_on_click(context)
                self.parent.redrawBezier(rmInfo)

            return {"RUNNING_MODAL"}

        # Refresh also in case of snapper events
        # except when digitsConfirmed (to give user opportunity to draw a straight line)
        # ~ if ((snapProc and not snapper.digitsConfirmed) \
        if snapProc or event.type == "MOUSEMOVE":
            bpy.context.window.cursor_set("DEFAULT")
            hdlPtIdxs = None
            if len(self.curvePts) > 0:
                if self.capture:
                    lastPt = self.curvePts[-1][:]
                    if self.grabRepos:
                        rtHandle = lastPt[2].copy()
                        xy2 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, lastPt[1])
                        xy1 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, rtHandle)
                        loc = snapper.get3dLocSnap(
                            rmInfo,
                            SnapParams(
                                snapper, xyDelta=[xy1[0] - xy2[0], xy1[1] - xy2[1]]
                            ),
                        )
                        delta = loc - lastPt[1]
                        self.moveBptElemByDelta("pt", delta)
                        self.moveBptElemByDelta("left", delta)
                        self.moveBptElemByDelta("right", delta)
                        if len(self.curvePts[-1]) > 5:
                            self.curvePts[-1][5] = snapper.lastFaceIdx
                    else:
                        loc = snapper.get3dLocSnap(rmInfo)
                        delta = loc - lastPt[1]
                        if self.dissociateHdl:
                            lastPt[3] = "FREE"
                            lastPt[4] = "FREE"
                        else:
                            lastPt[0] = lastPt[1] - delta
                            lastPt[3] = "ALIGNED"
                            lastPt[4] = "ALIGNED"
                        lastPt[2] = lastPt[1] + delta
                        self.setCurvePt(-1, lastPt)
                    hdlPtIdxs = {len(self.curvePts) - 1}
                else:
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc, snapper.lastFaceIdx)
                    hdlPtIdxs = {len(self.curvePts) - 2}

            self.parent.redrawBezier(rmInfo, hdlPtIdxs=hdlPtIdxs)
            return {"RUNNING_MODAL"}

        return {"PASS_THROUGH"} if not snapProc else {"RUNNING_MODAL"}

    def popLastSegment(self):
        if len(self.curvePts) <= 1:
            return
        active_pt = self.curvePts.pop()
        self.curvePts.pop()
        while len(self.curvePts) > 0 and len(self.curvePts[-1]) > 6 and self.curvePts[-1][6]:
            self.curvePts.pop()
        self.curvePts.append(active_pt)

    def resolve_segment_on_click(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        if params.snapOrient != 'SURFACE' or len(self.curvePts) < 3:
            return
            
        obj = context.active_object
        if obj is None or obj.type != 'MESH':
            return
            
        rmInfo = self.parent.rmInfo
        cached_bm = self.parent.get_cached_bmesh(obj)
        
        pt_active = self.curvePts.pop()
        
        end_idx = len(self.curvePts) - 1
        start_idx = end_idx - 1
        while start_idx > 0 and len(self.curvePts[start_idx]) > 6 and self.curvePts[start_idx][6]:
            start_idx -= 1
            
        pt_start = self.curvePts[start_idx]
        pt_end = self.curvePts[end_idx]
        
        # Resolve face indices to the current view only if not already set (keep mesh-relative face indices stable during view rotation)
        if len(pt_start) <= 5:
            pt_start.append(find_visible_face_at_loc(obj, rmInfo.region, rmInfo.rv3d, pt_start[1], None))
        elif pt_start[5] is None:
            pt_start[5] = find_visible_face_at_loc(obj, rmInfo.region, rmInfo.rv3d, pt_start[1], None)

        if len(pt_end) <= 5:
            pt_end.append(find_visible_face_at_loc(obj, rmInfo.region, rmInfo.rv3d, pt_end[1], None))
        elif pt_end[5] is None:
            pt_end[5] = find_visible_face_at_loc(obj, rmInfo.region, rmInfo.rv3d, pt_end[1], None)
        
        pt_start[3] = 'FREE'
        pt_start[4] = 'FREE'
        pt_end[3] = 'FREE'
        pt_end[4] = 'FREE'
        
        p_start = pt_start[1]
        h_start = pt_start[2]
        h_end = pt_end[0]
        p_end = pt_end[1]
        
        face_start = pt_start[5]
        face_end = pt_end[5]
        
        intersections, face_start, face_end = get_surface_intersection_points(
            obj, p_start, h_start, h_end, p_end, face_start, face_end,
            rmInfo.region, rmInfo.rv3d, cached_bm
        )
        
        if len(intersections) > 0:
            self.curvePts.pop()
            for p_int, face_prev, face_next, t_val in intersections:
                self.curvePts.append([p_int, p_int, p_int, 'FREE', 'FREE', face_next, True])
            if len(pt_end) > 5:
                pt_end[5] = intersections[-1][2]
            else:
                pt_end.append(intersections[-1][2])
            self.curvePts.append(pt_end)
            
            new_end_idx = len(self.curvePts) - 1
            t_vals = [0.0] + [item[3] for item in intersections] + [1.0]
            
            for sub_idx in range(len(intersections) + 1):
                idx_curr = start_idx + sub_idx
                pt_a = self.curvePts[idx_curr]
                pt_b = self.curvePts[idx_curr + 1]
                
                r0, r1, r2, r3 = get_bezier_sub_segment(
                    p_start, h_start, h_end, p_end,
                    t_vals[sub_idx], t_vals[sub_idx+1]
                )
                
                face_idx = face_start if sub_idx == 0 else pt_a[5]
                if face_idx is not None:
                    poly = obj.data.polygons[face_idx]
                    normal_world = obj.matrix_world.to_3x3() @ poly.normal
                    normal_world.normalize()
                    pt_a[2] = project_handle_to_face_plane(pt_a[1], r1, normal_world, rmInfo.region, rmInfo.rv3d)
                    pt_b[0] = project_handle_to_face_plane(pt_b[1], r2, normal_world, rmInfo.region, rmInfo.rv3d)
                else:
                    pt_a[2] = r1
                    pt_b[0] = r2
            
            if params.surfaceMode == 'SMOOTH':
                for idx in range(start_idx + 1, new_end_idx):
                    pt = self.curvePts[idx]
                    v_left = (pt[0] - pt[1]).normalized()
                    v_right = (pt[2] - pt[1]).normalized()
                    t_avg = (v_right - v_left).normalized()
                    l_left = (pt[0] - pt[1]).length
                    l_right = (pt[2] - pt[1]).length
                    pt[2] = pt[1] + t_avg * l_right
                    pt[0] = pt[1] - t_avg * l_left
        else:
            if pt_start[5] is not None and pt_end[5] is not None:
                face_idx = pt_start[5]
                poly = obj.data.polygons[face_idx]
                normal_world = obj.matrix_world.to_3x3() @ poly.normal
                normal_world.normalize()
                pt_start[2] = project_handle_to_face_plane(pt_start[1], pt_start[2], normal_world, rmInfo.region, rmInfo.rv3d)
                
                face_next_idx = pt_end[5]
                poly = obj.data.polygons[face_next_idx]
                normal_world = obj.matrix_world.to_3x3() @ poly.normal
                normal_world.normalize()
                pt_end[0] = project_handle_to_face_plane(pt_end[1], pt_end[0], normal_world, rmInfo.region, rmInfo.rv3d)
                
        self.curvePts.append(pt_active)

    def getRefLine(self):
        if len(self.curvePts) > 0:
            idx = 0
            if self.capture:
                if self.grabRepos and len(self.curvePts) > 1:
                    idx = -2
                else:
                    idx = -1
            # There should always be min 2 pts if not capture, check anyway
            elif len(self.curvePts) > 1:
                idx = -2
            if (len(self.curvePts) + (idx - 1)) >= 0:
                return [self.curvePts[idx - 1][1], self.curvePts[idx][1]]
            else:
                return [self.curvePts[idx][1]]
        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[-1] if len(refLine) > 0 else None


def get_bezier_sub_segment(p0, p1, p2, p3, t1, t2):
    def split_at(points, t):
        p0, p1, p2, p3 = points
        p01 = p0 + (p1 - p0) * t
        p12 = p1 + (p2 - p1) * t
        p23 = p2 + (p3 - p2) * t
        p012 = p01 + (p12 - p01) * t
        p123 = p12 + (p23 - p12) * t
        p0123 = p012 + (p123 - p012) * t
        return (p0, p01, p012, p0123), (p0123, p123, p23, p3)
        
    if t2 < 1.0:
        left, _ = split_at((p0, p1, p2, p3), t2)
    else:
        left = (p0, p1, p2, p3)
        
    if t1 > 0.0:
        t_scaled = t1 / t2 if t2 > 0.0 else 0.0
        _, right = split_at(left, t_scaled)
        return right
    else:
        return left


def find_visible_face_at_loc(obj, region, rv3d, p_co, fallback_face_idx):
    from bpy_extras.view3d_utils import location_3d_to_region_2d
    from ..utils.view_utils import getFaceUnderMouse
    from mathutils import Vector
    
    xy = location_3d_to_region_2d(region, rv3d, p_co)
    if xy is None:
        return fallback_face_idx
    _, _, face_idx = getFaceUnderMouse(obj, region, rv3d, xy, 1000000)
    if face_idx is not None:
        return face_idx
    # Try tiny offsets in pixels for floating-point tolerance on edges
    offsets = [
        (1, 0), (-1, 0), (0, 1), (0, -1)
    ]
    for dx, dy in offsets:
        xy_off = xy + Vector((dx, dy))
        _, _, face_idx = getFaceUnderMouse(obj, region, rv3d, xy_off, 1000000)
        if face_idx is not None:
            return face_idx
    return fallback_face_idx


def find_closest_face_locally(obj, bm, p_co, start_face_idx):
    if start_face_idx is None:
        return None
    try:
        bm.faces.ensure_lookup_table()
    except Exception:
        pass
    candidate_faces = {start_face_idx}
    try:
        face = bm.faces[start_face_idx]
        for edge in face.edges:
            for lf in edge.link_faces:
                candidate_faces.add(lf.index)
    except Exception:
        pass
    
    try:
        extended = list(candidate_faces)
        for f_idx in extended:
            face = bm.faces[f_idx]
            for edge in face.edges:
                for lf in edge.link_faces:
                    candidate_faces.add(lf.index)
    except Exception:
        pass
                
    best_face_idx = start_face_idx
    min_dist = 1e9
    import mathutils
    try:
        inv_mw = obj.matrix_world.inverted_safe()
        p_local = inv_mw @ p_co
        for f_idx in candidate_faces:
            face = bm.faces[f_idx]
            verts = [v.co for v in face.verts]
            face_dist = 1e9
            for idx in range(1, len(verts) - 1):
                v0, v1, v2 = verts[0], verts[idx], verts[idx+1]
                closest = mathutils.geometry.closest_point_on_tri(p_local, v0, v1, v2)
                dist = (p_local - closest).length
                if dist < face_dist:
                    face_dist = dist
            if face_dist < min_dist:
                min_dist = face_dist
                best_face_idx = f_idx
    except Exception:
        pass
    return best_face_idx


def get_surface_intersection_points(obj, p_start, h_start, h_end, p_end, face_start, face_end, region, rv3d, bm):
    from mathutils.geometry import intersect_line_line_2d, intersect_line_line
    from bpy_extras.view3d_utils import region_2d_to_origin_3d, region_2d_to_vector_3d
    
    print("\n--- INTERSECTION DEBUG ---")
    print(f"p_start: {p_start}, p_end: {p_end}")
    print(f"Original face_start: {face_start}, face_end: {face_end}")
    
    M = 100
    pts_3d = []
    pts_2d = []
    for i in range(M + 1):
        t = i / M
        p_t = (1-t)**3 * p_start + 3*(1-t)**2 * t * h_start + 3*(1-t)*t**2 * h_end + t**3 * p_end
        pts_3d.append(p_t)
        xy_t = getCoordFromLoc(region, rv3d, p_t)
        pts_2d.append(xy_t)
        
    print(f"Projected {len(pts_2d)} points. Start 2D: {pts_2d[0]}, End 2D: {pts_2d[-1]}")
    if any(pt.x == 9000 for pt in pts_2d):
        print("WARNING: Some projected curve points failed (returned 9000)!")
        
    # Resolve face_start/face_end locally to find the correct entry/exit faces
    if face_start is not None:
        face_start = find_closest_face_locally(obj, bm, pts_3d[1], face_start)
    if face_start is None:
        face_start = find_visible_face_at_loc(obj, region, rv3d, p_start, face_start)
        
    if face_end is not None:
        face_end = find_closest_face_locally(obj, bm, pts_3d[-2], face_end)
    if face_end is None:
        face_end = find_visible_face_at_loc(obj, region, rv3d, p_end, face_end)
        
    print(f"Resolved face_start: {face_start}, face_end: {face_end}")
    
    if face_start is not None:
        # Check if the curve immediately leaves the mesh surface / goes outside the silhouette
        def is_point_in_polygon_2d(x, y, poly):
            n = len(poly)
            inside = False
            p1x, p1y = poly[0].x, poly[0].y
            for i in range(n + 1):
                p2x, p2y = poly[i % n].x, poly[i % n].y
                if y > min(p1y, p2y):
                    if y <= max(p1y, p2y):
                        if x <= max(p1x, p2x):
                            if p1y != p2y:
                                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                            if p1x == p2x or x <= xints:
                                inside = not inside
                p1x, p1y = p2x, p2y
            return inside

        # Find the first point along the curve that is at least 1 pixel away in screen space
        p_check = None
        for k in range(1, len(pts_2d)):
            if (pts_2d[k] - pts_2d[0]).length > 1.0:
                p_check = pts_2d[k]
                break
        if p_check is None:
            p_check = pts_2d[1]

        # Get local neighborhood of faces
        candidate_faces = {face_start}
        try:
            for edge in bm.faces[face_start].edges:
                for lf in edge.link_faces:
                    candidate_faces.add(lf.index)
            extended = list(candidate_faces)
            for f_idx in extended:
                for edge in bm.faces[f_idx].edges:
                    for lf in edge.link_faces:
                        candidate_faces.add(lf.index)
        except Exception:
            pass

        ray_dir = region_2d_to_vector_3d(region, rv3d, p_check)
        found_face = None
        for f_idx in candidate_faces:
            poly = obj.data.polygons[f_idx]
            normal_world = obj.matrix_world.to_3x3() @ poly.normal
            normal_world.normalize()
            if normal_world.dot(ray_dir) < 0.0: # Front-facing
                f = bm.faces[f_idx]
                face_poly_2d = [getCoordFromLoc(region, rv3d, obj.matrix_world @ v.co) for v in f.verts]
                if is_point_in_polygon_2d(p_check.x, p_check.y, face_poly_2d):
                    found_face = f_idx
                    break

        if found_face is not None:
            face_start = found_face
        else:
            print("Start of segment immediately went outside the silhouette boundary!")
            face_start = None
            
    if face_start is None:
        if face_end is not None:
            # Start is outside, end is inside - trace backwards by reversing arguments recursively
            print("Reversing intersection trace because face_start is None and face_end is not None")
            rev_intersections, resolved_face_end, resolved_face_start = get_surface_intersection_points(
                obj, p_end, h_end, h_start, p_start, face_end, face_start, region, rv3d, bm
            )
            intersections = []
            for p_int, f_from, f_to, t_rev in reversed(rev_intersections):
                intersections.append((p_int, f_to, f_from, 1.0 - t_rev))
            return intersections, resolved_face_start, resolved_face_end
        else:
            print("Returning [] because both face_start and face_end are None")
            return [], face_start, face_end
        
    curr_face_idx = face_start
    t_min = 0.005
    intersections = []
    visited = {curr_face_idx}
    
    max_steps = 150
    for step in range(max_steps):
        if curr_face_idx == face_end or curr_face_idx is None:
            print(f"Ending trace loop: curr_face_idx={curr_face_idx}, face_end={face_end}")
            break
            
        bm_face = bm.faces[curr_face_idx]
        candidates = []
        
        print(f"Step {step}: curr_face_idx={curr_face_idx}")
        
        for bm_edge in bm_face.edges:
            v1 = obj.matrix_world @ bm_edge.verts[0].co
            v2 = obj.matrix_world @ bm_edge.verts[1].co
            xy1 = getCoordFromLoc(region, rv3d, v1)
            xy2 = getCoordFromLoc(region, rv3d, v2)
            print(f"  Edge {bm_edge.index} (verts {bm_edge.verts[0].index}-{bm_edge.verts[1].index}): 2D {xy1} to {xy2}")
            
            for j in range(M):
                t_A = j / M
                t_B = (j + 1) / M
                if t_B <= t_min:
                    continue
                
                res_2d = intersect_line_line_2d(pts_2d[j], pts_2d[j+1], xy1, xy2)
                if res_2d is not None:
                    eps = 1e-3
                    x_min1, x_max1 = min(pts_2d[j].x, pts_2d[j+1].x) - eps, max(pts_2d[j].x, pts_2d[j+1].x) + eps
                    y_min1, y_max1 = min(pts_2d[j].y, pts_2d[j+1].y) - eps, max(pts_2d[j].y, pts_2d[j+1].y) + eps
                    x_min2, x_max2 = min(xy1.x, xy2.x) - eps, max(xy1.x, xy2.x) + eps
                    y_min2, y_max2 = min(xy1.y, xy2.y) - eps, max(xy1.y, xy2.y) + eps
                    
                    if (x_min1 <= res_2d.x <= x_max1 and y_min1 <= res_2d.y <= y_max1 and
                        x_min2 <= res_2d.x <= x_max2 and y_min2 <= res_2d.y <= y_max2):
                        
                        seg_len = (pts_2d[j+1] - pts_2d[j]).length
                        if seg_len > 1e-6:
                            t_seg = (res_2d - pts_2d[j]).length / seg_len
                        else:
                            t_seg = 0.0
                        t_val = t_A + (t_B - t_A) * t_seg
                        
                        ray_origin = region_2d_to_origin_3d(region, rv3d, res_2d)
                        ray_dir = region_2d_to_vector_3d(region, rv3d, res_2d)
                        
                        res_3d = intersect_line_line(ray_origin, ray_origin + ray_dir, v1, v2)
                        if res_3d is not None:
                            p_ray, p_edge = res_3d
                            edge_vec = v2 - v1
                            edge_len = edge_vec.length
                            if edge_len > 1e-6:
                                factor = edge_vec.dot(p_edge - v1) / (edge_len * edge_len)
                                factor = max(0.0, min(1.0, factor))
                                p_edge = v1 + edge_vec * factor
                            
                            neighbor_face = None
                            for f in bm_edge.link_faces:
                                if f.index != curr_face_idx:
                                    neighbor_face = f.index
                                    break
                            if neighbor_face is not None:
                                poly_to = obj.data.polygons[neighbor_face]
                                normal_to_world = obj.matrix_world.to_3x3() @ poly_to.normal
                                normal_to_world.normalize()
                                if normal_to_world.dot(ray_dir) >= 0.0:
                                    # Neighbor face is back-facing (occluded), so we treat this edge as a silhouette boundary
                                    neighbor_face = None
                            candidates.append((t_val, p_edge, curr_face_idx, neighbor_face))
                            
        if len(candidates) > 0:
            valid_candidates = [c for c in candidates if c[0] > t_min + 1e-5]
            print(f"  Step {step}: Found {len(candidates)} total candidates, {len(valid_candidates)} valid candidates > t_min={t_min:.5f}")
            if not valid_candidates:
                print("  No valid candidates, breaking!")
                break
            best = min(valid_candidates, key=lambda x: x[0])
            t_int, p_int, f_from, f_to = best
            intersections.append((p_int, f_from, f_to, t_int))
            t_min = t_int
            print(f"  -> Moving to face {f_to} (t_min={t_min:.5f})")
            
            if f_to is None or f_to in visited:
                print(f"  -> Loop or boundary detected! f_to={f_to}, visited={visited}")
                break
            curr_face_idx = f_to
            visited.add(f_to)
        else:
            print("  No candidates found on this face, breaking!")
            break
            
    if curr_face_idx != face_end:
        face_end = None
    return intersections, face_start, face_end


def project_handle_to_face_plane(p_co, h_co, normal_world, region, rv3d):
    from bpy_extras.view3d_utils import region_2d_to_origin_3d, region_2d_to_vector_3d
    v = h_co - p_co
    v_len = v.length
    if v_len < 1e-4:
        # For zero or tiny handles, orthogonal projection is perfectly stable and correct
        proj = normal_world.dot(v) * normal_world
        return p_co + (v - proj)

    xy_h = getCoordFromLoc(region, rv3d, h_co)
    ray_org = region_2d_to_origin_3d(region, rv3d, xy_h)
    ray_dir = region_2d_to_vector_3d(region, rv3d, xy_h)
    denom = ray_dir.dot(normal_world)
    if abs(denom) > 0.2:
        u = (p_co - ray_org).dot(normal_world) / denom
        projected = ray_org + ray_dir * u
        if (projected - p_co).length > 3.0 * v_len:
            # Fallback if the view projection stretches the handle too much
            proj = normal_world.dot(v) * normal_world
            return p_co + (v - proj)
        return projected
    else:
        # Fallback to orthogonal projection
        proj = normal_world.dot(v) * normal_world
        return p_co + (v - proj)


def resolve_surface_crossovers(curvePts, obj, region, rv3d, surfaceMode, cached_bm=None):
    if len(curvePts) < 2:
        return curvePts

    is_cyclic = False
    if vectCmpWithMargin(curvePts[0][1], curvePts[-1][1]):
        is_cyclic = True

    import bmesh
    if cached_bm is not None:
        bm = cached_bm
    else:
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.faces.ensure_lookup_table()
        bm.edges.ensure_lookup_table()

    pts = [pt[:] for pt in curvePts]
    
    # Project original corner points visually on the surface
    from ..utils.view_utils import getFaceUnderMouse
    num_pts_to_proj = len(pts) - 1 if is_cyclic else len(pts)
    for idx in range(num_pts_to_proj):
        pt = pts[idx]
        p_co = pt[1]
        xy = getCoordFromLoc(region, rv3d, p_co)
        if xy is not None and xy.x != 9000:
            hit_loc, _, face_idx = getFaceUnderMouse(obj, region, rv3d, xy, 1000000)
            if face_idx is not None:
                orig_co = pt[1].copy()
                pt[1] = hit_loc.copy()
                if vectCmpWithMargin(pt[0], orig_co):
                    pt[0] = hit_loc.copy()
                if vectCmpWithMargin(pt[2], orig_co):
                    pt[2] = hit_loc.copy()
                if len(pt) > 5:
                    pt[5] = face_idx
                else:
                    pt.append(face_idx)
                    
    if is_cyclic:
        # Keep the left handle of pts[-1] (which is the left handle of the closing segment),
        # but synchronize the center, right handle, right handle type, and face index with pts[0].
        pts[-1][1] = pts[0][1]
        pts[-1][2] = pts[0][2]
        pts[-1][4] = pts[0][4]
        if len(pts[0]) > 5:
            if len(pts[-1]) > 5:
                pts[-1][5] = pts[0][5]
            else:
                pts[-1].append(pts[0][5])

    for pt in pts:
        if len(pt) > 4:
            if pt[3] != 'VECTOR':
                pt[3] = 'FREE'
            if pt[4] != 'VECTOR':
                pt[4] = 'FREE'

    if is_cyclic and not vectCmpWithMargin(pts[0][1], pts[-1][1]):
        pts.append(pts[0][:])

    new_pts = []
    num_segs = len(pts) - 1
    for i in range(num_segs):
        pt_curr = pts[i]
        pt_next = pts[i+1]
        
        p_start = pt_curr[1]
        h_start = pt_curr[2]
        h_end = pt_next[0]
        p_end = pt_next[1]
        
        face_start = pt_curr[5] if len(pt_curr) > 5 else None
        face_end = pt_next[5] if len(pt_next) > 5 else None
        
        intersections, face_start, face_end = get_surface_intersection_points(
            obj, p_start, h_start, h_end, p_end, face_start, face_end,
            region, rv3d, bm
        )
        
        start_idx = len(new_pts)
        new_pts.append(pt_curr)
        
        if len(intersections) > 0:
            for p_int, face_prev, face_next, t_val in intersections:
                new_pts.append([p_int, p_int, p_int, 'FREE', 'FREE', face_next])
            
            if len(pt_next) > 5:
                pt_next[5] = intersections[-1][2]
            else:
                pt_next.append(intersections[-1][2])
                
            t_vals = [0.0] + [item[3] for item in intersections] + [1.0]
            
            for sub_idx in range(len(intersections) + 1):
                idx_curr = start_idx + sub_idx
                pt_a = new_pts[idx_curr]
                pt_b = pt_next if sub_idx == len(intersections) else new_pts[idx_curr + 1]
                
                r0, r1, r2, r3 = get_bezier_sub_segment(
                    p_start, h_start, h_end, p_end,
                    t_vals[sub_idx], t_vals[sub_idx+1]
                )
                
                face_idx = face_start if sub_idx == 0 else pt_a[5]
                if face_idx is not None:
                    poly = obj.data.polygons[face_idx]
                    normal_world = obj.matrix_world.to_3x3() @ poly.normal
                    normal_world.normalize()
                    pt_a[2] = project_handle_to_face_plane(pt_a[1], r1, normal_world, region, rv3d)
                    pt_b[0] = project_handle_to_face_plane(pt_b[1], r2, normal_world, region, rv3d)
                else:
                    pt_a[2] = r1
                    pt_b[0] = r2
            
            if surfaceMode == 'SMOOTH':
                for idx in range(start_idx + 1, start_idx + len(intersections) + 1):
                    pt = new_pts[idx]
                    v_left = (pt[0] - pt[1]).normalized()
                    v_right = (pt[2] - pt[1]).normalized()
                    t_avg = (v_right - v_left).normalized()
                    l_left = (pt[0] - pt[1]).length
                    l_right = (pt[2] - pt[1]).length
                    pt[2] = pt[1] + t_avg * l_right
                    pt[0] = pt[1] - t_avg * l_left
        else:
            if len(pt_curr) > 5 and pt_curr[5] is not None and len(pt_next) > 5 and pt_next[5] is not None:
                face_idx = pt_curr[5]
                poly = obj.data.polygons[face_idx]
                normal_world = obj.matrix_world.to_3x3() @ poly.normal
                normal_world.normalize()
                pt_curr[2] = project_handle_to_face_plane(pt_curr[1], pt_curr[2], normal_world, region, rv3d)
                
                face_next_idx = pt_next[5]
                poly = obj.data.polygons[face_next_idx]
                normal_world = obj.matrix_world.to_3x3() @ poly.normal
                normal_world.normalize()
                pt_next[0] = project_handle_to_face_plane(pt_next[1], pt_next[0], normal_world, region, rv3d)
                
    new_pts.append(pts[-1])
    
    if is_cyclic:
        new_pts[0][0] = new_pts[-1][0]
        new_pts[0][3] = new_pts[-1][3]
        new_pts[0][4] = new_pts[-1][4]
        new_pts.pop()
        
    if cached_bm is None:
        bm.free()
        
    return new_pts
