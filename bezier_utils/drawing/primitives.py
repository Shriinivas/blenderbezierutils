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

    def newPoint(self, loc, ltype, rtype):
        self.curvePts.append([loc, loc, loc, ltype, rtype])

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

        if len(self.freeAxesN) > 2 and params.snapOrient in {
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
        import bpy
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
        from .math_fn import MathFnDraw

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

        snapOrigin = bpy.context.window_manager.bezierToolkitParams.snapOrigin
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
        metakeys = parent.snapper.getMetakeys()

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
        offset = params.drawStarOffset

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

        large_arc = 0
        rotation = 0

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

    def moveBezierPt(self, loc):
        if len(self.curvePts) > 0:
            pt = self.curvePts[-1][:]
            self.setCurvePt(-1, [loc, loc, loc, pt[3], pt[4]])

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
                        self.popCurvePt()

                # Because there is an extra point (the current one)
                if len(self.curvePts) <= 1:
                    self.initialize()
                    self.reset()
                else:
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc)
                self.capture = False
                self.parent.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        if not snapProc and (event.type == "RET" or event.type == "SPACE"):
            if event.value == "RELEASE":
                if snapper.digitsConfirmed:
                    self.reset()
                    snapper.digitsConfirmed = False
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.newPoint(loc, "ALIGNED", "ALIGNED")
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
                self.newPoint(loc, "ALIGNED", "ALIGNED")

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

                self.newPoint(loc, "ALIGNED", "ALIGNED")
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
                    self.moveBezierPt(loc)
                    hdlPtIdxs = {len(self.curvePts) - 2}

            self.parent.redrawBezier(rmInfo, hdlPtIdxs=hdlPtIdxs)
            return {"RUNNING_MODAL"}

        return {"PASS_THROUGH"} if not snapProc else {"RUNNING_MODAL"}

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
