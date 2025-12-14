import bpy
from mathutils import Vector, Matrix, kdtree, geometry
from ..core.props import FTProps
from bpy_extras.view3d_utils import region_2d_to_location_3d
from ..constants import (
    LARGE_NO,
    LARGE_VECT,
    EVT_NOT_CONS,
    EVT_CONS,
    DEF_ERR_MARGIN,
    unitMap,
    EVT_META_OR_SNAP,
)
from ..utils.math_utils import floatCmpWithMargin
from ..utils.view_utils import (
    getPtProjOnPlane,
    getLineTransMatrices,
    getSelFaceLoc,
    updateMetaBtns,
    showSnapToPlane,
    getUnit,
    getUnitScale,
    getCoordFromLoc,
    getViewDistRounding,
    roundedVect,
)
from .hotkeys import FTHotKeys

from math import sqrt, cos, sin, radians, tan, atan, atan2, degrees


class SnapDigits:
    digitMap = {
        "ZERO": "0",
        "ONE": "1",
        "TWO": "2",
        "THREE": "3",
        "FOUR": "4",
        "FIVE": "5",
        "SIX": "6",
        "SEVEN": "7",
        "EIGHT": "8",
        "NINE": "9",
        "PERIOD": ".",
    }
    numpadDigitMap = {
        "NUMPAD_0": "0",
        "NUMPAD_1": "1",
        "NUMPAD_2": "2",
        "NUMPAD_3": "3",
        "NUMPAD_4": "4",
        "NUMPAD_5": "5",
        "NUMPAD_6": "6",
        "NUMPAD_7": "7",
        "NUMPAD_8": "8",
        "NUMPAD_9": "9",
        "NUMPAD_PERIOD": ".",
    }

    def getValidFloat(sign, digits):
        delta = 0
        valid = True
        for i in range(len(digits), 0, -1):
            td = digits[0:i]
            try:
                delta = float(sign + "".join(td))
                return delta, valid
            except:
                valid = False
        return delta, valid

    def __init__(self, getFreeAxes, getEditCoPair):
        self.getFreeAxes = getFreeAxes
        self.getEditCoPair = getEditCoPair
        self.initialize()

    def initialize(self):
        self.axes = None  # [0, 1, 2] etc.
        self.deltaVec = None  # Always cartesian
        self.axisIdx = 0
        self.digitChars = []
        self.signChar = ""
        self.polar = False
        self.pDataIdx = 0  # Corresponding to axisIdx

    def hasVal(self):
        return self.axes is not None

    # Theta in degrees
    def getPolarCos(self):
        axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
        val0, val1 = self.deltaVec[axis0], self.deltaVec[axis1]
        if val0 == 0:
            if val1 == 0:
                return 0, 0
            theta = (val1 / abs(val1)) * 90
        else:
            theta = degrees(atan(abs(val1) / abs(val0)))
            if val0 < 0:
                theta = 180 - theta
            if val1 < 0:
                theta = -theta
        r = sqrt(val1 * val1 + val0 * val0)
        return [r, theta]

    # Theta in degrees
    def getDeltaFromPolar(self, polCos):
        r, theta = polCos
        return r * cos(radians(theta)), r * sin(radians(theta))

    def addToPolar(self, delta):
        polCos = list(self.getPolarCos())
        polCos[self.pDataIdx] += delta
        val0, val1 = self.getDeltaFromPolar(polCos)
        return val0, val1

    def digitsToVec(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        if not self.polar or self.pDataIdx == 0:
            delta /= getUnitScale()
        if self.polar:
            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
            self.deltaVec[axis0], self.deltaVec[axis1] = self.addToPolar(delta)
        else:
            self.deltaVec[self.axisIdx] += delta

    def vecToDigits(self):
        if self.polar:
            vals = self.getPolarCos()
            v = round(vals[self.pDataIdx], 4)
        else:
            v = round(self.deltaVec[self.axisIdx], 4)
        self.digitChars = list(str(abs(v))) if v != 0 else ""
        while len(self.digitChars) > 0 and self.digitChars[-1] == "0":
            self.digitChars.pop()
        if len(self.digitChars) > 0 and self.digitChars[-1] == ".":
            self.digitChars.pop()
        self.signChar = "-" if v < 0 else ""

    def updateDeltaFromEditCos(self):
        editCos = self.getEditCoPair()
        self.deltaVec = editCos[1] - editCos[0]

    def procEvent(self, context, event, metakeys):
        if FTHotKeys.isHotKey(FTHotKeys.hkTweakPos, event.type, metakeys):
            editCos = self.getEditCoPair()
            if len(editCos) == 2:
                if event.value == "RELEASE":
                    if self.axes is None:
                        self.axes = self.getFreeAxes()  # TODO: Always dynamic?
                        self.polar = False
                    # Only possible within a plane right now
                    elif len(self.getFreeAxes()) == 2:
                        self.polar = not self.polar
                        self.pDataIdx = 0
                    self.updateDeltaFromEditCos()
                    self.digitChars = []
                    self.axisIdx = self.axes[0]
            return True

        dmap = self.digitMap.copy()
        if FTProps.numpadEntry:
            dmap.update(self.numpadDigitMap)
        dval = dmap.get(event.type)
        if dval is not None:
            if event.value == "RELEASE":
                self.digitChars.append(dval)
                if self.axes is None:
                    self.axes = self.getFreeAxes()  # TODO: Always dynamic?
                    self.deltaVec = Vector()
                    self.axisIdx = self.axes[0]
            return True

        if not self.hasVal():  # No further processing if nothing to process
            return False

        retVal = True
        if event.type == "MINUS":
            if event.value == "RELEASE":
                self.signChar = "" if (self.signChar == "-") else "-"

        elif event.type == "ESC":
            if event.value == "RELEASE":
                self.initialize()

        elif event.type == "BACK_SPACE":
            if event.value == "RELEASE":
                if len(self.digitChars) > 0:
                    self.digitChars.pop()
                else:
                    if self.polar:
                        polarCos = self.getPolarCos()
                        if not floatCmpWithMargin(polarCos[self.pDataIdx], 0):
                            self.vecToDigits()

                            # TODO: Retain theta without this workaround
                            polarCos[self.pDataIdx] = 0.00001
                            delta = self.getDeltaFromPolar(polarCos)
                            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
                            self.deltaVec[axis0] = delta[0]
                            self.deltaVec[axis1] = delta[1]
                    else:
                        if self.deltaVec[self.axisIdx] != 0:
                            self.vecToDigits()
                            self.deltaVec[self.axisIdx] = 0

        elif event.type == "TAB":
            if event.value == "RELEASE":
                self.digitsToVec()
                if self.polar:
                    self.pDataIdx = (self.pDataIdx + 1) % 2
                else:
                    self.axisIdx = self.axes[
                        (self.axes.index(self.axisIdx) + 1) % len(self.axes)
                    ]

                self.digitChars = []
                self.signChar = ""
        else:
            retVal = False

        return retVal

    def getCurrDelta(self):
        if self.axes is None:
            return Vector()

        val = self.deltaVec.copy()
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        if not self.polar or self.pDataIdx == 0:
            delta /= getUnitScale()
        if self.polar:
            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
            val[axis0], val[axis1] = self.addToPolar(delta)
        else:
            val[self.axisIdx] += delta

        return val

    def getDeltaStrPolar(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        polCos = self.getPolarCos()
        polCos[0] *= getUnitScale()
        strs = ["["] * 2
        idx0 = self.pDataIdx
        idx1 = 1 - self.pDataIdx
        d = polCos[idx0]

        if d != 0:
            strs[idx0] += str(round(d, 4)) + ("+" if (self.signChar == "") else "")

        strs[idx0] += (
            self.signChar
            + "".join(self.digitChars)
            + "] = "
            + str(round((delta + d), 4))
            + ("" if valid else " <Invalid>")
        )

        strs[idx1] = "[" + str(round(polCos[1 - self.pDataIdx], 4)) + "]"

        return "r: " + strs[0] + " theta: " + strs[1]

    def getCurrDeltaStr(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        retStr = "["
        d = self.deltaVec[self.axisIdx] * getUnitScale()
        if d != 0:
            retStr += str(round(d, 4)) + ("+" if (self.signChar == "") else "")
        retStr += (
            self.signChar
            + "".join(self.digitChars)
            + "] = "
            + str(round((delta + d), 4))
            + ("" if valid else " <Invalid>")
        )

        return retStr


class CustomAxis:
    def __init__(self):
        # TODO: What's better?
        if bpy.data.scenes[0].get("btk_co1") is None:
            bpy.data.scenes[0]["btk_co1"] = LARGE_VECT
        if bpy.data.scenes[0].get("btk_co2") is None:
            bpy.data.scenes[0]["btk_co2"] = LARGE_VECT
        self.axisPts = [
            Vector(bpy.data.scenes[0]["btk_co1"]),
            Vector(bpy.data.scenes[0]["btk_co2"]),
        ]
        if bpy.data.scenes[0].get("btk_snapPtCnt") is None:
            bpy.data.scenes[0]["btk_snapPtCnt"] = 9  # 9 intervals = 11 points (10 divisions)
        self.snapCnt = bpy.data.scenes[0]["btk_snapPtCnt"]
        self.inDrawAxis = False  # User drawing the custom axis

    def length(self):
        pts = self.axisPts

        # Strange floating points!
        if all(pt < (LARGE_NO - 1000) for pt in pts[0] + pts[1]):
            return (pts[1] - pts[0]).length
        else:
            return 0

    def set(self, idx, co):
        self.axisPts[idx] = co
        if idx == 0:
            bpy.data.scenes[0]["btk_co1"] = [c for c in co]
        if idx == 1:
            bpy.data.scenes[0]["btk_co2"] = [c for c in co]
        bpy.data.scenes[0]["btk_snapPtCnt"] = self.snapCnt

    def getAngleFromGlobalX(self):
        """Get angle of custom axis from global X axis in degrees"""
        if self.length() == 0:
            return None, None

        vec = (self.axisPts[1] - self.axisPts[0]).normalized()

        # 2D angle (XY plane projection)
        vec_xy = Vector((vec.x, vec.y, 0))
        if vec_xy.length > 0.0001:
            vec_xy.normalize()
            angle_xy = atan2(vec_xy.y, vec_xy.x)
        else:
            angle_xy = 0

        # 3D angle from global X axis
        global_x = Vector((1, 0, 0))
        angle_3d = vec.angle(global_x)

        return angle_xy, angle_3d

    def getAngleString(self):
        """Get formatted angle string for display"""
        angle_xy, angle_3d = self.getAngleFromGlobalX()
        if angle_xy is None:
            return "Not defined"

        # Convert to degrees (atan2 returns radians, angle() returns radians)
        from math import degrees
        angle_xy_deg = degrees(angle_xy)
        angle_3d_deg = degrees(angle_3d)

        # Format: Show XY angle primarily, 3D angle in parentheses if different
        if abs(angle_3d_deg - abs(angle_xy_deg)) < 1.0:
            # Nearly planar in XY
            return f"{angle_xy_deg:.1f}°"
        else:
            # Significant Z component
            return f"{angle_xy_deg:.1f}° (3D: {angle_3d_deg:.1f}°)"

    def getSnapPts(self):  # ptCnt excluding end points
        pts = self.axisPts
        if self.length() == 0:
            return [pts[0], pts[1]]

        snapPts = [pts[0]]
        interval = self.snapCnt + 1
        incr = self.length() / interval
        diffV = pts[1] - pts[0]

        for i in range(1, interval):
            snapPts.append(pts[0] + (diffV * (incr * i) / self.length()))

        snapPts.append(pts[1])
        return snapPts

    def procDrawEvent(self, context, event, snapper, rmInfo):
        if event.type == "RIGHTMOUSE":
            params = bpy.context.window_manager.bezierToolkitParams
            snapOrigin = params.snapOrigin
            snapOrient = params.snapOrient
            if event.value == "RELEASE" and (snapOrigin == "AXIS" or snapOrient == "AXIS"):
                loc = snapper.get3dLocSnap(
                    rmInfo, SnapParams(snapper, snapToAxisLine=False)
                )
                if not self.inDrawAxis:
                    self.set(0, loc)
                else:
                    self.set(1, loc)
                self.inDrawAxis = not self.inDrawAxis
            return True

        if self.inDrawAxis:
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
                    return True
                elif event.type == "WHEELUPMOUSE" or event.type.endswith("PLUS"):
                    if self.snapCnt < 20:
                        self.snapCnt += 1
                elif event.type == "WHEELDOWNMOUSE" or event.type.endswith("MINUS"):
                    if self.snapCnt > 0:
                        self.snapCnt -= 1

            # Update axis point on:
            # - Mouse move
            # - Numeric input changes
            # - Snap key press/release (Ctrl for grid, Alt for vert, Shift for angle)
            needsUpdate = (
                event.type == "MOUSEMOVE"
                or snapper.snapDigits.hasVal()
                or event.type in {
                    "LEFT_CTRL", "RIGHT_CTRL",
                    "LEFT_ALT", "RIGHT_ALT",
                    "LEFT_SHIFT", "RIGHT_SHIFT"
                }
            )
            if needsUpdate:
                if snapper.snapDigits.hasVal():
                    # Numeric input: apply delta directly from first axis point
                    delta = snapper.snapDigits.getCurrDelta()
                    loc = self.axisPts[0] + delta
                else:
                    # Mouse move or snap key: use get3dLocSnap for snapping
                    loc = snapper.get3dLocSnap(
                        rmInfo, SnapParams(snapper, snapToAxisLine=False)
                    )
                self.set(1, loc)

            # Confirm numeric input with Return/Space - finalize axis
            if event.type in {"RET", "SPACE"} and snapper.snapDigits.hasVal():
                if event.value == "RELEASE":
                    # Apply the numeric input one final time before finalizing
                    delta = snapper.snapDigits.getCurrDelta()
                    loc = self.axisPts[0] + delta
                    self.set(1, loc)
                    # Finalize the axis
                    self.inDrawAxis = False
                    # Reset snapper for next operation
                    snapper.resetSnap()
                return True

            if event.type == "ESC":
                self.set(0, LARGE_VECT)
                self.set(1, LARGE_VECT)
                self.inDrawAxis = False
                snapper.resetSnap()  # Clear any numeric input
                return True

            return True

        return False


# TODO: Make independent of snapper
class SnapParams:
    def __init__(
        self,
        snapper,
        xyDelta=[0, 0],
        vec=None,
        refreshStatus=True,
        snapToAxisLine=True,
        lastCo1Axis=False,
        enableSnap=True,
        vertSnap=None,
        gridSnap=None,
        angleSnap=None,
        refLine=None,
        refLineOrig=None,
        selCo=None,
        inEdit=None,
        hasSel=None,
        transType=None,
        origType=None,
        offsetRef=None,
        axisScale=None,
        freeAxesN=None,
        dispAxes=True,
        snapToPlane=None,
    ):
        self.xyDelta = xyDelta
        self.vec = vec
        self.refreshStatus = refreshStatus
        self.snapToAxisLine = snapToAxisLine
        self.lastCo1Axis = lastCo1Axis

        if enableSnap:
            self.vertSnap = snapper.vertSnap if (vertSnap is None) else vertSnap
            self.gridSnap = snapper.gridSnap if (gridSnap is None) else gridSnap
            self.angleSnap = snapper.angleSnap if (angleSnap is None) else angleSnap
        else:
            self.vertSnap = False
            self.gridSnap = False
            self.angleSnap = False

        self.refLine = snapper.getRefLine() if (refLine is None) else refLine
        self.refLineOrig = (
            snapper.getRefLineOrig() if (refLineOrig is None) else refLineOrig
        )
        self.selCo = snapper.getSelCo() if (selCo is None) else selCo

        self.inEdit = snapper.isEditing() if (inEdit is None) else inEdit
        self.hasSel = snapper.hasSelection() if (hasSel is None) else hasSel

        params = bpy.context.window_manager.bezierToolkitParams
        self.transType = params.snapOrient if (transType is None) else transType
        self.origType = params.snapOrigin if (origType is None) else origType
        self.offsetRef = params.offsetRef if (offsetRef is None) else offsetRef
        self.axisScale = params.axisScale if (axisScale is None) else axisScale

        self.freeAxesN = (
            snapper.getFreeAxesNormalized() if (freeAxesN is None) else freeAxesN
        )

        self.dispAxes = dispAxes

        if snapToPlane is not None:
            self.snapToPlane = snapToPlane
        else:
            if showSnapToPlane(params):  # Only if Snap to Plane option is visible
                self.snapToPlane = params.snapToPlane
            else:
                self.snapToPlane = False


class Snapper:
    DEFAULT_ANGLE_SNAP_STEPS = 3
    MAX_SNAP_VERT_CNT = 1000
    MAX_SNAP_FACE_CNT = 1000

    def __init__(
        self,
        context,
        getSnapLocs,
        getRefLine,
        getRefLineOrig,
        getSelCo,
        getCurrLine,
        hasSelection,
        isEditing,
    ):
        self.getSnapLocs = getSnapLocs
        self.getRefLine = getRefLine
        self.getCurrLine = getCurrLine
        self.getRefLineOrig = getRefLineOrig
        self.getSelCo = getSelCo
        self.hasSelection = hasSelection
        self.isEditing = isEditing
        self.angleSnapSteps = Snapper.DEFAULT_ANGLE_SNAP_STEPS
        self.customAxis = CustomAxis()
        self.snapDigits = SnapDigits(self.getSnapParamsFreeAxes, self.getEditCoPair)
        self.initialize()
        self.updateSnapKeyMap()
        # ~ bpy.context.space_data.overlay.show_axis_x = False
        # ~ bpy.context.space_data.overlay.show_axis_y = False

    def resetMetakeys(self):
        # metakeys for all functions
        self.shift = False
        self.ctrl = False
        self.alt = False

    def resetSnapKeys(self):
        self.angleSnap = False
        self.gridSnap = False
        self.vertSnap = False

    def initialize(self):
        self.resetMetakeys()
        self.resetSnapKeys()

        self.tm = None
        self.orig = None
        self.snapCo = None
        self.freezeOrient = False

        self.lastSnapTypes = set()

        self.resetSnap()

    def updateSnapKeyMap(self):
        self.snapKeyMap = {}
        kIds = [FTHotKeys.hkSnapVert, FTHotKeys.hkSnapAngle, FTHotKeys.hkSnapGrid]
        varMap = {
            FTHotKeys.hkSnapVert: "vertSnap",
            FTHotKeys.hkSnapAngle: "angleSnap",
            FTHotKeys.hkSnapGrid: "gridSnap",
        }
        for kId in kIds:
            keys = FTHotKeys.getSnapHotKeys(kId)
            for key in keys:
                self.snapKeyMap[key] = varMap[kId]

    def resetSnap(self):  # Called even during isEditing
        self.freeAxes = []  # All axes free
        self.snapDigits.initialize()
        self.rmInfo = None
        self.snapParams = None

        # This variable lets caller know that return was pressed after digits were entered
        # Caller can reset snapper as per convenience
        self.digitsConfirmed = False
        self.lastSelCo = None

        # ~ self.snapStack = [] # TODO

    # Return user selection in header menu, [] for None
    # Values: "X","Y","Z" for axis; "shift-Z","shift-Y","shift-X" for planes
    def getFreeAxesGlobal(self):
        constrAxes = bpy.context.window_manager.bezierToolkitParams.constrAxes
        # Empty string can occur during dynamic enum transitions
        if not constrAxes or constrAxes == "NONE":
            return []
        idx = constrAxes.find("-")
        axis = constrAxes[idx + 1]
        freeAxes = [ord(axis) - ord("X")]
        if idx > 0:
            freeAxes = sorted(list({0, 1, 2} - set(freeAxes)))
        return freeAxes

    # Return actual free axes (menu + hotkey combined), [] for None
    def getFreeAxesCombined(self):
        constrAxes = self.getFreeAxesGlobal()
        if len(constrAxes) == 0:
            return self.freeAxes
        if len(self.freeAxes) == 0:
            return constrAxes
        freeAxes = sorted(list(set(constrAxes).intersection(set(self.freeAxes))))
        if len(freeAxes) == 0:
            return constrAxes
        else:
            return freeAxes

    # Return actual free axes (menu + hotkey combined), [0, 1, 2] for None
    def getFreeAxesNormalized(self):
        if len(self.getFreeAxesCombined()) == 0:
            return [0, 1, 2]
        else:
            return self.getFreeAxesCombined()

    def getSnapParamsFreeAxes(self):
        if self.snapParams is not None:
            return self.snapParams.freeAxesN
        else:
            return self.getFreeAxesNormalized()

    def getCurrOrig(self, rmInfo, obj, origType, refLineOrig, selCo):
        """Get the pivot point location based on origType (pivot point setting)."""
        if origType == "AXIS":
            if self.customAxis.length() != 0:
                return self.customAxis.axisPts[0]
        elif origType == "OBJECT" and obj is not None:
            return obj.location
        elif origType == "FACE" and rmInfo is not None:
            selObj, location, normal, faceIdx = getSelFaceLoc(
                rmInfo.region, rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT
            )
            if faceIdx is not None:
                return selObj.matrix_world @ selObj.data.polygons[faceIdx].center
        elif origType == "CURSOR":
            return bpy.context.scene.cursor.location
        # Default: GLOBAL origin
        return Vector((0, 0, 0))

    def getOffsetRefPoint(self, rmInfo, obj, origType, offsetRef, refLineOrig, selCo):
        """Get the offset reference point for numeric input and angle snapping.

        Args:
            offsetRef: 'PIVOT' to use pivot point, 'PREVIOUS' to use previous point
        Returns:
            Vector: The point from which offsets should be calculated
        """
        if offsetRef == "PREVIOUS":
            # Use previous point (refLineOrig) if available
            if refLineOrig is not None:
                return refLineOrig
            # Fall back to pivot if no previous point
        # Default: use pivot point
        return self.getCurrOrig(rmInfo, obj, origType, refLineOrig, selCo)

    def getTransMatsForOrient(self, rmInfo, obj, transType, axisScale):
        custAxis = self.customAxis
        if abs(custAxis.length()) <= DEF_ERR_MARGIN:
            custAxis = None

        refLine = self.getRefLine()
        currLine = self.getCurrLine()
        if refLine is not None and len(refLine) < 2:
            refLine = None
        if currLine is not None and len(currLine) < 2:
            currLine = None

        tmScale = Matrix()
        if custAxis is not None and axisScale == "AXIS":
            unitD = custAxis.length() / (custAxis.snapCnt + 1)
            tmScale = Matrix.Scale(1 / unitD, 4)

        elif refLine is not None and axisScale == "REFERENCE":
            unitD = (refLine[1] - refLine[0]).length / 10
            if unitD > DEF_ERR_MARGIN:
                tmScale = Matrix.Scale(1 / unitD, 4)

        tm = None
        if transType == "AXIS" and custAxis is not None:
            tm, invTm = getLineTransMatrices(custAxis.axisPts[0], custAxis.axisPts[1])

        elif transType == "REFERENCE" and refLine is not None:
            tm, invTm = getLineTransMatrices(refLine[0], refLine[1])

        elif transType == "CURR_POS" and currLine is not None:
            tm, invTm = getLineTransMatrices(currLine[0], currLine[1])

        elif transType == "VIEW":
            tm = rmInfo.rv3d.view_matrix

        elif obj is not None and transType == "OBJECT":
            tm = (obj.matrix_world).inverted_safe()

        elif transType == "FACE":
            selObj, location, normal, faceIdx = getSelFaceLoc(
                rmInfo.region, rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT
            )
            if faceIdx is not None:
                normal = selObj.data.polygons[faceIdx].normal
                quat = normal.to_track_quat("Z", "X").to_matrix().to_4x4()
                tm = (selObj.matrix_world @ quat).inverted_safe()

        if tm is not None:
            trans, quat, scale = tm.decompose()
            tm = quat.to_matrix().to_4x4() @ tmScale
        else:
            tm = tmScale

        return tm, tm.inverted_safe()

    def isLocked(self):
        return len(self.freeAxes) > 0 or (
            self.snapDigits.hasVal() and not self.digitsConfirmed
        )

    def getMetakeys(self):
        return [self.alt, self.ctrl, self.shift]

    # To be called in modal method of parent
    def procEvent(self, context, event):
        # update ctrl etc.
        retValMeta = updateMetaBtns(self, event)

        retValSnap = updateMetaBtns(self, event, self.snapKeyMap)

        if retValMeta or retValSnap:
            return EVT_META_OR_SNAP

        metakeys = self.getMetakeys()
        refLineOrig = (
            self.snapParams.refLineOrig
            if self.snapParams is not None
            else self.getRefLineOrig()
        )

        # Allow numeric input when:
        # 1. There's a reference line origin (normal bezier drawing), or
        # 2. Custom axis is being drawn (use axis start point as reference)
        if refLineOrig is not None or self.customAxis.inDrawAxis:
            snapDProc = self.snapDigits.procEvent(context, event, metakeys)
            if snapDProc:
                self.digitsConfirmed = (
                    False  # Always reset if there was any digit entered
                )
                return EVT_CONS

        if FTHotKeys.isHotKey(FTHotKeys.hkReorient, event.type, metakeys):
            if event.value == "RELEASE":
                self.freezeOrient = not self.freezeOrient
            return EVT_CONS

        if not self.ctrl and event.type in {"X", "Y", "Z"}:
            self.digitsConfirmed = False  # Always reset if there is any lock axis
            if event.value == "RELEASE":
                self.freeAxes = [ord(event.type) - ord("X")]
                if self.shift:
                    self.freeAxes = sorted(list({0, 1, 2} - set(self.freeAxes)))

                # if already part of global axes, don't store (no escape needed)
                if self.getFreeAxesCombined() == self.getFreeAxesGlobal():
                    self.freeAxes = []

            return EVT_CONS

        retVal = EVT_NOT_CONS
        # Consume escape or return / space only if there's something to process
        if self.isLocked():
            retVal = EVT_CONS
            if event.type == "RET" or event.type == "SPACE":
                if event.value == "RELEASE":
                    # ~ self.resetSnap() # This is the responsibility of the caller
                    self.digitsConfirmed = True  # Confirm first time
            elif event.type == "ESC":
                if event.value == "RELEASE":
                    self.resetSnap()
            else:
                retVal = EVT_NOT_CONS

        return retVal

    def getStatusStr(self, unit, invTm, refPt, newPt):
        manualEntry = self.snapDigits.hasVal()

        if manualEntry and self.snapDigits.polar:
            return self.snapDigits.getDeltaStrPolar()

        axes = self.getSnapParamsFreeAxes()
        diffV = self.snapDigits.getCurrDelta() if manualEntry else (newPt - refPt)

        diffV *= getUnitScale()

        diffVActual = invTm @ diffV

        retStr = ""
        transformed = invTm != Matrix()

        axisDeltaFormat = "D{axis}: {axisDelta}"
        axisDiffFormat = "{{{axisDiff}}}  "

        for i, d in enumerate(diffV):
            if i not in axes:
                continue
            v1 = chr(ord("x") + i)

            if manualEntry and i == self.snapDigits.axisIdx:
                v2 = self.snapDigits.getCurrDeltaStr()
            else:
                v2 = str(round(d, 4))

            retStr += axisDeltaFormat.format(axis=v1, axisDelta=v2)

            if transformed:
                v3 = str(round(diffVActual[i], 4))
                retStr += axisDiffFormat.format(axisDiff=v3)

            retStr += "  "

        unitT = ""
        unitA = ""
        if transformed:
            unitA = unit
        else:
            unitT = unit

        totalDeltaFormat = "({totalDelta}{unit})"
        diffVStr = str(round(diffV.length, 4))
        retStr += totalDeltaFormat.format(totalDelta=diffVStr, unit=unitT)

        totalDiffVFormat = "{{{totalDiffV}{unit}}}"
        if transformed:
            diffVStr = str(round(diffVActual.length, 4))
            retStr += totalDiffVFormat.format(totalDiffV=diffVStr, unit=unitA)

        return retStr

    def getAllSnapLocs(self, snapToAxisLine):
        snapLocs = self.getSnapLocs()
        snapLocs.append(bpy.context.scene.cursor.location)
        snapLocs.append(Vector((0, 0, 0)))

        if snapToAxisLine:
            snapLocs += self.customAxis.getSnapPts()

        vertCnt = 0
        aos = [bpy.context.object] if bpy.context.object is not None else []
        objs = bpy.context.selected_objects + aos
        for obj in objs:
            snapLocs.append(obj.location)
            if obj.type == "MESH":
                if vertCnt + len(obj.data.vertices) < self.MAX_SNAP_VERT_CNT:
                    snapLocs += [obj.matrix_world @ v.co for v in obj.data.vertices]
                    vertCnt = +len(obj.data.vertices)
                else:
                    break
        return snapLocs

    def getTMInfoAndOrig(
        self, rmInfo, transType, origType, freezeOrient, axisScale, refLineOrig, selCo
    ):
        obj = bpy.context.object

        if self.tm is not None and self.freezeOrient and transType == "FACE":
            tm, invTm = self.tm, self.tm.inverted_safe()
        else:
            tm, invTm = self.getTransMatsForOrient(rmInfo, obj, transType, axisScale)

        if self.orig is not None and freezeOrient and origType == "FACE":
            orig = self.orig
        else:
            orig = self.getCurrOrig(rmInfo, obj, origType, refLineOrig, selCo)

        return tm, invTm, orig

    def get3dLocSnap(self, rmInfo, snapParams=None):
        if snapParams is None:
            snapParams = SnapParams(self)

        self.rmInfo = rmInfo
        self.snapParams = snapParams

        refreshStatus = snapParams.refreshStatus
        snapToAxisLine = snapParams.snapToAxisLine
        lastCo1Axis = snapParams.lastCo1Axis

        vertSnap = snapParams.vertSnap
        gridSnap = snapParams.gridSnap
        angleSnap = snapParams.angleSnap

        refLine = snapParams.refLine
        refLineOrig = snapParams.refLineOrig
        selCo = snapParams.selCo

        inEdit = snapParams.inEdit
        hasSel = snapParams.hasSel

        transType = snapParams.transType
        origType = snapParams.origType
        offsetRef = snapParams.offsetRef
        axisScale = snapParams.axisScale

        vec = snapParams.vec
        freeAxesN = snapParams.freeAxesN
        snapToPlane = snapParams.snapToPlane

        self.snapCo = None
        region = rmInfo.region
        rv3d = rmInfo.rv3d

        xy = [
            rmInfo.xy[0] - snapParams.xyDelta[0],
            rmInfo.xy[1] - snapParams.xyDelta[1],
        ]

        loc = None

        tm, invTm, orig = self.getTMInfoAndOrig(
            rmInfo,
            transType,
            origType,
            self.freezeOrient,
            axisScale,
            refLineOrig,
            selCo,
        )

        # Must be done after the call to getTMInfoAndOrig
        if hasSel:
            self.freezeOrient = True

        vec = snapParams.vec if (snapParams.vec is not None) else orig

        self.lastSnapTypes = set()

        unit = unitMap.get(getUnit())
        if unit is None:
            unit = ""

        digitsValid = True
        # ~ freeAxesC = self.getFreeAxesCombined()
        # ~ freeAxesN = self.getFreeAxesNormalized()
        # ~ freeAxesG = self.getFreeAxesGlobal()

        if FTProps.dispSnapInd or vertSnap:
            # TODO: Called very frequently (store the tree [without duplicating data])
            snapLocs = self.getAllSnapLocs(
                (snapToAxisLine and "AXIS" in {transType, origType, axisScale})
            ) + [orig]

            kd = kdtree.KDTree(len(snapLocs))
            for i, l in enumerate(snapLocs):
                kd.insert(getCoordFromLoc(region, rv3d, l).to_3d(), i)
            kd.balance()

            coFind = Vector(xy).to_3d()
            searchResult = kd.find_range(coFind, FTProps.snapDist)

            if len(searchResult) != 0:
                co, idx, dist = min(searchResult, key=lambda x: x[2])
                self.snapCo = snapLocs[idx]

        if vertSnap:
            if self.snapCo is not None:
                loc = self.snapCo
            else:
                selObj, loc, normal, faceIdx = getSelFaceLoc(
                    region, rv3d, xy, self.MAX_SNAP_FACE_CNT, checkEdge=True
                )

        if loc is not None:
            loc = tm @ loc
            self.lastSnapTypes.add("loc")
        else:
            loc = region_2d_to_location_3d(region, rv3d, xy, vec)
            loc = tm @ loc

            # TODO: Get gridSnap and angleSnap out of this if
            # Apply constraint when: axis constraints active (freeAxesN < 3),
            # or other snap modes enabled
            if (
                (transType != "GLOBAL" and inEdit)
                or snapToPlane
                or gridSnap
                or self.snapDigits.hasVal()
                or len(freeAxesN) < 3  # Axis constraint - apply even for first point
                or (inEdit and angleSnap)
            ):
                # snapToPlane means global constrain axes selection is a plane
                # refCo is used for axis constraint calculations
                refCo = tm @ orig

                # offsetRefPt is used for numeric input delta calculations
                offsetRefPt = self.getOffsetRefPoint(
                    rmInfo, bpy.context.object, origType, offsetRef, refLineOrig, selCo
                )

                if self.snapDigits.hasVal():
                    delta = self.snapDigits.getCurrDelta()
                    # Apply delta from offset reference point, not pivot
                    loc = tm @ offsetRefPt + delta
                    self.lastSnapTypes.add("keyboard")
                else:
                    # Special condition for lock to single axis
                    # ~ if(len(freeAxesN) == 1 and refLineOrig != None):
                    # ~ refCo = tm @ refLineOrig
                    if len(freeAxesN) == 2:
                        constrAxes = freeAxesN
                        loc = refCo.copy()
                        # Any other two points on the plane
                        ppt1 = loc.copy()
                        ppt1[constrAxes[0]] = loc[constrAxes[0]] + 10
                        ppt2 = loc.copy()
                        ppt2[constrAxes[1]] = loc[constrAxes[1]] + 10
                        ppt3 = ppt2.copy()
                        ppt2[constrAxes[0]] = loc[constrAxes[0]] + 10

                        # Raycast from 2d point onto the plane
                        pt = getPtProjOnPlane(
                            region,
                            rv3d,
                            xy,
                            invTm @ loc,
                            invTm @ ppt1,
                            invTm @ ppt2,
                            invTm @ ppt3,
                        )

                        for axis in constrAxes:
                            # TODO: Better handling of boundary condition
                            if pt is None or pt[axis] > 1000:
                                loc[axis] = refCo[axis]
                            else:
                                loc[axis] = (tm @ pt)[axis]
                        self.lastSnapTypes.add("axis2")

                    if len(freeAxesN) == 1:
                        if lastCo1Axis:
                            refCo = self.lastSelCo  # TODO: More testing
                        axis = freeAxesN[0]
                        # Any one point on axis
                        ptOnAxis = refCo.copy()

                        # Convert everything to 2d
                        lastCo2d = getCoordFromLoc(region, rv3d, invTm @ refCo)

                        # Very small distance so that the point is on viewport
                        # TODO: This is not foolproof :(
                        ptOnAxis[axis] += 0.01
                        ptOnAxis2d = getCoordFromLoc(region, rv3d, invTm @ ptOnAxis)

                        # Find 2d projection (needed)
                        pt2d = geometry.intersect_point_line(xy, lastCo2d, ptOnAxis2d)[
                            0
                        ]
                        # Any other 2 points on the plane, on which the axis lies
                        ppt1 = refCo.copy()
                        ppt1[axis] += 10
                        newAxis = [i for i in range(0, 3) if i != axis][0]
                        ppt2 = refCo.copy()
                        ppt2[newAxis] += 10

                        # Raycast from 2d point onto the plane
                        pt = getPtProjOnPlane(
                            region,
                            rv3d,
                            pt2d,
                            invTm @ refCo,
                            invTm @ ppt1,
                            invTm @ ppt2,
                        )
                        loc = refCo.copy()
                        if pt is None or pt[axis] > 1000:
                            loc[axis] = refCo[axis]
                        else:
                            loc[axis] = (tm @ pt)[axis]
                        self.lastSnapTypes.add("axis1")

                if not self.snapDigits.hasVal() and gridSnap:
                    if axisScale in {"AXIS", "REFERENCE"}:
                        # Independent of view distance
                        diffV = loc - refCo
                        if diffV.length > DEF_ERR_MARGIN:  # Avoid division by zero
                            loc = refCo + round(diffV.length) * (diffV / diffV.length)
                    elif transType == "GLOBAL" and len(freeAxesN) == 2:
                        # Fixed 1.0 unit grid snapping for Global constraint plane
                        # This matches the visual grid spacing
                        gridSpacing = 1.0
                        locInTm = invTm @ loc
                        for axis in freeAxesN:
                            locInTm[axis] = round(locInTm[axis] / gridSpacing) * gridSpacing
                        loc = tm @ locInTm
                    else:
                        rounding = getViewDistRounding(rmInfo.space3d, rv3d)
                        loc = tm @ roundedVect(
                            rmInfo.space3d, invTm @ loc, rounding, freeAxesN
                        )
                    self.lastSnapTypes.add("grid")

                if not self.snapDigits.hasVal() and angleSnap and len(refLine) > 0:
                    # ~ freeAxesC = [0, 1, 2] if len(freeAxesC) == 0 else freeAxesC
                    # Use offset reference point for angle snapping
                    snapStart = tm @ offsetRefPt
                    actualLoc = loc.copy()

                    # First decide the main movement axis
                    diff = [abs(v) for v in (actualLoc - snapStart)]
                    maxDiff = max(diff)
                    axis = diff.index(maxDiff)

                    loc = snapStart.copy()
                    loc[axis] = actualLoc[axis]

                    snapIncr = 45 / self.angleSnapSteps
                    snapAngles = [
                        radians(snapIncr * a) for a in range(0, self.angleSnapSteps + 1)
                    ]

                    l1 = actualLoc[axis] - snapStart[axis]  # Main axis diff value

                    for i in range(0, 3):
                        if i != axis and (i in freeAxesN):
                            l2 = actualLoc[i] - snapStart[i]  # Minor axis value
                            angle = abs(atan(l2 / l1)) if l1 != 0 else 0
                            dirn = (l1 * l2) / abs(l1 * l2) if (l1 * l2) != 0 else 1
                            prevDiff = LARGE_NO
                            for j in range(0, len(snapAngles) + 1):
                                if j == len(snapAngles):
                                    loc[i] = snapStart[i] + dirn * l1 * tan(
                                        snapAngles[-1]
                                    )
                                    break
                                cmpAngle = snapAngles[j]
                                if abs(angle - cmpAngle) > prevDiff:
                                    loc[i] = snapStart[i] + dirn * l1 * tan(
                                        snapAngles[j - 1]
                                    )
                                    break
                                prevDiff = abs(angle - cmpAngle)

                    self.lastSnapTypes.add("angle")

        if refreshStatus and inEdit:
            refPt = tm @ orig
            newPt = loc.copy()
            text = self.getStatusStr(unit, invTm, refPt, newPt)
        else:
            text = None

        self.setStatus(rmInfo.area, text)
        loc = invTm @ loc
        self.lastSelCo = loc
        self.tm = tm
        self.orig = orig

        return loc

    def getEditCoPair(self):
        """Get the reference and current point pair for numeric input.

        Returns pair of (offset reference point, current point) in transformed coordinates.
        The offset reference point is where numeric deltas are calculated from.
        """
        # Support numeric input during custom axis drawing
        if self.customAxis.inDrawAxis and self.customAxis.length() > 0:
            # Use custom axis points: start point as origin, end point as current
            axisPts = self.customAxis.axisPts
            # Return in global coordinates (no transformation needed for custom axis)
            return (axisPts[0], axisPts[1])

        refLineOrig = self.getRefLineOrig()
        if self.lastSelCo is None or refLineOrig is None:
            return []

        params = bpy.context.window_manager.bezierToolkitParams
        if self.snapParams is None:
            origType = params.snapOrigin
            offsetRef = params.offsetRef
            refLineOrig = self.getRefLineOrig()
            selCo = self.getSelCo()
        else:
            origType = self.snapParams.origType
            offsetRef = self.snapParams.offsetRef
            refLineOrig = self.snapParams.refLineOrig
            selCo = self.snapParams.selCo

        # Use offset reference point for numeric input, not pivot
        offsetRefPt = self.getOffsetRefPoint(
            self.rmInfo, bpy.context.object, origType, offsetRef, refLineOrig, selCo
        )
        return (self.tm @ offsetRefPt, self.tm @ self.lastSelCo)

    def setStatus(self, area, text):  # TODO Global
        area.header_text_set(text)

    # Tightly coupled with get3dLocSnap
    def updateGuideBatches(self, bglDrawMgr):
        # Default values for resetting the bgl lines and points
        drawAxes = [0, 1, 2]
        axisLineCos = [[], [], []]
        axisLineCols = [[], [], []]

        snapLineCos = []
        snapLineCols = []

        custAxisLineCos = []
        custAxisLineCols = []
        custAxisYLineCos = []
        custAxisYLineCols = []
        custAxisZLineCos = []
        custAxisZLineCols = []
        custAxisGradStart = None
        custAxisGradEnd = None

        custAxisPtCos = []
        custAxisPtCols = []

        snapIndPtCos = []
        snapIndPtCols = []

        pivotPtCos = []
        pivotPtCols = []

        constraintPlaneCos = []
        constraintPlaneCols = []

        constraintGridCos = []
        constraintGridCols = []

        rmInfo = self.rmInfo

        if rmInfo is not None:  # self.snapParams is also not None
            snapParams = self.snapParams
            refLine = snapParams.refLine
            refLineOrig = snapParams.refLineOrig
            selCo = snapParams.selCo
            freeAxesN = snapParams.freeAxesN

            transType = snapParams.transType
            origType = snapParams.origType
            offsetRef = snapParams.offsetRef
            axisScale = snapParams.axisScale
            dispAxes = snapParams.dispAxes

            tm, invTm, orig = self.getTMInfoAndOrig(
                rmInfo,
                transType,
                origType,
                self.freezeOrient,
                axisScale,
                refLineOrig,
                selCo,
            )

            if (
                dispAxes
                and FTProps.showGuides
                and (
                    (
                        refLineOrig is not None
                        or transType == "VIEW"
                        or len(freeAxesN) == 1
                    )
                    or len(freeAxesN) > 0
                )
            ):
                colors = [(0.6, 0.2, 0.2, 1), (0.2, 0.6, 0.2, 1), (0.2, 0.4, 0.6, 1)]
                l = 2 * rmInfo.rv3d.view_distance

                if self.lastSelCo is not None and len(freeAxesN) == 1:
                    orig = self.lastSelCo

                refCo = tm @ orig
                for axis in freeAxesN[:2]:
                    axisLineCols[axis] = [colors[axis]]
                    pt1 = refCo.copy()
                    pt2 = refCo.copy()
                    pt1[axis] = l + refCo[axis]
                    pt2[axis] = -l + refCo[axis]
                    axisLineCos[axis] = [invTm @ pt1, invTm @ pt2]

            # Reference line highlighting
            # Show when: angle snap, polar input, or PREVIOUS offset mode with guides on
            isPreviousMode = offsetRef == "PREVIOUS"
            showRefLine = (
                refLineOrig is not None
                and self.lastSelCo is not None
                and (
                    self.angleSnap
                    or ("keyboard" in self.lastSnapTypes and self.snapDigits.polar)
                    or (FTProps.showGuides and isPreviousMode)
                )
            )
            if showRefLine:
                # Show line from offset reference point to current point
                offsetRefPt = self.getOffsetRefPoint(
                    rmInfo, bpy.context.object, origType, offsetRef, refLineOrig, selCo
                )
                snapLineCos = [offsetRefPt, self.lastSelCo]
                # Brighter color when in PREVIOUS mode
                if isPreviousMode and FTProps.showGuides:
                    snapLineCols = [(0.8, 0.5, 0.2, 1)]  # Orange-ish for reference highlight
                else:
                    snapLineCols = [(0.4, 0.4, 0.4, 1)]  # Default gray
                ptCol = (1, 1, 1, 1)

            # Custom axis visualization with full coordinate frame (X, Y, Z)
            # Drawn as separate line infos to avoid conflict with standard axes
            # Controlled by FTProps.showGuides preference
            if self.customAxis.length() != 0 and FTProps.showGuides and (
                self.customAxis.inDrawAxis or "AXIS" in {transType, origType, axisScale}
            ):
                apts = self.customAxis.axisPts

                # Get transformation matrix for custom axis to extract Y and Z axes
                tm_custom, invTm_custom = getLineTransMatrices(apts[0], apts[1])

                # Origin point for custom axis frame
                axis_orig = apts[0]
                axis_length = 2 * rmInfo.rv3d.view_distance  # Match standard axis length

                # X axis - along the custom axis (Red)
                custAxisLineCos = [apts[0], apts[1]]
                custAxisLineCols = [(0.8, 0.2, 0.2, 0.7)]  # Subtle red for X
                custAxisGradStart = 0.9
                custAxisGradEnd = 0.3

                # Y axis - perpendicular to custom axis (Green)
                # Drawn as separate line info with same gradient as X axis
                y_axis_start = axis_orig
                y_axis_end = invTm_custom @ (tm_custom @ axis_orig + Vector((0, axis_length, 0)))
                custAxisYLineCos = [y_axis_start, y_axis_end]
                custAxisYLineCols = [(0.2, 0.6, 0.2, 0.7)]  # Subtle green for Y

                # Z axis - perpendicular to custom axis (Blue)
                # Drawn as separate line info with same gradient as X axis
                z_axis_start = axis_orig
                z_axis_end = invTm_custom @ (tm_custom @ axis_orig + Vector((0, 0, axis_length)))
                custAxisZLineCos = [z_axis_start, z_axis_end]
                custAxisZLineCols = [(0.2, 0.2, 0.6, 0.7)]  # Subtle blue for Z

                # Snap division points on X axis
                custAxisPtCos = self.customAxis.getSnapPts()
                custAxisPtCols = [(1, 0.4, 0, 1)]

            # Pivot point marker - show origin when guides enabled
            if FTProps.showGuides and orig is not None:
                pivotPtCos = [orig]
                pivotPtCols = [(1.0, 0.6, 0.0, 1.0)]  # Orange color for pivot

            # Constraint plane visualization - show when 2 axes are free (plane constraint)
            if FTProps.showGuides and len(freeAxesN) == 2 and orig is not None:
                planeSize = rmInfo.rv3d.view_distance * 0.5  # Half the view distance
                refCo = tm @ orig
                # Get the two free axes
                ax1, ax2 = freeAxesN[0], freeAxesN[1]

                # For grid alignment: use integer-aligned origin in Global mode
                if transType == "GLOBAL":
                    # Align grid to integer world coordinates
                    gridOrigin = Vector([round(orig[i]) for i in range(3)])
                    gridRefCo = tm @ gridOrigin
                    # Calculate grid range for aligned plane sizing
                    gridSpacing = 1.0
                    numLines = int(planeSize / gridSpacing) + 1
                    actualPlaneSize = numLines * gridSpacing
                else:
                    gridRefCo = refCo
                    actualPlaneSize = planeSize

                # Create 4 corners of the plane
                corners = []
                for s1 in [-1, 1]:
                    for s2 in [-1, 1]:
                        corner = refCo.copy()
                        corner[ax1] += s1 * actualPlaneSize
                        corner[ax2] += s2 * actualPlaneSize
                        corners.append(invTm @ corner)
                # Draw as outline (4 lines forming a square)
                # Order: 0(-,-), 1(-,+), 2(+,-), 3(+,+)
                constraintPlaneCos = [
                    corners[0], corners[1],  # Left edge
                    corners[1], corners[3],  # Top edge
                    corners[3], corners[2],  # Right edge
                    corners[2], corners[0],  # Bottom edge
                ]
                constraintPlaneCols = [(0.5, 0.5, 0.5, 0.4)]  # Semi-transparent gray

                # Grid lines within the constraint plane
                gridSpacing = 1.0  # Fixed 1.0 unit spacing
                gridLines = []

                # Calculate grid range based on plane size and alignment
                # For aligned grids, extend range to include integer grid lines
                if transType == "GLOBAL":
                    # Calculate how many grid lines we need from aligned origin
                    numLines = int(planeSize / gridSpacing) + 1
                    gridMin = -numLines * gridSpacing
                    gridMax = numLines * gridSpacing
                else:
                    gridMin = -planeSize
                    gridMax = planeSize

                # Lines parallel to ax1 (varying along ax2)
                pos = gridMin
                while pos <= gridMax:
                    if abs(pos) < 0.01:  # Skip center line (near zero)
                        pos += gridSpacing
                        continue
                    pt1 = gridRefCo.copy()
                    pt2 = gridRefCo.copy()
                    pt1[ax1] = gridMin + gridRefCo[ax1]
                    pt1[ax2] = pos + gridRefCo[ax2]
                    pt2[ax1] = gridMax + gridRefCo[ax1]
                    pt2[ax2] = pos + gridRefCo[ax2]
                    gridLines.extend([invTm @ pt1, invTm @ pt2])
                    pos += gridSpacing

                # Lines parallel to ax2 (varying along ax1)
                pos = gridMin
                while pos <= gridMax:
                    if abs(pos) < 0.01:  # Skip center line (near zero)
                        pos += gridSpacing
                        continue
                    pt1 = gridRefCo.copy()
                    pt2 = gridRefCo.copy()
                    pt1[ax1] = pos + gridRefCo[ax1]
                    pt1[ax2] = gridMin + gridRefCo[ax2]
                    pt2[ax1] = pos + gridRefCo[ax1]
                    pt2[ax2] = gridMax + gridRefCo[ax2]
                    gridLines.extend([invTm @ pt1, invTm @ pt2])
                    pos += gridSpacing

                constraintGridCos = gridLines
                constraintGridCols = [(0.15, 0.15, 0.15, 0.6)]  # Dark gray, high opacity for visibility

            if FTProps.dispSnapInd and self.snapCo is not None:
                snapIndPtCos = [self.snapCo]
                snapIndPtCols = [FTProps.colHltTip]  # [(1, .4, 0, 1)]

        for i, axis in enumerate(drawAxes):
            axisGradStart = 0.2 if len(axisLineCos[i]) > 0 else None
            axisGradEnd = 0.9 if len(axisLineCos[i]) > 0 else None
            bglDrawMgr.addLineInfo(
                "snapAxis" + str(axis),
                FTProps.axisLineWidth,
                axisLineCols[i],
                axisLineCos[i],
                axisGradStart,
                axisGradEnd,
            )

        bglDrawMgr.addLineInfo(
            "SnapLine", FTProps.axisLineWidth, snapLineCols, snapLineCos
        )

        bglDrawMgr.addLineInfo(
            "CustAxisLine",
            FTProps.axisLineWidth,
            custAxisLineCols,
            custAxisLineCos,
            custAxisGradStart,
            custAxisGradEnd,
            mid=False,
        )

        # Custom axis Y and Z axes - drawn separately with same gradient as X
        bglDrawMgr.addLineInfo(
            "CustAxisYLine",
            FTProps.axisLineWidth,
            custAxisYLineCols,
            custAxisYLineCos,
            custAxisGradStart,
            custAxisGradEnd,
            mid=False,
        )

        bglDrawMgr.addLineInfo(
            "CustAxisZLine",
            FTProps.axisLineWidth,
            custAxisZLineCols,
            custAxisZLineCos,
            custAxisGradStart,
            custAxisGradEnd,
            mid=False,
        )

        bglDrawMgr.addPtInfo(
            "CustAxisPt", FTProps.snapPtSize, custAxisPtCols, custAxisPtCos
        )

        bglDrawMgr.addPtInfo("SnapPt", FTProps.snapPtSize, snapIndPtCols, snapIndPtCos)

        # Pivot point marker - larger size for visibility
        bglDrawMgr.addPtInfo("PivotPt", FTProps.snapPtSize * 1.5, pivotPtCols, pivotPtCos)

        # Constraint grid - drawn first so plane outline is on top
        bglDrawMgr.addLineInfo(
            "ConstraintGrid",
            FTProps.axisLineWidth * 0.4,  # Thin but visible lines for grid
            constraintGridCols,
            constraintGridCos,
        )

        # Constraint plane outline
        bglDrawMgr.addLineInfo(
            "ConstraintPlane",
            FTProps.axisLineWidth * 0.5,  # Thinner line for plane
            constraintPlaneCols,
            constraintPlaneCos,
        )


################################## Flexi Tool Classes ##################################
#
#                                         Operator
#                                            |
#                                     ModalBaseFlexiOp
#                                            |
#                          -------------------------------------
#                          |                                   |
#                 ModalDrawBezierOp                   ModalFlexiEditBezierOp
#                          |
#          ---------------------------------
#          |                               |
#    ModalFlexiDrawBezierOp       ModalFlexiDrawGreaseOp
