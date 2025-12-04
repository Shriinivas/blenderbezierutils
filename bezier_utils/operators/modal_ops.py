# bezier_utils/operators/modal_ops.py

import bpy
import blf
import gpu
from bpy.types import Operator
from bpy.app.handlers import persistent
from mathutils import Vector
import mathutils
import time
from ..core.props import FTProps
from ..core.hotkeys import FTHotKeys
from ..core.menus import FTMenu
from ..core.snap import Snapper
from ..drawing.primitives import *
from ..drawing.math_fn import MathFnDraw
from ..utils.view_utils import *
from ..utils.curve_utils import *
from ..utils.bezier_math import *
from ..utils.object_utils import *
from ..tools.workspace_tools import (
    FlexiDrawBezierTool,
    FlexiEditBezierTool,
    FlexiGreaseBezierTool,
)
from ..constants import *


class BptDisplayInfo:
    # handleNos: 0: seg1-left, 1: seg1-right
    # tipColors: leftHdl, pt, rightHdl
    # Caller to make sure there are no tips without handle
    def __init__(self, pt, tipColors, handleNos=None):
        self.pt = pt
        if len(tipColors) == 1:
            self.tipColors = [None, tipColors[0], None]
        elif len(tipColors) == 3:
            self.tipColors = tipColors
        if handleNos is None:
            self.handleNos = []
        else:
            self.handleNos = handleNos


class ModalBaseFlexiOp(Operator):
    running = False
    drawHdlrRef = None
    drawTxtHdlrRef = None
    drawFunc = None
    shader = None
    bglDrawMgr = None
    opObj = None

    pointSize = 4  # For Draw (Marker is of diff size)

    def drawKeyMap():
        regions = [
            r
            for area in bpy.context.screen.areas
            if area.type == "VIEW_3D"
            for r in area.regions
            if r.type == "WINDOW"
        ]
        maxArea = max(r.width * r.height for r in regions)
        currRegion = [r for r in bpy.context.area.regions if r.type == "WINDOW"][0]

        # Only display in window with max area
        if currRegion.width * currRegion.height < maxArea:
            return
        toolRegion = [r for r in bpy.context.area.regions if r.type == "TOOLS"][0]

        xOff1 = (
            (10 + toolRegion.width)
            if (FTProps.keyMapNextToTool)
            else FTProps.keyMapLocX
        )
        maxWidth = 0
        yStart = 10 if (FTProps.keyMapNextToTool) else FTProps.keyMapLocY
        yOff = yStart

        font_id = 0

        if ModalBaseFlexiOp.opObj is not None and FTProps.showKeyMap:
            descrCol = [0] + list(FTProps.colKeymapText)
            keyCol = [0] + list(FTProps.colKeymapKey)

            toolType = ModalBaseFlexiOp.opObj.getToolType()
            config, labels, keys = FTHotKeys.getHKDispLines(toolType)
            blf.size(font_id, FTProps.keyMapFontSize)

            maxLabelWidth = max(blf.dimensions(font_id, l + "XXXXX")[0] for l in labels)
            xOff2 = xOff1 + maxLabelWidth

            lineHeight = 1.2 * max(blf.dimensions(font_id, l)[1] for l in labels)

            blf.position(font_id, xOff1, yOff, 0)
            blf.color(*descrCol)
            blf.draw(font_id, "*")
            dim = blf.dimensions(font_id, "*")
            blf.position(font_id, xOff1 + dim[0], yOff, 0)
            blf.draw(font_id, " Indicates Configurable Hot Keys")
            yOff += 1.5 * lineHeight

            for i, label in enumerate(labels):
                blf.color(*descrCol)
                blf.position(font_id, xOff1, yOff, 0)
                blf.draw(font_id, label)
                if config[i]:
                    dim = blf.dimensions(font_id, label)
                    blf.position(font_id, xOff1 + dim[0], yOff, 0)
                    blf.draw(font_id, "*")
                blf.color(*keyCol)
                blf.position(font_id, xOff2, yOff, 0)
                blf.draw(font_id, keys[i])
                yOff += lineHeight

            maxWidth = maxLabelWidth + max(blf.dimensions(font_id, l)[0] for l in keys)
            header = toolType.title() + " Keymap"
            headerX = xOff1 + (maxWidth - blf.dimensions(font_id, header)[0]) / 2
            blf.position(font_id, headerX, yOff + lineHeight * 0.5, 0)
            blf.color(*descrCol)
            blf.draw(font_id, header)

        mathFnTxts = MathFnDraw.getMathFnTxts()

        mathFnCol = [0] + list(FTProps.colMathFnTxt)
        blf.size(font_id, FTProps.mathFnTxtFontSize)
        blf.color(*mathFnCol)
        lineHeight = blf.dimensions(font_id, "yX")[1]

        yOff = lineHeight / 2
        if mathFnTxts is not None:
            for t in mathFnTxts:
                mathFnPos = (currRegion.width - blf.dimensions(font_id, t)[0]) / 2
                # ~ mathFnPos = xOff1 + maxWidth + 10
                blf.position(font_id, mathFnPos, yOff, 0)
                blf.draw(font_id, t)
                yOff += 1.5 * lineHeight

        # Custom axis info overlay
        if ModalBaseFlexiOp.opObj is not None:
            from ..core.snap import CustomAxis
            try:
                customAxis = ModalBaseFlexiOp.opObj.snapper.customAxis
                params = bpy.context.window_manager.bezierToolkitParams

                # Show custom axis info only when AXIS mode is active (orient/origin/scale)
                if (params.snapOrient == 'AXIS' or params.snapOrigin == 'AXIS' or params.axisScale == 'AXIS'):
                    custAxisCol = [0, 1.0, 0.8, 0.2, 1.0]  # Yellow-orange
                    blf.size(font_id, FTProps.mathFnTxtFontSize)
                    blf.color(*custAxisCol)

                    # Build status text
                    snapper = ModalBaseFlexiOp.opObj.snapper
                    if customAxis.inDrawAxis:
                        # Currently drawing - show angle and length in real-time
                        angle_str = customAxis.getAngleString()
                        length = (customAxis.axisPts[1] - customAxis.axisPts[0]).length
                        statusText = f"Custom Axis: {angle_str} | Length: {length:.3f} | Snaps: {customAxis.snapCnt} (scroll)"

                        # Show numeric input if active
                        if snapper.snapDigits.hasVal():
                            if snapper.snapDigits.polar:
                                delta_str = snapper.snapDigits.getDeltaStrPolar()
                            else:
                                delta_str = snapper.snapDigits.getCurrDeltaStr()
                            statusText += f" | Input: {delta_str} (Enter to confirm)"
                    elif customAxis.length() > 0:
                        # Defined
                        angle_str = customAxis.getAngleString()
                        length = customAxis.length()
                        statusText = f"Custom Axis: {angle_str} | Length: {length:.3f} | Snaps: {customAxis.snapCnt}"
                    else:
                        statusText = "Custom Axis: Right-click to define"

                    # Draw centered at bottom of viewport
                    textWidth = blf.dimensions(font_id, statusText)[0]
                    textX = (currRegion.width - textWidth) / 2
                    textY = lineHeight  # Bottom of viewport
                    blf.position(font_id, textX, textY, 0)
                    blf.draw(font_id, statusText)
            except Exception:
                pass  # Silently ignore if custom axis info not available

    def addDrawHandlers(context):
        hdlr = ModalBaseFlexiOp.opObj.__class__.drawHandler
        ModalBaseFlexiOp.drawHdlrRef = bpy.types.SpaceView3D.draw_handler_add(
            hdlr, (), "WINDOW", "POST_VIEW"
        )
        ModalBaseFlexiOp.drawTxtHdlrRef = bpy.types.SpaceView3D.draw_handler_add(
            ModalBaseFlexiOp.drawKeyMap, (), "WINDOW", "POST_PIXEL"
        )

    def removeDrawHandlers():
        if ModalBaseFlexiOp.drawHdlrRef is not None:
            bpy.types.SpaceView3D.draw_handler_remove(
                ModalBaseFlexiOp.drawHdlrRef, "WINDOW"
            )
            ModalBaseFlexiOp.drawHdlrRef = None
        if ModalBaseFlexiOp.drawTxtHdlrRef is not None:
            bpy.types.SpaceView3D.draw_handler_remove(
                ModalBaseFlexiOp.drawTxtHdlrRef, "WINDOW"
            )
            ModalBaseFlexiOp.drawTxtHdlrRef = None

    def drawHandlerBase():
        if ModalBaseFlexiOp.shader is not None:
            ModalBaseFlexiOp.bglDrawMgr.redraw()

    def tagRedraw():
        areas = [a for a in bpy.context.screen.areas if a.type == "VIEW_3D"]
        for a in areas:
            a.tag_redraw()

    # Called back after add-on preferences are changed
    def propsChanged():
        if (
            ModalBaseFlexiOp.opObj is not None
            and ModalBaseFlexiOp.opObj.snapper is not None
        ):
            ModalBaseFlexiOp.opObj.snapper.updateSnapKeyMap()

    def resetDisplayBase():
        ModalBaseFlexiOp.bglDrawMgr.reset()
        ModalBaseFlexiOp.tagRedraw()

    def refreshDisplayBase(segDispInfos, bptDispInfos, snapper):
        areaRegionInfo = getAllAreaRegions()

        updateBezierBatches(
            ModalDrawBezierOp.bglDrawMgr, segDispInfos, bptDispInfos, areaRegionInfo
        )

        if snapper is not None:
            snapper.updateGuideBatches(ModalBaseFlexiOp.bglDrawMgr)

        ModalBaseFlexiOp.tagRedraw()

    @persistent
    def loadPostHandler(dummy):
        if ModalBaseFlexiOp.shader is not None:
            ModalBaseFlexiOp.resetDisplayBase()
        ModalBaseFlexiOp.running = False

    @persistent
    def loadPreHandler(dummy):
        ModalBaseFlexiOp.removeDrawHandlers()
        if ModalBaseFlexiOp.drawFunc is not None:
            bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc

    @classmethod
    def poll(cls, context):
        return not ModalBaseFlexiOp.running

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getToolType(self):
        raise NotImplementedError("Call to abstract method.")

    def preInvoke(self, context, event):
        pass  # placeholder

    def subInvoke(self, context, event):
        return {"RUNNING_MODAL"}  # placeholder

    def invoke(self, context, event):
        ModalBaseFlexiOp.opObj = self
        ModalBaseFlexiOp.running = True
        self.preInvoke(context, event)
        ModalBaseFlexiOp.addDrawHandlers(context)

        ModalBaseFlexiOp.drawFunc = bpy.types.VIEW3D_HT_tool_header.draw
        bpy.types.VIEW3D_HT_tool_header.draw = drawSettingsFT
        context.space_data.show_region_tool_header = True

        self.snapper = Snapper(
            context,
            self.getSnapLocs,
            self.getRefLine,
            self.getRefLineOrig,
            self.getSelCo,
            self.getCurrLine,
            self.hasSelection,
            self.isEditing,
        )

        self.rmInfo = None

        ModalBaseFlexiOp.shader = gpu.shader.from_builtin("FLAT_COLOR")
        ModalBaseFlexiOp.bglDrawMgr = BGLDrawMgr(ModalBaseFlexiOp.shader)

        # ~ ModalBaseFlexiOp.shader.bind()
        context.window_manager.modal_handler_add(self)

        ModalBaseFlexiOp.ColGreaseHltSeg = (0.3, 0.3, 0.3, 1)  # Not used

        FTProps.updateProps(None, context)
        FTHotKeys.updateHotkeys(None, context)
        FTHotKeys.initSnapMetaFromPref(context)

        self.clickT, self.pressT = None, None
        self.click, self.doubleClick = False, False

        return self.subInvoke(context, event)

    def modal(self, context, event):
        snapper = self.snapper
        if event.type == "WINDOW_DEACTIVATE" and event.value == "PRESS":
            snapper.initialize()
            return {"PASS_THROUGH"}

        if not self.isToolSelected(context):  # Subclass
            self.cancelOp(context)
            return {"CANCELLED"}

        self.click, self.doubleClick = False, False
        if event.type == "LEFTMOUSE":
            if event.value == "PRESS":
                self.pressT = time.time()
            elif event.value == "RELEASE":
                t = time.time()
                if self.clickT is not None and (t - self.clickT) < DBL_CLK_DURN:
                    self.clickT = None
                    self.doubleClick = True
                elif self.pressT is not None and (t - self.pressT) < SNGL_CLK_DURN:
                    self.clickT = t
                    self.click = True
                self.pressT = None

        snapProc = snapper.procEvent(context, event)
        metakeys = snapper.getMetakeys()

        if FTHotKeys.isHotKey(FTHotKeys.hkToggleKeyMap, event.type, metakeys):
            if event.value == "RELEASE":
                try:
                    addon_prefs = context.preferences.addons.get("bezier_utils")
                    prefs = addon_prefs.preferences if addon_prefs else None
                    prefs.showKeyMap = not prefs.showKeyMap
                except Exception as e:
                    print(e)
                    FTProps.showKeyMap = not FTProps.showKeyMap
            return {"RUNNING_MODAL"}

        # Special condition for case where user has configured a different snap key
        # In such case, pass through mouse clicks, if there's a meta key held down...
        # to allow e.g. alt + LMB to rotate view
        if (
            snapProc != EVT_CONS
            and any(metakeys)
            and not any([snapper.angleSnap, snapper.gridSnap, snapper.vertSnap])
            and (
                (event.type == "LEFTMOUSE" and event.value in {"PRESS", "RELEASE"})
                or event.type == "MOUSEMOVE"
            )
        ):
            return {"PASS_THROUGH"}

        rmInfo = RegionMouseXYInfo.getRegionMouseXYInfo(event, self.exclToolRegion())

        ret = FTMenu.procMenu(self, context, event, rmInfo is None)
        if ret:
            # Menu displayed on release, so retain metakeys till release
            if event.value == "RELEASE":
                snapper.resetMetakeys()
                snapper.resetSnapKeys()
            return {"RUNNING_MODAL"}

        if self.isEditing() and self.rmInfo != rmInfo:
            return {"RUNNING_MODAL"}
        if rmInfo is None:
            return {"PASS_THROUGH"}

        self.rmInfo = rmInfo

        ret = snapper.customAxis.procDrawEvent(context, event, snapper, rmInfo)
        evtCons = ret or snapProc == EVT_CONS

        # Ignore all PRESS events if consumed, since action is taken only on RELEASE...
        # ...except 1) wheelup / down where there is no release & 2) snap / meta where...
        # ...refresh is needed even on press
        # TODO: Simplify the condition (Maybe return EVT values from all proc methods)
        if (
            evtCons
            and event.value == "PRESS"
            and event.type != "MOUSEMOVE"
            and not event.type.startswith("WHEEL")
            and (snapProc != EVT_META_OR_SNAP)
        ):
            return {"RUNNING_MODAL"}

        if FTHotKeys.isHotKey(FTHotKeys.hkSwitchOut, event.type, metakeys):
            if event.value == "RELEASE":
                self.cancelOp(context)
                resetToolbarTool()
                return {"CANCELLED"}
            return {"RUNNING_MODAL"}

        evtCons = evtCons or (snapProc == EVT_META_OR_SNAP)
        return self.subModal(context, event, evtCons)

    def cancelOpBase(self):
        for a in bpy.context.screen.areas:
            if a.type == "VIEW_3D":
                a.header_text_set(None)

        ModalBaseFlexiOp.removeDrawHandlers()
        ModalBaseFlexiOp.running = False
        bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc
        self.snapper = None
        ModalBaseFlexiOp.opObj = None

    def getSnapLocs(self):
        return self.getSnapLocsImpl()


################################## Flexi Draw Bezier Curve ###############################
#
#                                         BaseDraw
#                                            |
#                                -------------------------------
#                                |                             |
#                           Primitive2DDraw                BezierDraw
#                                |
#                        -----------------------------
#                        |                           |
#               ClosedShapeDraw               MathFnDraw
#                        |
#        ---------------------------------
#        |               |               |
#  RectangleDraw    PolygonDraw     EllipseDraw


class ModalDrawBezierOp(ModalBaseFlexiOp):
    # Static members shared by flexi draw and flexi grease
    markerSize = 8
    h = False

    drawObjMap = {}

    # static method
    def drawHandler():
        ModalBaseFlexiOp.drawHandlerBase()

    def updateDrawType(dummy, context):
        opObj = ModalDrawBezierOp.opObj
        if opObj is not None:
            opObj.setDrawObj()
            opObj.initialize()
            opObj.resetDisplay()
            ModalDrawBezierOp.updateDrawSides(dummy, context)

    def updateDrawSides(dummy, context):
        params = bpy.context.window_manager.bezierToolkitParams
        opObj = ModalDrawBezierOp.opObj
        if opObj is not None and opObj.drawType != "BEZIER":
            opObj.drawObj.shapeSegCnt = params.drawSides

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getToolType(self):
        return TOOL_TYPE_FLEXI_DRAW

    def resetDisplay(self):
        ModalBaseFlexiOp.resetDisplayBase()

    def setDrawObj(self):
        params = bpy.context.window_manager.bezierToolkitParams
        self.drawObj = ModalDrawBezierOp.drawObjMap[params.drawObjType]
        self.drawType = params.drawObjType

    def preInvoke(self, context, event):
        self.bezierDrawObj = BezierDraw(self)
        self.rectangleDrawObj = RectangleDraw(self)
        self.ellipseDrawObj = EllipseDraw(self)
        self.polygonDrawObj = PolygonDraw(self)
        self.starDrawObj = PolygonDraw(self, star=True)
        self.mathDrawObj = MathFnDraw(self)

        ModalDrawBezierOp.drawObjMap = {
            "BEZIER": self.bezierDrawObj,
            "RECTANGLE": self.rectangleDrawObj,
            "ELLIPSE": self.ellipseDrawObj,
            "POLYGON": self.polygonDrawObj,
            "STAR": self.starDrawObj,
            "MATH": self.mathDrawObj,
        }
        self.setDrawObj()

    # This will be called multiple times not just at the beginning
    def initialize(self):
        self.drawObj.initialize()
        self.snapper.initialize()

    def subInvoke(self, context, event):
        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)

        try:
            ModalDrawBezierOp.markerSize = (
                context.preferences.addons.get("bezier_utils").preferences.markerSize
                if context.preferences.addons.get("bezier_utils")
                else 7
            )
        except Exception as e:
            # ~ print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
            ModalDrawBezierOp.markerSize = 8

        self.updateDrawType(context)
        self.updateDrawSides(context)
        return {"RUNNING_MODAL"}

    def cancelOp(self, context):
        self.resetDisplay()
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)
        return self.cancelOpBase()

    def postUndoRedo(self, scene, dummy=None):  # signature different in 2.8 and 2.81?
        self.updateSnapLocs()  # subclass method

    def confirm(self, context, event, location=None):
        metakeys = self.snapper.getMetakeys()
        shift = self.snapper.angleSnap  # Overloaded key op
        autoclose = (
            self.drawType == "BEZIER"
            and shift
            and (event.type == "SPACE" or event.type == "RET")
        )
        self.save(context, event, autoclose, location)
        self.resetDisplay()
        self.initialize()

    def exclToolRegion(self):
        return True

    def isEditing(self):
        return len(self.drawObj.curvePts) > 0

    def hasSelection(self):
        return self.isEditing()

    # Common subModal for Flexi Draw and Flexi Grease
    def baseSubModal(self, context, event, snapProc):
        return self.drawObj.procDrawEvent(context, event, snapProc)

    def refreshMarkerPos(self, rmInfo):
        colMap = self.getColorMap()
        colMarker = colMap["MARKER_COLOR"]
        markerLoc = self.snapper.get3dLocSnap(rmInfo)

        self.resetDisplay()
        self.bglDrawMgr.addPtInfo(
            "drawMarker", ModalDrawBezierOp.markerSize, [colMarker], [markerLoc]
        )

        ModalBaseFlexiOp.refreshDisplayBase(
            segDispInfos=[], bptDispInfos=[], snapper=self.snapper
        )

    def redrawBezier(self, rmInfo, lastSegOnly=False, hdlPtIdxs=None, hltEndSeg=True):
        curvePts = self.drawObj.curvePts

        ptCnt = len(curvePts)

        if ptCnt == 0:
            self.refreshMarkerPos(rmInfo)
            return

        self.bglDrawMgr.resetPtInfo("drawMarker")

        colMap = self.getColorMap()
        colSelSeg = colMap["SEL_SEG_COLOR"]
        colNonAdjSeg = colMap["NONADJ_SEG_COLOR"]
        colTip = colMap["TIP_COLOR"]
        colEndTip = colMap["ENDPT_TIP_COLOR"]

        segColor = colSelSeg
        tipColors = (
            [colTip, colEndTip, colTip]
            if (not ModalDrawBezierOp.h)
            else [None, colEndTip, None]
        )

        segDispInfos = []
        bptDispInfos = []

        if hdlPtIdxs is None:
            hdlPtIdxs = {ptCnt - 2}  # Default last but one
        elif len(hdlPtIdxs) == 0:
            hdlPtIdxs = range(ptCnt)

        for hdlPtIdx in hdlPtIdxs:
            bptDispInfos.append(
                BptDisplayInfo(
                    curvePts[hdlPtIdx],
                    tipColors,
                    handleNos=[0, 1] if (not ModalDrawBezierOp.h) else [],
                )
            )

        startIdx = 0
        if lastSegOnly and ptCnt > 1:
            startIdx = ptCnt - 2
        for i in range(startIdx, ptCnt - 1):
            if not hltEndSeg or i == ptCnt - 2:
                segColor = colSelSeg
            else:
                segColor = colNonAdjSeg
            segDispInfos.append(
                SegDisplayInfo([curvePts[i], curvePts[i + 1]], segColor)
            )

        ModalBaseFlexiOp.refreshDisplayBase(segDispInfos, bptDispInfos, self.snapper)

    def getRefLine(self):
        return self.drawObj.getRefLine()

    def getRefLineOrig(self):
        return self.drawObj.getRefLineOrig()

    def getSelCo(self):
        return self.getRefLineOrig()

    def getCurrLine(self):
        return self.getRefLine()


class ModalFlexiDrawBezierOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves in object mode"
    bl_idname = "wm.flexi_draw_bezier_curves"
    bl_label = "Flexi Draw Bezier Curves"
    bl_options = {"REGISTER", "UNDO"}

    # For some curve-changing ops (like reset rotation); possible in draw
    def updateAfterGeomChange(self, scene=None, dummy=None):  # 3 params in 2.81
        self.updateSnapLocs()

    def isToolSelected(self, context):
        if context.mode != "OBJECT":
            return False

        tool = context.workspace.tools.from_space_view3d_mode("OBJECT", create=False)

        if tool is None or tool.idname != FlexiDrawBezierTool.bl_idname:
            # if(tool == None or tool.idname != 'flexi_bezier.draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {
            "SEL_SEG_COLOR": FTProps.colDrawSelSeg,
            "NONADJ_SEG_COLOR": FTProps.colDrawNonHltSeg,
            "TIP_COLOR": FTProps.colHdlPtTip,
            "ENDPT_TIP_COLOR": FTProps.colBezPt,
            "MARKER_COLOR": FTProps.colDrawMarker,
        }

    def cancelOp(self, context):
        bpy.app.handlers.depsgraph_update_post.remove(self.updateAfterGeomChange)
        super(ModalFlexiDrawBezierOp, self).cancelOp(context)

    def preInvoke(self, context, event):
        super(ModalFlexiDrawBezierOp, self).preInvoke(context, event)
        # If the operator is invoked from context menu, enable the tool on toolbar
        if not self.isToolSelected(context) and context.mode == "OBJECT":
            bpy.ops.wm.tool_set_by_id(name=FlexiDrawBezierTool.bl_idname)
            # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')

        # Object name -> [spline index, [pts]]
        # Not used right now (maybe in case of large no of curves)
        self.snapInfos = {}
        self.updateSnapLocs()
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()

        if FTHotKeys.isHotKey(FTHotKeys.hkToggleDrwEd, event.type, metakeys):
            if event.value == "RELEASE":
                bpy.ops.wm.tool_set_by_id(name=FlexiEditBezierTool.bl_idname)
                # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.edit_tool')
            return {"RUNNING_MODAL"}

        return self.baseSubModal(context, event, snapProc)

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if len(self.drawObj.curvePts) > 0:
            locs += [pt[1] for pt in self.drawObj.curvePts[:-1]]
            # ~ locs += [self.curvePts[-1][0], self.curvePts[-1][1], self.curvePts[-1][2]]
        return locs

    def updateSnapLocs(self, objNames=None):
        updateCurveEndPtMap(self.snapInfos, addObjNames=objNames)

    def createCurveObj(
        self,
        context,
        startObj=None,
        startSplineIdx=None,
        endObj=None,
        endSplineIdx=None,
        autoclose=False,
    ):
        # First create the new curve
        collection = context.collection
        if collection is None:
            collection = context.scene.collection
        obj = createObjFromPts(self.drawObj.curvePts, "3D", collection, autoclose)

        # Undo stack in case the user does not want to join
        if endObj is not None or startObj is not None:
            obj.select_set(True)
            # ~ bpy.context.view_layer.objects.active = obj
            bpy.ops.ed.undo_push()
        else:
            return obj

        endObjs = []

        # Connect the end curve (if exists) first
        splineIdx = endSplineIdx

        if endObj is not None and startObj != endObj:
            # first separate splines
            endObjs, changeCnt = splitCurve([endObj], split="spline", newColl=False)

            # then join the selected spline from end curve with new obj
            obj = joinSegs(
                [endObjs[splineIdx], obj],
                optimized=True,
                straight=False,
                srcCurve=endObjs[splineIdx],
            )

            endObjs[splineIdx] = obj

            # Use this if there is no start curve
            objComps = endObjs

        if startObj is not None:
            # Repeat the above process for start curve
            startObjs, changeCnt = splitCurve([startObj], split="spline", newColl=False)

            obj = joinSegs(
                [startObjs[startSplineIdx], obj],
                optimized=True,
                straight=False,
                srcCurve=startObjs[startSplineIdx],
            )

            if startObj == endObj and startSplineIdx != endSplineIdx:
                # If startSplineIdx == endSplineIdx the join call above would take care
                # but if they are different they need to be joined with a separate call
                obj = joinSegs(
                    [startObjs[endSplineIdx], obj],
                    optimized=True,
                    straight=False,
                    srcCurve=startObjs[endSplineIdx],
                )

                # can't replace the elem with new obj as in case of end curve
                # (see the seq below)
                startObjs.pop(endSplineIdx)
                if endSplineIdx < startSplineIdx:
                    startSplineIdx -= 1

            # Won't break even if there were no endObjs
            objComps = (
                startObjs[:startSplineIdx]
                + endObjs[:splineIdx]
                + [obj]
                + endObjs[(splineIdx + 1) :]
                + startObjs[(startSplineIdx + 1) :]
            )

        obj = joinCurves(objComps)

        if any(p.co.z != 0 for s in obj.data.splines for p in s.bezier_points):
            obj.data.dimensions = "3D"

        return obj

    # TODO: At least store in map instead of linear search
    def getSnapObjs(self, context, locs):
        retVals = [[None, 0, 0]] * len(locs)
        foundVals = 0
        for obj in bpy.data.objects:
            if isBezier(obj):
                mw = obj.matrix_world
                for i, s in enumerate(obj.data.splines):
                    if s.use_cyclic_u or len(s.bezier_points) == 0:
                        continue
                    for j, loc in enumerate(locs):
                        p = s.bezier_points[0]
                        if vectCmpWithMargin(loc, mw @ p.co):
                            retVals[j] = [obj, i, 0]
                            foundVals += 1
                        p = s.bezier_points[-1]
                        if vectCmpWithMargin(loc, mw @ p.co):
                            retVals[j] = [obj, i, -1]
                            foundVals += 1
                        if foundVals == len(locs):
                            return retVals
        return retVals

    def save(self, context, event, autoclose, location, align=True):
        curvePts = self.drawObj.curvePts
        if len(curvePts) > 1:
            startObj, startSplineIdx, ptIdx2, endObj, endSplineIdx, ptIdx1 = [
                x
                for y in self.getSnapObjs(context, [curvePts[0][1], curvePts[-1][1]])
                for x in y
            ]

            ctrl = self.snapper.gridSnap  # Overloaded key op

            # ctrl pressed and there IS a snapped end obj,
            # so user does not want connection

            # (no option to only connect to starting curve when end object exists)
            if ctrl and endObj is not None:
                obj = self.createCurveObj(context, autoclose=False)
            else:
                startObjName = startObj.name if (startObj is not None) else ""
                endObjName = endObj.name if (endObj is not None) else ""

                obj = self.createCurveObj(
                    context, startObj, startSplineIdx, endObj, endSplineIdx, autoclose
                )

            if align and startObj is None and endObj is None:
                alignToNormal(obj)
                bpy.context.evaluated_depsgraph_get().update()
                if location is None:
                    location = getObjBBoxCenter(obj)

            if location is not None:
                shiftOrigin(obj, location)
                obj.location = location
                bpy.context.evaluated_depsgraph_get().update()

            params = bpy.context.window_manager.bezierToolkitParams
            copyProperties(params.copyPropsObj, obj)

            # TODO: Why try?
            try:
                obj.select_set(True)
                # ~ bpy.context.view_layer.objects.active = obj
                self.updateSnapLocs([obj.name, startObjName, endObjName])
            except Exception as e:
                pass
        bpy.ops.ed.undo_push()


################### Flexi Draw Grease Bezier ###################


class ModalFlexiDrawGreaseOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves as grease pencil strokes"
    bl_idname = "wm.flexi_draw_grease_bezier_curves"
    bl_label = "Flexi Draw Grease Bezier Curves"
    bl_options = {"REGISTER", "UNDO"}

    h = False

    def getToolType(self):
        return TOOL_TYPE_FLEXI_GREASE

    def isToolSelected(self, context):
        if context.mode != GP_CONTEXT_MODE:
            return False

        tool = context.workspace.tools.from_space_view3d_mode(
            GP_CONTEXT_MODE, create=False
        )

        if tool is None or tool.idname != FlexiGreaseBezierTool.bl_idname:
            # if(tool == None or tool.idname != 'flexi_bezier.grease_draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {
            "SEL_SEG_COLOR": FTProps.colGreaseSelSeg,
            "NONADJ_SEG_COLOR": ModalBaseFlexiOp.ColGreaseHltSeg,  # Not used
            "TIP_COLOR": FTProps.colHdlPtTip,
            "ENDPT_TIP_COLOR": FTProps.colGreaseBezPt,
            "MARKER_COLOR": FTProps.colGreaseMarker,
        }

    def preInvoke(self, context, event):
        super(ModalFlexiDrawGreaseOp, self).preInvoke(context, event)
        # If the operator is invoked from context menu, enable the tool on toolbar
        if not self.isToolSelected(context) and context.mode == GP_CONTEXT_MODE:
            bpy.ops.wm.tool_set_by_id(name=FlexiGreaseBezierTool.bl_idname)
            # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.grease_draw_tool')

        o = context.object
        if o is None or o.type != "GREASEPENCIL":
            d = bpy.data.grease_pencils.new("Grease Pencil Data")
            o = bpy.data.objects.new("Grease Pencil", d)
            context.scene.collection.objects.link(o)
        self.gpencil = o

        self.subdivCos = []
        self.interpPts = []
        self.subdivPerUnit = None

    # overridden
    def redrawBezier(self, rmInfo, hdlPtIdxs=None, hltEndSeg=True):
        curvePts = self.drawObj.curvePts
        ptCnt = len(curvePts)
        subdivCos = self.subdivCos if ptCnt > 1 else []
        self.bglDrawMgr.addLineInfo(
            "gpSubdivLines",
            FTProps.lineWidth,
            [FTProps.colGreaseNonHltSeg],
            getLinesFromPts(subdivCos),
        )
        if ModalFlexiDrawGreaseOp.h:
            self.bglDrawMgr.resetPtInfo("gpSubdivPts")
        else:
            self.bglDrawMgr.addPtInfo(
                "gpSubdivPts",
                FTProps.greaseSubdivPtSize,
                [FTProps.colGreaseSubdiv],
                subdivCos,
            )

        if self.drawType != "BEZIER" and len(curvePts) > 0:
            ModalBaseFlexiOp.refreshDisplayBase(
                segDispInfos=[], bptDispInfos=[], snapper=self.snapper
            )
        else:
            super(ModalFlexiDrawGreaseOp, self).redrawBezier(
                rmInfo, lastSegOnly=True, hdlPtIdxs={ptCnt - 1}, hltEndSeg=hltEndSeg
            )

    def initialize(self):
        super(ModalFlexiDrawGreaseOp, self).initialize()
        self.subdivCos = []
        self.interpPts = []
        self.updateSnapLocs()

    def subModal(self, context, event, snapProc):
        curvePts = self.drawObj.curvePts
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()

        if not metakeys[2]:
            cntIncr = 5  # if(self.isDrawShape) else 5

            if (
                self.drawType in {"BEZIER", "ELLIPSE"}
                and event.type
                in {
                    "WHEELDOWNMOUSE",
                    "WHEELUPMOUSE",
                    "NUMPAD_PLUS",
                    "NUMPAD_MINUS",
                    "PLUS",
                    "MINUS",
                }
                and len(curvePts) > 1
            ):
                if (
                    event.type in {"NUMPAD_PLUS", "NUMPAD_MINUS", "PLUS", "MINUS"}
                    and event.value == "PRESS"
                ):
                    return {"RUNNING_MODAL"}
                elif event.type == "WHEELUPMOUSE" or event.type.endswith("PLUS"):
                    self.subdivAdd(cntIncr)
                elif event.type == "WHEELDOWNMOUSE" or event.type.endswith("MINUS"):
                    self.subdivAdd(-cntIncr)

                self.redrawBezier(rmInfo)
                return {"RUNNING_MODAL"}

        if event.type == "H" or event.type == "h":
            if event.value == "RELEASE":
                ModalFlexiDrawGreaseOp.h = not ModalFlexiDrawGreaseOp.h
                self.redrawBezier(self.rmInfo)
            return {"RUNNING_MODAL"}

        ptCnt = len(curvePts)

        retVal = self.baseSubModal(context, event, snapProc)

        newPtCnt = len(curvePts)
        # ~ if(newPtCnt - ptCnt != 0):
        if len(curvePts) > 0:
            if self.subdivPerUnit is None:
                viewDist = context.space_data.region_3d.view_distance
                self.initSubdivPerUnit = (
                    5000.0 / viewDist
                )  # TODO: default configurable?
                self.subdivPerUnit = 0.02 * self.initSubdivPerUnit
                self.snapLocs.append(curvePts[0][1])
            if len(curvePts) > 1:
                slens = self.getCurveSegLens()
                self.updateInterpPts(slens)
                self.updateSubdivCos(sum(slens))
                self.redrawBezier(rmInfo)

        return retVal

    def getCurveSegLens(self):
        clen = []
        curvePts = self.drawObj.curvePts

        for i in range(1, len(curvePts) - 1):
            clen.append(
                getSegLen(
                    [
                        curvePts[i - 1][1],
                        curvePts[i - 1][2],
                        curvePts[i][0],
                        curvePts[i][1],
                    ]
                )
            )
        return clen

    def updateSubdivCos(self, clen=None):
        if self.drawType in {"POLYGON", "STAR", "RECTANGLE"}:
            self.subdivCos = [p[0] for p in self.drawObj.curvePts]
        elif self.interpPts != []:
            if clen is None:
                clen = sum(self.getCurveSegLens())
            cnt = round(self.subdivPerUnit * clen)
            if cnt > 0:
                self.subdivCos = getInterpolatedVertsCo(self.interpPts, cnt)  # [1:-1]
                return
            self.subdivCos = []

    def updateInterpPts(self, slens):
        curvePts = (
            self.drawObj.curvePts[:]
            if (self.drawType != "BEZIER")
            else self.drawObj.curvePts[:-1]
        )

        self.interpPts = getInterpBezierPts(curvePts, self.initSubdivPerUnit, slens)

        return self.interpPts

    def subdivAdd(self, addCnt):
        slens = self.getCurveSegLens()
        clen = sum(slens)
        if clen == 0:
            return
        cnt = self.subdivPerUnit * clen + addCnt
        if cnt < 1:
            cnt = 1

        self.subdivPerUnit = cnt / clen
        self.updateSubdivCos(clen)

    def getSnapLocsImpl(self):
        return self.snapLocs

    def updateSnapLocs(self):
        self.snapLocs = []
        gpencils = [o for o in bpy.data.objects if o.type == "GREASEPENCIL"]
        for gpencil in gpencils:
            mw = gpencil.matrix_world
            for layer in gpencil.data.layers:
                for f in layer.frames:
                    for s in f.drawing.strokes:
                        if len(s.points) > 0:  # Shouldn't be needed, but anyway...
                            self.snapLocs += [
                                mw @ s.points[0].position,
                                mw @ s.points[-1].position,
                            ]

    def save(self, context, event, autoclose, location):
        layer = self.gpencil.data.layers.active
        if layer is None:
            layer = self.gpencil.data.layers.new("GP_Layer", set_active=True)
        if len(layer.frames) == 0:
            layer.frames.new(0)
        frame = layer.frames[-1]

        invMw = self.gpencil.matrix_world.inverted_safe()
        point_count = len(self.subdivCos)
        if point_count > 0:
            brush = context.scene.tool_settings.gpencil_paint.brush
            lineWidth = brush.size
            strength = brush.gpencil_settings.pen_strength

            frame.drawing.add_strokes([1])
            stroke = frame.drawing.strokes[-1]
            # stroke.display_mode = '3DSPACE'
            stroke.add_points(count=point_count)
            for i in range(0, point_count):
                pt = self.subdivCos[i]
                stroke.points[i].position = (
                    self.gpencil.matrix_world.inverted_safe() @ pt
                )
                stroke.points[i].opacity = strength
                stroke.points[i].radius = lineWidth / 100
            # stroke.line_width = lineWidth
            self.snapLocs += [self.subdivCos[0][1], self.subdivCos[-1][1]]

            # For some reason an extra point is added at the end, so remove it (ver 4.3)
            stroke.remove_points(1)

            if autoclose:
                stroke.add_points(count=1)
                stroke.points[-1].position = stroke.points[0].position.copy()
                stroke.points[-1].opacity = strength

        bpy.ops.ed.undo_push()


################### Flexi Edit Bezier Curve ###################


class EditSegDisplayInfo(SegDisplayInfo):
    def __init__(self, segPts, segColor, subdivCos):
        super(EditSegDisplayInfo, self).__init__(segPts, segColor)
        self.subdivCos = subdivCos


def getSegPtsInSpline(wsData, splineIdx, ptIdx, cyclic):
    splinePts = wsData[splineIdx]
    if ptIdx < (len(splinePts) - 1):
        ptRange = [ptIdx, ptIdx + 1]
    elif ptIdx == (len(splinePts) - 1) and cyclic:
        ptRange = [-1, 0]
        if splinePts[-1][4] == "VECTOR":
            splinePts[-1][2] = splinePts[-1][1] + 1 / 3 * (
                splinePts[-1][1] - splinePts[0][1]
            )
        if splinePts[0][3] == "VECTOR":
            splinePts[0][0] = splinePts[0][1] + 1 / 3 * (
                splinePts[0][1] - splinePts[-1][1]
            )
    else:
        return []

    return [[splinePts[x][i] for i in range(5)] for x in ptRange]


def getInterpSegPts(wsData, splineIdx, ptIdx, cyclic, res, maxRes):
    segPts = getSegPtsInSpline(wsData, splineIdx, ptIdx, cyclic)
    areaRegionInfo = getAllAreaRegions()  # TODO: To be passed from caller

    return getPtsAlongBezier2D(segPts, areaRegionInfo, res, getCoordFromLoc, maxRes)


# Wrapper for spatial search within segment
def getClosestPt2dWithinSeg(
    region, rv3d, coFind, selObj, selSplineIdx, selSegIdx, withHandles, withBezPts
):
    infos = {selObj: {selSplineIdx: [[selSegIdx], []]}}

    # set selObj in objs for CurveBezPts
    return getClosestPt2d(
        region,
        rv3d,
        coFind,
        [selObj],
        infos,
        withHandles,
        withBezPts,
        withObjs=False,
        maxSelObjRes=MAX_SEL_CURVE_RES,
    )


def getClosestPt2d(
    region,
    rv3d,
    coFind,
    objs,
    selObjInfos,
    withHandles=True,
    withBezPts=True,
    withObjs=True,
    maxSelObjRes=MAX_NONSEL_CURVE_RES,
    withShapeKey=True,
):
    objLocMap = {}

    objLocList = []  # For mapping after search returns
    objInterpLocs = []
    objInterpCounts = []
    objBezPtCounts = []

    objSplineEndPts = []

    for obj in objs:
        # TODO: Check of shape key bounding box
        if obj.active_shape_key is None and not isPtIn2dBBox(
            obj, region, rv3d, coFind, FTProps.snapDist
        ):
            continue

        wsDataSK = None
        # Curve data with shape key value applied (if shape key exists)
        wsData = getBptData(obj, fromMix=True, updateDeps=True)
        if withShapeKey and obj.active_shape_key is not None:
            # active shape key data with value = 1
            wsDataSK = getBptData(obj, fromMix=False)

        for i, spline in enumerate(obj.data.splines):
            for j, pt in enumerate(spline.bezier_points):
                objLocList.append([obj, i, j])
                if withObjs:
                    interpLocs = getInterpSegPts(
                        wsData,
                        i,
                        j,
                        spline.use_cyclic_u,
                        res=SEARCH_CURVE_RES,
                        maxRes=MAX_NONSEL_CURVE_RES,
                    )[1:-1]
                    if wsDataSK is not None:
                        interpLocs += getInterpSegPts(
                            wsDataSK,
                            i,
                            j,
                            spline.use_cyclic_u,
                            res=SEARCH_CURVE_RES,
                            maxRes=MAX_NONSEL_CURVE_RES,
                        )[1:-1]

                    objInterpLocs += interpLocs
                    objInterpCounts.append(len(interpLocs))

                if withBezPts:
                    cnt = 1
                    objSplineEndPts.append(wsData[i][j][1])  # mw @ pt.co)
                    if wsDataSK is not None:
                        objSplineEndPts.append(wsDataSK[i][j][1])  # mw @ pt.co)
                        cnt += 1
                    objBezPtCounts.append(cnt)

    selObjLocList = []  # For mapping after search returns
    selObjHdlList = []  # Better to create a new one, even if some redundancy

    segInterpLocs = []
    selObjInterpCounts = []
    selObjHdlCounts = []

    hdls = []

    for selObj in selObjInfos.keys():
        wsDataSK = None
        # Curve data with shape key value applied (if shape key exists)
        wsData = getBptData(selObj, fromMix=True, updateDeps=True)
        if withShapeKey and selObj.active_shape_key is not None:
            # active shape key data with value = 1
            wsDataSK = getBptData(selObj, fromMix=False)

        info = selObjInfos[selObj]
        for splineIdx in info.keys():
            cyclic = selObj.data.splines[splineIdx].use_cyclic_u
            segIdxs = info[splineIdx][0]
            for segIdx in segIdxs:
                selObjLocList.append([selObj, splineIdx, segIdx])
                interpLocs = getInterpSegPts(
                    wsData,
                    splineIdx,
                    segIdx,
                    cyclic,
                    res=SEARCH_CURVE_RES * 5,
                    maxRes=maxSelObjRes,
                )[1:-1]
                if wsDataSK is not None:
                    interpLocs += getInterpSegPts(
                        wsDataSK,
                        splineIdx,
                        segIdx,
                        cyclic,
                        res=SEARCH_CURVE_RES * 5,
                        maxRes=maxSelObjRes,
                    )[1:-1]
                segInterpLocs += interpLocs
                selObjInterpCounts.append(len(interpLocs))

            if withHandles:
                ptIdxs = info[splineIdx][1]
                for ptIdx in ptIdxs:
                    selObjHdlList.append([selObj, splineIdx, ptIdx])
                    hdlCnt = 2
                    if wsDataSK is not None:
                        pt = wsDataSK[splineIdx][ptIdx]
                        hdls += [pt[0], pt[2]]
                    else:
                        pt = wsData[splineIdx][ptIdx]
                        hdls += [pt[0], pt[2]]

                    selObjHdlCounts.append(hdlCnt)

    searchPtsList = [[], [], [], [], [], []]
    retStr = [[], [], [], [], [], []]

    #'SelHandles', 'SegLoc', 'CurveBezPt', 'CurveLoc'

    searchPtsList[0], retStr[0] = hdls, "SelHandles"
    searchPtsList[1], retStr[1] = objSplineEndPts, "CurveBezPt"
    searchPtsList[2], retStr[2] = segInterpLocs, "SegLoc"
    searchPtsList[3], retStr[3] = objInterpLocs, "CurveLoc"

    # TODO: Remove duplicates before sending for search?
    searchPtsList = [
        [getCoordFromLoc(region, rv3d, pt).to_3d() for pt in pts]
        for pts in searchPtsList
    ]

    srs = NestedListSearch(searchPtsList).findInLists(
        coFind, searchRange=FTProps.snapDist
    )

    if len(srs) == 0:
        return None

    sr = min(srs, key=lambda x: x[3])

    if sr[0] > 1:
        # If seg loc then first priority to the nearby handle, end pt (even if farther)
        sr = min(srs, key=lambda x: (x[0], x[3]))

    idx = sr[1]
    retId = retStr[sr[0]]

    if sr[0] == 0:  # SelHandles
        listIdx = NestedListSearch.findListIdx(selObjHdlCounts, idx)
        obj, splineIdx, ptIdx = selObjHdlList[listIdx]
        # ~ obj, splineIdx, ptIdx = selObjHdlList[int(idx / 2)]
        return retId, obj, splineIdx, ptIdx, 2 * (idx % 2)

    elif sr[0] == 1:  # CurveBezPt
        listIdx = NestedListSearch.findListIdx(objBezPtCounts, idx)
        obj, splineIdx, ptIdx = objLocList[listIdx]
        # ~ obj, splineIdx, ptIdx = objLocList[int(idx / ptIdxCnt)]
        return retId, obj, splineIdx, ptIdx, 1  # otherInfo = segIdx

    elif sr[0] == 2:  # SegLoc
        listIdx = NestedListSearch.findListIdx(selObjInterpCounts, idx)
        obj, splineIdx, segIdx = selObjLocList[listIdx]
        return retId, obj, splineIdx, segIdx, segInterpLocs[idx]

    else:  # CurveLoc
        listIdx = NestedListSearch.findListIdx(objInterpCounts, idx)
        obj, splineIdx, segIdx = objLocList[listIdx]
        return retId, obj, splineIdx, segIdx, objInterpLocs[idx]


class NestedListSearch:
    # Find the list element containing the given idx from flattened list
    # return the index of the list element containing the idx
    def findListIdx(counts, idx):
        cumulCnt = 0
        cntIdx = 0
        while idx >= cumulCnt:
            cumulCnt += counts[cntIdx]  # cntIdx can never be >= len(counts)
            cntIdx += 1
        return cntIdx - 1

    def __init__(self, ptsList):
        self.ptsList = ptsList
        self.kd = mathutils.kdtree.KDTree(sum(len(pts) for pts in ptsList))
        idx = 0
        self.counts = []
        for i, pts in enumerate(ptsList):
            self.counts.append(len(pts))
            for j, pt in enumerate(pts):
                self.kd.insert(pt, idx)
                idx += 1
        self.kd.balance()

    def findInLists(self, coFind, searchRange):
        if searchRange is None:
            foundVals = [self.kd.find(coFind)]
        else:
            foundVals = self.kd.find_range(coFind, searchRange)
            foundVals = sorted(foundVals, key=lambda x: x[2])

        searchResults = []
        for co, idx, dist in foundVals:
            listIdx = NestedListSearch.findListIdx(self.counts, idx)
            ptIdxInList = idx - sum(len(self.ptsList[i]) for i in range(0, listIdx))
            searchResults.append([listIdx, ptIdxInList, co, dist])

        return searchResults


class SelectCurveInfo:
    def __init__(self, obj, splineIdx):
        self.obj = obj
        self.splineIdx = splineIdx
        self.updateWSData()

        # User Selection (mouse click); format ptIdx: set(sel)...
        # where sel: -1->seg, 0->left hdl, 1->bezier pt, 2->right hdl
        self.ptSels = {}

        # Highlighted point (mouse move)
        # 'ptIdx': ptIdx, 'hltIdx':hltIdx {-1, 0, 1, 2} (just as in sel above)
        self.hltInfo = {}

        # obj.name gives exception if obj is not in bpy.data.objects collection,
        # so keep a copy
        self.objName = obj.name
        self.interpPts = {}

        # Format 'ptIdx': segIdx, 'hdlIdx': hdlIdx, 'loc':loc, 't':t
        # hdlIdx - {-1, 0, 1, 2} similar to sel in ptSels
        self.clickInfo = {}

    def __hash__(self):
        return hash((self.objName, self.splineIdx))

    def updateWSData(self):
        self.hasShapeKey = self.obj.active_shape_key is not None
        self.shapeKeyIdx = self.obj.active_shape_key_index if self.hasShapeKey else -1

        # for shape keys
        self.keyStartIdx = sum(
            len(self.obj.data.splines[i].bezier_points) for i in range(self.splineIdx)
        )

        # WS Data of the shape key (if exists)
        self.wsData = getBptData(self.obj, fromMix=False)[self.splineIdx]

    # For convenience
    def getAdjIdx(self, ptIdx, offset=1):
        return getAdjIdx(self.obj, self.splineIdx, ptIdx, offset)

    def getBezierPt(self, ptIdx):
        return self.obj.data.splines[self.splineIdx].bezier_points[ptIdx]

    def getShapeKeyData(self, ptIdx, keyIdx=None):
        if not self.hasShapeKey:
            return None
        if keyIdx is None:
            keyIdx = self.obj.active_shape_key_index
        keydata = self.obj.data.shape_keys.key_blocks[keyIdx].data
        keyIdx = self.keyStartIdx + ptIdx
        return keydata[keyIdx] if (keyIdx < len(keydata)) else None

    def getAllShapeKeysData(self, ptIdx):
        if not self.hasShapeKey:
            return None
        pts = []
        for keyIdx in range(len(self.obj.data.shape_keys.key_blocks)):
            pts.append(self.getShapeKeyData(ptIdx, keyIdx))
        return pts

    def getSegPtsInfo(self, ptIdx):
        nextIdx = self.getAdjIdx(ptIdx)
        pt0 = self.wsData[ptIdx]
        pt1 = self.wsData[nextIdx]
        return nextIdx, pt0, pt1

    def getSegPts(self, ptIdx):
        nextIdx, pt0, pt1 = self.getSegPtsInfo(ptIdx)
        return (pt0, pt1)

    # All selected points which have handles displayed (for example for snapping)
    def getAllPtsWithHdls(self):
        ptIdxs = set(self.ptSels.keys())
        nextIdxs = set(self.getSegPtsInfo(p)[0] for p in ptIdxs if -1 in self.ptSels[p])
        ptIdxs = ptIdxs.union(nextIdxs)
        return sorted(self.wsData[p] for p in ptIdxs)

    def setClickInfo(self, ptIdx, hdlIdx, clickLoc, lowerT=0.001, higherT=0.999):
        self.clickInfo = None
        if clickLoc is not None:
            nextIdx, pt0, pt1 = self.getSegPtsInfo(ptIdx)
            t = getTForPt([pt0[1], pt0[2], pt1[0], pt1[1]], clickLoc)

            if t is not None and (t < lowerT or t > higherT):
                hdlIdx = 1
                if t > higherT:
                    ptIdx = nextIdx
            else:
                self.clickInfo = {
                    "ptIdx": ptIdx,
                    "hdlIdx": hdlIdx,
                    "loc": clickLoc,
                    "t": t,
                }

        if self.clickInfo is None:
            self.clickInfo = {"ptIdx": ptIdx, "hdlIdx": hdlIdx}

    def addSel(self, ptIdx, sel, toggle=False):
        self.addSels(ptIdx, set([sel]), toggle)

    def addSels(self, ptIdx, sels, toggle=False):
        # TODO: Check this condition at other places
        if -1 in sels and self.getAdjIdx(ptIdx) is None:
            sels.remove(-1)
        if len(sels) == 0:
            return

        currSels = self.ptSels.get(ptIdx)
        if currSels is None:
            currSels = set()
        modSels = currSels.union(sels)
        if toggle:
            modSels -= currSels.intersection(sels)
        self.ptSels[ptIdx] = modSels
        if len(self.ptSels[ptIdx]) == 0:
            self.ptSels.pop(ptIdx)

    def removeSel(self, ptIdx, sel):
        self.removeSels(ptIdx, {sel})

    def removeSels(self, ptIdx, sels):
        for sel in sels:
            currSels = self.ptSels.get(ptIdx)
            if currSels is not None and sel in currSels:
                currSels.remove(sel)
                if len(currSels) == 0:
                    self.ptSels.pop(ptIdx)

    def resetClickInfo(self):
        self.clickInfo = {}

    def resetPtSel(self):
        self.ptSels = {}

    def resetHltInfo(self):
        self.hltInfo = {}

    def getHltInfo(self):
        return self.hltInfo

    def setHltInfo(self, ptIdx, hltIdx):
        self.hltInfo = {"ptIdx": ptIdx, "hltIdx": hltIdx}

    def getClickLoc(self):
        return self.clickInfo.get("loc")

    def getSelCo(self):
        if len(self.clickInfo) > 0:
            hdlIdx = self.clickInfo["hdlIdx"]
            if hdlIdx == -1:
                return self.clickInfo["loc"]
            else:
                ptIdx = self.clickInfo["ptIdx"]
                pt0 = self.wsData[ptIdx]
                return pt0[hdlIdx]

        return None

    def subdivSeg(self, subdivCnt):
        if self.hasShapeKey:
            return False
        if subdivCnt > 1:
            invMw = self.obj.matrix_world.inverted_safe()
            ts = []
            addCnt = 0
            for ptIdx in sorted(self.ptSels.keys()):
                if -1 in self.ptSels[ptIdx]:
                    vertCos = getInterpolatedVertsCo(self.interpPts[ptIdx], subdivCnt)[
                        1:-1
                    ]
                    changedIdx = ptIdx + addCnt
                    insertBezierPts(
                        self.obj,
                        self.splineIdx,
                        changedIdx,
                        [invMw @ v for v in vertCos],
                        "FREE",
                    )
                    addCnt += len(vertCos)
        return addCnt > 0

    def bevelPts(self, bevelCnt, deltaPos):
        if self.hasShapeKey:
            return False
        pts, ptSels = self.getBevelPts(bevelCnt, self.wsData, deltaPos)
        spline = self.obj.data.splines[self.splineIdx]
        newPtCnt = len(pts) - len(self.wsData)
        spline.bezier_points.add(newPtCnt)
        for pt in spline.bezier_points:
            pt.handle_left_type = "FREE"
            pt.handle_right_type = "FREE"

        invMw = self.obj.matrix_world.inverted_safe()
        for i, pt in enumerate(spline.bezier_points):
            pt.handle_left = invMw @ pts[i][0]
            pt.co = invMw @ pts[i][1]
            pt.handle_right = invMw @ pts[i][2]
            pt.handle_left_type = pts[i][3]
            pt.handle_right_type = pts[i][4]

        self.updateWSData()
        self.ptSels = ptSels
        return True

    def initSubdivMode(self, rv3d):
        if self.hasShapeKey:
            return False
        changed = False
        for ptIdx in self.ptSels.keys():
            if -1 in self.ptSels[ptIdx]:
                self.interpPts[ptIdx] = getPtsAlongBezier3D(
                    self.getSegPts(ptIdx), rv3d, curveRes=1000, minRes=1000
                )
                changed = True
        return changed

    def isBevelabel(self, rv3d):
        if self.hasShapeKey:
            return False
        changed = False
        for ptIdx in self.ptSels.keys():
            ptIdxs = [ptIdx]
            if -1 in self.ptSels[ptIdx]:
                ptIdxs.append(self.getAdjIdx(ptIdx))
            elif 1 not in self.ptSels[ptIdx]:
                continue  # only pt and seg selection
            for idx in ptIdxs:
                prevIdx = self.getAdjIdx(idx, -1)
                nextIdx = self.getAdjIdx(idx)
                if (
                    nextIdx is not None
                    and prevIdx is not None
                    and not hasAlignedHandles(self.wsData[idx])
                ):
                    changed = True
                    break
        return changed

    def getLastSegIdx(self):
        return getLastSegIdx(self.obj, self.splineIdx)

    # Remove all selected segments
    # Returns map with spline index and seg index change after every seg removal
    def removeSegs(self):
        changedSelMap = {}
        if self.hasShapeKey:
            return changedSelMap
        segSels = [p for p in self.ptSels if -1 in self.ptSels[p]]
        cumulSegIdxIncr = 0
        changedSplineIdx = self.splineIdx
        segIdxIncr = 0

        for segIdx in sorted(segSels):
            changedSegIdx = segIdx + cumulSegIdxIncr
            splineIdxIncr, segIdxIncr = removeBezierSeg(
                self.obj, changedSplineIdx, changedSegIdx
            )
            changedSplineIdx += splineIdxIncr
            cumulSegIdxIncr += segIdxIncr
            changedSelMap[segIdx] = [splineIdxIncr, segIdxIncr]
        return changedSelMap

    def straightenHandle(self, ptIdx, hdlIdx, allShapekeys=False):
        bpt = self.getBezierPt(ptIdx)
        prevIdx = self.getAdjIdx(ptIdx, -1)
        nextIdx = self.getAdjIdx(ptIdx)
        prevPts = None
        nextPts = None
        if self.hasShapeKey:
            if allShapekeys:
                pts = self.getAllShapeKeysData(ptIdx)
                if prevIdx is not None:
                    prevPts = self.getAllShapeKeysData(prevIdx)
                if nextIdx is not None:
                    nextPts = self.getAllShapeKeysData(nextIdx)
            else:
                pts = [self.getShapeKeyData(ptIdx)]
                if prevIdx is not None:
                    prevPts = [self.getShapeKeyData(prevIdx)]
                if nextIdx is not None:
                    nextPts = [self.getShapeKeyData(nextIdx)]

        else:
            pts = [bpt]
            if prevIdx is not None:
                prevPts = [self.getBezierPt(prevIdx)]
            if nextIdx is not None:
                nextPts = [self.getBezierPt(nextIdx)]

        if hdlIdx == 0:
            if bpt.handle_left_type != "VECTOR":
                bpt.handle_left_type = "FREE"
            for i in range(len(pts)):
                pt = pts[i]
                if prevPts is not None:
                    diffV = pt.co - prevPts[i].co
                else:
                    diffV = nextPts[i].co - pt.co
                pt.handle_left = pt.co - diffV / 3
        elif hdlIdx == 2:
            if bpt.handle_right_type != "VECTOR":
                bpt.handle_right_type = "FREE"
            for i in range(len(pts)):
                pt = pts[i]
                if nextPts is not None:
                    diffV = nextPts[i].co - pt.co
                else:
                    diffV = pt.co - prevPts[i].co
                pt.handle_right = pt.co + diffV / 3

    def straightenSelHandles(self):
        changed = False
        for ptIdx in self.ptSels:
            for hdlIdx in self.ptSels[ptIdx]:
                self.straightenHandle(ptIdx, hdlIdx)
                changed = True
        return changed

    def alignHandle(self, ptIdx, hdlIdx, allShapekeys=False):
        if hdlIdx == -1:
            return False
        oppIdx = 2 - hdlIdx
        if self.hasShapeKey:
            if allShapekeys:
                pts = self.getAllShapeKeysData(ptIdx)
            else:
                pts = [self.getShapeKeyData(ptIdx)]
        else:
            pts = [self.getBezierPt(ptIdx)]
        bpt = self.getBezierPt(ptIdx)
        if hdlIdx == 0 and bpt.handle_left_type != "ALIGNED":
            bpt.handle_left_type = "FREE"
        if hdlIdx == 2 and bpt.handle_right_type != "ALIGNED":
            bpt.handle_right_type = "FREE"

        for pt in pts:
            if hdlIdx == 0:
                pt.handle_left = (
                    pt.co
                    - (pt.co - pt.handle_left).length
                    * (pt.handle_right - pt.co).normalized()
                )
            else:
                pt.handle_right = (
                    pt.co
                    + (pt.co - pt.handle_right).length
                    * (pt.co - pt.handle_left).normalized()
                )
        return True

    def alignSelHandles(self):
        changed = False
        invMw = self.obj.matrix_world.inverted_safe()
        for ptIdx in self.ptSels:
            sels = self.ptSels[ptIdx]
            for hdlIdx in sels:
                changed = self.alignHandle(ptIdx, hdlIdx) or changed
        return changed

    def insertNode(self, handleType, select=True):
        if self.hasShapeKey:
            return False
        invMw = self.obj.matrix_world.inverted_safe()
        insertBezierPts(
            self.obj,
            self.splineIdx,
            self.clickInfo["ptIdx"],
            [invMw @ self.clickInfo["loc"]],
            handleType,
        )
        return True

    def removeNode(self):
        if self.hasShapeKey:
            return False
        changed = False

        toRemove = set()  # Bezier points to remove from object
        toRemoveSel = set()  # Selection entry to remove from ptSels

        nodeSels = [p for p in self.ptSels if 1 in self.ptSels[p]]

        for ptIdx in nodeSels:
            self.ptSels.pop(ptIdx)

        if len(nodeSels) > 0:
            removeBezierPts(self.obj, self.splineIdx, nodeSels)
            changed = True

        if changed:
            selIdxs = sorted(self.ptSels.keys())
            cnt = 0
            for ptIdx in nodeSels:
                cIdxs = [i for i in selIdxs if i >= (ptIdx - cnt)]
                for idx in cIdxs:
                    sels = self.ptSels.pop(idx - cnt)
                    self.ptSels[idx - cnt - 1] = sels
                cnt += 1

        return changed

    def getBevelPts(self, bevelCnt, pts, deltaPos):
        deltaLen = deltaPos.length
        if floatCmpWithMargin(deltaLen, DEF_ERR_MARGIN):
            return pts, self.ptSels

        # http://launchpadlibrarian.net/12692602/rcp.svg
        kFact = (bevelCnt / 3) * (sqrt(2) - 1)
        maxT = 0.5

        # Deep copy
        pts = [[p if isinstance(p, str) else p.copy() for p in pt] for pt in pts]

        newPts = []
        newSelPtIdxs = []

        # Add both points of the selected segments in selection
        bevelPtIdxs = set()
        for ptIdx in self.ptSels.keys():
            if 1 in self.ptSels[ptIdx]:
                bevelPtIdxs.add(ptIdx)
            if -1 in self.ptSels[ptIdx]:
                adjIdx = self.getAdjIdx(ptIdx)
                bevelPtIdxs.add(ptIdx)
                bevelPtIdxs.add(adjIdx)

        ptSels = {k: self.ptSels[k].copy() for k in self.ptSels.keys()}

        # Extra loop because next points need to be determined beforehand
        for ptIdx in bevelPtIdxs:
            if ptSels.get(ptIdx) is None:
                ptSels[ptIdx] = {1}
            else:
                ptSels[ptIdx].add(1)
            nextIdx = self.getAdjIdx(ptIdx)
            prevIdx = self.getAdjIdx(ptIdx, -1)
            if (
                prevIdx is not None
                and nextIdx is not None
                and not hasAlignedHandles(pts[ptIdx])
            ):
                newSelPtIdxs.append(ptIdx)

        for ptIdx, pt in enumerate(pts):
            if ptIdx in newSelPtIdxs:
                nextIdx = self.getAdjIdx(ptIdx)
                prevIdx = self.getAdjIdx(ptIdx, -1)
                prevPt = pts[prevIdx][:]
                diffV = pt[1] - prevPt[1]
                segLen = diffV.length
                if segLen < 0.001:
                    newPts.append(pt)
                else:
                    t = deltaLen / segLen
                    if t > maxT and (prevIdx in newSelPtIdxs):
                        t = maxT
                        k = kFact * (segLen / 2)
                    elif t > 1:
                        t = 1
                        k = kFact * segLen
                    else:
                        k = kFact * deltaLen

                    seg = [pts[prevIdx][1], pts[prevIdx][2], pt[0], pt[1]]
                    partialSeg = getPartialSeg(seg, t0=0, t1=1 - t)
                    newPt = partialSeg[3]

                    tangent0 = getTangentAtT(
                        pts[prevIdx][1], pts[prevIdx][2], pt[0], pt[1], 1 - t
                    )
                    pt0_2 = newPt + k * (tangent0.normalized())

                    pt0 = [partialSeg[2], newPt, pt0_2, "FREE", "FREE"]
                    newPts.append(pt0)

                    prevPt[2] = partialSeg[1]

                nextPt = pts[nextIdx][:]
                diffV = nextPt[1] - pt[1]
                segLen = diffV.length
                if floatCmpWithMargin(segLen, 0):
                    newPts.append(pt)
                else:
                    t = deltaLen / segLen
                    if t > maxT and (nextIdx in newSelPtIdxs):
                        t = maxT
                        k = kFact * (segLen / 2)
                    elif t > 1:
                        t = 1
                        k = kFact * segLen
                    else:
                        k = kFact * deltaLen

                    seg = [pt[1], pt[2], pts[nextIdx][0], pts[nextIdx][1]]
                    partialSeg = getPartialSeg(seg, t0=t, t1=1)
                    newPt = partialSeg[0]

                    tangent1 = getTangentAtT(
                        pt[1], pt[2], pts[nextIdx][0], pts[nextIdx][1], t
                    )

                    pt1_0 = newPt - k * (tangent1.normalized())
                    pt1 = [pt1_0, newPt, partialSeg[1], "FREE", "FREE"]

                    newPts.append(pt1)

                    nextPt[0] = partialSeg[2]
            else:
                newPts.append(pt)

        newPtSels = {}
        cnt = 0
        for ptIdx in sorted(ptSels.keys()):
            newPtSels[ptIdx + cnt] = ptSels[ptIdx].copy()
            if ptIdx in newSelPtIdxs:
                newPtSels[ptIdx + cnt].add(-1)
                newPtSels[ptIdx + cnt + 1] = {1}
                cnt += 1

        return newPts, newPtSels

    def getDisplayInfos(
        self, hideHdls=False, subdivCnt=0, bevelCnt=0, newPos=None, deltaPos=None
    ):
        # Making long short
        cHltTip = FTProps.colHltTip
        cBezPt = FTProps.colBezPt
        cHdlPt = FTProps.colHdlPtTip
        cAdjBezTip = FTProps.colAdjBezTip
        cNonHltSeg = FTProps.colDrawNonHltSeg

        segDispInfos = []
        bptDispInfos = []

        pts = self.wsData[:]
        ptSels = self.ptSels
        if newPos is not None:
            # TODO: This method is in EditCurveInfo
            nPtIdxs, nPts = self.getOffsetSegPts(newPos)
            # Update list with new position (editing)
            for i, ptIdx in enumerate(nPtIdxs):
                pts[ptIdx] = nPts[i]
        elif deltaPos is not None:
            pts, ptSels = self.getBevelPts(bevelCnt, pts, deltaPos)

        # Default display of spline
        for i, pt in enumerate(pts):
            bptDispInfos.append(BptDisplayInfo(pt, [cAdjBezTip]))
            if i > 0:
                segDispInfos.append(SegDisplayInfo([pts[i - 1], pt], cNonHltSeg))
        lastIdx = self.getAdjIdx(len(self.wsData) - 1)  # In case cyclic...
        if lastIdx is not None:
            segDispInfos.append(SegDisplayInfo([pts[-1], pts[0]], cNonHltSeg))

        hltInfo = self.getHltInfo()
        hltPtIdx = hltInfo.get("ptIdx")
        hltIdx = hltInfo.get("hltIdx")

        # Process highlighted segments before selected ones because...
        # selected segments take priority over highlighted
        if hltIdx == -1:
            segDispInfos[hltPtIdx].segColor = FTProps.colDrawHltSeg
            bptDispInfos[hltPtIdx].tipColors[1] = cBezPt
            nextIdx = self.getAdjIdx(hltPtIdx)
            bptDispInfos[nextIdx].tipColors[1] = cBezPt

        # Process selections
        for ptIdx in sorted(ptSels.keys()):
            sels = ptSels[ptIdx]

            if hideHdls:
                tipColors = [None, cBezPt, None]
                handleNos = []
            else:
                tipColors = [cHdlPt, cBezPt, cHdlPt]
                handleNos = [0, 1]

            bptDispInfos[ptIdx].tipColors = tipColors[:]
            bptDispInfos[ptIdx].handleNos = handleNos

            for hdlIdx in sorted(sels):  # Start with seg selection i. e. -1
                if hdlIdx == -1:
                    nextIdx = getAdjIdx(self.obj, self.splineIdx, ptIdx, ptCnt=len(pts))
                    segPts = [pts[ptIdx], pts[nextIdx]]

                    # process next only if there are no selection pts with that idx
                    if nextIdx not in ptSels.keys():
                        bptDispInfos[nextIdx].tipColors = tipColors[:]
                        bptDispInfos[nextIdx].handleNos = handleNos

                    vertCos = []
                    if subdivCnt > 1:
                        vertCos = getInterpolatedVertsCo(
                            self.interpPts[ptIdx], subdivCnt
                        )[1:-1]

                    selSegDispInfo = EditSegDisplayInfo(
                        segPts, FTProps.colDrawSelSeg, vertCos
                    )
                    segDispInfos[ptIdx] = selSegDispInfo
                elif hdlIdx == 1 or not hideHdls:
                    bptDispInfos[ptIdx].tipColors[hdlIdx] = FTProps.colSelTip

        # Process highlighted points after selected ones because...
        # highlighted points take priority over selected
        if hltIdx in {0, 1, 2}:
            bptDispInfos[hltPtIdx].tipColors[hltIdx] = cHltTip

        return [segDispInfos, bptDispInfos]


class EditCurveInfo(SelectCurveInfo):
    def __init__(self, obj, splineIdx, ptSels=None):
        super(EditCurveInfo, self).__init__(obj, splineIdx)
        if ptSels is not None:
            self.ptSels = ptSels

    def syncAlignedHdl(self, pt, ctrlPLoc, hdlIdx):
        typeIdx = 3 if hdlIdx == 0 else 4
        if pt[typeIdx] == "ALIGNED":
            oppTypeIdx = 4 if hdlIdx == 0 else 3
            if pt[oppTypeIdx] in {"VECTOR", "ALIGNED"}:
                oppHdlIdx = 2 if hdlIdx == 0 else 0
                oppHdlV = ctrlPLoc - pt[oppHdlIdx]
                if oppHdlV.length != 0:
                    currL = (ctrlPLoc - pt[hdlIdx]).length
                    pt[hdlIdx] = ctrlPLoc + currL * oppHdlV / oppHdlV.length

    def setAlignedHdlsCo(self, pt, hdlIdx, ctrlPLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        if pt[typeIdx] == "ALIGNED":
            oppTypeIdx = 4 if hdlIdx == 0 else 3
            if pt[oppTypeIdx] != "VECTOR":
                pt[hdlIdx] += ctrlPLoc - pt[1]
            else:
                self.syncAlignedHdl(pt, ctrlPLoc, hdlIdx)

    def setFreeHdlsCo(self, pt, hdlIdx, newLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        if pt[typeIdx] == "FREE":
            pt[hdlIdx] += newLoc - pt[1]

    def setVectHdlsCo(self, pt, newLoc, hdlIdx, prevPt, nextPt):
        typeIdx = 3 if hdlIdx == 0 else 4
        if pt[typeIdx] == "VECTOR":
            typeIdx = 3 if hdlIdx == 0 else 4
            pts = [prevPt, nextPt] if hdlIdx == 0 else [nextPt, prevPt]
            diffV = None
            if pts[0] is not None:
                diffV = pts[0][1] - newLoc
            if diffV is None and pts[1] is not None:
                diffV = newLoc - pts[1][1]
            if diffV is None:
                pt[hdlIdx] = newLoc
            else:
                pt[hdlIdx] = newLoc + diffV * 1 / 3

    # Calculate both handle and adjacent pt handles in case of Vector type
    # TODO: AUTO has a separate logic set to ALIGNED for now
    def syncCtrlPtHdls(self, ptIdx, newLoc):
        wsData = getBptData(self.obj, fromMix=False)
        pt = wsData[self.splineIdx][ptIdx]
        prevIdx = self.getAdjIdx(ptIdx, -1)
        prevPt = None if prevIdx is None else wsData[self.splineIdx][prevIdx]
        nextIdx = self.getAdjIdx(ptIdx)
        nextPt = None if nextIdx is None else wsData[self.splineIdx][nextIdx]

        ptIdxs = [ptIdx]
        pts = [pt]

        for typeIdx in [3, 4]:
            if pt[typeIdx] == "AUTO":
                pt[typeIdx] = "ALIGNED"
        for hdlIdx in [0, 2]:
            self.setVectHdlsCo(pt, newLoc, hdlIdx, prevPt, nextPt)
        for hdlIdx in [0, 2]:
            self.setFreeHdlsCo(pt, hdlIdx, newLoc)
        for hdlIdx in [0, 2]:
            self.setAlignedHdlsCo(pt, hdlIdx, newLoc)

        pt[1] = newLoc

        if prevPt is not None and prevPt[4] == "VECTOR":
            pPrevIdx = self.getAdjIdx(prevIdx, -1)
            pPrevPt = None if pPrevIdx is None else wsData[self.splineIdx][pPrevIdx]
            self.setVectHdlsCo(prevPt, prevPt[1], 2, pPrevPt, pt)
            self.setAlignedHdlsCo(prevPt, 0, prevPt[1])

            ptIdxs.append(prevIdx)
            pts.append(prevPt)

        if nextPt is not None and nextPt[3] == "VECTOR":
            nNextIdx = self.getAdjIdx(nextIdx)
            nNextPt = None if nNextIdx is None else wsData[self.splineIdx][nNextIdx]
            self.setVectHdlsCo(nextPt, nextPt[1], 0, pt, nNextPt)
            self.setAlignedHdlsCo(nextPt, 2, nextPt[1])

            ptIdxs.append(nextIdx)
            pts.append(nextPt)

        return ptIdxs, pts

    # Calculate the opposite handle values in case of ALIGNED and AUTO handles
    # Also set the type(s) of current (opposite) handle(s)
    def syncHdls(self, pt, hdlIdx, newLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        oppTypeIdx = 4 if hdlIdx == 0 else 3

        if pt[typeIdx] == "VECTOR":
            pt[typeIdx] = "FREE"
        if pt[typeIdx] == "AUTO":
            pt[typeIdx] = "ALIGNED"
        if pt[oppTypeIdx] == "AUTO" and pt[typeIdx] != "FREE":
            pt[oppTypeIdx] = "ALIGNED"

        pt[hdlIdx] = newLoc

        self.syncAlignedHdl(pt, pt[1], 2 - hdlIdx)  # First opposite

        self.syncAlignedHdl(pt, pt[1], hdlIdx)

    # Get seg points after change in position of handles or drag curve
    # The only function called on all 3 events: grab curve pt, grab handle, grab Bezier pt
    def getOffsetSegPts(self, newLoc):
        inf = self.clickInfo
        ptIdx = inf["ptIdx"]
        hdlIdx = inf["hdlIdx"]
        wsData = getBptData(self.obj, fromMix=False)
        pt = wsData[self.splineIdx][ptIdx]

        if hdlIdx == -1:  # Grab point on curve
            adjIdx = self.getAdjIdx(ptIdx)
            adjPt = wsData[self.splineIdx][adjIdx]

            ptIdxs = [ptIdx, adjIdx]
            pts = [pt, adjPt]

            delta = newLoc - inf["loc"]
            if delta == 0:
                return ptIdxs, pts
            t = inf["t"]

            # ****************************************************************
            # Magic Bezier Drag Equations (Courtesy: Inkscape)             #*
            # ****************************************************************
            # *
            if t <= 1.0 / 6.0:  # *
                weight = 0  # *
            elif t <= 0.5:  # *
                weight = (pow((6 * t - 1) / 2.0, 3)) / 2  # *
            elif t <= 5.0 / 6.0:  # *
                weight = (1 - pow((6 * (1 - t) - 1) / 2.0, 3)) / 2 + 0.5  # *
            else:  # *
                weight = 1  # *
                # *
            offset0 = ((1 - weight) / (3 * t * (1 - t) * (1 - t))) * delta  # *
            offset1 = (weight / (3 * t * t * (1 - t))) * delta  # *
            # *
            # ****************************************************************

            # If the segment is edited, the 1st pt right handle...
            pts[0][2] += offset0

            if pts[0][4] == "VECTOR":
                pts[0][4] = "FREE"
            if pts[0][4] == "AUTO":
                pts[0][4] = "ALIGNED"
            # opposite handle must be changed if this is not FREE
            if pts[0][3] == "VECTOR" and pts[0][4] != "FREE":
                pts[0][3] = "FREE"
            if pts[0][3] == "AUTO" and pts[0][4] != "FREE":
                pts[0][3] = "ALIGNED"

            self.syncAlignedHdl(pts[0], pts[0][1], hdlIdx=0)

            # ...and 2nd pt left handle impacted
            pts[1][0] += offset1

            if pts[1][3] == "VECTOR":
                pts[1][3] = "FREE"
            if pts[1][3] == "AUTO":
                pts[1][3] = "ALIGNED"
            # opposite handle must be changed if this is not FREE
            if pts[1][4] == "VECTOR" and pts[1][3] != "FREE":
                pts[1][4] = "FREE"
            if pts[1][4] == "AUTO" and pts[1][3] != "FREE":
                pts[1][4] = "ALIGNED"

            self.syncAlignedHdl(pts[1], pts[1][1], hdlIdx=2)
            return ptIdxs, pts

        elif hdlIdx in {0, 2}:  # Grab one of the handles
            self.syncHdls(pt, hdlIdx, newLoc)
            return [ptIdx], [pt]

        else:  # Grab the Bezier point
            return self.syncCtrlPtHdls(ptIdx, newLoc)

    def moveSeg(self, newPos):
        ptIdxs, pts = self.getOffsetSegPts(newPos)

        invMw = self.obj.matrix_world.inverted_safe()
        spline = self.obj.data.splines[self.splineIdx]
        bpts = [spline.bezier_points[idx] for idx in ptIdxs]

        for i, bpt in enumerate(bpts):
            bpt.handle_right_type = "FREE"
            bpt.handle_left_type = "FREE"

        if self.hasShapeKey:
            for i, ptIdx in enumerate(ptIdxs):
                keydata = self.getShapeKeyData(ptIdx)
                keydata.handle_left = invMw @ pts[i][0]
                keydata.co = invMw @ pts[i][1]
                keydata.handle_right = invMw @ pts[i][2]
                if pts[i][3] == "AUTO":
                    pts[i][3] = "ALIGNED"
                if pts[i][4] == "AUTO":
                    pts[i][4] = "ALIGNED"
                impIdxs = [ptIdx, self.getAdjIdx(ptIdx, -1), self.getAdjIdx(ptIdx, 1)]
                for idx in impIdxs:
                    if idx is None:
                        continue
                    if spline.bezier_points[idx].handle_left_type == "AUTO":
                        spline.bezier_points[idx].handle_left_type = "ALIGNED"
                    if spline.bezier_points[idx].handle_right_type == "AUTO":
                        spline.bezier_points[idx].handle_right_type = "ALIGNED"
        else:
            for i, bpt in enumerate(bpts):
                bpt.handle_left = invMw @ pts[i][0]
                bpt.co = invMw @ pts[i][1]
                bpt.handle_right = invMw @ pts[i][2]

        for i, bpt in enumerate(bpts):
            bpt.handle_left_type = pts[i][3]
            bpt.handle_right_type = pts[i][4]

        self.updateWSData()


class ModalFlexiEditBezierOp(ModalBaseFlexiOp):
    bl_description = "Flexi editing of Bezier curves in object mode"
    bl_idname = "wm.modal_flexi_edit_bezier"
    bl_label = "Flexi Edit Curve"
    bl_options = {"REGISTER", "UNDO"}

    h = False

    def drawHandler():
        ModalBaseFlexiOp.drawHandlerBase()

    def resetDisplay():
        ModalBaseFlexiOp.resetDisplayBase()

    # static method
    def refreshDisplay(segDispInfos, bptDispInfos, locOnCurve=None, snapper=None):
        ptCos = [
            co
            for d in segDispInfos
            if isinstance(d, EditSegDisplayInfo)
            for co in d.subdivCos
        ]

        # ~ if(locOnCurve != None): ptCos.append(locOnCurve) # For debugging

        ModalBaseFlexiOp.bglDrawMgr.addPtInfo(
            "editSubdiv", FTProps.editSubdivPtSize, [FTProps.colEditSubdiv], ptCos
        )

        ModalBaseFlexiOp.refreshDisplayBase(segDispInfos, bptDispInfos, snapper)

    def getToolType(self):
        return TOOL_TYPE_FLEXI_EDIT

    # Refresh display with existing curves (nonstatic)
    def refreshDisplaySelCurves(
        self,
        hltSegDispInfos=None,
        hltBptDispInfos=None,
        locOnCurve=None,
        refreshPos=False,
    ):
        if self.rmInfo is None:
            return  # Possible in updateAfterGeomChange

        newPos = None
        if FTProps.liveUpdate and self.editCurveInfo is not None:
            newPos = self.getNewPos(refreshStatus=True)
            self.editCurveInfo.moveSeg(newPos)
            clickInfo = self.editCurveInfo.clickInfo
            if clickInfo["hdlIdx"] == -1:
                self.editCurveInfo.setClickInfo(
                    clickInfo["ptIdx"], clickInfo["hdlIdx"], newPos
                )
            self.xyPress = self.rmInfo.xy[:]

        segDispInfos = []
        bptDispInfos = []
        # ~ curveInfos = self.selectCurveInfos.copy()
        # ~ if(self.editCurveInfo != None):
        # ~ curveInfos.add(self.editCurveInfo)
        if self.bevelMode:
            deltaPos = self.getNewDeltaPos(refreshStatus=True)
        else:
            deltaPos = None
        for c in self.selectCurveInfos:
            if refreshPos and c == self.editCurveInfo and newPos is None:
                newPos = self.getNewPos(refreshStatus=True)
            else:
                newPos = None
            info1, info2 = c.getDisplayInfos(
                hideHdls=ModalFlexiEditBezierOp.h,
                subdivCnt=self.subdivCnt,
                bevelCnt=self.bevelCnt,
                newPos=newPos,
                deltaPos=deltaPos,
            )
            segDispInfos += info1
            bptDispInfos += info2

        # Highlighted at the top
        if hltSegDispInfos is not None:
            segDispInfos += hltSegDispInfos
        if hltBptDispInfos is not None:
            bptDispInfos += hltBptDispInfos

        ModalFlexiEditBezierOp.refreshDisplay(
            segDispInfos, bptDispInfos, locOnCurve, self.snapper
        )

    def reset(self):
        self.editCurveInfo = None
        self.selectCurveInfos = set()
        # TODO: freezeOrient logic should be internal to Snapper
        if self.snapper is not None:
            self.snapper.freezeOrient = False
        ModalFlexiEditBezierOp.resetDisplay()

    def postUndoRedo(self, scene, dummy=None):  # signature different in 2.8 and 2.81?
        # ~ self.snapper.customAxis.reload()
        self.updateAfterGeomChange()
        for ci in self.selectCurveInfos:
            ci.resetPtSel()

    def cancelOp(self, context):
        self.reset()
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)
        bpy.app.handlers.depsgraph_update_post.remove(self.updateAfterGeomChange)
        return self.cancelOpBase()

    def isToolSelected(self, context):
        if context.mode != "OBJECT":
            return False

        tool = context.workspace.tools.from_space_view3d_mode("OBJECT", create=False)
        if tool is None or tool.idname != FlexiEditBezierTool.bl_idname:
            # if(tool == None or tool.idname != 'flexi_bezier.edit_tool'):
            return False
        return True

    # Will be called after the curve is changed (by the tool or externally)
    # So handle all possible conditions
    def updateAfterGeomChange(self, scene=None, dummy=None):  # 3 params in 2.81
        ciRemoveList = []

        removeObjNames = set()  # For snaplocs
        addObjNames = set()
        self.htlCurveInfo = None

        # TODO: check if self.editCurveInfo is to be set to None
        if not FTProps.liveUpdate:
            self.editCurveInfo = None  # Reset if editing (capture == True)

        for ci in self.selectCurveInfos:
            if bpy.data.objects.get(ci.objName) is not None:
                ci.obj = bpy.data.objects.get(ci.objName)  # refresh anyway
                splines = ci.obj.data.splines
                if ci.splineIdx >= len(ci.obj.data.splines):
                    ciRemoveList.append(ci)
                    continue
                spline = splines[ci.splineIdx]
                bpts = spline.bezier_points
                bptsCnt = len(bpts)
                # Don't keep a point object / spline
                if bptsCnt <= 1:
                    if len(splines) == 1:
                        ciRemoveList.append(ci)
                    else:
                        splines.remove(spline)
                        if ci.splineIdx >= (len(splines)):
                            ci.splineIdx = len(splines) - 1
                        ci.resetPtSel()
                else:
                    # If any of the current selections is / are...
                    # greater than last idx (pt) or last but one (seg)...
                    # move it / them to last idx (pt) or last but one idx (seg)
                    changeSegSels = set()
                    changePtSels = set()
                    lastIdx = bptsCnt - 1
                    lastSegIdx = ci.getLastSegIdx()
                    for ptIdx in ci.ptSels.keys():
                        sels = ci.ptSels[ptIdx]
                        if -1 in sels and ptIdx > lastSegIdx:
                            changeSegSels.add(ptIdx)
                            sels.remove(-1)
                        if ptIdx > lastIdx and len(sels) > 0:
                            changePtSels.add(ptIdx)
                    for ptIdx in changePtSels:
                        sels = ci.ptSels.pop(ptIdx)
                        ci.addSels(lastSegIdx, sels)
                    for pt in changeSegSels:
                        ci.addSels((lastIdx - 1), set([-1]))

                addObjNames.add(ci.objName)
                ci.updateWSData()
            else:
                ciRemoveList.append(ci)
                removeObjNames.add(ci.objName)

        if len(ciRemoveList) > 0:
            for c in ciRemoveList:
                self.selectCurveInfos.remove(c)

        if self.editCurveInfo is None:  # exclude live update condition
            self.updateSnapLocs(addObjNames, removeObjNames)
            self.refreshDisplaySelCurves()

    def subInvoke(self, context, event):
        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

        self.editCurveInfo = None
        self.htlCurveInfo = None
        self.selectCurveInfos = set()
        self.subdivCnt = 0
        self.bevelCnt = 4
        self.bevelMode = False

        # For double click (TODO: remove; same as editCurveInfo == None?)
        self.capture = False
        self.xyPress = None  # ...to avoid jerky movement at the beginning
        self.xyLoc = None  # for bevel

        self.snapInfos = {}
        self.updateSnapLocs()

        return {"RUNNING_MODAL"}

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if not ModalFlexiEditBezierOp.h:
            for ci in self.selectCurveInfos:
                pts = ci.getAllPtsWithHdls()
                for pt in pts:  # Already world space
                    locs.append(pt[0])
                    locs.append(pt[2])
        return locs

    def updateSnapLocs(self, addObjNames=None, removeObjNames=None):
        updateCurveEndPtMap(self.snapInfos, addObjNames, removeObjNames)

    def getRefLine(self):
        if self.editCurveInfo is not None:
            ei = self.editCurveInfo
            ptIdx = ei.clickInfo["ptIdx"]
            hdlIdx = ei.clickInfo["hdlIdx"]
            pt0 = ei.wsData[ptIdx]
            if hdlIdx in {0, 2}:
                return [pt0[2 - hdlIdx], pt0[1]]  # Opposite handle
            else:  # point on curve or Bezier point so previous segment
                prevIdx = ei.getAdjIdx(ptIdx, -1)
                pPrevIdx = ei.getAdjIdx(ptIdx, -2)
                if prevIdx is not None and pPrevIdx is not None:
                    return [ei.wsData[pPrevIdx][1], ei.wsData[prevIdx][1]]
                else:
                    nextIdx = ei.getAdjIdx(ptIdx, 1)
                    nNextIdx = ei.getAdjIdx(ptIdx, 2)
                    if nextIdx is not None and nNextIdx is not None:
                        return [ei.wsData[nNextIdx][1], ei.wsData[nextIdx][1]]
        return self.getCurrLine()

    def getCurrLine(self):
        ei = self.editCurveInfo
        if ei is not None:
            ptIdx = ei.clickInfo["ptIdx"]
            hdlIdx = ei.clickInfo["hdlIdx"]
            pt0 = ei.wsData[ptIdx]
            clickLoc = ei.getClickLoc()
            if clickLoc is not None:
                return [pt0[1], clickLoc]
            if hdlIdx in {0, 2}:
                return [pt0[1], pt0[hdlIdx]]  # Current handle
            elif hdlIdx == 1:
                adjIdx = ei.getAdjIdx(ptIdx, -1)
                if adjIdx is None:
                    adjIdx = ei.getAdjIdx(ptIdx, 1)
                if adjIdx is None:
                    return [pt0[1]]
                else:
                    return [ei.wsData[adjIdx][1], pt0[1]]
        return []

    def getRefLineOrig(self):
        ei = self.editCurveInfo
        refLine = self.getRefLine()
        if ei is not None and len(refLine) > 0:
            return refLine[-1]
        return None

    def getSelCo(self):
        if self.editCurveInfo is not None:
            return self.editCurveInfo.getSelCo()
        return None

    def getEditableCurveObjs(self):
        return [
            b
            for b in bpy.data.objects
            if isBezier(b)
            and b.visible_get()
            and not b.hide_select
            and len(b.data.splines[0].bezier_points) > 1
        ]

    def getSearchQueryInfo(self):  # TODO: Simplify if possible
        queryInfo = {}
        for ci in self.selectCurveInfos:
            info = queryInfo.get(ci.obj)
            if info is None:
                info = {}
                queryInfo[ci.obj] = info

            segPtIdxs = info.get(ci.splineIdx)
            if segPtIdxs is None:
                # First is for seg search, second for handles
                segPtIdxs = [[], []]
                info[ci.splineIdx] = segPtIdxs

            for ptIdx in ci.ptSels.keys():
                sels = ci.ptSels[ptIdx]
                if -1 in sels:
                    # TODO: Could be duplicate
                    adjIdx = ci.getAdjIdx(ptIdx)
                    segPtIdxs[0].append(ptIdx)
                    segPtIdxs[1].append(ptIdx)
                    segPtIdxs[1].append(adjIdx)
                else:
                    segPtIdxs[1].append(ptIdx)

        return queryInfo

    def getSelInfoObj(self, obj, splineIdx):
        for ci in self.selectCurveInfos:
            if ci.obj == obj and ci.splineIdx == splineIdx:
                return ci
        return None

    # Delete selected segments and synchronize remaining selections
    # TODO: Way too complicated, maybe there exists a much simpler way to do this
    def delSelSegs(self):
        changed = False
        curveInfoList = sorted(
            self.selectCurveInfos, key=lambda x: (x.objName, x.splineIdx)
        )

        # Process one spline at a time
        for cIdx, c in enumerate(curveInfoList):
            c.resetHltInfo()

            spline = c.obj.data.splines[c.splineIdx]
            wasCyclic = spline.use_cyclic_u
            oldPtCnt = len(spline.bezier_points)

            changedSelMap = c.removeSegs()

            if len(changedSelMap) == 0:
                continue
            changed = True

            # Shift all the splineIdxs after the changed one by spline incr count
            totalSplineIdxIncr = sum(x[0] for x in changedSelMap.values())

            # Order doesn't matter (different curveInfo)
            for i in range(cIdx + 1, len(curveInfoList)):
                if curveInfoList[i].objName != c.objName:
                    break
                curveInfoList[i].splineIdx += totalSplineIdxIncr

            # TODO: Remove the try after sufficient testing (or better replacement)
            # Exception means no selected points will be deleted
            try:
                # Copy old selections as they will change
                oIdxs = sorted(c.ptSels.keys())
                ptSelsCopy = c.ptSels.copy()
                newSplineIdx = c.splineIdx

                c.resetPtSel()
                currCurveInfo = c

                # Reflects new selection after every seg removal
                modifiedSegIdxs = {idx: idx for idx in oIdxs}

                # First get the segment selections out of the way
                for i, segIdx in enumerate(sorted(changedSelMap.keys())):
                    ptSelsCopy[segIdx].remove(-1)
                    if len(ptSelsCopy[segIdx]) == 0:
                        ptSelsCopy.pop(segIdx)

                # Each of the 'if' blocks in the loop iterate over all the selections
                # and update them iteratively for each seg removal from spline
                # The updated selected seg idx for each iteration is in modifiedSegIdxs

                # segIdx and oIdx don't change throughout
                # they always refer to the selections that were there before removal
                for i, segIdx in enumerate(sorted(changedSelMap.keys())):
                    splineIdxIncr = changedSelMap[segIdx][0]
                    segIdxIncr = changedSelMap[segIdx][1]

                    # This will be executed only once (if at all), at first iteration
                    if wasCyclic and i == 0:
                        for j, oIdx in enumerate(oIdxs):
                            # First iteration, so no need to refer to modifiedSegIdxs
                            newSegIdx = oIdx + segIdxIncr
                            if newSegIdx < 0:
                                newSegIdx += oldPtCnt
                            modifiedSegIdxs[oIdx] = newSegIdx
                            if ptSelsCopy.get(oIdx) is not None:
                                currCurveInfo.ptSels[newSegIdx] = ptSelsCopy[
                                    oIdx
                                ].copy()

                    # 'removed' segment at one of the either ends
                    elif splineIdxIncr == 0:
                        ptCnt = len(c.obj.data.splines[newSplineIdx].bezier_points)
                        for j, oIdx in enumerate(oIdxs):
                            prevIdx = modifiedSegIdxs[oIdx]
                            # segIdxIncr: only two values possible: 0, -1
                            newSegIdx = prevIdx + segIdxIncr
                            # ~ if(currCurveInfo.ptSels.get(prevIdx) != None):
                            # ~ currCurveInfo.ptSels.pop(prevIdx)
                            if (
                                ptSelsCopy.get(oIdx) is not None
                                and newSegIdx >= 0
                                and newSegIdx < ptCnt
                            ):
                                currCurveInfo.ptSels[newSegIdx] = ptSelsCopy[
                                    oIdx
                                ].copy()
                            modifiedSegIdxs[oIdx] = newSegIdx

                    # Most likely condition
                    elif splineIdxIncr > 0:
                        splineCnt = len(c.obj.data.splines)
                        prevCurveInfo = currCurveInfo
                        newSplineIdx += 1
                        # No overwriting since the higher splineIdxs already moved above
                        # But it's possible this spline was removed in subsequent
                        # iterations by removeSegs, so check...
                        if newSplineIdx < splineCnt:
                            currCurveInfo = SelectCurveInfo(c.obj, newSplineIdx)
                            self.selectCurveInfos.add(currCurveInfo)
                        # idxs and prevCurve have to be updated so continue
                        else:
                            currCurveInfo = None  # Fail fast

                        for oIdx in oIdxs:
                            prevIdx = modifiedSegIdxs[oIdx]
                            # If prevIdx itself is negative, this is previous to previous
                            # So won't change
                            if prevIdx < 0:
                                continue
                            newSegIdx = prevIdx + segIdxIncr
                            # newSegIdx negative... first part of the split spline
                            if newSegIdx < 0 and ptSelsCopy.get(oIdx) is not None:
                                prevCurveInfo.ptSels[prevIdx] = ptSelsCopy[oIdx].copy()
                            # newSegIdx positive... second part of the split spline
                            elif ptSelsCopy.get(oIdx) is not None and newSegIdx >= 0:
                                if newSplineIdx < splineCnt:
                                    currCurveInfo.ptSels[newSegIdx] = ptSelsCopy[
                                        oIdx
                                    ].copy()
                                if prevCurveInfo.ptSels.get(prevIdx) is not None:
                                    prevCurveInfo.ptSels.pop(prevIdx)
                            modifiedSegIdxs[oIdx] = newSegIdx

                    elif splineIdxIncr < 0:
                        # This is not the same as c
                        # (could be a new spline added in between)
                        toRemList = [
                            x
                            for x in self.selectCurveInfos
                            if x.splineIdx == newSplineIdx
                        ]
                        if len(toRemList) > 0:
                            self.selectCurveInfos.remove(toRemList[0])

            except Exception as e:
                c.resetPtSel()

        return changed

    def mnSelect(self, opt):
        if opt[0] == "miSelAllSplines":
            curves = self.getEditableCurveObjs()
            allCurveInfos = [
                SelectCurveInfo(curve, i)
                for curve in curves
                for i in range(len(curve.data.splines))
            ]
            for c in allCurveInfos:
                if c not in self.selectCurveInfos:
                    self.selectCurveInfos.add(c)
            # ~ for c in allCurveInfos:
            # ~ for ptIdx in range(len(c.wsData)):
            # ~ c.addSel(ptIdx, 1)
        else:
            h = ModalFlexiEditBezierOp.h
            self.selHltInfo(makeActive=True)
            for i, c in enumerate(self.selectCurveInfos):
                if opt[0] == "miSelObj":
                    c.obj.select_set(True)
                    if (
                        self.htlCurveInfo is None
                        and i == len(self.selectCurveInfos) - 1
                    ):
                        bpy.context.view_layer.objects.active = c.obj
                else:
                    for ptIdx in range(len(c.wsData)):
                        if opt[0] == "miSelSegs":
                            c.addSel(ptIdx, -1)
                        if opt[0] == "miSelBezPts":
                            c.addSel(ptIdx, 1)
                        if opt[0] == "miSelHdls" and not h:
                            c.addSels(ptIdx, {0, 2})
                        if opt[0] == "miSelAll":
                            c.addSels(ptIdx, {-1, 1}.union({0, 2} if not h else set()))
        self.htlCurveInfo = None

    def mnDeselect(self, opt):
        h = ModalFlexiEditBezierOp.h
        if opt[0] == "miDeselObj":
            if self.htlCurveInfo is not None:
                self.htlCurveInfo.obj.select_set(False)
                self.htlCurveInfo = None
        for c in self.selectCurveInfos:
            if opt[0] == "miDeselObj":
                c.obj.select_set(False)
            else:
                for ptIdx in range(len(c.wsData)):
                    if opt[0] == "miDeselSegs":
                        c.removeSel(ptIdx, -1)
                    if opt[0] == "miDeselBezPts":
                        c.removeSel(ptIdx, 1)
                    if opt[0] == "miDeselHdls" and not h:
                        c.removeSels(ptIdx, {0, 2})
                    if opt[0] == "miDeselInvert":
                        c.addSels(
                            ptIdx,
                            {-1, 1}.union({0, 2} if not h else set()),
                            toggle=True,
                        )

    def mnSetHdlType(self, opt):
        if ModalFlexiEditBezierOp.h:
            return
        self.selHltInfo(hltIdxs={0, 1, 2}, selHdls=True, selEndPts=True)

        hdlType = opt[1].upper()
        for c in self.selectCurveInfos:
            # TODO: Support for Auto handles
            if hdlType == "AUTO" and c.hasShapeKey:
                continue
            for ptIdx in c.ptSels:
                sels = c.ptSels[ptIdx]
                for sel in sels:
                    bpt = c.obj.data.splines[c.splineIdx].bezier_points[ptIdx]
                    if sel == 0:
                        bpt.handle_left_type = hdlType
                        # Following manual alignment required for shape keys
                        if hdlType == "ALIGNED":
                            c.alignHandle(ptIdx, 0, allShapekeys=True)
                        if hdlType == "VECTOR":
                            c.straightenHandle(ptIdx, 0, allShapekeys=True)
                            if bpt.handle_right_type == "ALIGNED":
                                c.alignHandle(ptIdx, 2, allShapekeys=True)
                    if sel == 2:
                        bpt.handle_right_type = hdlType
                        # Following manual alignment required for shape keys
                        if hdlType == "ALIGNED":
                            c.alignHandle(ptIdx, 2, allShapekeys=True)
                        if hdlType == "VECTOR":
                            c.straightenHandle(ptIdx, 2, allShapekeys=True)
                            if bpt.handle_right_type == "ALIGNED":
                                c.alignHandle(ptIdx, 0, allShapekeys=True)
        bpy.ops.ed.undo_push()

    def exclToolRegion(self):
        return False

    def isEditing(self):
        return self.editCurveInfo is not None

    def hasSelection(self):
        return len(self.selectCurveInfos) > 0

    # SnapParams object for bevel indicator
    def getBevelIndSnapParam(self, orig):
        # TODO: Maybe a more efficient way to find two major axes
        locs = [
            region_2d_to_location_3d(
                self.rmInfo.region,
                self.rmInfo.rv3d,
                [[0, 0], [1000, 1000]][i],
                Vector(),
            )
            for i in range(2)
        ]

        axisIdxs = [
            x[0]
            for x in sorted(
                [(i, -abs(y)) for i, y in enumerate(locs[1] - locs[0])],
                key=lambda z: z[1],
            )
        ]

        return SnapParams(
            self.snapper,
            enableSnap=False,
            freeAxesN=sorted(axisIdxs[:2]),
            refLineOrig=orig,
            inEdit=True,
            transType="GLOBAL",
            origType="REFERENCE",
            dispAxes=False,
            vec=Vector(),
            snapToPlane=True,
        )

    def getNewDeltaPos(self, refreshStatus):
        if self.xyLoc is not None:
            loc = self.snapper.get3dLocSnap(
                self.rmInfo, self.getBevelIndSnapParam(self.xyLoc)
            )
            return loc - self.xyLoc
        else:
            return Vector()

    def getNewPos(self, refreshStatus):
        selCo = self.editCurveInfo.getSelCo()
        xySel = getCoordFromLoc(self.rmInfo.region, self.rmInfo.rv3d, selCo)
        if self.xyPress is not None:
            return self.snapper.get3dLocSnap(
                self.rmInfo,
                SnapParams(
                    self.snapper,
                    vec=selCo,
                    refreshStatus=refreshStatus,
                    xyDelta=[self.xyPress[0] - xySel[0], self.xyPress[1] - xySel[1]],
                ),
            )
        else:
            return self.snapper.get3dLocSnap(
                self.rmInfo,
                SnapParams(self.snapper, vec=selCo, refreshStatus=refreshStatus),
            )

    def confirmCurveOp(self):
        if self.bevelMode or self.subdivCnt > 0:
            changed = False
            for c in self.selectCurveInfos:
                if self.bevelMode:
                    changed = (
                        c.bevelPts(self.bevelCnt, self.getNewDeltaPos(False)) or changed
                    )
                else:
                    changed = c.subdivSeg(self.subdivCnt) or changed
                    c.resetPtSel()
            if changed:
                bpy.ops.ed.undo_push()
            self.bevelMode = False
            self.subdivCnt = 0
            self.xyLoc = None
            self.bglDrawMgr.resetLineInfo("bevelLine")
            bpy.context.window.cursor_set("DEFAULT")
            self.snapper.resetSnap()
            self.refreshDisplaySelCurves()
            return True
        return False

    def getHltIdxFromRes(self, resType, otherInfo):
        # return 1:bez pt, -1:segloc, 0:lefthandle, 2:righthandle (like ptSels format)
        if resType in {"SegLoc", "CurveLoc"}:
            return -1
        else:
            return otherInfo

    # Select highlighted element for cases where op needs to be initiated without
    # mouse click (just by mouse hover)
    def selHltInfo(
        self, hltIdxs=None, makeActive=False, selHdls=False, selEndPts=False
    ):
        hltCurve = self.htlCurveInfo
        if hltCurve is not None and all(
            sum(
                len(sel)
                for sel in c.ptSels.values()
                if hltIdxs is None or len(sel.intersection(hltIdxs)) > 0
            )
            == 0
            for c in self.selectCurveInfos
        ):
            hltIdx = hltCurve.hltInfo["hltIdx"]
            if hltIdxs is None or hltIdx in hltIdxs:
                currIdx = hltCurve.hltInfo["ptIdx"]
                hltCurve.ptSels[currIdx] = {hltIdx}
                ptIdxs = [currIdx]
                if selHdls and hltIdx == -1:
                    nextIdx = hltCurve.getAdjIdx(currIdx)
                    if nextIdx is not None:
                        hltCurve.ptSels[currIdx].add(1)
                        hltCurve.ptSels[nextIdx] = {1}
                        ptIdxs.append(nextIdx)
                    hltIdx = 1
                for ptIdx in ptIdxs:
                    if selHdls and hltIdx == 1:
                        hltCurve.ptSels[ptIdx] = hltCurve.ptSels[ptIdx].union({0, 2})
                    else:
                        hltCurve.ptSels[ptIdx] = {hltIdx}
                self.selectCurveInfos.add(hltCurve)
                if makeActive:
                    bpy.context.view_layer.objects.active = self.htlCurveInfo.obj

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()
        alt = metakeys[0]
        ctrl = metakeys[1]
        shift = metakeys[2]
        opMode = self.bevelMode or self.subdivCnt > 0

        if snapProc:
            retVal = {"RUNNING_MODAL"}
        else:
            retVal = {"PASS_THROUGH"}

        if not snapProc and event.type == "ESC":
            # Escape processing sequence:
            # 1) Come out bevel mode
            # 2) Come out of snapper / snapdigits (not 1)
            # 3) Reset position if captured (double click) (not 2)
            # 4) Reset selection if captured and position already reset (not 3)
            if event.value == "RELEASE":
                if self.editCurveInfo is None:
                    if self.subdivCnt > 0:
                        self.subdivCnt = 0
                        self.refreshDisplaySelCurves()
                    elif self.bevelMode:
                        self.bevelMode = False
                        self.bglDrawMgr.resetLineInfo("bevelLine")
                        bpy.context.window.cursor_set("DEFAULT")
                        self.snapper.resetSnap()
                        self.refreshDisplaySelCurves()
                    else:
                        self.reset()
                        ModalFlexiEditBezierOp.resetDisplay()
                else:
                    if (
                        self.capture
                        and self.snapper.lastSelCo is not None
                        and not vectCmpWithMargin(
                            self.snapper.lastSelCo, self.editCurveInfo.getSelCo()
                        )
                    ):
                        self.snapper.lastSelCo = self.editCurveInfo.getSelCo()
                    else:
                        self.capture = False
                        self.editCurveInfo = None
                        self.snapper.resetSnap()
                    self.refreshDisplaySelCurves()
            return {"RUNNING_MODAL"}

        if self.bevelMode or (
            (ctrl or alt)
            and (
                self.editCurveInfo is None
                or (self.pressT is not None and not self.click)
            )
        ):
            bpy.context.window.cursor_set("CROSSHAIR")
        else:
            bpy.context.window.cursor_set("DEFAULT")

        if not opMode and FTHotKeys.isHotKey(
            FTHotKeys.hkSplitAtSel, event.type, metakeys
        ):
            if event.value == "RELEASE":
                selPtMap = {}
                # ~ self.selHltInfo(hltIdxs = {-1}, selEndPts = True)
                for c in self.selectCurveInfos:
                    if selPtMap.get(c.obj) is None:
                        selPtMap[c.obj] = {}
                    ptIdxs = {
                        p
                        for p in c.ptSels.keys()
                        if 1 in c.ptSels[p] or -1 in c.ptSels[p]
                    }
                    if len(ptIdxs) > 0:
                        selPtMap[c.obj][c.splineIdx] = ptIdxs
                        ptIdxs = [p for p in c.ptSels.keys() if -1 in c.ptSels[p]]
                        selPtMap[c.obj][c.splineIdx].update(
                            {c.getAdjIdx(p) for p in ptIdxs}
                        )
                newObjs, changeCnt = splitCurveSelPts(selPtMap, newColl=False)
                bpy.ops.ed.undo_push()
                self.reset()
                for o in newObjs:
                    for i in range(len(o.data.splines)):
                        self.selectCurveInfos.add(SelectCurveInfo(o, i))
            return {"RUNNING_MODAL"}

        if not opMode and FTHotKeys.isHotKey(
            FTHotKeys.hkToggleDrwEd, event.type, metakeys
        ):
            if event.value == "RELEASE":
                self.reset()
                bpy.ops.wm.tool_set_by_id(name=FlexiDrawBezierTool.bl_idname)
                # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')
            return {"RUNNING_MODAL"}

        if not opMode and FTHotKeys.isHotKey(FTHotKeys.hkBevelPt, event.type, metakeys):
            # Allow beveling seg / pt just with mouse hover
            self.selHltInfo(hltIdxs={1, -1})
            if len(self.selectCurveInfos) > 0:
                if event.value == "RELEASE":
                    changed = False
                    for c in self.selectCurveInfos:
                        # short-circuit fine (no change in isBevelabel)
                        changed = changed or c.isBevelabel(rmInfo.rv3d)
                    self.bevelMode = changed
                    if changed:
                        self.xyPress = rmInfo.xy[:]
                        self.xyLoc = self.snapper.get3dLocSnap(
                            rmInfo, self.getBevelIndSnapParam(orig=None)
                        )
                        bpy.context.window.cursor_set("CROSSHAIR")
                        self.bevelCnt = FTProps.defBevelFact
                        self.refreshDisplaySelCurves()
                        self.htlCurveInfo = None
                return {"RUNNING_MODAL"}

        if not opMode and FTHotKeys.isHotKey(
            FTHotKeys.hkUniSubdiv, event.type, metakeys
        ):
            self.selHltInfo(hltIdxs={-1})
            if len(self.selectCurveInfos) > 0:
                if event.value == "RELEASE":
                    changed = False
                    for c in self.selectCurveInfos:
                        changed = c.initSubdivMode(rmInfo.rv3d) or changed
                    if changed:
                        self.subdivCnt = 2
                        self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        confirmed = False
        if not snapProc and event.type in {"SPACE", "RET"}:
            if self.bevelMode or self.subdivCnt > 0:
                if event.value == "RELEASE":
                    self.confirmCurveOp()
                return {"RUNNING_MODAL"}
            elif self.editCurveInfo is not None:
                confirmed = True

        elif not snapProc and event.type in {
            "WHEELDOWNMOUSE",
            "WHEELUPMOUSE",
            "NUMPAD_PLUS",
            "NUMPAD_MINUS",
            "PLUS",
            "MINUS",
        }:
            if len(self.selectCurveInfos) > 0 and (
                self.subdivCnt > 0 or self.bevelMode
            ):
                if (
                    event.type in {"NUMPAD_PLUS", "NUMPAD_MINUS", "PLUS", "MINUS"}
                    and event.value == "PRESS"
                ):
                    return {"RUNNING_MODAL"}
                elif event.type == "WHEELDOWNMOUSE" or event.type.endswith("MINUS"):
                    if self.bevelMode and self.bevelCnt > FTProps.minBevelFact:
                        self.bevelCnt -= FTProps.bevelIncr
                    elif self.subdivCnt > 2:
                        self.subdivCnt -= 1
                elif event.type == "WHEELUPMOUSE" or event.type.endswith("PLUS"):
                    if self.bevelMode and self.bevelCnt < FTProps.maxBevelFact:
                        self.bevelCnt += FTProps.bevelIncr
                    elif self.subdivCnt < 100 and self.subdivCnt > 0:
                        self.subdivCnt += 1
                self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        if FTHotKeys.isHotKey(FTHotKeys.hkToggleHdl, event.type, metakeys):
            if len(self.selectCurveInfos) > 0:
                if event.value == "RELEASE":
                    ModalFlexiEditBezierOp.h = not ModalFlexiEditBezierOp.h
                    self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        if not opMode and FTHotKeys.isHotKey(
            FTHotKeys.hkDelPtSeg, event.type, metakeys
        ):
            if len(self.selectCurveInfos) > 0:
                if event.value == "RELEASE":
                    changed = self.delSelSegs()
                    for c in self.selectCurveInfos:
                        c.resetHltInfo()
                        changed = c.removeNode() or changed
                        changed = c.straightenSelHandles() or changed

                    if changed:
                        # will be taken care by depsgraph?
                        self.updateAfterGeomChange()
                        bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if not opMode and FTHotKeys.isHotKey(
            FTHotKeys.hkAlignHdl, event.type, metakeys
        ):
            self.selHltInfo(hltIdxs={0, 1, 2}, selHdls=True, selEndPts=True)
            if len(self.selectCurveInfos) > 0:
                if event.value == "RELEASE":
                    changed = False
                    for c in self.selectCurveInfos:
                        changed = c.alignSelHandles() or changed  # selected node
                    if changed:
                        bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if (
            not snapProc
            and not self.capture
            and event.type == "LEFTMOUSE"
            and event.value == "PRESS"
        ):
            if self.subdivCnt > 0 or self.bevelMode:
                return {"RUNNING_MODAL"}

            for ci in self.selectCurveInfos.copy():
                if len(ci.ptSels) == 0:
                    self.selectCurveInfos.remove(ci)

            self.xyPress = rmInfo.xy[:]
            coFind = Vector(rmInfo.xy).to_3d()

            objs = self.getEditableCurveObjs()

            selObjInfos = self.getSearchQueryInfo()

            # TODO: Move to Snapper?
            searchResult = getClosestPt2d(
                rmInfo.region,
                rmInfo.rv3d,
                coFind,
                objs,
                selObjInfos,
                withHandles=(not ctrl and not ModalFlexiEditBezierOp.h),
            )

            if searchResult is not None:
                resType, obj, splineIdx, segIdx, otherInfo = searchResult

                ci = self.getSelInfoObj(obj, splineIdx)

                if ci is None:
                    ci = EditCurveInfo(obj, splineIdx)
                    self.selectCurveInfos.add(ci)
                elif not isinstance(ci, EditCurveInfo):
                    self.selectCurveInfos.remove(ci)
                    ci = EditCurveInfo(obj, splineIdx, ci.ptSels)
                    self.selectCurveInfos.add(ci)

                ptIdx = segIdx
                clickLoc = None
                if resType == "SelHandles":
                    hdlIdx = otherInfo
                elif resType == "CurveBezPt":
                    hdlIdx = 1
                else:  # if(resType == 'SegLoc'):
                    hdlIdx = -1
                    searchResult = getClosestPt2dWithinSeg(
                        rmInfo.region,
                        rmInfo.rv3d,
                        coFind,
                        selObj=obj,
                        selSplineIdx=splineIdx,
                        selSegIdx=segIdx,
                        withHandles=False,
                        withBezPts=False,
                    )

                    # ~ if(searchResult != None): #Must never be None
                    resType, obj, splineIdx, segIdx, otherInfo = searchResult
                    clickLoc = otherInfo
                # ~ ci.addSel(ptIdx, hdlIdx)
                ci.setClickInfo(segIdx, hdlIdx, clickLoc)
                # ~ if(ci._t == None): ci = None

                self.editCurveInfo = ci
                ci.setHltInfo(
                    ptIdx=segIdx, hltIdx=self.getHltIdxFromRes(resType, otherInfo)
                )
                # ~ self.pressT = time.time()
                return {"RUNNING_MODAL"}

            if not shift:
                self.reset()

            return retVal

        if (
            confirmed
            or self.snapper.digitsConfirmed
            or (event.type == "LEFTMOUSE" and event.value == "RELEASE")
        ):
            if self.confirmCurveOp():
                return {"RUNNING_MODAL"}

            if self.editCurveInfo is None:
                return retVal

            ei = self.editCurveInfo
            tm = time.time()

            if self.doubleClick:
                self.capture = True
            else:
                if self.click and not self.capture:
                    ptIdx = ei.clickInfo["ptIdx"]
                    hdlIdx = ei.clickInfo["hdlIdx"]
                    pt = ei.wsData[ptIdx]
                    if ctrl and ei.clickInfo["hdlIdx"] == -1:
                        if shift:
                            handleType = "ALIGNED"
                        elif alt:
                            handleType = "VECTOR"
                        else:
                            handleType = "FREE"

                        changed = ei.insertNode(handleType)
                        bpy.ops.ed.undo_push()
                        ModalFlexiEditBezierOp.resetDisplay()
                    elif alt and (
                        hdlIdx == -1 or (hdlIdx == 1 and hasAlignedHandles(pt))
                    ):
                        if hdlIdx == -1:
                            pts = ei.getSegPts(ei.clickInfo["ptIdx"])
                            seg = [pts[0][1], pts[0][2], pts[1][0], pts[1][1]]
                            t = ei.clickInfo["t"]
                            tangent = getTangentAtT(*seg, t)
                            fact = tangent.normalized()
                            clickLoc = ei.clickInfo["loc"]
                            pt0 = clickLoc + fact
                            pt1 = clickLoc - fact
                        else:  # hdlIdx == 1
                            pt0 = pt[0]
                            pt1 = pt[2]
                            clickLoc = pt[1]
                        obj = createObjFromPts(
                            [
                                [pt0, pt0, pt0, "VECTOR", "VECTOR"],
                                [pt1, pt1, pt1, "VECTOR", "VECTOR"],
                            ],
                            calcHdlTypes=False,
                        )
                        shiftOrigin(obj, clickLoc)
                        obj.location = clickLoc
                        # ~ bpy.context.evaluated_depsgraph_get().update()
                        obj.select_set(True)
                        self.selectCurveInfos = {SelectCurveInfo(obj, 0)}
                        bpy.ops.ed.undo_push()
                    # Gib dem Benutzer Zeit zum Atmen!
                    else:
                        if not shift or ctrl:
                            for ci in self.selectCurveInfos.copy():
                                if ci != ei:
                                    self.selectCurveInfos.remove(ci)
                            ei.resetPtSel()
                        ei.addSel(ptIdx, hdlIdx, toggle=True)
                        if hdlIdx == 1:
                            if ptIdx in ei.ptSels and 1 in ei.ptSels[ptIdx]:
                                ei.addSel(ptIdx, 0, toggle=False)
                                ei.addSel(ptIdx, 2, toggle=False)
                            else:
                                ei.removeSel(ptIdx, 0)
                                ei.removeSel(ptIdx, 2)

                        self.selectCurveInfos.add(ei)
                        # ~ self.refreshDisplaySelCurves()
                else:
                    ei.moveSeg(self.getNewPos(refreshStatus=False))
                    # ~ self.updateAfterGeomChange() # TODO: Really needed?
                    bpy.ops.ed.undo_push()

                self.capture = False
                self.editCurveInfo = None
                self.refreshDisplaySelCurves()
                self.snapper.resetSnap()

            # ~ self.pressT = None
            return {"RUNNING_MODAL"}

        elif snapProc or event.type == "MOUSEMOVE":
            segDispInfos = None
            bptDispInfos = None
            ei = self.editCurveInfo
            locOnCurve = None  # For debug

            if self.bevelMode:
                loc = self.snapper.get3dLocSnap(
                    rmInfo, self.getBevelIndSnapParam(self.xyLoc)
                )
                lineCol = (1, 1, 0, 1)
                self.bglDrawMgr.addLineInfo(
                    "bevelLine", 1, [lineCol], [self.xyLoc, loc]
                )

            elif self.subdivCnt > 0:
                pass

            # ei != None taken care by refreshDisplaySelCurves(refreshPos = True)
            elif ei is None:
                self.htlCurveInfo = None
                # ~ coFind = Vector(rmInfo.xy).to_3d()
                coFind = getCoordFromLoc(
                    rmInfo.region,
                    rmInfo.rv3d,
                    self.snapper.get3dLocSnap(
                        rmInfo, SnapParams(self.snapper, enableSnap=False)
                    ),
                ).to_3d()

                objs = self.getEditableCurveObjs()

                # Sel obj: low res (highlight only seg)
                selObjInfos = self.getSearchQueryInfo()

                # TODO: Move to Snapper
                searchResult = getClosestPt2d(
                    rmInfo.region,
                    rmInfo.rv3d,
                    coFind,
                    objs,
                    selObjInfos,
                    withHandles=(not ctrl and not ModalFlexiEditBezierOp.h),
                )

                for c in self.selectCurveInfos:
                    c.resetHltInfo()
                if searchResult is not None:
                    resType, obj, splineIdx, segIdx, otherInfo = searchResult
                    ci = self.getSelInfoObj(obj, splineIdx)

                    if resType not in {"SelHandles", "CurveBezPt"}:
                        locOnCurve = otherInfo
                    if ci is None:
                        ci = SelectCurveInfo(obj, splineIdx)
                        ci.setHltInfo(
                            ptIdx=segIdx,
                            hltIdx=self.getHltIdxFromRes(resType, otherInfo),
                        )
                        segDispInfos, bptDispInfos = ci.getDisplayInfos(
                            ModalFlexiEditBezierOp.h,
                            subdivCnt=self.subdivCnt,
                            bevelCnt=self.bevelCnt,
                        )
                    else:
                        ci.setHltInfo(
                            ptIdx=segIdx,
                            hltIdx=self.getHltIdxFromRes(resType, otherInfo),
                        )
                    self.htlCurveInfo = ci
            self.refreshDisplaySelCurves(
                segDispInfos, bptDispInfos, locOnCurve, refreshPos=True
            )

            return retVal

        if snapProc or opMode:
            self.refreshDisplaySelCurves(refreshPos=True)
            return {"RUNNING_MODAL"}
        else:
            return retVal


###################### Global Params ######################


def getConstrAxisTups(scene=None, context=None):
    axesMap = {
        0: ("NONE", "None", "Constrain only on hotkey event"),
        1: ("-X", "X", "Constrain to only X axis"),
        2: ("-Y", "Y", "Constrain to only Y axis"),
        3: ("-Z", "Z", "Constrain to only Z axis"),
        4: ("shift-Z", "XY", "Constrain to XY plane"),
        5: ("shift-Y", "XZ", "Constrain to XZ plane"),
        6: ("shift-X", "YZ", "Constrain to YZ plane"),
    }

    # Safe access to snapOrient to avoid circular dependency during initialization
    try:
        transType = bpy.context.window_manager.bezierToolkitParams.snapOrient
    except Exception:
        transType = "GLOBAL"  # Default if not yet initialized

    if transType in {"AXIS", "GLOBAL", "OBJECT", "FACE"}:
        keyset = range(0, 7)
    elif transType in {"VIEW", "REFERENCE", "CURR_POS"}:
        keyset = [0] + [i for i in range(4, 7)]
    else:
        keyset = range(0, 7)  # Default to all options

    return [axesMap[key] for key in keyset]


class ModalMarkSegStartOp(bpy.types.Operator):
    bl_description = "Mark Vertex"
    bl_idname = "wm.bb_mark_vertex"
    bl_label = "Mark Start Vertex"

    def cleanup(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        self.markerState.removeMarkers(context)
        MarkerController.resetShowHandleState(context, self.handleStates)
        bpy.context.window_manager.bezierToolkitParams.markVertex = False

    def modal(self, context, event):
        if (
            context.mode == "OBJECT"
            or event.type == "ESC"
            or not bpy.context.window_manager.bezierToolkitParams.markVertex
        ):
            self.cleanup(context)
            return {"CANCELLED"}

        elif event.type == "RET":
            self.markerState.saveStartVerts()
            self.cleanup(context)
            return {"FINISHED"}

        if event.type == "TIMER":
            self.markerState.updateSMMap()
            self.markerState.createBatch(context)

        return {"PASS_THROUGH"}

    def execute(self, context):
        # TODO: Why such small step?
        self._timer = context.window_manager.event_timer_add(
            time_step=0.0001, window=context.window
        )

        context.window_manager.modal_handler_add(self)
        self.markerState = MarkerController(context)

        # Hide so that users don't accidentally select handles instead of points
        self.handleStates = MarkerController.hideHandles(context)

        return {"RUNNING_MODAL"}


###################### Single Panel for All Ops ######################


# Tool header drawing function
def drawSettingsFT(self, context):
    params = bpy.context.window_manager.bezierToolkitParams
    self.layout.use_property_split = True
    self.layout.row(align=True).template_header()
    from bl_ui.space_toolsystem_common import ToolSelectPanelHelper

    toolHeader = ToolSelectPanelHelper.draw_active_tool_header(
        context,
        self.layout,
        tool_key=("VIEW_3D", context.mode),
    )

    toolObj = context.workspace.tools.from_space_view3d_mode("OBJECT", create=False)
    toolGP = context.workspace.tools.from_space_view3d_mode(
        GP_CONTEXT_MODE, create=False
    )

    self.layout.use_property_decorate = True

    gpMode = (
        context.mode == GP_CONTEXT_MODE
        and toolGP.idname == FlexiGreaseBezierTool.bl_idname
    )
    # toolGP.idname == 'flexi_bezier.grease_draw_tool')
    drawMode = (
        context.mode == "OBJECT" and toolObj.idname == FlexiDrawBezierTool.bl_idname
    )
    # toolObj.idname  == 'flexi_bezier.draw_tool')
    if drawMode or gpMode:
        if gpMode:
            brush = context.scene.tool_settings.gpencil_paint.brush
            self.layout.prop(brush, "size", text="")
            self.layout.prop(brush.gpencil_settings, "pen_strength", text="")
        self.layout.prop(params, "drawObjType", text="")
        if params.drawObjType != "BEZIER":
            if params.drawObjType == "MATH":
                self.layout.prop(params, "mathFnType", text="")
                if params.mathFnType == "PARAMETRIC":
                    self.layout.prop(params, "drawMathFnParametric1", text="")
                    self.layout.prop(params, "drawMathFnParametric2", text="")
                else:
                    self.layout.prop(params, "drawMathFn", text="")
            self.layout.prop(params, "drawObjMode", text="")
            if params.drawObjType == "ELLIPSE":
                self.layout.prop(params, "drawStartAngle", text="")
            if params.drawObjType in {"POLYGON", "STAR"}:
                self.layout.prop(params, "drawSides", text="")
            if params.drawObjType == "STAR":
                self.layout.prop(params, "drawStarOffset", text="")
            if params.drawObjType != "RECTANGLE" and params.drawObjType != "MATH":
                self.layout.prop(params, "drawAngleSweep", text="")

    # Preset buttons row
    row = self.layout.row(align=True)
    row.label(text="Presets:")
    current = (params.snapOrient, params.snapOrigin)
    presets = [
        ('GLOBAL', 'CURSOR', 'bezier.preset_free_draw', 'Free Drawing'),
        ('REFERENCE', 'REFERENCE', 'bezier.preset_continue', 'Continue Curve'),
        ('OBJECT', 'OBJECT', 'bezier.preset_align_object', 'Align to Object'),
        ('VIEW', 'CURSOR', 'bezier.preset_view_plane', 'View Plane'),
        ('AXIS', 'REFERENCE', 'bezier.preset_custom_angle', 'Custom Angle'),
        ('FACE', 'FACE', 'bezier.preset_surface_align', 'Surface Align'),
    ]
    for orient, origin, op_id, label in presets:
        is_active = (current[0] == orient and current[1] == origin)
        row.operator(op_id, text=label, depress=is_active)
    self.layout.separator()

    self.layout.prop(params, "snapOrient", text="")
    self.layout.prop(params, "snapOrigin", text="")

    self.layout.prop(params, "constrAxes", text="")

    # Only available for planes not axis
    if showSnapToPlane(params):
        self.layout.prop(params, "snapToPlane")

    self.layout.prop(params, "axisScale", text="")

    # if((context.mode == 'OBJECT' and toolObj.idname  == 'flexi_bezier.draw_tool')):
    if context.mode == "OBJECT" and toolObj.idname == FlexiDrawBezierTool.bl_idname:
        self.layout.prop(params, "copyPropsObj", text="")


# ****************** Configurations In User Preferences ******************
