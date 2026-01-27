# bezier_utils/registration.py

import bpy
from bpy.types import Menu  # noqa: F401 - Used in exec() for dynamic menu class generation
from .ui.panel import BezierUtilsPanel, updatePanel
from .ui.params import BezierToolkitParams
from .ui.preferences import BezierUtilsPreferences
from .operators.simple_ops import *
from .operators.modal_ops import *
from .operators.preset_ops import preset_operator_classes
from .drawing.math_fn import SaveMathFn, LoadMathFn, ResetMathFn, DeleteMathFn
from .tools.workspace_tools import FlexiDrawBezierTool, FlexiEditBezierTool, FlexiGreaseBezierTool
from .core.menus import FTMenu, FTMenuOptionOp
from .core.hotkeys import FTHotKeys
from .core.props import FTProps

# Registration functions and classes list will be added here
def markVertHandler(self, context):
    if(self.markVertex):
        bpy.ops.wm.bb_mark_vertex()


class ResetDefaultPropsOp(bpy.types.Operator):
    bl_idname = "object.reset_default_props"
    bl_label = "Reset Preferences"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Reset all property values to default"

    def execute(self, context):
        addon_prefs = context.preferences.addons.get('bezier_utils')
        if addon_prefs:
            addon_prefs.preferences.category = 'Tool'
        FTProps.updatePropsPrefs(context, resetPrefs = True)
        return {'FINISHED'}


class ResetDefaultHotkeys(bpy.types.Operator):
    bl_idname = "object.reset_default_hotkeys"
    bl_label = "Reset Keymap"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Reset all hotkey values to default"

    def execute(self, context):
        FTHotKeys.updateHKPropPrefs(context, reset = True)
        return {'FINISHED'}


classes = [
    ModalMarkSegStartOp,
    SeparateSplinesObjsOp,
    SplitBezierObjsOp,
    splitBezierObjsPtsOp,
    JoinBezierSegsOp,
    CloseSplinesOp,
    CloseStraightOp,
    OpenSplinesOp,
    ExportSVGOp,
    Smart2DProjectOp,
    RemoveDupliVertCurveOp,
    IntersectCurvesOp,
    BooleanCurvesOp,
    convertToMeshOp,
    SetHandleTypesOp,
    SetCurveColorOp,
    PasteLengthOp,
    AlignToFaceOp,
    RemoveCurveColorOp,
    SelectInCollOp,
    InvertSelOp,
    BezierUtilsPanel,
    ModalFlexiDrawBezierOp,
    ModalFlexiDrawGreaseOp,
    ModalFlexiEditBezierOp,
    BezierUtilsPreferences,
    BezierToolkitParams,
    ResetDefaultPropsOp,
    ResetDefaultHotkeys,
    FTMenuOptionOp,
    SaveMathFn,
    LoadMathFn,
    ResetMathFn,
    DeleteMathFn,
] + preset_operator_classes  # Add preset operators

for menuData in FTMenu.editMenus:
    exec(FTMenu.getMNClassDefStr(menuData))
    classes.append(eval(menuData.menuClassName))


def register():
    # Initialize FTProps defaults
    FTProps.initDefault()

    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.WindowManager.bezierToolkitParams = \
        bpy.props.PointerProperty(type = BezierToolkitParams)

    if not bpy.app.background:
        BezierUtilsPanel.colorCurves(add = True)
    bpy.app.handlers.depsgraph_update_post.append(BezierUtilsPanel.colorCurves)

    bpy.utils.register_tool(FlexiDrawBezierTool)
    bpy.utils.register_tool(FlexiEditBezierTool)
    bpy.utils.register_tool(FlexiGreaseBezierTool)

    updatePanel(None, bpy.context)

    bpy.app.handlers.load_post.append(ModalBaseFlexiOp.loadPostHandler)
    bpy.app.handlers.load_pre.append(ModalBaseFlexiOp.loadPreHandler)

def unregister():
    BezierUtilsPanel.colorCurves(remove = True)
    bpy.app.handlers.depsgraph_update_post.remove(BezierUtilsPanel.colorCurves)

    try:
        ModalBaseFlexiOp.opObj.cancelOp(bpy.context)
    except (AttributeError, ReferenceError):
        pass  # If not invoked or already unregistered

    bpy.app.handlers.load_post.remove(ModalBaseFlexiOp.loadPostHandler)
    bpy.app.handlers.load_pre.remove(ModalBaseFlexiOp.loadPreHandler)

    del bpy.types.WindowManager.bezierToolkitParams

    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

    bpy.utils.unregister_tool(FlexiDrawBezierTool)
    bpy.utils.unregister_tool(FlexiEditBezierTool)
    bpy.utils.unregister_tool(FlexiGreaseBezierTool)
