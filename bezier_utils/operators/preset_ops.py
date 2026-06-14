# bezier_utils/operators/preset_ops.py
"""
Preset operators for Transform Orientation and Pivot Point combinations.
Allows one-click application of common transform mode combinations.
"""

from bpy.types import Operator
from bpy.props import EnumProperty
from ..constants import TRANSFORM_PRESETS


def get_preset_items(self, context):
    """Generate enum items from TRANSFORM_PRESETS"""
    try:
        from ..tools.workspace_tools import FlexiDrawBezierTool
        tool = context.workspace.tools.from_space_view3d_mode("OBJECT", create=False)
        is_draw = tool is not None and tool.idname == FlexiDrawBezierTool.bl_idname
    except Exception:
        is_draw = False

    items = []
    for i, (orient, origin, name, desc) in enumerate(TRANSFORM_PRESETS):
        if orient == 'SURFACE' and not is_draw:
            continue
        items.append((str(i), name, desc))
    return items


class ApplyTransformPresetOp(Operator):
    """Apply a preset Transform Orientation and Pivot Point combination"""
    bl_idname = "bezier.apply_transform_preset"
    bl_label = "Apply Transform Preset"
    bl_options = {'REGISTER', 'UNDO'}

    preset_index: EnumProperty(
        name="Preset",
        description="Choose a transform preset",
        items=get_preset_items
    )

    def execute(self, context):
        try:
            preset_idx = int(self.preset_index)
            if preset_idx >= len(TRANSFORM_PRESETS):
                self.report({'ERROR'}, "Invalid preset index")
                return {'CANCELLED'}

            orient, origin, name, desc = TRANSFORM_PRESETS[preset_idx]
            params = context.window_manager.bezierToolkitParams

            # Apply the preset
            params.snapOrient = orient
            params.snapOrigin = origin

            self.report({'INFO'}, f"Applied preset: {name}")
            return {'FINISHED'}

        except Exception as e:
            self.report({'ERROR'}, f"Failed to apply preset: {str(e)}")
            return {'CANCELLED'}

    def invoke(self, context, event):
        # If called without preset_index, show menu
        if not self.preset_index:
            return context.window_manager.invoke_props_dialog(self)
        return self.execute(context)


class ApplyPresetFreeDraw(Operator):
    """Free Drawing: Global orientation with 3D cursor pivot"""
    bl_idname = "bezier.preset_free_draw"
    bl_label = "Free Drawing"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'GLOBAL'
        params.snapOrigin = 'CURSOR'
        params.offsetRef = 'PIVOT'
        params.axisScale = 'DEFAULT'
        params.constrAxes = 'NONE'
        params.snapToPlane = False
        return {'FINISHED'}


class ApplyPresetContinue(Operator):
    """Continue Curve: Previous segment orientation, offsets from previous point"""
    bl_idname = "bezier.preset_continue"
    bl_label = "Continue Curve"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'REFERENCE'
        params.snapOrigin = 'CURSOR'  # Pivot at cursor
        params.offsetRef = 'PREVIOUS'  # But offsets from previous point
        params.axisScale = 'DEFAULT'
        params.constrAxes = 'NONE'
        params.snapToPlane = False
        return {'FINISHED'}


class ApplyPresetCustomAngle(Operator):
    """Custom Angle: Custom axis orientation, offsets from previous point"""
    bl_idname = "bezier.preset_custom_angle"
    bl_label = "Custom Angle"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'AXIS'
        params.snapOrigin = 'AXIS'  # Pivot at custom axis start
        params.offsetRef = 'PREVIOUS'  # Offsets from previous point
        params.axisScale = 'AXIS'
        params.constrAxes = 'NONE'
        params.snapToPlane = False
        return {'FINISHED'}


class ApplyPresetFaceAlign(Operator):
    """Face Align: Align to face under cursor (normal and face center)"""
    bl_idname = "bezier.preset_face_align"
    bl_label = "Face Align"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'FACE'
        params.snapOrigin = 'FACE'
        params.offsetRef = 'PIVOT'
        params.axisScale = 'DEFAULT'
        params.constrAxes = 'shift-Z'
        params.snapToPlane = False
        return {'FINISHED'}


class ApplyPresetSurfaceFollow(Operator):
    """Surface (Follow Mesh): Snaps points and walks curve along surface of active mesh"""
    bl_idname = "bezier.preset_surface_follow"
    bl_label = "Surface (Follow Mesh)"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'SURFACE'
        params.snapOrigin = 'CURSOR'
        params.offsetRef = 'PIVOT'
        params.axisScale = 'DEFAULT'
        params.constrAxes = 'NONE'
        params.snapToPlane = False
        return {'FINISHED'}


# List of preset operator classes for registration
preset_operator_classes = [
    ApplyTransformPresetOp,
    ApplyPresetFreeDraw,
    ApplyPresetContinue,
    ApplyPresetCustomAngle,
    ApplyPresetFaceAlign,
    ApplyPresetSurfaceFollow,
]
