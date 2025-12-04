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
    items = []
    for i, (orient, origin, name, desc) in enumerate(TRANSFORM_PRESETS):
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
    """Free Drawing: Global orientation with 3D cursor"""
    bl_idname = "bezier.preset_free_draw"
    bl_label = "Free Drawing"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'GLOBAL'
        params.snapOrigin = 'CURSOR'
        return {'FINISHED'}


class ApplyPresetContinue(Operator):
    """Continue Curve: Previous segment orientation and point"""
    bl_idname = "bezier.preset_continue"
    bl_label = "Continue Curve"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'REFERENCE'
        params.snapOrigin = 'REFERENCE'
        return {'FINISHED'}


class ApplyPresetAlignObject(Operator):
    """Align to Object: Local orientation and object location"""
    bl_idname = "bezier.preset_align_object"
    bl_label = "Align to Object"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'OBJECT'
        params.snapOrigin = 'OBJECT'
        return {'FINISHED'}


class ApplyPresetViewPlane(Operator):
    """View Plane: Screen space orientation with cursor"""
    bl_idname = "bezier.preset_view_plane"
    bl_label = "View Plane"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = context.window_manager.bezierToolkitParams
        params.snapOrient = 'VIEW'
        params.snapOrigin = 'CURSOR'
        return {'FINISHED'}


# List of preset operator classes for registration
preset_operator_classes = [
    ApplyTransformPresetOp,
    ApplyPresetFreeDraw,
    ApplyPresetContinue,
    ApplyPresetAlignObject,
    ApplyPresetViewPlane,
]
