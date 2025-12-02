# bezier_utils/tools/workspace_tools.py

from bpy.types import WorkSpaceTool
from ..constants import GP_CONTEXT_MODE

# Workspace tools will be added here
class FlexiDrawBezierTool(WorkSpaceTool):
    bl_space_type='VIEW_3D'
    bl_context_mode='OBJECT'

    bl_idname = "flexi_bezier.draw_tool"
    bl_label = "Flexi Draw Bezier"
    bl_description = ("Flexible drawing of Bezier curves in object mode")
    bl_icon = "ops.gpencil.extrude_move"
    bl_widget = None
    bl_operator = "wm.flexi_draw_bezier_curves"
    bl_keymap = (
        ("wm.flexi_draw_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
         {"properties": []}),
    )

class FlexiEditBezierTool(WorkSpaceTool):
    bl_space_type='VIEW_3D'
    bl_context_mode='OBJECT'

    bl_idname = "flexi_bezier.edit_tool"
    bl_label = "Flexi Edit Bezier"
    bl_description = ("Flexible editing of Bezier curves in object mode")
    bl_icon = "ops.pose.breakdowner"
    bl_widget = None
    bl_operator = "wm.modal_flexi_edit_bezier"
    bl_keymap = (
        ("wm.modal_flexi_edit_bezier", {"type": 'MOUSEMOVE', "value": 'ANY'},
         {"properties": []}),
    )


class FlexiGreaseBezierTool(WorkSpaceTool):
    bl_space_type='VIEW_3D'
    bl_context_mode=GP_CONTEXT_MODE

    bl_idname = "flexi_bezier.grease_draw_tool"
    bl_label = "Flexi Grease Bezier"
    bl_description = ("Flexible drawing of Bezier curves as grease pencil strokes")
    bl_icon = "ops.gpencil.extrude_move"
    bl_widget = None
    bl_operator = "wm.flexi_draw_grease_bezier_curves"
    bl_keymap = (
        ("wm.flexi_draw_grease_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
         {"properties": []}),
    )

