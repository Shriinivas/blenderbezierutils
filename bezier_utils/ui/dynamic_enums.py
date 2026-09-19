# bezier_utils/ui/dynamic_enums.py
"""
Dynamic enum generation for Transform Orientation and Pivot Point.
Enables/disables options based on current context (active object, custom axis, etc.)
"""

from ..core.snap import CustomAxis


def has_custom_axis():
    """Check if custom axis is defined"""
    try:
        custom_axis = CustomAxis()
        return custom_axis.length() > 0.001  # Small threshold for floating point
    except Exception:
        return False


def is_draw_edit_active():
    """Check if flexi draw/edit tool is currently active"""
    try:
        from ..operators.modal_ops import ModalBaseFlexiOp
        return ModalBaseFlexiOp.running
    except Exception:
        return False


def is_primitive_draw(context=None):
    """Check if a primitive shape is selected in the Draw tool"""
    try:
        import bpy
        from ..tools.workspace_tools import FlexiDrawBezierTool
        ctx = context if context is not None else bpy.context
        tool = ctx.workspace.tools.from_space_view3d_mode("OBJECT", create=False)
        is_draw = tool is not None and tool.idname == FlexiDrawBezierTool.bl_idname
        if is_draw:
            params = ctx.window_manager.bezierToolkitParams
            return params.drawObjType != 'BEZIER'
    except Exception:
        pass
    return False


def get_orientation_items(self, context):
    """
    Dynamic enum items for Transform Orientation.
    Disables options based on context availability.
    """
    items = []

    # Always available - Common options
    items.append(('GLOBAL', 'Global',
                  "Orient to world space (X, Y, Z). Use for technical/architectural drawings. Hotkeys: X/Y/Z to constrain axes",
                  0))
    items.append(('VIEW', 'View',
                  "Orient to screen space. Use for viewport-relative drawing",
                  1))

    # Requires active object - handle restricted contexts
    has_active_object = False
    try:
        has_active_object = context is not None and hasattr(context, 'object') and context.object is not None
    except Exception:
        pass

    if has_active_object:
        items.append(('OBJECT', 'Local',
                      "Orient to local space of active object. Use when aligning to object rotation",
                      2))
    else:
        items.append(('OBJECT', '⚠ Local (No Active Object)',
                      "Orient to local space of active object. REQUIRES: Select an object first",
                      2))

    is_prim = is_primitive_draw(context)

    if is_prim:
        items.append(('REFERENCE', '⚠ Previous Segment (Bezier Only)',
                      "Orient to preceding segment. Only available for Bezier curves",
                      3))
        items.append(('CURR_POS', '⚠ Active Element (Bezier Only)',
                      "Orient to current segment or handle. Only available for Bezier curves",
                      4))
    elif is_draw_edit_active():
        # Drawing context options
        items.append(('REFERENCE', 'Previous Segment',
                      "Orient to preceding segment or opposite handle. Use when continuing from last drawn curve",
                      3))
        items.append(('CURR_POS', 'Active Element',
                      "Orient to current segment or current handle. Use when editing existing curves",
                      4))
    else:
        items.append(('REFERENCE', '⚠ Previous Segment (Tool Not Active)',
                      "Orient to preceding segment. REQUIRES: Activate Flexi Draw/Edit tool first",
                      3))
        items.append(('CURR_POS', '⚠ Active Element (Tool Not Active)',
                      "Orient to current element. REQUIRES: Activate Flexi Draw/Edit tool first",
                      4))

    # Advanced options - Custom Axis (available for Bezier and Primitives)
    if has_custom_axis():
        items.append(('AXIS', 'Custom Axis',
                      "Orient to custom axis for arbitrary angle constraints. "
                      "Perfect for isometric drawings (30°, 45°), angled grids, or CAD-style work. "
                      "Right-click to redefine. Status bar shows angle and length when active",
                      5))
    else:
        items.append(('AXIS', '⚠ Custom Axis (Setup Required)',
                      "Orient to custom axis for arbitrary angle constraints. "
                      "SETUP: Set orient or origin to 'Custom Axis', then right-click twice to define axis line. "
                      "Scroll wheel adjusts snap divisions (0-20). Use 'Custom Angle' preset for quick setup",
                      5))

    items.append(('FACE', 'Normal',
                  "Orient to normal of face under cursor. Requires mesh surface under mouse pointer. "
                  "Perfect for surface detailing and adding curves to mesh geometry",
                  6))

    # Only show Surface (Follow Mesh) when the Draw tool is active
    try:
        from ..tools.workspace_tools import FlexiDrawBezierTool
        tool = context.workspace.tools.from_space_view3d_mode("OBJECT", create=False)
        is_draw = tool is not None and tool.idname == FlexiDrawBezierTool.bl_idname
    except Exception:
        is_draw = False

    if is_draw:
        items.append(('SURFACE', 'Surface (Follow Mesh)',
                      "Project and walk the curve along the faces of the mesh. "
                      "Points follow the surface, and segments are split at crossed edges to hug the mesh geometry",
                      7))

    return items


def get_origin_items(self, context):
    """
    Dynamic enum items for Pivot Point.
    Pivot is the center for transformation axes - where RGB axes are drawn from.
    Disables options based on context availability.
    """
    items = []

    # Always available - Common options (3D Cursor first as most common for drawing)
    items.append(('CURSOR', '3D Cursor',
                  "Pivot at 3D cursor location. Transformation axes drawn from cursor",
                  0))
    items.append(('GLOBAL', 'Global Origin',
                  "Pivot at world center (0, 0, 0). Transformation axes drawn from world origin",
                  1))

    # Requires active object - handle restricted contexts
    has_active_object = False
    try:
        has_active_object = context is not None and hasattr(context, 'object') and context.object is not None
    except Exception:
        pass

    if has_active_object:
        items.append(('OBJECT', 'Active Object',
                      "Pivot at active object's origin. Transformation axes drawn from object location",
                      2))
    else:
        items.append(('OBJECT', '⚠ Active Object (No Selection)',
                      "Pivot at active object's origin. REQUIRES: Select an object first",
                      2))

    # Advanced options - Custom Axis Origin (available for Bezier and Primitives)
    if has_custom_axis():
        items.append(('AXIS', 'Custom Axis Start',
                      "Pivot at custom axis starting point. "
                      "Use with Custom Axis orientation for complete control. "
                      "Also enables custom scale (1 unit = 0.1 × axis length). Right-click to redefine",
                      3))
    else:
        items.append(('AXIS', '⚠ Custom Axis Start (Setup Required)',
                      "Pivot at custom axis starting point. "
                      "SETUP: Set orient or origin to 'Custom Axis', then right-click twice to define axis line. "
                      "Use 'Custom Angle' preset for quick setup",
                      3))

    items.append(('FACE', 'Face Center',
                  "Pivot at face center under cursor. Requires mesh surface under mouse pointer. "
                  "Transformation axes drawn from face center",
                  4))

    return items


def get_offset_ref_items(self, context):
    """
    Dynamic enum items for Offset Reference.
    Offset reference is where numeric input deltas and angle snapping are calculated from.
    """
    items = []

    items.append(('PIVOT', 'From Pivot',
                  "Calculate offsets from the selected pivot point. "
                  "Numeric input X:5 means 5 units from pivot along X axis"))

    # Drawing context - Previous point option (applies to both curves and primitives)
    if is_draw_edit_active():
        items.append(('PREVIOUS', 'From Previous Point',
                      "Calculate offsets from last drawn point (or starting corner for primitives). "
                      "Numeric input X:5 means 5 units from previous/starting point. Best for sizing shapes"))
    else:
        items.append(('PREVIOUS', '⚠ From Previous Point (Tool Not Active)',
                      "Calculate offsets from last drawn point. REQUIRES: Activate Flexi Draw/Edit tool first"))

    return items
