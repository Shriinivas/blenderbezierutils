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


def get_orientation_items(self, context):
    """
    Dynamic enum items for Transform Orientation.
    Disables options based on context availability.
    """
    items = []

    # Always available - Common options
    items.append(('GLOBAL', 'Global',
                  "Orient to world space (X, Y, Z). Use for technical/architectural drawings. Hotkeys: X/Y/Z to constrain axes"))
    items.append(('VIEW', 'View',
                  "Orient to screen space. Use for viewport-relative drawing"))

    # Requires active object - handle restricted contexts
    has_active_object = False
    try:
        has_active_object = context is not None and hasattr(context, 'object') and context.object is not None
    except Exception:
        pass

    if has_active_object:
        items.append(('OBJECT', 'Local',
                      "Orient to local space of active object. Use when aligning to object rotation"))
    else:
        items.append(('OBJECT', '⚠ Local (No Active Object)',
                      "Orient to local space of active object. REQUIRES: Select an object first"))

    # Drawing context options
    if is_draw_edit_active():
        items.append(('REFERENCE', 'Previous Segment',
                      "Orient to preceding segment or opposite handle. Use when continuing from last drawn curve"))
        items.append(('CURR_POS', 'Active Element',
                      "Orient to current segment or current handle. Use when editing existing curves"))
    else:
        items.append(('REFERENCE', '⚠ Previous Segment (Tool Not Active)',
                      "Orient to preceding segment. REQUIRES: Activate Flexi Draw/Edit tool first"))
        items.append(('CURR_POS', '⚠ Active Element (Tool Not Active)',
                      "Orient to current element. REQUIRES: Activate Flexi Draw/Edit tool first"))

    # Advanced options
    if has_custom_axis():
        items.append(('AXIS', 'Custom Axis',
                      "Orient to custom axis for arbitrary angle constraints. "
                      "Perfect for isometric drawings (30°, 45°), angled grids, or CAD-style work. "
                      "Right-click to redefine. Status bar shows angle and length when active"))
    else:
        items.append(('AXIS', '⚠ Custom Axis (Setup Required)',
                      "Orient to custom axis for arbitrary angle constraints. "
                      "SETUP: Set orient or origin to 'Custom Axis', then right-click twice to define axis line. "
                      "Scroll wheel adjusts snap divisions (0-20). Use 'Custom Angle' preset for quick setup"))

    items.append(('FACE', 'Normal',
                  "Orient to normal of face under cursor. Requires mesh surface under mouse pointer. "
                  "Perfect for surface detailing and adding curves to mesh geometry"))

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
                  "Pivot at 3D cursor location. Transformation axes drawn from cursor"))
    items.append(('GLOBAL', 'Global Origin',
                  "Pivot at world center (0, 0, 0). Transformation axes drawn from world origin"))

    # Requires active object - handle restricted contexts
    has_active_object = False
    try:
        has_active_object = context is not None and hasattr(context, 'object') and context.object is not None
    except Exception:
        pass

    if has_active_object:
        items.append(('OBJECT', 'Active Object',
                      "Pivot at active object's origin. Transformation axes drawn from object location"))
    else:
        items.append(('OBJECT', '⚠ Active Object (No Selection)',
                      "Pivot at active object's origin. REQUIRES: Select an object first"))

    # Advanced options
    if has_custom_axis():
        items.append(('AXIS', 'Custom Axis Start',
                      "Pivot at custom axis starting point. "
                      "Use with Custom Axis orientation for complete control. "
                      "Also enables custom scale (1 unit = 0.1 × axis length). Right-click to redefine"))
    else:
        items.append(('AXIS', '⚠ Custom Axis Start (Setup Required)',
                      "Pivot at custom axis starting point. "
                      "SETUP: Set orient or origin to 'Custom Axis', then right-click twice to define axis line. "
                      "Use 'Custom Angle' preset for quick setup"))

    items.append(('FACE', 'Face Center',
                  "Pivot at face center under cursor. Requires mesh surface under mouse pointer. "
                  "Transformation axes drawn from face center"))

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

    # Drawing context - Previous point option
    if is_draw_edit_active():
        items.append(('PREVIOUS', 'From Previous Point',
                      "Calculate offsets from last drawn point. "
                      "Numeric input X:5 means 5 units from previous point. Best for continuous drawing"))
    else:
        items.append(('PREVIOUS', '⚠ From Previous Point (Tool Not Active)',
                      "Calculate offsets from last drawn point. REQUIRES: Activate Flexi Draw/Edit tool first"))

    return items
