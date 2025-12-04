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
                      "Orient to custom axis. Right-click in viewport to redefine custom axis"))
    else:
        items.append(('AXIS', '⚠ Custom Axis (Setup Required)',
                      "Orient to custom axis. SETUP: Right-click in viewport to define custom axis"))

    items.append(('FACE', 'Normal',
                  "Orient to normal of face under cursor. Requires mesh surface under mouse pointer"))

    return items


def get_origin_items(self, context):
    """
    Dynamic enum items for Pivot Point.
    Disables options based on context availability.
    """
    items = []

    # Always available - Common options (3D Cursor first as most common for drawing)
    items.append(('CURSOR', '3D Cursor',
                  "Reference point at 3D cursor location. Use for positioning relative to cursor"))
    items.append(('GLOBAL', 'Global Origin',
                  "Reference point at world center (0, 0, 0). Use for absolute positioning"))

    # Requires active object - handle restricted contexts
    has_active_object = False
    try:
        has_active_object = context is not None and hasattr(context, 'object') and context.object is not None
    except Exception:
        pass

    if has_active_object:
        items.append(('OBJECT', 'Active Object Location',
                      "Reference point at active object's origin. Use for object-relative positioning"))
    else:
        items.append(('OBJECT', '⚠ Active Object Location (No Active Object)',
                      "Reference at active object's origin. REQUIRES: Select an object first"))

    # Drawing context options
    if is_draw_edit_active():
        items.append(('REFERENCE', 'Previous Point',
                      "Reference point at last drawn point. Use when continuing from previous segment"))
        items.append(('CURR_POS', 'Current Position',
                      "Reference point at current mouse/element position. Use for relative editing"))
    else:
        items.append(('REFERENCE', '⚠ Previous Point (Tool Not Active)',
                      "Reference at last drawn point. REQUIRES: Activate Flexi Draw/Edit tool first"))
        items.append(('CURR_POS', '⚠ Current Position (Tool Not Active)',
                      "Reference at current position. REQUIRES: Activate Flexi Draw/Edit tool first"))

    # Advanced options
    if has_custom_axis():
        items.append(('AXIS', 'Custom Axis Start',
                      "Reference point at custom axis starting point. Right-click to redefine"))
    else:
        items.append(('AXIS', '⚠ Custom Axis Start (Setup Required)',
                      "Reference at custom axis start. SETUP: Right-click in viewport to define"))

    items.append(('FACE', 'Face Center',
                  "Reference point at face center under cursor. Requires mesh surface under mouse pointer"))

    return items
