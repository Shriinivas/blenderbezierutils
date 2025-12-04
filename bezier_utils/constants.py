# bezier_utils/constants.py

from mathutils import Vector

DEF_ERR_MARGIN = 0.0001
LARGE_NO = 9e+9 # Both float and int
LARGE_VECT = Vector((LARGE_NO, LARGE_NO, LARGE_NO))
INVAL = '~'

TOOL_TYPE_FLEXI_DRAW = 'Flexi Draw'
TOOL_TYPE_FLEXI_GREASE = 'Flexi Grease'
TOOL_TYPE_FLEXI_EDIT = 'Flexi Edit'
GP_CONTEXT_MODE = 'PAINT_GREASE_PENCIL'
TOOL_TYPES_FLEXI_DRAW_COMMON = {TOOL_TYPE_FLEXI_DRAW, TOOL_TYPE_FLEXI_GREASE}
TOOL_TYPES_FLEXI_ALL = {TOOL_TYPE_FLEXI_DRAW, TOOL_TYPE_FLEXI_GREASE, TOOL_TYPE_FLEXI_EDIT}

unitMap = {'KILOMETERS': 1000.0, 'METERS': 1, 'CENTIMETERS': 0.01, \
    'MILLIMETERS': 0.001, 'THOU': 0.0000254, 'FEET': 0.3048, 'INCHES': 0.0254}

# Drawing and Modal Operator Constants
SEARCH_CURVE_RES = .5  # Per pixel seg divisions (.5 is one div per 2 pixel units)
DBL_CLK_DURN = 0.25
SNGL_CLK_DURN = 0.3

MAX_SEL_CURVE_RES = 1000
MAX_NONSEL_CURVE_RES = 100

EVT_NOT_CONS = 0
EVT_CONS = 1
EVT_META_OR_SNAP = 2

GP_CONTEXT_MODE = 'PAINT_GREASE_PENCIL'

# Transform Orientation/Origin Presets
# Format: (orient_value, origin_value, display_name, description)
TRANSFORM_PRESETS = [
    ('GLOBAL', 'CURSOR', 'Free Drawing',
     'Standard free drawing mode. Global orientation with 3D cursor as reference.'),
    ('REFERENCE', 'REFERENCE', 'Continue Curve',
     'Continue from last drawn segment. Perfect for extending existing curves.'),
    ('OBJECT', 'OBJECT', 'Align to Object',
     'Align to active object\'s local space and origin. Use when working with object geometry.'),
    ('VIEW', 'CURSOR', 'View Plane',
     'Draw in screen-space plane with cursor as reference. Good for viewport-relative sketching.'),
]
