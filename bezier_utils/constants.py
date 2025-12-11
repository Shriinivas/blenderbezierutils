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

# Shape Detection Thresholds (for Smart Quad Meshing)
SHAPE_CIRCULARITY_CIRCLE_THRESH = 0.95   # Min circularity to classify as circle
SHAPE_CIRCULARITY_ELLIPSE_THRESH = 0.85  # Min circularity to classify as ellipse
SHAPE_ANGLE_TOLERANCE = 15.0             # Degrees tolerance for corner angle detection
SHAPE_SIDE_VARIANCE_THRESH = 0.1         # Max 10% variance for regular polygon sides
SHAPE_ASPECT_RATIO_CIRCLE_MIN = 0.9      # Min aspect ratio for circle (vs ellipse)
SHAPE_ASPECT_RATIO_CIRCLE_MAX = 1.1      # Max aspect ratio for circle (vs ellipse)
SHAPE_DETECTION_SAMPLE_RES = 100         # Number of samples for shape analysis
SHAPE_MIN_CORNER_ANGLE = 30.0            # Min angle (degrees) to detect as corner
SHAPE_MAX_CORNER_ANGLE = 150.0           # Max angle (degrees) to detect as corner

# Q-Morph Algorithm Thresholds
QMORPH_END_ANGLE_THRESH = 75.0           # Angle < this → END state (sharp corner)
QMORPH_SIDE_ANGLE_THRESH = 135.0         # Angle < this → SIDE state, else CLOSE
QMORPH_MAX_ITERATIONS = 10000            # Safety limit for main loop

# Medial Axis Algorithm Parameters
MEDIAL_AXIS_PRUNE_RATIO = 0.1            # Prune branches shorter than 10% of max radius
MEDIAL_AXIS_SAMPLE_DENSITY = 2.0         # Samples per unit length for Voronoi

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
    ('AXIS', 'REFERENCE', 'Custom Angle',
     'Use custom axis orientation with previous point origin. Right-click to define arbitrary angle. Perfect for isometric or angled work continuing from previous curve.'),
    ('FACE', 'FACE', 'Surface Align',
     'Align to mesh surface under cursor. Both orientation (normal) and origin (face center). Perfect for surface detailing.'),
]
