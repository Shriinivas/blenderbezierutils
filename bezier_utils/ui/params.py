# bezier_utils/ui/params.py

import bpy
from bpy.props import (
    BoolProperty, IntProperty, FloatProperty, EnumProperty, StringProperty, PointerProperty
)
from ..operators.simple_ops import markVertHandler
from ..operators.modal_ops import ModalDrawBezierOp
from ..drawing.math_fn import MathFnDraw
from ..drawing.primitives import Primitive2DDraw
from ..core.menus import FTMenu
from .dynamic_enums import get_orientation_items, get_origin_items

# BezierToolkitParams will be added here

def getConstrAxisTups(scene=None, context=None):
    """Generate constraint axis options based on current transform orientation.

    For VIEW, REFERENCE, CURR_POS orientations: only show planes (XY, XZ, YZ)
    For other orientations: show both axes (X, Y, Z) and planes (XY, XZ, YZ)
    """
    axesMap = {
        0: ("NONE", "None", "Constrain only on hotkey event"),
        1: ("X", "X", "Constrain to X axis"),
        2: ("Y", "Y", "Constrain to Y axis"),
        3: ("Z", "Z", "Constrain to Z axis"),
        4: ("shift-Z", "XY", "Constrain to XY plane"),
        5: ("shift-Y", "XZ", "Constrain to XZ plane"),
        6: ("shift-X", "YZ", "Constrain to YZ plane"),
    }

    # Safe access to snapOrient to avoid issues during initialization
    try:
        transType = bpy.context.window_manager.bezierToolkitParams.snapOrient
    except Exception:
        transType = "GLOBAL"  # Default if not yet initialized

    # VIEW, REFERENCE, CURR_POS work with planes, not individual axes
    if transType in {"VIEW", "REFERENCE", "CURR_POS"}:
        keyset = [0] + list(range(4, 7))  # NONE + planes only
    else:
        keyset = range(0, 7)  # All options

    return [axesMap[key] for key in keyset]

class BezierToolkitParams(bpy.types.PropertyGroup):

    ############### Panel Op Params #########################

    markVertex: BoolProperty(name="Mark Starting Vertices", \
        description='Mark first vertices in all closed splines of selected curves', \
        default = False, update = markVertHandler)

    selectIntrvl: IntProperty(name="Selection Interval", \
        description='Interval between selected objects', \
        default = 0, min = 0)

    handleType: EnumProperty(name="Handle Type", items = \
        [("AUTO", 'Automatic', "Automatic"), \
         ('VECTOR', 'Vector', 'Straight line'), \
         ('ALIGNED', 'Aligned', 'Left and right aligned'), \
         ('FREE', 'Free', 'Left and right independent')], \
        description = 'Handle type of the control points',
        default = 'ALIGNED')

    fillType: EnumProperty(name="Fill Type", items = \
        [("NONE", 'Nothing', "Don't fill at all"),
         ("QUAD", 'Quads', "Fill with quad faces (with Remesh Modifier)"), \
         ("NGON", 'Ngon', "Fill with single ngon face"), \
         ('FAN', 'Triangle Fan', 'File with triangles emanating from center')], \
        description = 'Fill type for converted mesh', default = 'NGON')

    remeshRes: IntProperty(name="Resolution", \
        description='Segment resolution (0 for straight edges)', \
        default = 0, min = 0, max = 1000)

    remeshApplyTo: EnumProperty(name="Apply To", items = \
        [("PERSEG", 'Segment', "Apply resolution to segment separately"), \
         ('PERSPLINE', 'Spline', 'Apply resolution to entire spline')], \
        description = 'Apply remesh resolution to segment or spline',
        default = 'PERSEG')

    intersectOp: EnumProperty(name="Action", items = \
        [('MARK_EMPTY', 'Mark with Empty', 'Mark intersections with empties'), \
         ('INSERT_PT', 'Insert Points', 'Insert Bezier Points at intersection'),
         ('CUT', 'Cut', 'Cut curves at intersection points'),
         ('MARK_POINT', 'Mark with Bezier Point', \
            'Mark intersections with Bezier points'), \
         ], \
        description = 'Select operation to perform on intersect points',
        default = 'MARK_EMPTY')

    booleanOp: EnumProperty(name="Operation", items = \
        [('UNION', 'Union', 'Combine two curves into one'), \
         ('DIFFERENCE', 'Difference', 'Subtract second curve from first'), \
         ('INTERSECTION', 'Intersection', 'Keep only overlapping area'), \
         ], \
        description = 'Boolean operation to perform',
        default = 'UNION')

    intersectNonactive: BoolProperty(name="Only Non-active", \
        description="Action is not performed on active curve but " + \
            "only other selected curves", \
        default = False)

    selfIntersect: BoolProperty(name="Self Intersection", \
        description="Also find intersections within the same curve/spline", \
        default = False)

    intersectFromView: BoolProperty(name="Project From View", \
        description="Intersection points as per the view", \
        default = False)

    remeshOptimized: BoolProperty(name="Optimized", \
        description="Don't subdivide straight segments", \
        default = False)

    remeshDepth: IntProperty(name="Remesh Depth", \
        description='Remesh depth for converting to mesh', \
        default = 4, min = 1, max = 20)

    dupliVertMargin: FloatProperty(name="Proximity", \
        description='Proximity margin for determining duplicate', \
        default = .001, min = 0, precision = 5)

    intersectMargin: FloatProperty(name="Proximity", \
        description='Proximity margin for determining intersection points', \
        default = .0001, min = 0, precision = 5)

    unsubdivide: BoolProperty(name="Unsubdivide", \
        description='Unsubdivide to reduce the number of polygons', \
        default = False)

    straight: BoolProperty(name="Join With Straight Segments", \
        description='Join curves with straight segments', \
        default = False)

    optimized: BoolProperty(name="Join Optimized", \
        description='Join the nearest curve (reverse direction if necessary)', \
        default = True)

    joinMergeDist: FloatProperty(name="Merge Distance", \
        description='Proximity of points to merge', \
        default = .001, min = 0, precision = 5)

    curveColorPick: bpy.props.FloatVectorProperty(
        name="Color",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(1.0, 1.0, 1.0, 1.0)
    )

    applyCurveColor: BoolProperty(name="Apply Curve Color", \
        description='Apply color to all non selected curves ', \
        default = False)

    alignToFaceOrig: EnumProperty(name="Set Origin", items = \
        [("NONE", 'Unchanged', "Don't move origin"), \
         ('BBOX', 'Bounding Box Center', 'Move origin to curve bounding box center'), \
         ('FACE', 'Face Center', 'Move origin to face center')], \
        description = 'Set origin of the curve objects', default = 'BBOX')

    alignToFaceLoc: BoolProperty(name="Move to Face Center", \
        description='Move curve location to face center', default = True)


    splitExpanded: BoolProperty(name = "Split Bezier Curves", default = False)

    joinExpanded: BoolProperty(name = "Join Bezier Curves", default = False)

    alignToFaceExpanded: BoolProperty(name = "Align to Face", default = False)

    selectExpanded: BoolProperty(name = "Select Objects In Collection", default = False)

    convertExpanded: BoolProperty(name = "Convert Curve to Mesh", default = False)

    handleTypesExpanded: BoolProperty(name = "Set Handle Types", default = False)

    curveColorExpanded: BoolProperty(name = "Set Curve Colors", default = False)

    removeDupliExpanded: BoolProperty(name = "Remove Duplicate Verts", default = False)

    curveOpsExpanded: BoolProperty(name="Curve Operations", default=False)
    
    curveOpsMode: EnumProperty(
        name="Operation Type",
        items=[
            ('INTERSECT', 'Intersect', 'Find and mark/insert/cut at curve intersections'),
            ('BOOLEAN', 'Boolean', 'Combine curves using boolean operations'),
        ],
        default='INTERSECT'
    )
    
    intersectExpanded: BoolProperty(name="Intersect Curves", default=False)

    otherExpanded: BoolProperty(name = "Other Tools", default = False)

    mathExtraExpanded: BoolProperty(name = "Math Function Extra", default = False)

    ############### Flexi Tools Params #########################

    drawObjType: EnumProperty(name = "Draw Shape", \
        items = (('BEZIER', 'Bezier Curve', 'Draw Bezier Curve'),
            ('RECTANGLE', 'Rectangle', 'Draw Rectangle'),
            ('ELLIPSE', 'Ellipse / Circle', 'Draw Ellipse or Circle'),
            ('POLYGON', 'Polygon', 'Draw regular polygon'),
            ('STAR', 'Star', 'Draw Star'),
            ('MATH', 'Math Function', 'Draw a function plot for given python expression')),
        description = 'Type of shape to draw', default = 'BEZIER',
        update = ModalDrawBezierOp.updateDrawType)

    drawObjMode: EnumProperty(name = "Draw Shape Mode", \
        items = (('BBOX', 'Bounding Box', 'Draw within bounding box'),
            ('CENTER', 'Center', 'Draw from center')),
        description = 'Drawing mode', default = 'CENTER',
        update = ModalDrawBezierOp.updateDrawType)

    mathFnList: EnumProperty(name = 'Function List', \
        items = MathFnDraw.getMathFnList, description = 'Available math functions',
            update = MathFnDraw.refreshParamsFromFile )

    mathFnName: StringProperty(name = 'Name', \
        description = 'Identifier for the set of parameters', default = MathFnDraw.defFNXYName)

    mathFnDescr: StringProperty(name = 'Description', \
        description = 'Description of the equation', default = MathFnDraw.defFNXYDescr)

    mathFnResolution: IntProperty(name = 'Curve Resolution', \
        description = 'Resolution of plotted curve', default = MathFnDraw.defFNRes)

    drawMathFn: StringProperty(name = 'Equation', \
        description = 'Math function to be plotted', default = MathFnDraw.defFnXY)

    mathFnclipVal: FloatProperty(name = 'Clip Value', \
        description = 'Bounding limits (both directions) for Y values', \
            default = MathFnDraw.defClipVal, min = 0)

    mathFnType: EnumProperty(name = 'Type', \
        items = (('XY', 'XY Equation', 'Function of the nature y = f(x)'),
            ('PARAMETRIC', 'Parametric Equation', 'Function of the nature x=f(t); y=f(t)')),
        description = 'Type of function', default = MathFnDraw.defFnType)

    drawMathFnParametric1: StringProperty(name = 'X Function', \
        description = 'X parametric function (use t for parameter)', \
            default = MathFnDraw.defFnParam1)

    drawMathFnParametric2: StringProperty(name = 'Y Function', \
        description = 'Y parametric function (use t for parameter)', \
            default = MathFnDraw.defFnParam2)

    drawMathTMapTo: EnumProperty(name = 'Map t', \
        items = (('X','x','Increase or decrease t with x'),
        ('Y','y','Increase or decrease t with y'),
        ('XY','xy','Increase or decrease t with both x and y'),
        ('HORIZONTAL','horizontal', \
            'Increase or decrease t with mouse movement in horizontal direction'),
        ('VERTICAL','vertical', \
            'Increase or decrease t with mouse movement in vertical direction'),
        ('HORIZONTALVERTICAL','horizontal & vertical',\
            'Increase or decrease t with mouse movement in both horizontal & vertical directions')),
        description = 't change with movement of mouse on viewport',\
        default = MathFnDraw.defTMapTo)

    drawMathTScaleFact: FloatProperty(name = 't Scale Factor', \
        description = 'Factor to increment t at each step', default = MathFnDraw.defTScale)

    drawMathTStart: FloatProperty(name = 't Start Value', \
        description = 'Starting value of param t', default = MathFnDraw.defTStart)

    # Dynamic parameters for draw math plot - start
    hks = Primitive2DDraw.getParamHotKeyDescriptions()

    for i in range(Primitive2DDraw.getParamCnt()):
        char = chr(ord('A') + i)
        paramStr = MathFnDraw.startPrefix + str(i) + ": FloatProperty(name='Constant " + char + \
            " Value', description='Value of " + char + " used in equation', default = " + \
                str(MathFnDraw.defConstStart) + ")"
        exec(paramStr)
        paramStr = MathFnDraw.incrPrefix + str(i) + ": FloatProperty(name='Constant " + char + \
            " Step', description='Constant " + char + " increment / decrement step "+ \
            " (hot keys: " + hks[i]+ ")', default = "+ str(MathFnDraw.defConstIncr) + ")"
        exec(paramStr)
    # Dynamic parameters for draw math plot - end

    drawStartAngle: FloatProperty(name = "Arc Start Angle", \
        description = 'Start angle in degrees', default = 90, max = 360, min = -360)

    drawSides: IntProperty(name = "Polygon / Star Sides", description = 'Sides of polygon', \
        default = 4, max = 100, min = 3, update = ModalDrawBezierOp.updateDrawSides)

    drawAngleSweep: FloatProperty(name = "Arc Sweep", \
        description = 'Arc sweep in degrees', default = 360, max = 360, min = -360)

    drawStarOffset: FloatProperty(name = "Offset", \
        description = 'Offset of star sides', default = .3)

    snapOrient: EnumProperty(
        name = 'Transform Orientation',
        items = get_orientation_items,
        description='Transform orientation for Draw / Edit operations. Context-aware: unavailable options show ⚠ warning')

    snapOrigin: EnumProperty(
        name = 'Pivot Point',
        items = get_origin_items,
        description='Pivot point for Draw / Edit operations. Context-aware: unavailable options show ⚠ warning')

    constrAxes: EnumProperty(name = 'Constrain Axis', #"Constrain axis for draw and edit ops",
        items = getConstrAxisTups,
        description='Constrain Draw / Edit Axes')

    snapToPlane: BoolProperty(name="Snap to Plane",
        description='During draw / edit snap the point to the selected plane', \
                    default = False)

    axisScale: EnumProperty(name="Scale", \
        items = (('DEFAULT', 'Default Scale', 'Use default scale'),
                 ('REFERENCE', 'Reference Line Scale', \
                  'Use Reference Line scale (1 Unit = 0.1 x Reference Line Length)'), \
                 ('AXIS', 'Custom Axis Scale', \
                  'Use Custom Axis scale (1 Unit = 0.1 x Custom Axis Length)')),
        description='Scale to use for grid snap and transform values entered', \
                    default = 'DEFAULT')

    customAxisSnapCnt: IntProperty(default = 3, min = 0)

    copyPropsObj : PointerProperty(
            name = 'Copy Object Properties',
            description = "Copy properties (Material, Bevel Depth etc.) from object",
            type = bpy.types.Object)



    ############################ Menu ###############################

    for menudata in FTMenu.editMenus:
        exec(FTMenu.getMNPropDefStr(menudata))


