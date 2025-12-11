import bpy
from bpy.types import Panel
from bpy.app.handlers import persistent
from gpu_extras.batch import batch_for_shader
import gpu
from ..operators.simple_ops import *
from ..utils.curve_utils import getBezierDataForSeg
from ..utils.bezier_math import getPtsAlongBezier2D, getLinesFromPts
from ..utils.view_utils import getAllAreaRegions, getCoordFromLoc
from ..core.props import FTProps
from ..constants import MAX_NONSEL_CURVE_RES
from ..tools.workspace_tools import FlexiDrawBezierTool
from ..drawing.primitives import Primitive2DDraw
from ..drawing.math_fn import MathFnDraw

# BezierUtilsPanel will be added here
class BezierUtilsPanel(Panel):
    bl_label = "Bezier Utilities"
    bl_idname = "CURVE_PT_bezierutils"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Tool'

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_CURVE'}

    def draw(self, context):
        params = bpy.context.window_manager.bezierToolkitParams

        layout = self.layout
        layout.use_property_decorate = False

        if(context.mode == 'OBJECT'):


            # Combined Curve Operations Section
            row = layout.row()
            row.prop(params, "curveOpsExpanded",
                icon="TRIA_DOWN" if params.curveOpsExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Curve Operations', icon='MOD_BOOLEAN')
            
            if params.curveOpsExpanded:
                box = layout.box()
                
                # Shared proximity parameter
                col = box.column().split()
                row = col.row()
                row.prop(params, 'intersectMargin', text='Proximity')
                
                # Operation mode selector
                col = box.column().split()
                row = col.row()
                row.prop(params, 'curveOpsMode', expand=True)
                
                # Intersect-specific options
                if params.curveOpsMode == 'INTERSECT':
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'intersectOp', text='Action')
                    
                    if params.intersectOp in {'INSERT_PT', 'CUT'}:
                        row = col.row()
                        row.prop(params, 'intersectNonactive', text='Only Non-active')
                    
                    row = col.row()
                    row.prop(params, 'selfIntersect', text='Self Intersection')
                    
                    col = box.column().split()
                    col.operator('object.intersect_curves')
                
                # Boolean-specific options
                elif params.curveOpsMode == 'BOOLEAN':
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'booleanOp', text='Operation')
                    
                    col = box.column().split()
                    col.operator('object.boolean_curves')


            row = layout.row()
            row.prop(params, "splitExpanded",
                icon="TRIA_DOWN" if params.splitExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False)
            row.label(text="Split Curves", icon = 'UNLINKED')

            if params.splitExpanded:
                box = layout.box()
                col = box.column().split()
                col.operator('object.separate_splines')
                col = box.column().split()
                col.operator('object.separate_segments')
                col = box.column().split()
                col.operator('object.separate_points')

            row = layout.row()
            row.prop(params, "joinExpanded",
                icon="TRIA_DOWN" if params.joinExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text="Join Curves", icon = 'LINKED')

            if params.joinExpanded:
                box = layout.box()
                col = box.column().split()
                col.prop(params, 'straight')
                col = box.column().split()
                col.prop(params, 'optimized')
                col = box.column().split()
                col.prop(params, 'joinMergeDist')
                col = box.column().split()
                col.operator('object.join_curves')

            row = layout.row()
            row.prop(params, "alignToFaceExpanded",
                icon="TRIA_DOWN" if params.alignToFaceExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text="Align to Face", icon = 'FCURVE')

            if params.alignToFaceExpanded:
                box = layout.box()
                col = box.column().split()
                col.prop(params, 'alignToFaceOrig')
                col = box.column().split()
                col.prop(params, 'alignToFaceLoc')
                col = box.column().split()
                col.operator('object.align_to_face')

            row = layout.row()
            row.prop(params, "selectExpanded",
                icon="TRIA_DOWN" if params.selectExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Select Objects In Collection', icon='RESTRICT_SELECT_OFF')
            if params.selectExpanded:
                box = layout.box()
                col = box.column().split()
                row = col.row()
                row.prop(params, 'selectIntrvl')
                row.operator('object.select_in_collection')
                col = box.column().split()
                col.operator('object.invert_sel_in_collection')

            row = layout.row()
            row.prop(params, "convertExpanded",
                icon="TRIA_DOWN" if params.convertExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Convert Curve to Mesh', icon='MESH_DATA')

            if params.convertExpanded:
                box = layout.box()
                col = box.column().split()
                col.prop(params, 'fillType')
                if params.fillType == 'QUAD':
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'remeshDepth')
                    row.prop(params, 'unsubdivide')
                
                elif params.fillType in {'GRID', 'OFFSET'}:
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'fillDetail')
                    
                    if params.fillType == 'OFFSET':
                        row.prop(params, 'offsetSize')
                        
                elif params.fillType == 'QUADRIFLOW':
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'quadriflowFaces')
                    row.prop(params, 'quadriflowSeed')

                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'quadriflowPreserveSharp')
                    row.prop(params, 'quadriflowPreserveBoundary')

                elif params.fillType in {'SMART', 'MEDIAL'}:
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'fillDetail')
                    row.prop(params, 'offsetSize')

                # Common Params (Curve Resolution) - Valid for all or most
                # Assuming 'remeshRes' controls the initial curve sampling
                col = box.column().split()
                row = col.row()
                row.prop(params, 'remeshRes') 
                # row.prop(params, 'remeshApplyTo') # Make optional/advanced? Keeping for now if relevant
                col = box.column().split()
                col.operator('object.convert_2d_mesh')

            row = layout.row()
            row.prop(params, "handleTypesExpanded",
                icon="TRIA_DOWN" if params.handleTypesExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Set Handle Type', icon='MOD_CURVE')

            if params.handleTypesExpanded:
                box = layout.box()
                col = box.column().split()
                row = col.row()
                col = row.column()
                col.prop(params, 'handleType')
                col = box.column().split()
                col.operator('object.set_handle_types')

            row = layout.row()
            row.prop(params, "removeDupliExpanded",
                icon="TRIA_DOWN" if params.removeDupliExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Remove Duplicate Vertices', icon='X')
            if params.removeDupliExpanded:
                box = layout.box()
                col = box.column().split()
                row = col.row()
                row.prop(params, 'dupliVertMargin', text = 'Proximity')
                col = box.column().split()
                col.operator('object.remove_dupli_vert_curve')

            ######## Curve Color #########

            row = layout.row()
            row.prop(params, "curveColorExpanded",
                icon="TRIA_DOWN" if params.curveColorExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Set Curve Colors', icon='MATERIAL')

            if params.curveColorExpanded:
                box = layout.box()
                col = box.column().split()
                row = col.row()
                row.prop(params, "curveColorPick", text = 'Curve Color')
                row.operator('object.set_curve_color')
                row.operator('object.remove_curve_color')
                col = box.column().split()
                row = col.row()
                row.prop(params, 'applyCurveColor', toggle = True)

            ######## Other Tools #########

            row = layout.row()
            row.prop(params, "otherExpanded",
                icon="TRIA_DOWN" if params.otherExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Other Tools', icon='TOOL_SETTINGS')

            if params.otherExpanded:
                box = layout.box()
                col = box.column().split()
                col.operator('object.export_svg')
                col = box.column().split()
                col.operator('object.paste_length')
                col = box.column().split()
                col.operator('object.close_splines')
                col = box.column().split()
                col.operator('object.close_straight')
                col = box.column().split()
                col.operator('object.open_splines')

            tool = context.workspace.tools.from_space_view3d_mode('OBJECT', \
                create = False)

            if(tool.idname == FlexiDrawBezierTool.bl_idname and \
                params.drawObjType == 'MATH'):
                row = layout.row()
                row.prop(params, "mathExtraExpanded",
                    icon="TRIA_DOWN" if params.mathExtraExpanded else "TRIA_RIGHT",
                    icon_only=True, emboss=False
                )
                row.label(text='Flexi Draw Math Function', icon='GRAPH')
                if params.mathExtraExpanded:
                    # Equation params are duplicated in the toolbar...
                    # any changes here should also reflect there
                    box = layout.box()

                    col = box.column().split()
                    col.prop(params, "mathFnList")

                    col = box.column().split()
                    col.prop(params, "mathFnName")
                    col = box.column().split()
                    col.prop(params, "mathFnDescr")

                    col = box.column().split()
                    col.prop(params, "mathFnType")
                    col = box.column().split()
                    col.prop(params, "mathFnResolution")

                    if(params.mathFnType == 'PARAMETRIC'):
                        col = box.column().split()
                        col.prop(params, "drawMathFnParametric1")
                        col = box.column().split()
                        col.prop(params, "drawMathFnParametric2")
                        col = box.column().split()
                        col.prop(params, "drawMathTMapTo")
                        col = box.column().split()
                        col.prop(params, "drawMathTScaleFact")
                        col = box.column().split()
                        col.prop(params, "drawMathTStart")
                    else:
                        col = box.column().split()
                        col.prop(params, "drawMathFn")
                        col = box.column().split()
                        col.prop(params, "mathFnclipVal")

                    paramCol = box.column()
                    for i in range(Primitive2DDraw.getParamCnt()):
                        char = chr(ord('A') + i)
                        innerBox = paramCol.box()
                        col = innerBox.column().split()
                        row = col.row()
                        row.label(text = char)
                        row.prop(params, MathFnDraw.startPrefix + str(i), \
                            text = '') # Value
                        row.prop(params, MathFnDraw.incrPrefix + str(i), \
                            text = '') # Step

                    col = box.column().split()
                    row = col.row()
                    row.operator('object.save_math_fn')
                    row.operator('object.reset_math_fn')
                    row.operator('object.load_math_fn', text = 'Import')
                    row.operator('object.delete_math_fn', text = 'Delete')
        else:
            col = layout.column()
            col.operator('object.separate_segments', text = 'Split At Selected Points')
            col = layout.column()
            col.prop(params, 'markVertex', toggle = True)


    ################ Stand-alone handler for changing curve colors #################

    drawHandlerRef = None
    shader = None
    lineBatch = None
    lineWidth = 1.5

    @persistent
    def colorCurves(scene = None, add = False, remove = False):
        def ccDrawHandler():
            if(bpy.context.window_manager.bezierToolkitParams.applyCurveColor):
                # bgl.glLineWidth(BezierUtilsPanel.lineWidth)
                gpu.state.line_width_set(BezierUtilsPanel.lineWidth)
                if(BezierUtilsPanel.lineBatch is not None):
                    BezierUtilsPanel.lineBatch.draw(BezierUtilsPanel.shader)

        if(add and BezierUtilsPanel.drawHandlerRef is None):
            try:
                addon_prefs = bpy.context.preferences.addons.get('bezier_utils')
                if addon_prefs and hasattr(addon_prefs.preferences, 'lineWidth'):
                    lineWidth = addon_prefs.preferences.lineWidth
                else:
                    lineWidth = 1.5  # Default value
            except Exception as e:
                print("BezierUtils: Error fetching line width in ColorCurves: ", e)
                lineWidth = 1.5  # Default value
            BezierUtilsPanel.lineWidth = lineWidth

            BezierUtilsPanel.drawHandlerRef = \
                bpy.types.SpaceView3D.draw_handler_add(ccDrawHandler, \
                    (), "WINDOW", "POST_VIEW")
            BezierUtilsPanel.shader = gpu.shader.from_builtin('FLAT_COLOR')
            return

        elif(remove):
            if(BezierUtilsPanel.drawHandlerRef is not None):
                bpy.types.SpaceView3D.draw_handler_remove(BezierUtilsPanel.drawHandlerRef, \
                    "WINDOW")
                BezierUtilsPanel.drawHandlerRef = None
                return

        if(bpy.context.screen is None):
            return

        if(bpy.context.window_manager.bezierToolkitParams.applyCurveColor):
            objs = [o for o in bpy.context.scene.objects if(isBezier(o) and \
                o.visible_get() and len(o.modifiers) == 0 and not o.select_get())]

            lineCos = []
            lineColors = []
            for o in objs:
                colorVal = o.data.get('curveColor')
                if(colorVal is not None):
                    for i, spline in enumerate(o.data.splines):
                        for j in range(0, len(spline.bezier_points)):
                            segPts = getBezierDataForSeg(o, i, j, withShapeKey = True, \
                                shapeKeyIdx = None, fromMix = True)
                            if(segPts is None):
                                continue
                            pts = getPtsAlongBezier2D(segPts, getAllAreaRegions(), \
                                FTProps.dispCurveRes, getCoordFromLoc, maxRes = MAX_NONSEL_CURVE_RES)
                            linePts = getLinesFromPts(pts)
                            lineCos += linePts
                            lineColors += [colorVal for i in range(0, len(linePts))]
            BezierUtilsPanel.lineBatch = batch_for_shader(BezierUtilsPanel.shader, \
                "LINES", {"pos": lineCos, "color": lineColors})
        # ~ else:
            # ~ BezierUtilsPanel.lineBatch = batch_for_shader(BezierUtilsPanel.shader, \
                # ~ "LINES", {"pos": [], "color": []})

            areas = [a for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
            for a in areas:
                a.tag_redraw()


################### Common Bezier Functions & Classes ###################


def updatePanel(self, context):
    try:
        panel = BezierUtilsPanel
        if "bl_rna" in panel.__dict__:
            bpy.utils.unregister_class(panel)
    except Exception as e:
        print("BezierUtils: Unregistering Panel has failed", e)
        return

    try:
        # Use the root package name instead of __name__
        addon_prefs = context.preferences.addons.get('bezier_utils')
        if addon_prefs and hasattr(addon_prefs.preferences, 'category'):
            panel.bl_category = addon_prefs.preferences.category
        else:
            # Default to 'Tool' category if preferences aren't available yet
            panel.bl_category = 'Tool'
        bpy.utils.register_class(panel)

    except Exception as e:
        print("BezierUtils: Updating Panel locations has failed", e)
        # Fallback to default category
        panel.bl_category = 'Tool'
        try:
            bpy.utils.register_class(panel)
        except:
            pass  # Panel might already be registered


