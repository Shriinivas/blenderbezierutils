# bezier_utils/operators/simple_ops.py

import bpy
from bpy.types import Operator
from bpy.props import StringProperty, BoolProperty, EnumProperty, FloatProperty
from ..utils.curve_utils import (
    splitCurve, isBezier, safeRemoveObj, intersectCurves, booleanCurves
)

# Import operators will be added here after extraction
class SeparateSplinesObjsOp(Operator):

    bl_idname = "object.separate_splines"
    bl_label = "Separate Splines"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate splines of selected Bezier curves as new objects"

    def execute(self, context):

        selObjs = bpy.context.selected_objects
        newObjs, changeCnt = splitCurve(selObjs, split = 'spline')

        if(changeCnt > 0):
            self.report({'INFO'}, "Separated "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + \
                    " into " +str(len(newObjs)) + " new ones")

        return {'FINISHED'}


class SplitBezierObjsOp(Operator):

    bl_idname = "object.separate_segments"
    bl_label = "Separate Segments"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate segments of selected Bezier curves as new objects"

    def execute(self, context):
        selObjs = [o for o in bpy.context.selected_objects if(isBezier(o))]
        if(context.mode == 'EDIT_CURVE'):
            selPtMap = {}

            for o in selObjs:
                selPtMap[o] = {}
                for i, s in enumerate(o.data.splines):
                    pts = s.bezier_points
                    ptIdxs = [x for x in range(0, len(pts)) if pts[x].select_control_point]
                    if(len(ptIdxs) > 0): selPtMap[o][i] = ptIdxs

            newObjs, changeCnt = splitCurveSelPts(selPtMap)
        else:
            newObjs, changeCnt = splitCurve(selObjs, split = 'seg')

        if(changeCnt > 0):
            bpy.context.view_layer.objects.active = newObjs[-1]
            self.report({'INFO'}, "Split "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + \
                    " into " + str(len(newObjs)) + " new objects")

        return {'FINISHED'}


class splitBezierObjsPtsOp(Operator):

    bl_idname = "object.separate_points"
    bl_label = "Separate Points"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate bezier points of selected curves as new objects"

    def execute(self, context):
        selObjs = bpy.context.selected_objects
        newObjs, changeCnt = splitCurve(selObjs, split = 'point')

        if(changeCnt > 0):
            bpy.context.view_layer.objects.active = newObjs[-1]
            self.report({'INFO'}, "Split "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + \
                    " into " + str(len(newObjs)) + " new objects")

        return {'FINISHED'}


class JoinBezierSegsOp(Operator):
    bl_idname = "object.join_curves"
    bl_label = "Join"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Join selected curves (with new segments if required)"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        straight = bpy.context.window_manager.bezierToolkitParams.straight
        optimized = bpy.context.window_manager.bezierToolkitParams.optimized
        mergeDist = bpy.context.window_manager.bezierToolkitParams.joinMergeDist

        newCurve = joinSegs(curves, optimized = optimized, straight = straight, \
            margin = mergeDist)

        # ~ removeShapeKeys(newCurve)
        bpy.context.view_layer.objects.active = newCurve

        return {'FINISHED'}


class InvertSelOp(Operator):
    bl_idname = "object.invert_sel_in_collection"
    bl_label = "Invert Selection"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Invert selection within collection of active object"

    def execute(self, context):
        if(bpy.context.active_object is not None):
            collections = bpy.context.active_object.users_collection

            for collection in collections:
                for o in collection.objects:
                    o.select_set(not o.select_get())

        return {'FINISHED'}


class SelectInCollOp(Operator):
    bl_idname = "object.select_in_collection"
    bl_label = "Select"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Select objects within collection of active object"

    def execute(self, context):
        if(bpy.context.active_object is not None):
            selectIntrvl = bpy.context.window_manager.bezierToolkitParams.selectIntrvl

            for obj in bpy.context.selected_objects:
                collections = obj.users_collection
                for collection in collections:
                    objs = [o for o in collection.objects]
                    idx = objs.index(obj)
                    for i, o in enumerate(objs[idx:]):
                        if( i % (selectIntrvl + 1) == 0): o.select_set(True)
                        else: o.select_set(False)
                    for i, o in enumerate(reversed(objs[:idx+1])):
                        if( i % (selectIntrvl + 1) == 0): o.select_set(True)
                        else: o.select_set(False)

        return {'FINISHED'}


class CloseStraightOp(Operator):
    bl_idname = "object.close_straight"
    bl_label = "Close Splines With Straight Segment"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Close all splines in selected curves with straight segmennt"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        for curve in curves:
            for spline in curve.data.splines:
                spline.bezier_points[0].handle_right_type = 'FREE'
                spline.bezier_points[0].handle_left_type = 'VECTOR'
                spline.bezier_points[-1].handle_left_type = 'FREE'
                spline.bezier_points[-1].handle_right_type = 'VECTOR'
                spline.use_cyclic_u = True

        return {'FINISHED'}


class CloseSplinesOp(Operator):
    bl_idname = "object.close_splines"
    bl_label = "Close Splines"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Close all splines in selected curves"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        for curve in curves:
            for spline in curve.data.splines:
                # ~ spline.bezier_points[0].handle_left_type = 'ALIGNED'
                # ~ spline.bezier_points[-1].handle_right_type = 'ALIGNED'
                spline.use_cyclic_u = True

        return {'FINISHED'}


class OpenSplinesOp(Operator):
    bl_idname = "object.open_splines"
    bl_label = "Open up Splines"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Open up all splines in selected curves"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        for curve in curves:
            for spline in curve.data.splines:
                spline.use_cyclic_u = False

        return {'FINISHED'}


class SetHandleTypesOp(Operator):
    bl_idname = "object.set_handle_types"
    bl_label = "Set Handle Type of All Points"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Set the handle type of all the points of the selected curves"

    def execute(self, context):
        ht = bpy.context.window_manager.bezierToolkitParams.handleType
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        for curve in curves:
            for spline in curve.data.splines:
                for pt in spline.bezier_points:
                    pt.handle_left_type = ht
                    pt.handle_right_type = ht

        return {'FINISHED'}


class RemoveDupliVertCurveOp(Operator):
    bl_idname = "object.remove_dupli_vert_curve"
    bl_label = "Remove Duplicate Curve Vertices"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Remove duplicate vertices and mark splines as cyclic if applicable"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        for curve in curves:
            removeDupliVert(curve, \
                bpy.context.window_manager.bezierToolkitParams.dupliVertMargin)

        return {'FINISHED'}


class convertToMeshOp(Operator):
    bl_idname = "object.convert_2d_mesh"
    bl_label = "Convert"
    bl_description = "Convert 2D curve to mesh with quad faces"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]
        params = bpy.context.window_manager.bezierToolkitParams
        remeshDepth = params.remeshDepth
        unsubdivide = params.unsubdivide
        fillType = params.fillType
        optimized = params.remeshOptimized
        remeshRes = params.remeshRes
        perSeg = (params.remeshApplyTo == 'PERSEG')

        for curve in curves:
            center, normal = None, None
            if(fillType == 'QUAD'):
                for spline in curve.data.splines:
                    spline.use_cyclic_u = True
                curve.data.dimensions = '2D'
                curve.data.fill_mode = 'BOTH'
                meshObj = convertToMesh(curve)

                applyMeshModifiers(meshObj, remeshDepth)

                if(unsubdivide):
                    unsubdivideObj(meshObj)
            else:
                meshObj, center, normal = \
                    convertToFace(curve, remeshRes, perSeg, fillType, optimized)

            meshObj.matrix_world = curve.matrix_world.copy()
            meshObj.select_set(True)
            if(center is not None and normal is not None):
                newOrig = meshObj.matrix_world @ center
                shiftOrigin(meshObj, newOrig)
                quatMat = normal.to_track_quat('Z', 'X').to_matrix().to_4x4()
                tm = meshObj.matrix_world @ quatMat
                shiftMatrixWorld(meshObj, tm)
                meshObj.location = newOrig # Location not from tm

            safeRemoveObj(curve)

        return {'FINISHED'}


class AlignToFaceOp(Operator):
    bl_idname = "object.align_to_face"
    bl_label = "Align to Face"
    bl_description = "Align all points of selected curves with " + \
                     "nearest face of selected meshes"

    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):

        alignLoc = bpy.context.window_manager.bezierToolkitParams.alignToFaceLoc
        alignOrig = bpy.context.window_manager.bezierToolkitParams.alignToFaceOrig

        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o) and o.visible_get()]

        if(len(curves) == 0): return {'FINISHED'}

        mesheObjs = [o for o in bpy.data.objects if o in bpy.context.selected_objects \
            and  o.type == 'MESH' and o.visible_get()]

        if(len(mesheObjs) == 0): return {'FINISHED'}

        depsgraph = bpy.context.evaluated_depsgraph_get()
        medianLists = []

        for o in mesheObjs:
            medianLists.append([o.matrix_world @ f.center for f in o.data.polygons])

        searchTree = NestedListSearch(medianLists)

        for curve in curves:
            center = getObjBBoxCenter(curve)
            oLoc = curve.location.copy()

            srs = searchTree.findInLists(center, searchRange = None)
            if(len(srs) != 1): continue

            objIdx, faceIdx, median, dist = srs[0]
            meshObj = mesheObjs[objIdx]
            faceCenter = meshObj.matrix_world @ meshObj.data.polygons[faceIdx].center
            normal = meshObj.data.polygons[faceIdx].normal
            quatMat = normal.to_track_quat('Z', 'X').to_matrix().to_4x4()
            tm = meshObj.matrix_world @ quatMat

            shiftMatrixWorld(curve, tm)

            invTm = tm.inverted_safe()
            centerLocal =  invTm @ center

            for spline in curve.data.splines:
                for bpt in spline.bezier_points:
                    lht = bpt.handle_left_type
                    rht = bpt.handle_right_type
                    bpt.handle_left_type = 'FREE'
                    bpt.handle_right_type = 'FREE'
                    bpt.co[2] = centerLocal[2]
                    bpt.handle_right[2] = centerLocal[2]
                    bpt.handle_left[2] = centerLocal[2]
                    bpt.handle_left_type = lht
                    bpt.handle_right_type = rht

            if(alignOrig == 'FACE'):
                shiftOrigin(curve, faceCenter)
            elif(alignOrig == 'BBOX'):
                depsgraph.update()
                center = getObjBBoxCenter(curve) # Recalculate after alignment
                shiftOrigin(curve, center)
            else:
                shiftOrigin(curve, oLoc)

            if(alignLoc): curve.location = faceCenter

        return {'FINISHED'}


class SetCurveColorOp(bpy.types.Operator):
    bl_idname = "object.set_curve_color"
    bl_label = "Set Color"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Set color of selected curves"

    def execute(self, context):
        curves = bpy.context.selected_objects

        for curve in curves:
            curve.data['curveColor'] = \
                bpy.context.window_manager.bezierToolkitParams.curveColorPick

        return {'FINISHED'}


class RemoveCurveColorOp(bpy.types.Operator):
    bl_idname = "object.remove_curve_color"
    bl_label = "Remove Color"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Remove color of selected curves"

    def execute(self, context):
        curves = bpy.context.selected_objects

        for curve in curves:
            if(curve.data.get('curveColor')):
                del curve.data['curveColor']

        return {'FINISHED'}


class PasteLengthOp(Operator):
    bl_idname = "object.paste_length"
    bl_label = "Paste Length"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Make the selected curves the same length as the active one"

    def execute(self, context):
        src = bpy.context.object
        if(src is not None and isBezier(src)):
            dests = [o for o in bpy.context.selected_objects if(isBezier(o) and o != src)]
            if(len(dests) > 0):
                pasteLength(src, dests)
        return {'FINISHED'}

class IntersectCurvesOp(Operator):
    bl_idname = "object.intersect_curves"
    bl_label = "Intersect Curves"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Intersect the selected curves"

    def execute(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        rounding = 2
        actCurve = bpy.context.active_object

        if(params.intersectNonactive):
            if(not isBezier(actCurve)):
                self.report({'WARNING'}, "Active Object Not A Curve")
                return {'FINISHED'}

        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o) and o != actCurve]

        if(actCurve is not None and isBezier(actCurve)):
            curves = [actCurve] + curves

        if(len(curves) < 2):
            self.report({'INFO'}, "Please select at least two Bezier curve objects")
        else:
            intersectCurves(curves, params.intersectOp, params.intersectNonactive, \
            params.intersectMargin, rounding)

        return {'FINISHED'}


class BooleanCurvesOp(Operator):
    bl_idname = "object.boolean_curves"
    bl_label = "Boolean Curves"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Perform boolean operation on two closed curves"

    def execute(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        rounding = 2
        actCurve = bpy.context.active_object

        if not isBezier(actCurve):
            self.report({'WARNING'}, "Active object is not a Bezier curve")
            return {'FINISHED'}

        curves = [o for o in bpy.context.selected_objects 
                  if isBezier(o) and o != actCurve]

        if len(curves) < 1:
            self.report({'INFO'}, "Select two closed Bezier curves")
            return {'FINISHED'}

        # Active curve first, then other selected
        curves = [actCurve, curves[0]]

        # Check all splines are closed
        for c in curves:
            if not all(s.use_cyclic_u for s in c.data.splines):
                self.report({'WARNING'}, "Boolean requires closed curves")
                return {'FINISHED'}

        result = booleanCurves(curves, params.booleanOp, 
                               params.intersectMargin, rounding)

        if result:
            for obj in result:
                obj.select_set(True)
            if result:
                bpy.context.view_layer.objects.active = result[0]
            self.report({'INFO'}, f"Boolean {params.booleanOp} completed")
        else:
            self.report({'WARNING'}, "Boolean operation produced no result")

        return {'FINISHED'}

class ExportSVGOp(Operator):

    bl_idname = "object.export_svg"
    bl_label = "Export to SVG"
    bl_options = {'REGISTER', 'UNDO'}

    def getExportViewList(scene = None, context = None):
        cameras = [o for o in bpy.data.objects if o.type == 'CAMERA']
        vlist = [('ACTIVE_VIEW', 'Viewport View', "Export Viewport View")]
        for c in cameras:
             vlist.append((c.name, c.name, 'Export view from ' + c.name))
        return vlist

    filepath : StringProperty(subtype='FILE_PATH')

    #User input
    clipView : BoolProperty(name="Clip View", \
        description = "Clip objects to view boundary", \
            default = True)

    exportView: EnumProperty(name = 'Export View',
        items = getExportViewList,
        description='View to export')

    lineWidth: FloatProperty(name="Line Width", \
        description='Line width in exported SVG', default = 3, min = 0)

    lineColorOpts: EnumProperty(name = 'Line Color',
        items = (('RANDOM', 'Random', 'Use random color for curves'),
                 ('PICK', 'Pick', 'Pick color'),
        ),
        description='Color to draw curve lines')

    lineColor: bpy.props.FloatVectorProperty(
        name="Line Color",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(0.5, 0.5, 0.5, 1.0)
    )

    fillColorOpts: EnumProperty(name = 'Fill Color',
        items = (('RANDOM', 'Random', 'Use random fill color'),
                 ('PICK', 'Pick', 'Pick color'),
        ),
        description='Color to fill solid curves')

    fillColor: bpy.props.FloatVectorProperty(
        name="Fill Color",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(0, 0.3, 0.5, 1.0)
    )

    def execute(self, context):
        exportSVG(context, self.filepath, self.exportView, self.clipView, self.lineWidth, \
            self.lineColorOpts, self.lineColor, self.fillColorOpts, self.fillColor)
        return {'FINISHED'}

    def draw(self, context):
        layout = self.layout
        col = layout.column()
        row = col.row()
        row.prop(self, "exportView")
        col = layout.column()
        row = col.row()
        row.prop(self, "clipView")
        col = layout.column()
        row = col.row()
        row.prop(self, "lineWidth")
        col = layout.column()
        row = col.row()
        row.prop(self, "lineColorOpts")
        if(self.lineColorOpts == 'PICK'):
            row = row.split()
            row.prop(self, "lineColor", text = '')
        col = layout.column()
        row = col.row()
        row.prop(self, "fillColorOpts")
        if(self.fillColorOpts == 'PICK'):
            row = row.split()
            row.prop(self, "fillColor", text = '')

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


def markVertHandler(self, context):
    if(self.markVertex):
        bpy.ops.wm.bb_mark_vertex()
