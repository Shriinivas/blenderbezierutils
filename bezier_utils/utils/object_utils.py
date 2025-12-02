# bezier_utils/utils/object_utils.py

import bpy
from mathutils import Vector, Matrix, geometry

def isBezier(bObj):
    return bObj.type == 'CURVE' and len(bObj.data.splines) > 0 \
        and bObj.data.splines[0].type == 'BEZIER' and  \
            len(bObj.data.splines[0].bezier_points) > 0

def safeRemoveObj(obj):
    try:
        collections = obj.users_collection

        for c in collections:
            c.objects.unlink(obj)

        if(obj.data.users == 1):
            if(obj.type == 'MESH'):
                bpy.data.meshes.remove(obj.data)
            elif(obj.type == 'CURVE'):
                bpy.data.curves.remove(obj.data)
            #else? TODO
        bpy.data.objects.remove(obj)
    except:
        pass

#TODO: Fix this hack if possible
def copyObjAttr(src, dest, invDestMW = Matrix(), mw = Matrix()):
    for att in dir(src):
        try:
            if(att not in ['co', 'handle_left', 'handle_right', \
                'handle_left_type', 'handle_right_type']):
                setattr(dest, att, getattr(src, att))
        except Exception as e:
            pass
    try:
        lt = src.handle_left_type
        rt = src.handle_right_type
        dest.handle_left_type = 'FREE'
        dest.handle_right_type = 'FREE'
        dest.co = invDestMW @ (mw @ src.co)
        dest.handle_left = invDestMW @ (mw @ src.handle_left)
        dest.handle_right = invDestMW @ (mw @ src.handle_right)
        dest.handle_left_type = lt
        dest.handle_right_type = rt
        pass
    except Exception as e:
        pass

def getObjBBoxCenter(obj):
    bbox = obj.bound_box
    return obj.matrix_world @ Vector(((bbox[0][0] + bbox[4][0]) / 2, \
        (bbox[0][1] + bbox[3][1]) / 2, (bbox[0][2] + bbox[1][2]) / 2))

# Only mesh and Bezier curve
def shiftOrigin(obj, origin):
    oLoc = obj.location.copy()
    mw = obj.matrix_world
    invMw = mw.inverted_safe()
    if(obj.type == 'MESH'):
        for vert in obj.data.vertices:
            vert.co += invMw @ oLoc - invMw @ origin
    elif(obj.type == 'CURVE'):
        for s in obj.data.splines:
            bpts = s.bezier_points
            for bpt in bpts:
                lht = bpt.handle_left_type
                rht = bpt.handle_right_type
                bpt.handle_left_type = 'FREE'
                bpt.handle_right_type = 'FREE'

                bpt.co += invMw @ oLoc - invMw @ origin
                bpt.handle_left += invMw @ oLoc - invMw @ origin
                bpt.handle_right += invMw @ oLoc - invMw @ origin

                bpt.handle_left_type = lht
                bpt.handle_right_type = rht

    obj.location = origin

# Only mesh and Bezier curve; depsgraph not updated
def shiftMatrixWorld(obj, mw):
    invMw = mw.inverted_safe()
    omw = obj.matrix_world

    if(obj.type == 'MESH'):
        for vert in obj.data.vertices:
            vert.co = invMw @ (omw @ vert.co)
    elif(obj.type == 'CURVE'):
        for s in obj.data.splines:
            bpts = s.bezier_points
            for bpt in bpts:
                lht = bpt.handle_left_type
                rht = bpt.handle_right_type
                bpt.handle_left_type = 'FREE'
                bpt.handle_right_type = 'FREE'
                bpt.co = invMw @ (omw @ bpt.co)
                bpt.handle_left = invMw @ (omw @ bpt.handle_left)
                bpt.handle_right = invMw @ (omw @ bpt.handle_right)

                bpt.handle_left_type = lht
                bpt.handle_right_type = rht

    obj.matrix_world = mw

# Also shifts origin; depsgraph not updated
def alignToNormal(curve):
    depsgraph = bpy.context.evaluated_depsgraph_get()
    depsgraph.update()
    mw = curve.matrix_world.copy()
    normals = []
    for spline in curve.data.splines:
        bpts = spline.bezier_points
        bptCnt = len(bpts)
        if(bptCnt > 2):
            normals.append(geometry.normal(mw @ bpts[i].co for i in range(bptCnt)))
    cnt = len(normals)
    if(cnt > 0):
        normal =  Vector([sum(normals[i][j] for i in range(cnt)) \
            for j in range(3)]) / cnt
        quatMat = normal.to_track_quat('Z', 'X').to_matrix().to_4x4()
        shiftMatrixWorld(curve, quatMat)

def copyProperties(srcObj, destCurve):
    if(srcObj is None or destCurve is None):
        return

    destData = destCurve.data
    srcData = srcObj.data

    # If object is bezier curve copy curve properties and material
    if(isBezier(srcObj)):
        # Copying just a few attributes
        destData.dimensions = srcData.dimensions

        destData.resolution_u = srcData.resolution_u
        destData.render_resolution_u = srcData.render_resolution_u
        destData.fill_mode = srcData.fill_mode

        destData.use_fill_deform = srcData.use_fill_deform
        destData.use_radius = srcData.use_radius
        destData.use_stretch = srcData.use_stretch
        destData.use_deform_bounds = srcData.use_deform_bounds

        destData.twist_smooth = srcData.twist_smooth
        destData.twist_mode = srcData.twist_mode

        destData.offset = srcData.offset
        destData.extrude = srcData.extrude
        destData.bevel_depth = srcData.bevel_depth
        destData.bevel_resolution = srcData.bevel_resolution
        destData.bevel_object = srcData.bevel_object
        destData.taper_object = srcData.taper_object

        destData.use_fill_caps = srcData.use_fill_caps

    if(hasattr(srcData, 'materials') and len(srcData.materials) > 0):
        mat = srcData.materials[srcObj.active_material_index]
        if(len(destData.materials) == 0 or mat.name not in destData.materials):
            destData.materials.append(mat)
            activeIdx = -1 #Last
        else:
            activeIdx = destData.materials.find(mat.name)

        destCurve.active_material_index = activeIdx

def reverseCurve(curve):
    cp = curve.data.copy()
    curve.data.splines.clear()
    for s in reversed(cp.splines):
        ns = curve.data.splines.new('BEZIER')
        copyObjAttr(s, ns)
        ns.bezier_points.add(len(s.bezier_points) - 1)
        for i, p in enumerate(reversed(s.bezier_points)):
            copyObjAttr(p, ns.bezier_points[i])
            ns.bezier_points[i].handle_left_type = 'FREE'
            ns.bezier_points[i].handle_right_type = 'FREE'
            ns.bezier_points[i].handle_left = p.handle_right
            ns.bezier_points[i].handle_right = p.handle_left

            ns.bezier_points[i].handle_left_type = p.handle_right_type
            ns.bezier_points[i].handle_right_type = p.handle_left_type
    bpy.data.curves.remove(cp)
