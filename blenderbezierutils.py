#
#
# Blender add-on with tools to draw and edit Bezier curves along with other utility ops
#
# Supported Blender Version: 2.8x
#
# Copyright (C) 2019  Shrinivas Kulkarni

# License: GPL (https://github.com/Shriinivas/blenderbezierutils/blob/master/LICENSE)
#

import os, bpy, bmesh, bgl, blf, gpu
from bpy.props import BoolProperty, IntProperty, EnumProperty, \
FloatProperty, StringProperty, CollectionProperty, FloatVectorProperty, PointerProperty
from bpy.types import Panel, Operator, WorkSpaceTool, AddonPreferences, Menu
from mathutils import Vector, Matrix, geometry, kdtree
from math import e, pi, log, sin, cos, tan, radians, degrees, sqrt, asin, acos, atan, floor, \
ceil, pow, exp
from bpy_extras.view3d_utils import region_2d_to_vector_3d, region_2d_to_location_3d, \
region_2d_to_origin_3d
from bpy_extras.view3d_utils import location_3d_to_region_2d
from bpy_extras.object_utils import world_to_camera_view
from gpu_extras.batch import batch_for_shader
import random, time, math
from bpy.app.handlers import persistent
from xml.dom import minidom
from shutil import copyfile

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "version": (0, 9, 96),
    "location": "Properties > Active Tool and Workspace Settings > Bezier Utilities",
    "description": "Collection of Bezier curve utility ops",
    "category": "Object",
    "wiki_url": "https://github.com/Shriinivas/blenderbezierutils/blob/master/README.md",
    "blender": (2, 80, 0),
}

DEF_ERR_MARGIN = 0.0001

# Markers for invalid data
LARGE_NO = 9e+9 # Both float and int
LARGE_VECT = Vector((LARGE_NO, LARGE_NO, LARGE_NO))
INVAL = '~'

###################### Common functions ######################

def floatCmpWithMargin(float1, float2, margin = DEF_ERR_MARGIN):
    return abs(float1 - float2) < margin

def vectCmpWithMargin(v1, v2, margin = DEF_ERR_MARGIN):
    return all(floatCmpWithMargin(v1[i], v2[i], margin) for i in range(0, len(v1)))

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

#TODO combine with copyObjAttr
def copyBezierPt(src, target, freeHandles = None, srcMw = Matrix(), invDestMW = Matrix()):
    target.handle_left_type = 'FREE'
    target.handle_right_type = 'FREE'

    target.co = invDestMW @ (srcMw @ src.co)
    target.handle_left = invDestMW @ (srcMw @ src.handle_left)
    target.handle_right = invDestMW @ (srcMw @ src.handle_right)

    if(freeHandles == None or freeHandles[0] == False):
        target.handle_left_type = src.handle_left_type
    if(freeHandles == None or freeHandles[1] == False):
        target.handle_right_type =  src.handle_right_type

def createSplineForSeg(curveData, bezierPts):
    spline = curveData.splines.new('BEZIER')
    spline.bezier_points.add(len(bezierPts)-1)
    spline.use_cyclic_u = False

    for i, pt in enumerate(bezierPts):
        if(i == 0): freeHandles = [False, True]
        elif(i == len(bezierPts) - 1): freeHandles = [True, False]
        else: freeHandles = None
        copyBezierPt(pt, spline.bezier_points[i], freeHandles = freeHandles)

    return spline

def createSpline(curveData, srcSpline, excludePtIdxs = {}):
    spline = curveData.splines.new('BEZIER')
    spline.bezier_points.add(len(srcSpline.bezier_points) - len(excludePtIdxs) - 1)
    spline.use_cyclic_u = srcSpline.use_cyclic_u

    ptIdx = 0
    for i in range(0, len(srcSpline.bezier_points)):
        if(i not in excludePtIdxs):
            copyBezierPt(srcSpline.bezier_points[i], \
                spline.bezier_points[ptIdx], freeHandles = None)
            ptIdx += 1

    return spline

def createSkeletalCurve(obj, collections):
    objCopy = obj.copy()
    objCopy.name = obj.name
    dataCopy = obj.data.copy()
    dataCopy.splines.clear()
    objCopy.data = dataCopy

    for coll in collections:
        coll.objects.link(objCopy)

    return objCopy

def createObjFromPts(curvePts, dimensions = '3D', collection = None, \
    closed = False, calcHdlTypes = True):
    data = bpy.data.curves.new('BezierCurve', 'CURVE')
    data.dimensions = dimensions
    obj = bpy.data.objects.new('BezierCurve', data)
    # ~ collection = context.collection
    if(collection == None):
        collection = bpy.context.scene.collection
    collection.objects.link(obj)
    # ~ obj.location = context.scene.cursor.location

    # ~ depsgraph = context.evaluated_depsgraph_get()
    # ~ depsgraph.update()

    # ~ invM = obj.matrix_world.inverted_safe()

    spline = data.splines.new('BEZIER')
    spline.use_cyclic_u = False

    if(vectCmpWithMargin(curvePts[0][1], curvePts[-1][1])):
        curvePts[0][0] = curvePts[-1][0]
        spline.use_cyclic_u = True
        curvePts.pop()

    if(closed): spline.use_cyclic_u = True

    spline.bezier_points.add(len(curvePts) - 1)
    prevPt = None
    for i, pt in enumerate(curvePts):
        currPt = spline.bezier_points[i]
        currPt.co = pt[1]
        currPt.handle_right = pt[2]
        if(not calcHdlTypes and len(pt) > 3):
            currPt.handle_right_type = pt[3]
            currPt.handle_left_type = pt[4]
        elif(prevPt != None and prevPt.handle_right == prevPt.co \
            and pt[0] == pt[1] and currPt.co != prevPt.co): # straight line
                if(prevPt.handle_left_type != 'VECTOR'):
                    prevPt.handle_left_type = 'FREE'
                prevPt.handle_right_type = 'VECTOR'
                currPt.handle_right_type = 'FREE'
                currPt.handle_left_type = 'VECTOR'
        else:
            currPt.handle_left_type = 'FREE'
            currPt.handle_right_type = 'FREE'
            currPt.handle_left = pt[0]
            ldiffV = (pt[1] - pt[0])
            rdiffV = (pt[2] - pt[1])
            if(vectCmpWithMargin(ldiffV, rdiffV) and \
                not floatCmpWithMargin(ldiffV.length, 0)):
                currPt.handle_left_type = 'ALIGNED'
                currPt.handle_right_type = 'ALIGNED'
        prevPt = currPt

    bpts = spline.bezier_points
    if(spline.use_cyclic_u and vectCmpWithMargin(bpts[-1].handle_right, bpts[-1].co) \
        and vectCmpWithMargin(bpts[0].handle_left, bpts[0].co)):
            if(bpts[-1].handle_left_type != 'VECTOR'):
                bpts[-1].handle_left_type = 'FREE'
            bpts[-1].handle_right_type = 'VECTOR'
            if(bpts[0].handle_right_type != 'VECTOR'):
                bpts[0].handle_right_type = 'FREE'
            bpts[0].handle_left_type = 'VECTOR'

    return obj

def removeShapeKeys(obj):
    if(obj.data.shape_keys == None):
        return

    keyblocks = reversed(obj.data.shape_keys.key_blocks)
    for sk in keyblocks:
        obj.shape_key_remove(sk)

def getShapeKeyInfo(obj):
    keyData = []
    keyNames = []

    if(obj.data.shape_keys != None):
        keyblocks = obj.data.shape_keys.key_blocks
        for key in keyblocks:
            keyData.append([[d.handle_left.copy(), d.co.copy(), d.handle_right.copy()] \
                for d in key.data])
            keyNames.append(key.name)

    return keyNames, keyData

def updateShapeKeyData(obj, keyData, keyNames, startIdx, cnt = None, add = False):
    if(obj.data.shape_keys == None and not add):
        return

    currIdx = obj.active_shape_key_index
    if(not add): removeShapeKeys(obj)
    if(cnt == None): cnt = len(keyData[0])

    for i, name in enumerate(keyNames):
        key = obj.shape_key_add(name = name)
        for j in range(0, cnt):
            keyIdx = j + startIdx
            key.data[j].handle_left = keyData[i][keyIdx][0].copy()
            key.data[j].co = keyData[i][keyIdx][1].copy()
            key.data[j].handle_right = keyData[i][keyIdx][2].copy()

    obj.active_shape_key_index = currIdx

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

def getLastSegIdx(obj, splineIdx):
    spline = obj.data.splines[splineIdx]
    ptCnt = len(spline.bezier_points)
    return ptCnt - 1 if(spline.use_cyclic_u) else ptCnt - 2

def addLastSeg(spline):
    if(spline.use_cyclic_u):
        lt0 = spline.bezier_points[0].handle_left_type
        rt0 = spline.bezier_points[0].handle_right_type
        pt = spline.bezier_points[0]
        pt.handle_left_type = 'FREE'
        pt.handle_right_type = 'FREE'
        spline.use_cyclic_u = False
        spline.bezier_points.add(1)
        copyObjAttr(spline.bezier_points[0], spline.bezier_points[-1])
        spline.bezier_points[0].handle_left_type = lt0
        spline.bezier_points[-1].handle_right_type = rt0

def moveSplineStart(obj, splineIdx, idx):
    pts = obj.data.splines[splineIdx].bezier_points
    cnt = len(pts)

    ptCopy = [[p.co.copy(), p.handle_right.copy(), \
        p.handle_left.copy(), p.handle_right_type, \
            p.handle_left_type] for p in pts]

    for i, pt in enumerate(pts):
        srcIdx = (idx + i) % cnt
        p = ptCopy[srcIdx]

        pt.handle_left_type = 'FREE'
        pt.handle_right_type = 'FREE'
        pt.co = p[0]
        pt.handle_right = p[1]
        pt.handle_left = p[2]
        pt.handle_right_type = p[3]
        pt.handle_left_type = p[4]

def joinCurves(curves):
    obj = curves[0]
    invMW = obj.matrix_world.inverted_safe()
    for curve in curves[1:]:
        mw = curve.matrix_world
        for spline in curve.data.splines:
            newSpline = obj.data.splines.new('BEZIER')
            copyObjAttr(spline, newSpline)
            newSpline.bezier_points.add(len(spline.bezier_points)-1)
            for i, pt in enumerate(spline.bezier_points):
                copyObjAttr(pt, newSpline.bezier_points[i], \
                    invDestMW = invMW, mw = mw)
        safeRemoveObj(curve)
    return obj

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
    if(srcObj == None or destCurve == None):
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

# Insert spline at location insertIdx, duplicated from existing spline at
# location srcSplineIdx and remove points with indices in removePtIdxs from new spline
def insertSpline(obj, srcSplineIdx, insertIdx, removePtIdxs):
    srcSpline = obj.data.splines[srcSplineIdx]
    # Appended at end
    createSpline(obj.data, srcSpline, removePtIdxs)
    splineCnt = len(obj.data.splines)
    nextIdx = insertIdx
    for idx in range(nextIdx, splineCnt - 1):
        srcSpline = obj.data.splines[nextIdx]
        createSpline(obj.data, srcSpline)
        obj.data.splines.remove(srcSpline)

def removeBezierPts(obj, splineIdx, removePtIdxs):
    oldSpline = obj.data.splines[splineIdx]
    bpts = oldSpline.bezier_points
    if(min(removePtIdxs) >= len(bpts)):
        return

    if(len(set(range(len(bpts))) - set(removePtIdxs)) == 0) :
        obj.data.splines.remove(oldSpline)
        if(len(obj.data.splines) == 0):
            safeRemoveObj(obj)
        return

    insertSpline(obj, splineIdx, splineIdx, removePtIdxs)
    obj.data.splines.remove(obj.data.splines[splineIdx + 1])

# Returns a tuple with first value indicating change in spline index (-1, 0, 1)
# and second indicating shift in seg index (negative) due to removal
def removeBezierSeg(obj, splineIdx, segIdx):
    nextIdx = getAdjIdx(obj, splineIdx, segIdx)
    if(nextIdx == None): return
    spline = obj.data.splines[splineIdx]
    bpts = spline.bezier_points
    ptCnt = len(bpts)
    lastSegIdx = getLastSegIdx(obj, splineIdx)
    splineIdxIncr = 0
    segIdxIncr = 0
    if(ptCnt <= 2):
        removeBezierPts(obj, splineIdx, {segIdx, nextIdx})
        # Spline removed by above call
        splineIdxIncr = -1
    else:
        bpt = obj.data.splines[splineIdx].bezier_points[segIdx]
        bpt.handle_right_type = 'FREE'
        bpt.handle_left_type = 'FREE'
        nextIdx = getAdjIdx(obj, splineIdx, segIdx)
        bpt = obj.data.splines[splineIdx].bezier_points[nextIdx]
        bpt.handle_right_type = 'FREE'
        bpt.handle_left_type = 'FREE'
        if(spline.use_cyclic_u):
            spline.use_cyclic_u = False
            if(segIdx != lastSegIdx):
                moveSplineStart(obj, splineIdx, getAdjIdx(obj, splineIdx, segIdx))
                segIdxIncr = - (segIdx + 1)
        else:
            if(segIdx == lastSegIdx):
                removeBezierPts(obj, splineIdx, {lastSegIdx + 1})
            elif(segIdx == 0):
                removeBezierPts(obj, splineIdx, {0})
                segIdxIncr = -1
            else:
                insertSpline(obj, splineIdx, splineIdx, set(range(segIdx + 1, ptCnt)))
                removeBezierPts(obj, splineIdx + 1, range(segIdx + 1))
                splineIdxIncr = 1
                segIdxIncr = - (segIdx + 1)
    return splineIdxIncr, segIdxIncr

def insertBezierPts(obj, splineIdx, startIdx, cos, handleType, margin = DEF_ERR_MARGIN):

    spline = obj.data.splines[splineIdx]
    bpts = spline.bezier_points

    nextIdx = getAdjIdx(obj, splineIdx, startIdx)

    firstPt = bpts[startIdx]
    nextPt = bpts[nextIdx]

    if(firstPt.handle_right_type == 'AUTO'):
        firstPt.handle_left_type = 'ALIGNED'
        firstPt.handle_right_type = 'ALIGNED'
    if(nextPt.handle_left_type == 'AUTO'):
        nextPt.handle_left_type = 'ALIGNED'
        nextPt.handle_right_type = 'ALIGNED'

    fhdl = firstPt.handle_right_type
    nhdl = nextPt.handle_left_type

    firstPt.handle_right_type = 'FREE'
    nextPt.handle_left_type = 'FREE'

    ptCnt = len(bpts)
    addCnt = len(cos)

    bpts.add(addCnt)
    nextIdx = startIdx + 1

    for i in range(0, (ptCnt - nextIdx)):
        idx = ptCnt - i - 1# reversed
        offsetIdx = idx + addCnt
        copyObjAttr(bpts[idx], bpts[offsetIdx])

    endIdx = getAdjIdx(obj, splineIdx, nextIdx, addCnt)
    firstPt = bpts[startIdx]
    nextPt = bpts[endIdx]

    prevPt = firstPt
    for i, pt in enumerate(bpts[nextIdx:nextIdx + addCnt]):
        pt.handle_left_type = 'FREE'
        pt.handle_right_type = 'FREE'

        co = cos[i]
        seg = [prevPt.co, prevPt.handle_right, nextPt.handle_left, nextPt.co]
        t = getTForPt(seg, co, margin)
        ctrlPts0 = getPartialSeg(seg, 0, t)
        ctrlPts1 = getPartialSeg(seg, t, 1)

        segPt = [ctrlPts0[2], ctrlPts1[0], ctrlPts1[1]]

        prevRight = ctrlPts0[1]
        nextLeft =  ctrlPts1[2]

        pt.handle_left = segPt[0]
        pt.co = segPt[1]
        pt.handle_right = segPt[2]
        pt.handle_left_type = handleType
        pt.handle_right_type = handleType
        prevPt.handle_right = prevRight

        prevPt = pt
        nextPt.handle_left = nextLeft

    firstPt.handle_right_type = fhdl
    nextPt.handle_left_type = nhdl

# https://devtalk.blender.org/t/get-hex-gamma-corrected-color/2422/2
def toHexStr(rgba):
    ch = []
    for c in rgba[:3]:
        if c < 0.0031308:
            cc = 0.0 if c < 0.0 else c * 12.92
        else:
            cc = 1.055 * pow(c, 1.0 / 2.4) - 0.055
        ch.append(hex(max(min(int(cc * 255 + 0.5), 255), 0))[2:])
    return ''.join(ch), str(rgba[-1])

# Change position of bezier points according to new matrix_world
def changeMW(obj, newMW):
    invMW = newMW.inverted_safe()
    for spline in obj.data.splines:
        for pt in spline.bezier_points:
            pt.co = invMW @ (obj.mw @ pt.co)
            pt.handle_left = invMW @ (obj.mw @ pt.handle_left)
            pt.handle_right = invMW @ (obj.mw @ pt.handle_right)
    obj.matrix_world = newMW

# Return map in the form of objName->[splineIdx, [startPt, endPt]]
# Remove the invalid keys (if any) from it.
def updateCurveEndPtMap(endPtMap, addObjNames = None, removeObjNames = None):
    invalOs = set()
    if(addObjNames == None):
        addObjNames = [o.name for o in bpy.context.scene.objects]
        invalOs = endPtMap.keys() - set(addObjNames) # In case of redo

    if(removeObjNames != None):
        invalOs.union(set(removeObjNames))

    for o in invalOs:
        del endPtMap[o]

    for objName in addObjNames:
        obj = bpy.context.scene.objects.get(objName)
        if(obj != None and isBezier(obj) and obj.visible_get()):
            endPtMap[objName] = []
            mw = obj.matrix_world
            for i, s in enumerate(obj.data.splines):
                pts = [mw @ pt.co for pt in s.bezier_points]
                endPtMap[objName].append([i, pts])
        elif(endPtMap.get(objName) != None):
            del endPtMap[objName]

    return endPtMap

#Round to logarithmic scale .1, 0, 10, 100 etc.
#(47.538, -1) -> 47.5; (47.538, 0) -> 48.0; (47.538, 1) -> 50.0; (47.538, 2) -> 0,
# TODO: Rework after grid subdiv is enabled (in a version later than 2.8)
def roundedVect(space3d, vect, rounding, axes):
    rounding += 1
    subdiv = getGridSubdiv(space3d)
    # TODO: Separate logic for 1
    if(subdiv == 1): subdiv = 10
    fact = ((subdiv ** rounding) / subdiv) / getUnitScale()
    retVect = vect.copy()
    # ~ Vector([round(vect[i] / fact) * fact for i in axes])
    for i in axes: retVect[i] = round(vect[i] / fact) * fact
    return retVect

###################### Screen functions ######################

def getGridSubdiv(space3d):
    return space3d.overlay.grid_subdivisions

def getUnit():
    return bpy.context.scene.unit_settings.length_unit

def getUnitSystem():
    return bpy.context.scene.unit_settings.system

def getUnitScale():
    fact = 3.28084 if(getUnitSystem() == 'IMPERIAL') else 1
    return fact * bpy.context.scene.unit_settings.scale_length

def get3dLoc(region, rv3d, xy, vec = None):
    if(vec == None):
        vec = region_2d_to_vector_3d(region, rv3d, xy)
    return region_2d_to_location_3d(region, rv3d, xy, vec)

# TODO: Rework after grid subdiv is enabled (in a version later than 2.8)
def  getViewDistRounding(space3d, rv3d):
    viewDist = rv3d.view_distance * getUnitScale()
    gridDiv = getGridSubdiv(space3d)
    subFact = 1
    # TODO: Separate logic for 1
    if(gridDiv == 1): gridDiv = 10
    elif(gridDiv == 2): subFact = 5
    elif(viewDist < 0.5): subFact = 2
    return int(log(viewDist, gridDiv)) - subFact

# Return axis-indices (x:0, y:1, z:2) of plane with closest orientation
# to view
def getClosestPlaneToView(rv3d):
    viewmat = rv3d.view_matrix
    trans, quat, scale = viewmat.decompose()
    tm = quat.to_matrix().to_4x4()
    viewnormal = tm.inverted() @ Vector((0, 0, 1))

    xynormal = [Vector((0, 0, 1)), [0, 1]]
    yznormal = [Vector((1, 0, 0)), [1, 2]]
    xznormal = [Vector((0, 1, 0)), [0, 2]]
    normals = [xynormal, yznormal, xznormal]
    minAngle = pi
    minIdx = None
    for idx, normal in enumerate(normals):
        rotDiff = viewnormal.rotation_difference(normal[0]).angle
        if(rotDiff > pi / 2):
            rotDiff = pi - rotDiff
        if(rotDiff < minAngle):
            minIdx = idx
            minAngle = rotDiff
    return normals[minIdx][1]

def getCoordFromLoc(region, rv3d, loc):
    coord = location_3d_to_region_2d(region, rv3d, loc)
    # return a unlocatable pt if None to avoid errors
    return coord if(coord != None) else Vector((9000, 9000))

# To be called only from 3d view
def getCurrAreaRegion(context):
    a, r = [(a, r) for a in bpy.context.screen.areas if a.type == 'VIEW_3D' \
        for r in a.regions if(r == context.region)][0]
    return a, r

def isOutside(context, event, exclInRgns = True):
    x = event.mouse_region_x
    y = event.mouse_region_y
    region = context.region

    if(x < 0 or x > region.width or y < 0 or y > region.height):
        return True

    elif(not exclInRgns):
        return False

    area, r = getCurrAreaRegion(context)

    for r in area.regions:
        if(r == region):
            continue
        xR = r.x - region.x
        yR = r.y - region.y
        if(x >= xR and y >= yR and x <= (xR + r.width) and y <= (yR + r.height)):
            return True

    return False

def getPtProjOnPlane(region, rv3d, xy, p1, p2, p3, p4 = None):
    vec = region_2d_to_vector_3d(region, rv3d, xy)
    orig = region_2d_to_origin_3d(region, rv3d, xy)

    pt = geometry.intersect_ray_tri(p1, p2, p3, vec, orig, False)#p4 != None)
    # ~ if(not pt and p4):
        # ~ pt = geometry.intersect_ray_tri(p2, p4, p3, vec, orig, True)
    return pt

# find the location on 3d line p1-p2 if xy is already on 2d projection (in rv3d) of p1-p2
def getPtProjOnLine(region, rv3d, xy, p1, p2):
    # Just find a non-linear point (TODO: simpler way)
    pd1 = p2 - p1
    pd2 = Vector(sorted(pd1, key=lambda x: abs(x), reverse = True))
    maxIdx0 = [i for i in range(3) if abs(pd1[i]) == abs(pd2[0])][0]
    maxIdx1 = [i for i in range(3) if abs(pd1[i]) == abs(pd2[1])][0]
    pd = Vector()
    pd[maxIdx0] = -pd2[1]
    pd[maxIdx1] = pd2[0]
    p3 = p2 + pd

    # Raycast from 2d point onto the plane
    return getPtProjOnPlane(region, rv3d, xy[:2], p1, p2, p3)

def getLineTransMatrices(pt0, pt1):
    diffV = (pt1 - pt0)
    invTm = diffV.to_track_quat('X', 'Z').to_matrix().to_4x4()
    tm = invTm.inverted_safe()
    return tm, invTm

def getWindowRegionIdx(area, regionIdx): # For finding quad view index
    idx = 0
    for j, r in enumerate(area.regions):
        if(j == regionIdx): return idx
        if(r.type == 'WINDOW'): idx += 1
    return None

def getAreaRegionIdxs(xy, exclInRgns = True):
    x, y = xy
    areas = [a for a in bpy.context.screen.areas]
    idxs = None
    for i, a in enumerate(areas):
        if(a.type != 'VIEW_3D'): continue
        regions = [r for r in a.regions]
        for j, r in enumerate(regions):
            if(x > r.x and x < r.x + r.width and y > r.y and y < r.y + r.height):
                if(r.type == 'WINDOW'):
                    if(not exclInRgns):
                        return [i, j]
                    idxs = [i, j]
                elif(exclInRgns):
                    return None
    return idxs

# ~ def getMinViewDistRegion():
    # ~ viewDist = 9e+99
    # ~ rv3d = None
    # ~ area = None
    # ~ region = None
    # ~ areas = [a for a in bpy.context.screen.areas if(a.type == 'VIEW_3D')]
    # ~ for a in areas:
        # ~ regions = [r for r in a.regions if r.type == 'WINDOW']
        # ~ if(len(a.spaces[0].region_quadviews) > 0):
            # ~ for i, r in enumerate(a.spaces[0].region_quadviews):
                # ~ if(r.view_distance < viewDist):
                    # ~ viewDist = r.view_distance
                    # ~ rv3d = r
                    # ~ area = a
                    # ~ region = regions[i]
        # ~ else:
            # ~ r = a.spaces[0].region_3d
            # ~ if(r.view_distance < viewDist):
                # ~ viewDist = r.view_distance
                # ~ rv3d = r
                # ~ area = a
                # ~ region = regions[0]
    # ~ return a, region, rv3d#, viewDist

def getAllAreaRegions():
    info = []
    areas = []
    i = 0

    areas = [a for a in bpy.context.screen.areas if(a.type == 'VIEW_3D')]

    # bpy.context.screen doesn't work in case of Add-on Config window
    while(len(areas) == 0 and i < len(bpy.data.screens)):
        areas = [a for a in bpy.data.screens[i].areas if(a.type == 'VIEW_3D')]
        i += 1

    for a in areas:
        regions = [r for r in a.regions if r.type == 'WINDOW']
        if(len(a.spaces[0].region_quadviews) > 0):
            for i, r in enumerate(a.spaces[0].region_quadviews):
                info.append([a, regions[i], r])
        else:
            r = a.spaces[0].region_3d
            info.append([a, regions[0], r])
    return info

def getResetBatch(shader, btype): # "LINES" or "POINTS"
    return batch_for_shader(shader, btype, {"pos": [], "color": []})

# From python template
def getFaceUnderMouse(obj, region, rv3d, xy, maxFaceCnt):
    if(obj == None or obj.type != 'MESH' \
        or len(obj.data.polygons) > maxFaceCnt):
        return None, None, None
    viewVect = region_2d_to_vector_3d(region, rv3d, xy)
    rayOrig = region_2d_to_origin_3d(region, rv3d, xy)
    mw = obj.matrix_world
    invMw = mw.inverted_safe()

    rayTarget = rayOrig + viewVect
    rayOrigObj = invMw @ rayOrig
    rayTargetObj = invMw @ rayTarget
    rayDirObj = rayTargetObj - rayOrigObj

    success, location, normal, faceIdx = obj.ray_cast(rayOrigObj, rayDirObj)
    if(success):
        return mw @ location, normal, faceIdx
    else:
        return None, None, None

def getSnappableObjs(region, rv3d, xy):
    objs = bpy.context.selected_objects
    if(bpy.context.object != None):
        objs.append(bpy.context.object)
    return [o for o in objs if(o.type == 'MESH' and len(o.modifiers) == 0 and \
            isPtIn2dBBox(o, region, rv3d, xy))]

# precise can be pretty expensive with large vert count
def get2dBBox(obj, region, rv3d, precise = False):
    mw = obj.matrix_world
    if(precise):
        co2ds = [getCoordFromLoc(region, rv3d, mw @ Vector(v.co)) \
            for v in obj.data.vertices]
    else:
        co2ds = [getCoordFromLoc(region, rv3d, mw @ Vector(b)) for b in obj.bound_box]
    minX = min(c[0] for c in co2ds)
    maxX = max(c[0] for c in co2ds)
    minY = min(c[1] for c in co2ds)
    maxY = max(c[1] for c in co2ds)

    return minX, minY, maxX, maxY

def isPtIn2dBBox(obj, region, rv3d, xy, extendBy = 0, precise = False):
    minX, minY, maxX, maxY = get2dBBox(obj, region, rv3d, precise)
    if(xy[0] > (minX - extendBy) and xy[0] < (maxX + extendBy) \
        and xy[1] > (minY - extendBy) and xy[1] < (maxY + extendBy)):
        return True
    else: return False

# ~ def isLocIn2dBBox(obj, region, rv3d, loc, extendBy = 0, precise = False):
    # ~ xy = getCoordFromLoc(region, rv3d, loc)
    # ~ return isPtIn2dBBox(obj, region, rv3d, xy, extendBy, precise)

def getClosestEdgeLoc2d(obj, region, rv3d, xy, faceIdx = None):
    mw = obj.matrix_world
    minDist = LARGE_NO
    closestLoc = None
    pt = Vector(xy).to_3d()
    edgeWSCos = None
    edgeIdx = None
    closestIntersect = None
    vertPairs = obj.data.polygons[faceIdx].edge_keys if faceIdx != None \
        else [e.vertices for e in obj.data.edges]
    for i, vertPair in enumerate(vertPairs):
        co0 = mw @ obj.data.vertices[vertPair[0]].co
        co1 = mw @ obj.data.vertices[vertPair[1]].co

        pt0 = getCoordFromLoc(region, rv3d, co0).to_3d()
        pt1 = getCoordFromLoc(region, rv3d, co1).to_3d()
        intersect, percDist = geometry.intersect_point_line(pt, pt0, pt1)
        if(percDist < 0):
            intersect = pt0
            percDist = 0
        elif(percDist > 1):
            intersect = pt1
            percDist = 1
        dist = (intersect - pt).length
        if(dist < minDist):
            minDist = dist
            closestIntersect = intersect
            edgeWSCos = [co0, co1]
            edgeIdx = i
    if(edgeWSCos != None):
        closestLoc = getPtProjOnLine(region, rv3d, closestIntersect, \
            edgeWSCos[0], edgeWSCos[1])

    return edgeIdx, edgeWSCos, closestLoc, minDist

# TODO: Fix the signature
def getSelFaceLoc(region, rv3d, xy, maxFaceCnt, objs = None, checkEdge = False):
    if(objs == None): objs = getSnappableObjs(region, rv3d, xy)
    if(len(objs) > maxFaceCnt): return None, None, None, None
    for obj in objs:
        loc, normal, faceIdx = getFaceUnderMouse(obj, region, rv3d, xy, maxFaceCnt)
        if(loc != None):
            if(checkEdge):
                edgeIdx, edgeWSCos, closestLoc, minDist = \
                    getClosestEdgeLoc2d(obj, region, rv3d, xy, faceIdx)
                if(closestLoc != None and minDist < FTProps.snapDist):
                    return obj, closestLoc, normal, faceIdx
            return obj, loc, normal, faceIdx
    return None, None, None, None

# TODO: Fix the signature
# ~ def getSelEdgeLoc(region, rv3d, xy, maxFaceCnt, objs = None):

    # ~ if(objs == None): objs = getSnappableObjs()
    # ~ if(len(objs) > maxFaceCnt): return None, None, None, None

    # ~ minDist = LARGE_NO
    # ~ closestLoc = None
    # ~ closestEdgeIdx = None
    # ~ closestObj = None
    # ~ objCnt = 0
    # ~ objs = bpy.context.selected_objects
    # ~ if(bpy.context.object != None): objs.append(bpy.context.object)
    # ~ for obj in objs:
        # ~ edgeIdx, edgeWSCos, loc, dist = \
            # ~ getClosestEdgeLoc2d(obj, region, rv3d, xy)
        # ~ if(dist < minDist):
            # ~ minDist = dist
            # ~ closestLoc = loc
            # ~ closestEdgeIdx = edgeIdx
            # ~ closestObj = obj
    # ~ return closestObj, closestLoc, minDist, closestEdgeIdx

###################### Op Specific functions ######################

def closeSplines(curve, htype = None):
    for spline in curve.data.splines:
        if(htype != None):
            spline.bezier_points[0].handle_left_type = htype
            spline.bezier_points[-1].handle_right_type = htype
        spline.use_cyclic_u = True

# TODO: Update shapekey (not working due to moving of start pt in cyclic)
def splitCurveSelPts(selPtMap, newColl = True):
    changeCnt = 0
    newObjs = []

    if(len(selPtMap) == 0): return newObjs, changeCnt

    for obj in selPtMap.keys():
        splinePtMap = selPtMap.get(obj)

        if((len(obj.data.splines) == 1 and \
            len(obj.data.splines[0].bezier_points) <= 2 and \
                not obj.data.splines[0].use_cyclic_u) or len(splinePtMap) == 0):
            continue

        keyNames, keyData = getShapeKeyInfo(obj)
        collections = obj.users_collection

        if(newColl):
            objGrp = bpy.data.collections.new(obj.name)
            parentColls = [objGrp]
        else:
            parentColls = collections

        splineCnt = len(obj.data.splines)


        endSplineIdx = splineCnt- 1
        if(endSplineIdx not in splinePtMap.keys()):
            splinePtMap[endSplineIdx] = \
                [len(obj.data.splines[endSplineIdx].bezier_points) - 1]

        splineIdxs = sorted(splinePtMap.keys())

        lastSplineIdx = -1
        objCopy = createSkeletalCurve(obj, parentColls)
        newObjs.append(objCopy)
        for i in splineIdxs:
            for j in range(lastSplineIdx + 1, i):
                srcSpline = obj.data.splines[j]
                createSpline(objCopy.data, srcSpline)
                # ~ updateShapeKeyData(objCopy, keyData, keyNames, skStart, ptCnt)
            srcSpline = obj.data.splines[i]
            selPtIdxs = sorted(splinePtMap[i])

            if(len(selPtIdxs) == 0):
                createSpline(objCopy.data, srcSpline)
            else:
                bpts = srcSpline.bezier_points
                cyclic = srcSpline.use_cyclic_u
                if(cyclic):
                    firstIdx = selPtIdxs[0]
                    moveSplineStart(obj, i, firstIdx)
                    selPtIdxs = [getAdjIdx(obj, i, s, -firstIdx) for s in selPtIdxs]
                    addLastSeg(srcSpline)
                if(len(selPtIdxs) > 0 and selPtIdxs[0] == 0):
                    selPtIdxs.pop(0)
                if(len(selPtIdxs) > 0 and selPtIdxs[-1] == len(bpts) - 1):
                    selPtIdxs.pop(-1)
                bpts = srcSpline.bezier_points

                if(len(selPtIdxs) == 0):
                    segBpts = bpts[:len(bpts)]
                    createSplineForSeg(objCopy.data, segBpts)
                else:
                    lastSegIdx = 0
                    bpts = srcSpline.bezier_points
                    for j in selPtIdxs:
                        segBpts = bpts[lastSegIdx:j + 1]
                        createSplineForSeg(objCopy.data, segBpts)
                        # ~ updateShapeKeyData(objCopy, keyData, keyNames, \
                            # ~ len(newObjs), 2)
                        objCopy = createSkeletalCurve(obj, parentColls)
                        newObjs.append(objCopy)
                        lastSegIdx = j
                    if(j != len(bpts) - 1): createSplineForSeg(objCopy.data, bpts[j:])

            lastSplineIdx = i

        if(len(objCopy.data.splines) == 0):
            newObjs.remove(objCopy)
            safeRemoveObj(objCopy)

        if(newColl):
            for collection in collections:
                collection.children.link(objGrp)

        safeRemoveObj(obj)
        changeCnt += 1

    for obj in newObjs:
        obj.data.splines.active = obj.data.splines[0]

    return newObjs, changeCnt

#split value is one of {'spline', 'seg', 'point'} (TODO: Enum)
def splitCurve(selObjs, split, newColl = True):
    changeCnt = 0
    newObjs = []

    if(len(selObjs) == 0):
        return newObjs, changeCnt

    for obj in selObjs:

        if(not isBezier(obj) or len(obj.data.splines) == 0):
            continue

        if(len(obj.data.splines) == 1):
            if(split == 'spline'):
                newObjs.append(obj)
                continue
            if(split == 'seg' and len(obj.data.splines[0].bezier_points) <= 2):
                newObjs.append(obj)
                continue
            if(split == 'point' and len(obj.data.splines[0].bezier_points) == 1):
                newObjs.append(obj)
                continue

        keyNames, keyData = getShapeKeyInfo(obj)
        collections = obj.users_collection

        if(newColl):
            objGrp = bpy.data.collections.new(obj.name)
            parentColls = [objGrp]
        else:
            parentColls = collections

        segCnt = 0

        for i, spline in enumerate(obj.data.splines):
            if(split == 'seg' or split == 'point'):
                ptLen = len(spline.bezier_points)
                if(split == 'seg'):
                    ptLen -= 1

                for j in range(0, ptLen):
                    objCopy = createSkeletalCurve(obj, parentColls)
                    if(split == 'seg'):
                        createSplineForSeg(objCopy.data, \
                            spline.bezier_points[j:j+2])
                        updateShapeKeyData(objCopy, keyData, keyNames, \
                            len(newObjs), 2)
                    else: #(split == 'point')
                        mw = obj.matrix_world.copy()
                        newSpline = objCopy.data.splines.new('BEZIER')
                        newPtCo = mw @ spline.bezier_points[j].co.copy()
                        newWM = Matrix()
                        newWM.translation = newPtCo
                        objCopy.matrix_world = newWM
                        copyObjAttr(spline.bezier_points[j], \
                            newSpline.bezier_points[0], newWM.inverted_safe(), mw)

                        # No point having shapekeys (pun intended :)
                        removeShapeKeys(objCopy)

                    newObjs.append(objCopy)

                if(split == 'seg' and spline.use_cyclic_u):
                    objCopy = createSkeletalCurve(obj, parentColls)
                    createSplineForSeg(objCopy.data, \
                        [spline.bezier_points[-1], spline.bezier_points[0]])
                    updateShapeKeyData(objCopy, keyData, keyNames, -1, 2)
                    newObjs.append(objCopy)

            else: #split == 'spline'
                objCopy = createSkeletalCurve(obj, parentColls)
                createSpline(objCopy.data, spline)
                currSegCnt = len(objCopy.data.splines[0].bezier_points)
                updateShapeKeyData(objCopy, keyData, keyNames, segCnt, currSegCnt)
                newObjs.append(objCopy)
                segCnt += currSegCnt

        if(newColl):
            for collection in collections:
                collection.children.link(objGrp)

        safeRemoveObj(obj)
        changeCnt += 1

    for obj in newObjs:
        obj.data.splines.active = obj.data.splines[0]

    return newObjs, changeCnt

def getClosestCurve(srcMW, pt, curves, minDist = 9e+99):
    closestCurve = None
    for i, curve in enumerate(curves):
        mw = curve.matrix_world
        addLastSeg(curve.data.splines[0])

        start = curve.data.splines[0].bezier_points[0]
        end = curve.data.splines[-1].bezier_points[-1]
        dist = ((mw @ start.co) - (srcMW @ pt)).length
        if(dist < minDist):
            minDist = dist
            closestCurve = curve
        dist = ((mw @ end.co) - (srcMW @ pt)).length
        if(dist < minDist):
            minDist = dist
            reverseCurve(curve)
            closestCurve = curve
    return closestCurve, minDist

def getCurvesArrangedByDist(curves):
    idMap = {c.name:c for c in curves}
    orderedCurves = [curves[0].name]
    nextCurve = curves[0]
    remainingCurves = curves[1:]

    #Arrange in order
    while(len(remainingCurves) > 0):
        addLastSeg(nextCurve.data.splines[-1])

        srcMW = nextCurve.matrix_world

        ncEnd = nextCurve.data.splines[-1].bezier_points[-1]
        closestCurve, dist = getClosestCurve(srcMW, ncEnd.co, remainingCurves)

        #Check the start also for the first curve
        if(len(orderedCurves) == 1):
            ncStart = nextCurve.data.splines[0].bezier_points[0]
            closestCurve2, dist2 = getClosestCurve(srcMW, ncStart.co, \
                remainingCurves, dist)
            if(closestCurve2 != None):
                reverseCurve(nextCurve)
                closestCurve = closestCurve2

        orderedCurves.append(closestCurve.name)
        nextCurve = closestCurve
        remainingCurves.remove(closestCurve)
    return [idMap[cn] for cn in orderedCurves]

def joinSegs(curves, optimized, straight, srcCurve = None, margin = DEF_ERR_MARGIN):
    if(len(curves) == 0):
        return None
    if(len(curves) == 1):
        return curves[0]

    if(optimized):
        curves = getCurvesArrangedByDist(curves)

    firstCurve = curves[0]

    if(srcCurve == None):
        srcCurve = firstCurve

    elif(srcCurve != firstCurve):
        srcCurveData = srcCurve.data.copy()
        changeMW(firstCurve, srcCurve.matrix_world)
        srcCurve.data = firstCurve.data

    srcMW = srcCurve.matrix_world
    invSrcMW = srcMW.inverted_safe()
    newCurveData = srcCurve.data

    for curve in curves[1:]:
        if(curve == srcCurve):
            curveData = srcCurveData
        else:
            curveData = curve.data

        mw = curve.matrix_world

        currSpline = newCurveData.splines[-1]
        nextSpline = curveData.splines[0]

        addLastSeg(currSpline)
        addLastSeg(nextSpline)

        currBezierPt = currSpline.bezier_points[-1]
        nextBezierPt = nextSpline.bezier_points[0]

        #Don't add new point if the last one and the current one are the 'same'
        if(vectCmpWithMargin(srcMW @ currBezierPt.co, mw @ nextBezierPt.co, margin)):
            currBezierPt.handle_right_type = nextBezierPt.handle_right_type
            if(currBezierPt.handle_right_type != 'VECTOR'):
                currBezierPt.handle_right_type = 'FREE'
            currBezierPt.handle_right = invSrcMW @ (mw @ nextBezierPt.handle_right)
            ptIdx = 1
        else:
            ptIdx = 0

        if(straight and ptIdx == 0):
            currBezierPt.handle_left_type = 'FREE'
            currBezierPt.handle_right_type = 'VECTOR'
            # ~ currBezierPt.handle_right = currBezierPt.co

        for i in range(ptIdx, len(nextSpline.bezier_points)):
            if((i == len(nextSpline.bezier_points) - 1) and
                vectCmpWithMargin(mw @ nextSpline.bezier_points[i].co, \
                    srcMW @ currSpline.bezier_points[0].co, margin)):

                    currSpline.bezier_points[0].handle_left_type = 'FREE'
                    currSpline.bezier_points[0].handle_left = \
                        invSrcMW @ (mw @ nextSpline.bezier_points[i].handle_left)

                    currSpline.use_cyclic_u = True
                    break

            currSpline.bezier_points.add(1)
            currBezierPt = currSpline.bezier_points[-1]
            copyObjAttr(nextSpline.bezier_points[i], currBezierPt, invSrcMW, mw)

            if(straight and i ==  0):
                currBezierPt.handle_right_type = 'FREE'
                currBezierPt.handle_left_type = 'VECTOR'
                # ~ currBezierPt.handle_left = currBezierPt.co

        #Simply add the remaining splines
        for spline in curveData.splines[1:]:
            newSpline = newCurveData.splines.new('BEZIER')
            copyObjAttr(spline, newSpline)
            for i, pt in enumerate(spline.bezier_points):
                if(i > 0):
                    newSpline.bezier_points.add(1)
                copyObjAttr(pt, newSpline.bezier_points[-1], invSrcMW, mw)

        if(curve != srcCurve):
            safeRemoveObj(curve)

    if(firstCurve != srcCurve):
        safeRemoveObj(firstCurve)

    return srcCurve

def removeDupliVert(curve, margin):
    newCurveData = curve.data.copy()
    newCurveData.splines.clear()
    dupliFound = False
    for spline in curve.data.splines:
        newCurveData.splines.new('BEZIER')
        currSpline = newCurveData.splines[-1]
        copyObjAttr(spline, currSpline)

        if(len(spline.bezier_points) == 1):
            copyObjAttr(spline.bezier_points[0], currSpline.bezier_points[0])
            continue

        cmpPts = spline.bezier_points[:]
        pt0 = cmpPts[0]
        while(vectCmpWithMargin(cmpPts[-1].co, pt0.co, margin) and
            len(cmpPts) > 1):
            endPt = cmpPts.pop()
            pt0.handle_left_type = 'FREE'
            pt0.handle_right_type = 'FREE'
            pt0.handle_left = endPt.handle_left
            currSpline.use_cyclic_u = True
            dupliFound = True

        prevPt = None
        for pt in cmpPts:
            if(prevPt != None and vectCmpWithMargin(prevPt.co, pt.co, margin)):
                copyObjAttr(pt, currSpline.bezier_points[-1])
                dupliFound = True
            else:
                if(prevPt != None): currSpline.bezier_points.add(1)
                copyObjAttr(pt, currSpline.bezier_points[-1])
            prevPt = pt

    if(dupliFound):
        curve.data = newCurveData
    else:
        bpy.data.curves.remove(newCurveData)

def convertToFace(curve, remeshRes, perSeg, fillType, optimized):
    bptData = getBptData(curve, local = True)
    bm = bmesh.new()
    splineLens = [spline.calc_length() for spline in curve.data.splines]
    maxSplineLen = max(splineLens)
    centers = []
    normals = []

    for splineIdx, spline in enumerate(curve.data.splines):
        bpts = spline.bezier_points
        verts = []
        addLastVert = spline.use_cyclic_u or fillType != 'NONE'
        if(not perSeg and remeshRes > 0):
            segPts = [bptData[splineIdx][x] for x in range(len(bpts))]
            if(addLastVert): segPts.append(segPts[0])
            numSegs = int(remeshRes * splineLens[splineIdx] / maxSplineLen)
            if(numSegs <= 2):
                vertCos = [bpts[0].co, bpts[-1].co]
            else:
                pts = getInterpBezierPts(segPts, subdivPerUnit = 100, segLens = None)
                vertCos = getInterpolatedVertsCo(pts, numSegs)
            for co in vertCos:
                verts.append(bm.verts.new(co))
        else:
            if(remeshRes > 0):
                segPtPairs = [getSegPtsInSpline(bptData, splineIdx, \
                    ptIdx, addLastVert) for ptIdx in range(len(bpts))]
                if(not addLastVert): segPtPairs.pop()
                segs = [[segPts[0][1], segPts[0][2], segPts[1][0], segPts[1][1]] \
                    for segPts in segPtPairs]
                segLens = [getSegLen(seg) for seg in segs]
                maxLen = max(segLens)
                for i, segPts in enumerate(segPtPairs):
                    if(optimized and isStraightSeg(segPts)):
                        verts.append(bm.verts.new(segPts[0][1]))
                        verts.append(bm.verts.new(segPts[1][1]))
                    else:
                        pts = getInterpBezierPts(segPts, subdivPerUnit = 100, \
                            segLens = [segLens[i]])
                        numSegs = ceil(remeshRes * segLens[i] / maxLen)
                        vertCos = getInterpolatedVertsCo(pts, numSegs)
                        for j, co in enumerate(vertCos):
                            if(j == 0 and i > 0): continue
                            verts.append(bm.verts.new(co))
            else:
                for ptIdx, bpt in enumerate(bpts):
                    verts.append(bm.verts.new(bpts[ptIdx].co))
                if(addLastVert):
                    verts.append(bm.verts.new(bpts[0].co))
        if(len(verts) < 2):
            pass
        elif(len(verts) == 2):
            bm.edges.new(verts)
        else:
            if(spline.use_cyclic_u):
                bm.verts.remove(verts[-1])
                verts.pop()

            vertCos = [v.co for v in verts]
            center =  Vector([sum(vertCos[i][j] for i in range(len(vertCos))) \
                for j in range(3)]) / len(vertCos)
            normal = geometry.normal(vertCos)

            centers.append(center)
            normals.append(normal)

            if(fillType == 'NGON'):
                bm.faces.new(verts)

            elif(fillType == 'NONE'):
                for i in range(1, len(verts)):
                    bm.edges.new([verts[i-1], verts[i]])
                if(spline.use_cyclic_u): bm.edges.new([verts[-1], verts[0]])
            elif(fillType == 'FAN'):
                centerVert = bm.verts.new(center)
                for i in range(1, len(verts)):
                    bm.faces.new([centerVert, verts[i-1], verts[i]])
                if(spline.use_cyclic_u): bm.faces.new([centerVert, verts[-1], verts[0]])
    cnt = len(centers)
    if(cnt > 0):
        center =  Vector([sum(centers[i][j] for i in range(cnt)) \
            for j in range(3)]) / cnt
        normal =  Vector([sum(normals[i][j] for i in range(cnt)) \
            for j in range(3)]) / cnt
    else:
        center =  None
        normal =  None
    m = bpy.data.meshes.new(curve.data.name)
    bm.to_mesh(m)
    meshObj = bpy.data.objects.new(curve.name, m)
    collections = curve.users_collection
    for c in collections:
        c.objects.link(meshObj)

    return meshObj, center, normal

def convertToMesh(curve):
    mt = curve.to_mesh()#Can't be used directly
    bm = bmesh.new()
    bm.from_mesh(mt)
    m = bpy.data.meshes.new(curve.data.name)
    bm.to_mesh(m)
    meshObj = bpy.data.objects.new(curve.name, m)
    collections = curve.users_collection
    for c in collections:
        c.objects.link(meshObj)
    return meshObj

def applyMeshModifiers(meshObj, remeshDepth):
    bpy.context.view_layer.objects.active = meshObj
    normal = geometry.normal([v.co for v in meshObj.data.vertices])
    normal = Vector([round(c, 5) for c in normal])
    if(vectCmpWithMargin(normal, Vector())):
        return meshObj

    planeVert = Vector([round(c, 5) for c in meshObj.data.vertices[0].co])
    mod = meshObj.modifiers.new('mod', type='SOLIDIFY')
    mod.thickness = 20
    bpy.ops.object.modifier_apply(modifier = mod.name)

    mod = meshObj.modifiers.new('mod', type='REMESH')
    mod.octree_depth = remeshDepth
    mod.use_remove_disconnected = False

    bpy.ops.object.modifier_apply(modifier = mod.name)

    bm = bmesh.new()
    bm.from_mesh(meshObj.data)
    bm.verts.ensure_lookup_table()
    toRemove = []
    for i, v in enumerate(bm.verts):
        co = Vector([round(c, 5) for c in v.co])
        if (abs(geometry.distance_point_to_plane(co, planeVert, normal)) > DEF_ERR_MARGIN):
            toRemove.append(v)
    for v in toRemove:
        bm.verts.remove(v)
    bm.to_mesh(meshObj.data)

def unsubdivideObj(meshObj):
    bm = bmesh.new()
    bm.from_object(meshObj, bpy.context.evaluated_depsgraph_get())
    bmesh.ops.unsubdivide(bm, verts = bm.verts)
    bm.to_mesh(meshObj.data)

def pasteLength(src, dests):
    tmp = bpy.data.curves.new('t', 'CURVE')
    ts = tmp.splines.new('BEZIER')
    ts.bezier_points.add(1)
    mw = src.matrix_world
    srcLen = sum(getSplineLenTmpObj(ts, s, mw) for s in src.data.splines)

    for c in dests:
        mw = c.matrix_world
        destLen = sum(getSplineLenTmpObj(ts, s, mw) for s in c.data.splines)
        fact = (srcLen / destLen)
        for s in c.data.splines:
            lts = []
            rts = []
            for pt in s.bezier_points:
                lts.append(pt.handle_left_type)
                rts.append(pt.handle_right_type)
                pt.handle_left_type = 'FREE'
                pt.handle_right_type = 'FREE'
            for pt in s.bezier_points:
                pt.co = fact * pt.co
                pt.handle_left = fact * pt.handle_left
                pt.handle_right = fact * pt.handle_right
            for i, pt in enumerate(s.bezier_points):
                pt.handle_left_type = lts[i]
                pt.handle_right_type = rts[i]
    bpy.data.curves.remove(tmp)

def intersectCurves(curves, action, firstActive, margin, rounding):

    allIntersectCos, intersectMap = getCurveIntersectPts(curves, firstActive, margin, \
        rounding)

    if(action in {'MARK_EMPTY', 'MARK_POINT'}):
        newObjs = []
        collection = bpy.data.collections.new('Intersect Markers')
        bpy.context.scene.collection.children.link(collection)
        objName = 'Marker'
        for co in allIntersectCos:
            if(action  == 'MARK_EMPTY'):
                obj = bpy.data.objects.new(objName, None)
            elif(action  == 'MARK_POINT'):
                curveData = bpy.data.curves.new(objName, 'CURVE')
                curveData.splines.new('BEZIER')
                obj = bpy.data.objects.new(objName, curveData)
            obj.location = co
            newObjs.append(obj)
            collection.objects.link(obj)
        for obj in newObjs:
            obj.select_set(True)

    if(action in {'INSERT_PT', 'CUT'}):
        mapKeys = sorted(intersectMap.keys(), key = lambda x:(x[0], x[1], x[2]))
        selPtMap = {}
        prevCnt = 0
        prevCurveIdx = None
        prevSplineIdx = None
        for i, key in enumerate(mapKeys):
            intersectPts = intersectMap[key]

            curveIdx, splineIdx, segIdx = key
            if(prevCurveIdx == curveIdx and prevSplineIdx == splineIdx):
                segIdx += prevCnt
            else:
                prevCnt = 0

            curve = curves[curveIdx]

            if(firstActive and curve == curves[0]):
                continue

            mw = curve.matrix_world
            pts = curve.data.splines[splineIdx].bezier_points
            nextIdx = getAdjIdx(curve, splineIdx, segIdx)

            seg = [mw @ pts[segIdx].co, mw @ pts[segIdx].handle_right, \
                mw @ pts[nextIdx].handle_left, mw @ pts[nextIdx].co]

            sortedCos = getCosSortedByT(seg, intersectPts, margin)
            sortedCos = removeDupliCos(sortedCos, margin)
            insertCos = sortedCos.copy()
            start = seg[0].freeze()
            end = seg[1].freeze()
            startIdxIncr = 1
            endIdxIncr = 1
            if(start in insertCos):
                startIdxIncr = 0 # Keep in split list
                insertCos.remove(start) # Remove from insert list
            if(end in insertCos):
                insertCos.remove(end) # Remove from insert list
                endIdxIncr = 2 # Keep in split list
            for j, co in enumerate(insertCos):
                insertCos[j] = mw.inverted_safe() @ co
            insertBezierPts(curve, splineIdx, segIdx, insertCos, 'FREE', margin)
            prevCurveIdx = curveIdx
            prevSplineIdx = splineIdx
            prevCnt += len(insertCos)

            if(action == 'CUT'):
                if(selPtMap.get(curve) == None):
                    selPtMap[curve] = {}
                if(selPtMap[curve].get(splineIdx) == None):
                    selPtMap[curve][splineIdx] = []
                selPtMap[curve][splineIdx] += \
                    list(range(segIdx + startIdxIncr, \
                        segIdx + len(sortedCos) + endIdxIncr))

        if(action == 'CUT'):
            for curve in list(selPtMap.keys()):
                splineIdxs = sorted(selPtMap[curve].keys())
                if(len(curve.data.splines) > 1):
                    newObjs, changeCnt = splitCurve([curve], 'spline', \
                        curve.users_collection)
                    splineIdxs = list(selPtMap[curve].keys())
                    for idx in splineIdxs:
                        newCurve = newObjs[idx]
                        selPtMap[newCurve] = {}
                        selPtMap[newCurve][0] = selPtMap[curve][idx]
                        selPtMap[curve].pop(idx)
                        if(len(selPtMap[curve]) == 0):
                            selPtMap.pop(curve)
            splitCurveSelPts(selPtMap)

def getSVGPt(co, docW, docH, camera = None, region = None, rv3d = None):
    if(camera != None):
        scene = bpy.context.scene
        xy = world_to_camera_view(scene, camera, co)
        return complex(xy[0] * docW, docH - (xy[1] * docH))
    elif(region != None and rv3d != None):
        xy = getCoordFromLoc(region, rv3d, co)
        return complex(xy[0], docH - xy[1])

def getPathD(path):
    curve = ''

    for i, part in enumerate(path):
        comps = []
        for j, segment in enumerate(part):
            if(j == 0):
                comps.append('M {},{} C'.format(segment[0].real, segment[0].imag))
            args = (segment[1].real, segment[1].imag,
                    segment[2].real, segment[2].imag,
                    segment[3].real, segment[3].imag)
            comps.append('{},{} {},{} {},{}'.format(*args))
        curve += ' ' .join(comps)

    return curve

def getPathBBox(path):
    minX, minY, maxX, maxY = [None, None, None, None]
    for part in path:
        for seg in part:
            seg3d = [(seg[i].real, seg[i].imag, 0) for i in range(len(seg))]
            leftBotFront, rgtTopBack = getBBox(seg3d)
            if(minX == None or leftBotFront[0] < minX):
                minX = leftBotFront[0]
            if(minY == None or leftBotFront[1] < minY):
                minY = leftBotFront[1]
            if(maxX == None or rgtTopBack[0] > maxX):
                maxX = rgtTopBack[0]
            if(maxY == None or rgtTopBack[1] > maxY):
                maxY = rgtTopBack[1]
    return minX, minY, maxX, maxY

def createClipElem(doc, svgElem, docW, docH, clipElemId):
    elem = doc.createElement('defs')
    svgElem.appendChild(elem)
    clipElem = doc.createElement('clipPath')
    clipElem.setAttribute('clipPathUnits', 'userSpaceOnUse')
    clipElem.setAttribute('id', clipElemId)
    elem.appendChild(clipElem)
    rectElem = doc.createElement('rect')
    rectElem.setAttribute('x', '0')
    rectElem.setAttribute('y', '0')
    rectElem.setAttribute('width', str(docW))
    rectElem.setAttribute('height', str(docH))
    clipElem.appendChild(rectElem)

def getSVGPathElem(doc, docW, docH, path, idx, lineWidth, lineCol, lineAlpha, \
    fillCol, fillAlpha, clipView, clipElemId):

    idPrefix = 'id'
    style= {'opacity':'1', 'stroke':'#000000', 'stroke-width':'1', \
        'fill':'none', 'stroke-linecap':'round', 'stroke-linejoin':'miter', \
        'stroke-miterlimit':'4'}

    clipped = False
    if(clipView):
        minX, minY, maxX, maxY = getPathBBox(path)
        if(maxX < 0 or maxY < 0 or minX > docW or minY > docH):
            return None

        if(minX < 0 or minY < 0 or maxX > docW or maxY > docH):
            clipped = True

    elem = doc.createElement('path')
    elem.setAttribute('id', idPrefix + str(idx).zfill(3))
    elem.setAttribute('d', getPathD(path))
    style['stroke-width'] =  str(lineWidth)
    style['stroke'] =  '#' + lineCol
    style['opacity'] =  lineAlpha
    if(fillCol != None):
        style['fill'] =  '#' + fillCol
        style['opacity'] =  fillAlpha # Overwrite
    styleStr = ';'.join([k + ':' + style[k] for k in style])
    elem.setAttribute('style', styleStr)

    if(clipped):
        elem.setAttribute('clip-path', 'url(#' + clipElemId + ')')

    return elem

def exportSVG(context, filepath, exportView, clipView, lineWidth, lineColorOpts, \
    lineColor, fillColorOpts, fillColor):

    svgXML = '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
    clipElemId = 'BBoxClipElem'

    if(lineColorOpts == 'PICK'):
        lineCol, lineAlpha = toHexStr(lineColor)

    if(fillColorOpts == 'PICK'):
        fillCol, fillAlpha = toHexStr(fillColor)

    if exportView == 'ACTIVE_VIEW':
        area = context.area
        if(area.type != 'VIEW_3D'):
            area = [a for a in bpy.context.screen.areas if  a.type == 'VIEW_3D'][0]
        region = [r for r in area.regions if r.type == 'WINDOW'][0]
        space3d = area.spaces[0]
        if(len(space3d.region_quadviews) > 0):
            rv3d = space3d.region_quadviews[3]
        else:
            rv3d = space3d.region_3d
        camera = None
        docW = region.width
        docH = region.height
    else:
        region = None
        rv3d = None
        camera = bpy.data.objects[exportView]
        docW = bpy.context.scene.render.resolution_x
        docH = bpy.context.scene.render.resolution_y

    doc = minidom.parseString(svgXML)
    svgElem = doc.documentElement

    svgElem.setAttribute('width', str(docW))
    svgElem.setAttribute('height', str(docH))

    if(clipView):
        createClipElem(doc, svgElem, docW, docH, clipElemId)

    idx = 0
    for o in bpy.context.scene.objects:
        mw = o.matrix_world
        if(isBezier(o) and o.visible_get()):
            path = []
            filledPath = []
            for spline in o.data.splines:
                part = []
                bpts = spline.bezier_points
                for i in range(1, len(bpts)):
                    prevBezierPt = bpts[i-1]
                    pt = bpts[i]
                    seg = [prevBezierPt.co, prevBezierPt.handle_right, pt.handle_left, pt.co]
                    part.append([getSVGPt(mw @ co, docW, docH, camera, region, rv3d) for co in seg])

                if(spline.use_cyclic_u):
                    seg = [bpts[-1].co, bpts[-1].handle_right, bpts[0].handle_left, bpts[0].co]
                    part.append([getSVGPt(mw @ co, docW, docH, camera, region, rv3d) for co in seg])

                if(len(part) > 0):
                    if (spline.use_cyclic_u and o.data.dimensions == '2D' \
                        and o.data.fill_mode != 'NONE'):
                        filledPath.append(part)
                    else:
                        path.append(part)

            for p in [path, filledPath]:

                if(len(p) == 0): continue

                if(lineColorOpts == 'RANDOM'):
                    lineColor = [random.random() for i in range(3)] + [1]
                    lineCol, lineAlpha = toHexStr(lineColor)

                if(p == path):
                    fc, fa = None, None
                elif(fillColorOpts == 'RANDOM'):
                    fillColor = [random.random() for i in range(3)] + [1]
                    fc, fa = toHexStr(fillColor)
                else:
                    fc, fa = fillCol, fillAlpha

                svgPathElem = getSVGPathElem(doc, docW, docH, p, idx, lineWidth, \
                    lineCol, lineAlpha, fc, fa, clipView, clipElemId)
                if(svgPathElem != None):
                    svgElem.appendChild(svgPathElem)
                    idx += 1

    doc.writexml(open(filepath,"w"))


###################### Operators ######################

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
        if(bpy.context.active_object != None):
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
        if(bpy.context.active_object != None):
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
            if(center != None and normal != None):
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
        if(src != None and isBezier(src)):
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

        if(actCurve != None and isBezier(actCurve)):
            curves = [actCurve] + curves

        if(len(curves) < 2):
            self.report({'INFO'}, "Please select at least two Bezier curve objects")
        else:
            intersectCurves(curves, params.intersectOp, params.intersectNonactive, \
            params.intersectMargin, rounding)

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


class MarkerController:
    drawHandlerRef = None
    defPointSize = 6
    ptColor = (0, .8, .8, 1)

    def createSMMap(self, context):
        objs = context.selected_objects
        smMap = {}
        for curve in objs:
            if(not isBezier(curve)):
                continue

            smMap[curve.name] = {}
            mw = curve.matrix_world
            for splineIdx, spline in enumerate(curve.data.splines):
                if(not spline.use_cyclic_u):
                    continue

                #initialize to the curr start vert co and idx
                smMap[curve.name][splineIdx] = \
                    [mw @ curve.data.splines[splineIdx].bezier_points[0].co, 0]

                for pt in spline.bezier_points:
                    pt.select_control_point = False

            if(len(smMap[curve.name]) == 0):
                del smMap[curve.name]

        return smMap

    def createBatch(self, context):
        positions = [s[0] for cn in self.smMap.values() for s in cn.values()]
        colors = [MarkerController.ptColor for i in range(0, len(positions))]

        self.batch = batch_for_shader(self.shader, \
            "POINTS", {"pos": positions, "color": colors})

        if context.area:
            context.area.tag_redraw()

    def drawHandler(self):
        bgl.glPointSize(MarkerController.defPointSize)
        self.batch.draw(self.shader)

    def removeMarkers(self, context):
        if(MarkerController.drawHandlerRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(MarkerController.drawHandlerRef, \
                "WINDOW")

            if(context.area and hasattr(context.space_data, 'region_3d')):
                context.area.tag_redraw()

            MarkerController.drawHandlerRef = None
        self.deselectAll()

    def __init__(self, context):
        self.smMap = self.createSMMap(context)
        self.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
        # ~ self.shader.bind()

        try:
            MarkerController.defPointSize = \
                context.preferences.addons[__name__].preferences.markerSize
        except Exception as e:
            # ~ print("BezierUtils: Fetching marker size", e)
            MarkerController.defPointSize = 6

        MarkerController.drawHandlerRef = \
            bpy.types.SpaceView3D.draw_handler_add(self.drawHandler, \
                (), "WINDOW", "POST_VIEW")

        self.createBatch(context)

    def saveStartVerts(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            spMap = self.smMap[curveName]

            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                if(markerInfo[1] != 0):
                    loc, idx = markerInfo[0], markerInfo[1]
                    moveSplineStart(curve, splineIdx, idx)

    def updateSMMap(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            spMap = self.smMap[curveName]
            mw = curve.matrix_world

            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                loc, idx = markerInfo[0], markerInfo[1]
                pts = curve.data.splines[splineIdx].bezier_points

                selIdxs = [x for x in range(0, len(pts)) \
                    if pts[x].select_control_point == True]

                selIdx = selIdxs[0] if(len(selIdxs) > 0 ) else idx
                co = mw @ pts[selIdx].co
                self.smMap[curveName][splineIdx] = [co, selIdx]

    def deselectAll(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            for spline in curve.data.splines:
                for pt in spline.bezier_points:
                    pt.select_control_point = False

    def getSpaces3D(context):
        areas3d  = [area for area in context.window.screen.areas \
            if area.type == 'VIEW_3D']

        return [s for a in areas3d for s in a.spaces if s.type == 'VIEW_3D']

    def hideHandles(context):
        states = []
        spaces = MarkerController.getSpaces3D(context)
        for s in spaces:
            if(hasattr(s.overlay, 'show_curve_handles')):
                states.append(s.overlay.show_curve_handles)
                s.overlay.show_curve_handles = False
            elif(hasattr(s.overlay, 'display_handle')): # 2.90
                states.append(s.overlay.display_handle)
                s.overlay.display_handle = 'NONE'
        return states

    def resetShowHandleState(context, handleStates):
        spaces = MarkerController.getSpaces3D(context)
        for i, s in enumerate(spaces):
            if(hasattr(s.overlay, 'show_curve_handles')):
                s.overlay.show_curve_handles = handleStates[i]
            elif(hasattr(s.overlay, 'display_handle')): # 2.90
                s.overlay.display_handle = handleStates[i]



class ModalMarkSegStartOp(bpy.types.Operator):
    bl_description = "Mark Vertex"
    bl_idname = "wm.bb_mark_vertex"
    bl_label = "Mark Start Vertex"

    def cleanup(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        self.markerState.removeMarkers(context)
        MarkerController.resetShowHandleState(context, self.handleStates)
        bpy.context.window_manager.bezierToolkitParams.markVertex = False

    def modal (self, context, event):

        if(context.mode  == 'OBJECT' or event.type == "ESC" or \
            not bpy.context.window_manager.bezierToolkitParams.markVertex):
            self.cleanup(context)
            return {'CANCELLED'}

        elif(event.type == "RET"):
            self.markerState.saveStartVerts()
            self.cleanup(context)
            return {'FINISHED'}

        if(event.type == 'TIMER'):
            self.markerState.updateSMMap()
            self.markerState.createBatch(context)

        return {"PASS_THROUGH"}

    def execute(self, context):
        #TODO: Why such small step?
        self._timer = context.window_manager.event_timer_add(time_step = 0.0001, \
            window = context.window)

        context.window_manager.modal_handler_add(self)
        self.markerState = MarkerController(context)

        #Hide so that users don't accidentally select handles instead of points
        self.handleStates = MarkerController.hideHandles(context)

        return {"RUNNING_MODAL"}

###################### Single Panel for All Ops ######################

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

            row = layout.row()
            row.prop(params, "intersectExpanded",
                icon="TRIA_DOWN" if params.intersectExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Intersect Curves', icon='GRAPH')
            if params.intersectExpanded:
                box = layout.box()
                col = box.column().split()
                row = col.row()
                row.prop(params, 'intersectMargin', text = 'Proximity')
                col = box.column().split()
                row = col.row()
                row.prop(params, 'intersectOp', text = 'Action')
                if(params.intersectOp in {'INSERT_PT', 'CUT'}):
                    row = col.row()
                    row.prop(params, 'intersectNonactive', text = 'Only Non-active')
                col = box.column().split()
                col.operator('object.intersect_curves')

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
                if(params.fillType == 'QUAD'):
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'remeshDepth')
                    row.prop(params, 'unsubdivide')
                else:
                    col = box.column().split()
                    row = col.row()
                    row.prop(params, 'remeshRes')
                    row.prop(params, 'remeshApplyTo')
                    if(params.remeshApplyTo == 'PERSEG'):
                        col = box.column().split()
                        row = col.row()
                        row.prop(params, 'remeshOptimized')
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
                bgl.glLineWidth(BezierUtilsPanel.lineWidth)
                if(BezierUtilsPanel.lineBatch != None):
                    BezierUtilsPanel.lineBatch.draw(BezierUtilsPanel.shader)

        if(add and BezierUtilsPanel.drawHandlerRef == None):
            try:
                BezierUtilsPanel.lineWidth = \
                    bpy.context.preferences.addons[__name__].preferences.lineWidth
            except Exception as e:
                print("BezierUtils: Error fetching line width in ColorCurves: ", e)
                BezierUtilsPanel.lineWidth = 1.5

            BezierUtilsPanel.drawHandlerRef = \
                bpy.types.SpaceView3D.draw_handler_add(ccDrawHandler, \
                    (), "WINDOW", "POST_VIEW")
            BezierUtilsPanel.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
            return

        elif(remove):
            if(BezierUtilsPanel.drawHandlerRef != None):
                bpy.types.SpaceView3D.draw_handler_remove(BezierUtilsPanel.drawHandlerRef, \
                    "WINDOW")
                BezierUtilsPanel.drawHandlerRef = None
                return

        if(bpy.context.screen == None):
            return

        if(bpy.context.window_manager.bezierToolkitParams.applyCurveColor):
            objs = [o for o in bpy.context.scene.objects if(isBezier(o) and \
                o.visible_get() and len(o.modifiers) == 0 and not o.select_get())]

            lineCos = []
            lineColors = []
            for o in objs:
                colorVal = o.data.get('curveColor')
                if(colorVal != None):
                    for i, spline in enumerate(o.data.splines):
                        for j in range(0, len(spline.bezier_points)):
                            segPts = getBezierDataForSeg(o, i, j, withShapeKey = True, \
                                shapeKeyIdx = None, fromMix = True)
                            if(segPts == None):
                                continue
                            pts = getPtsAlongBezier2D(segPts, getAllAreaRegions(), \
                                FTProps.dispCurveRes, maxRes = MAX_NONSEL_CURVE_RES)
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

def getPtFromT(p0, p1, p2, p3, t):
    c = (1 - t)
    pt = (c ** 3) * p0 + 3 * (c ** 2) * t * p1 + \
        3 * c * (t ** 2) * p2 + (t ** 3) * p3
    return pt

def getTangentAtT(p0, p1, p2, p3, t):
    c = (1 - t)
    tangent = -3 * (c * c) * p0 + 3 * c * c * p1 - 6 * t * c * p1 - \
        3 * t * t * p2 + 6 * t * c * p2 + 3 * t * t * p3
    return tangent

# iterative brute force, not optimized, some iterations maybe redundant
def getTsForPt(p0, p1, p2, p3, co, coIdx, tolerance = 0.000001, maxItr = 1000):
    ts = set()
    # check t from start to end and end to start
    for T in [1., 0.]:
        # check clockwise as well as anticlockwise
        for dirn in [1, -1]:
            t = T
            t2 = 1
            rhs = getPtFromT(p0, p1, p2, p3, t)[coIdx]
            error = rhs - co
            i = 0

            while(abs(error) > tolerance and i < maxItr):
                t2 /= 2
                if(dirn * error < 0):
                    t += t2
                else:
                    t -= t2
                rhs = getPtFromT(p0, p1, p2, p3, t)[coIdx]
                error = rhs - co

                i += 1

            if(i < maxItr and t >= 0 and t <= 1):
                ts.add(round(t, 3))
    return ts

#TODO: There may be a more efficient approach, but this seems foolproof
def getTForPt(curve, testPt, tolerance = .000001):
    minLen = LARGE_NO
    retT = None
    for coIdx in range(0, 3):
        ts = getTsForPt(curve[0], curve[1], curve[2], curve[3], \
            testPt[coIdx], coIdx, tolerance)
        for t in ts:
            pt = getPtFromT(curve[0], curve[1], curve[2], curve[3], t)
            pLen = (testPt - pt).length
            if(pLen < minLen):
                minLen = pLen
                retT = t
    return retT

def getCosSortedByT(seg, cos, margin):
    coInfo = set()

    for co in cos:
        t = getTForPt(seg, co, margin)
        # ~ if(all(abs(co[i] - seg[3][i]) < margin for i in range(3)) or t >= 1):
        if(t >= 1):
            coInfo.add(((seg[3]).freeze(), 1))
        # ~ elif(all(abs(co[i] - seg[0][i]) < margin for i in range(3)) or t <= 0):
        elif(t <= 0):
            coInfo.add(((seg[0]).freeze(), 0))
        else:
            coInfo.add((co.freeze(), t))

    return [inf[0] for inf in sorted(coInfo, key = lambda x: x[1])]

# TODO: Check for initial and end points of the segment (here or in getCosSortedByT)
# Currently inserting points more than once if intersect is invoked repeatedly
def removeDupliCos(sortedCos, margin):
    prevCo = sortedCos[0]
    newCos = [prevCo]
    for i in range(1, len(sortedCos)):
        co = sortedCos[i]
        if(not vectCmpWithMargin(co, prevCo, margin)):
            newCos.append(co)
        prevCo = co
    return newCos

# https://stackoverflow.com/questions/24809978/calculating-the-bounding-box-of-cubic-bezier-curve
#(3 D - 9 C + 9 B - 3 A) t^2 + (6 A - 12 B + 6 C) t + 3 (B - A)
def getBBox(seg):
    A = seg[0]
    B = seg[1]
    C = seg[2]
    D = seg[3]

    leftBotFront = Vector([min([A[i], D[i]]) for i in range(0, 3)])
    rgtTopBack = Vector([max([A[i], D[i]]) for i in range(0, 3)])

    a = [3 * D[i] - 9 * C[i] + 9 * B[i] - 3 * A[i] for i in range(0, 3)]
    b = [6 * A[i] - 12 * B[i] + 6 * C[i] for i in range(0, 3)]
    c = [3 * (B[i] - A[i]) for i in range(0, 3)]

    solnsxyz = []
    for i in range(0, 3):
        solns = []
        if(a[i] == 0):
            if(b[i] == 0):
                solns.append(0)#Independent of t so lets take the starting pt
            else:
                solns.append(c[i] / b[i])
        else:
            rootFact = b[i] * b[i] - 4 * a[i] * c[i]
            if(rootFact >=0 ):
                #Two solutions with + and - sqrt
                solns.append((-b[i] + sqrt(rootFact)) / (2 * a[i]))
                solns.append((-b[i] - sqrt(rootFact)) / (2 * a[i]))
        solnsxyz.append(solns)

    for i, soln in enumerate(solnsxyz):
        for j, t in enumerate(soln):
            if(t <= 1 and t >= 0):
                co = getPtFromT(A[i], B[i], C[i], D[i], t)
                if(co < leftBotFront[i]): leftBotFront[i] = co
                if(co > rgtTopBack[i]): rgtTopBack[i] = co

    return leftBotFront, rgtTopBack

def getBBoxOverlapInfo(seg0, seg1):
    bbox0 = getBBox(seg0)
    max0 = [max(bbox0[i][axis] for i in range(2)) for axis in range(3)]
    min0 = [min(bbox0[i][axis] for i in range(2)) for axis in range(3)]
    bbox1 = getBBox(seg1)

    overlap = True
    if(any(all(bbox1[i][j] < min0[j] for i in range(2)) for j in range(3)) or \
        any(all(bbox1[i][j] > max0[j] for i in range(2)) for j in range(3))):
        overlap = False

    return overlap, bbox0, bbox1

def getBBoxCenter(bbox): # bbox -> [leftBotFront, rightTopBack]
    return Vector(((bbox[0][i] + bbox[1][i]) / 2 for i in range(3)))

def getIntersectPts(seg0, seg1, soln, solnRounded, recurs, margin, rounding, \
    maxRecurs = 100):
    overlap, bbox0, bbox1 = getBBoxOverlapInfo(seg0, seg1)
    if(overlap):
        if(recurs == maxRecurs):
            print('Maximum recursions in getIntersectPts!')
            return False

        if(all(abs(bbox0[0][i] - bbox0[1][i]) < margin for i in range(3))):
            center = getBBoxCenter(bbox0)
            roundedVect = Vector([round(x, rounding) for x in center]).freeze()
            if(roundedVect not in solnRounded):
                soln.append(center)
                solnRounded.add(roundedVect)
            return True
        elif(all(abs(bbox1[0][i] - bbox1[1][i]) < margin for i in range(3))):
            center = getBBoxCenter(bbox1)
            roundedVect = Vector([round(x, rounding) for x in center]).freeze()
            if(roundedVect not in solnRounded):
                soln.append(center)
                solnRounded.add(roundedVect)
            return True
        else:
            seg01 = getPartialSeg(seg0, t0 = 0, t1 = 0.5)
            seg02 = getPartialSeg(seg0, t0 = 0.5, t1 = 1)

            seg11 = getPartialSeg(seg1, t0 = 0, t1 = 0.5)
            seg12 = getPartialSeg(seg1, t0 = 0.5, t1 = 1)

            r0 = getIntersectPts(seg01, seg11, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r1 = getIntersectPts(seg01, seg12, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r2 = getIntersectPts(seg02, seg11, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r3 = getIntersectPts(seg02, seg12, soln, solnRounded, \
                recurs + 1, margin, rounding)

            return any((r0, r1, r2, r3))
    return False

# splineInfos: [(curveIdx0, splineIdx0), (curveIdx0, splineIdx1),...]
# First active means intersections only with the first curve (otherwise all combinations)
def getSplineIntersectPts(curves, splineInfos, firstActive, margin, rounding):
    segPairMap = {}
    areas = [a for a in bpy.context.screen.areas if  a.type == 'VIEW_3D']
    for i, c0Info in enumerate(splineInfos):
        idxCurve0, idxSpline0 = c0Info
        if(firstActive and idxCurve0 != splineInfos[0][0]):
            break
        c0 = curves[idxCurve0]
        spline0 = c0.data.splines[idxSpline0]
        for j in range(i + 1, len(splineInfos)):
            idxCurve1, idxSpline1 = splineInfos[j]
            c1 = curves[idxCurve1]
            spline1 = c1.data.splines[idxSpline1]

            segPairMap[((idxCurve0, idxSpline0), (idxCurve1, idxSpline1))] = \
                [((idxSeg0, s0), (idxSeg1, s1)) for idxSeg0, s0 in \
                    enumerate(getRoundedSplineSegs(c0.matrix_world, spline0)) for \
                        idxSeg1, s1 in enumerate(getRoundedSplineSegs(c1.matrix_world, \
                            spline1))]

    intersectMap = {}
    allIntersectCos = []
    for key in segPairMap:
        (idxCurve0, idxSpline0), (idxCurve1, idxSpline1) = key
        segPairInfo = segPairMap[key]
        for info in segPairInfo:
            solnRounded = set()
            soln = []
            (idxSeg0, seg0), (idxSeg1, seg1) = info
            ret = getIntersectPts(seg0, seg1, soln, solnRounded, recurs = 0, \
                margin = margin, rounding = rounding)
            if(ret):
                extKey = (idxCurve0, idxSpline0, idxSeg0)
                if(intersectMap.get(extKey) == None):
                    intersectMap[extKey] = []
                intersectMap[extKey] += soln

                extKey = (idxCurve1, idxSpline1, idxSeg1)
                if(intersectMap.get(extKey) == None):
                    intersectMap[extKey] = []
                intersectMap[extKey] += soln

                allIntersectCos += soln
    return allIntersectCos, intersectMap

def getCurveIntersectPts(curves, firstActive, margin, rounding):
    splineInfos = [(x, y) for x in range(len(curves)) \
        for y in range(len(curves[x].data.splines))]

    return getSplineIntersectPts(curves, splineInfos, firstActive, margin, rounding)

def getRoundedSplineSegs(mw, spline, reverse = False, rounding = 5):

    def getRoundedVect(mw, co, rounding):
        return Vector((round(x, rounding) for x in (mw @ co)))

    bpts = spline.bezier_points
    if(reverse): bpts = reversed(bpts)
    segs = []
    for i in range(1, len(bpts)):
        segPts = [bpts[i-1].co, bpts[i-1].handle_right, bpts[i].handle_left, bpts[i].co]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    if(spline.use_cyclic_u):
        segPts = [bpts[-1].co, bpts[-1].handle_right, bpts[0].handle_left, bpts[0].co]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    return segs

# Because there is some discrepancy between this and getSegLen
# This seems to be more accurate
def getSegLenTmpObj(tmpSpline, bpts, mw = Matrix()):
    tmpSpline.bezier_points[0].co = mw @ bpts[0].co
    tmpSpline.bezier_points[0].handle_right = mw @ bpts[0].handle_right
    tmpSpline.bezier_points[1].handle_left = mw @ bpts[1].handle_left
    tmpSpline.bezier_points[1].co = mw @ bpts[1].co
    return tmpSpline.calc_length()

def getSplineLenTmpObj(tmpSpline, spline, mw):
    l = 0
    bpts = spline.bezier_points
    l += sum(getSegLenTmpObj(tmpSpline, bpts[i:i+2], mw) for i in range(len(bpts) -1))
    if(spline.use_cyclic_u): l += getSegLenTmpObj(tmpSpline, [bpts[-1], bpts[0]], mw)
    return l

def getSegLen(pts, error = DEF_ERR_MARGIN, start = None, end = None, t1 = 0, t2 = 1):
    if(start == None): start = pts[0]
    if(end == None): end = pts[-1]

    t1_5 = (t1 + t2)/2
    mid = getPtFromT(*pts, t1_5)
    l = (end - start).length
    l2 = (mid - start).length + (end - mid).length
    if (l2 - l > error):
        return (getSegLen(pts, error, start, mid, t1, t1_5) +
                getSegLen(pts, error, mid, end, t1_5, t2))
    return l2

def hasAlignedHandles(pt):
    if(len(pt) == 5 and 'ALIGNED' in {pt[3], pt[4]} and 'FREE' not in {pt[3], pt[4]}):
        return True
    diffV1 = pt[1] - pt[0]
    diffV2 = pt[2] - pt[1]
    if(vectCmpWithMargin(diffV1.normalized(), diffV2.normalized())):
        return True
    return False

def isStraightSeg(segPts):
    if(len(segPts) != 2): return False
    if((len(segPts[0]) == 5 or len(segPts[1]) == 5) and \
        segPts[0][4] == 'VECTOR' and segPts[1][3] == 'VECTOR'):
            return True
    if(vectCmpWithMargin((segPts[0][2]-segPts[0][1]).normalized(), \
        (segPts[1][1] - segPts[1][0]).normalized())):
        return True
    return False

# Get pt coords along curve defined by the four control pts (segPts)
# subdivPerUnit: No of subdivisions per unit length
# (which is the same as no of pts excluding the end pts)
def getInterpBezierPts(segPts, subdivPerUnit, segLens = None, maxRes = None):
    if(len(segPts) < 2):
        return []

    curvePts = []
    for i in range(1, len(segPts)):
        seg = [segPts[i-1][1], segPts[i-1][2], segPts[i][0], segPts[i][1]]
        if(segLens != None and len(segLens) > (i-1)):
            res = int(segLens[i-1] * subdivPerUnit)
        else:
            res = int(getSegLen(seg) * subdivPerUnit)
        if(res < 2): res = 2
        if(maxRes != None and res > maxRes): res = maxRes
        curvePts += geometry.interpolate_bezier(*seg, res)

    return curvePts

# Used in functions where actual locs of pts on curve matter (like subdiv Bezier)
# (... kind of expensive)
def getPtsAlongBezier3D(segPts, rv3d, curveRes, minRes = 200):

    viewDist = rv3d.view_distance

    # The smaller the view dist (higher zoom level),
    # the higher the num of subdivisions
    curveRes = curveRes / viewDist

    if(curveRes < minRes): curveRes = minRes

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes)

# Used in functions where only visual resolution of curve matters (like draw Bezier)
# (... not so expensive)
# TODO: Calculate maxRes dynamically
def getPtsAlongBezier2D(segPts, areaRegionInfo, curveRes, maxRes = None):
    segLens = []
    for i in range(1, len(segPts)):
        seg = [segPts[i-1][1], segPts[i-1][2], segPts[i][0], segPts[i][1]]

        #TODO: A more optimized solution... (Called very frequently)
        segLen = 0
        for info in areaRegionInfo:
            seg2D = [getCoordFromLoc(info[1], info[2], loc) for loc in seg]
            sl = getSegLen(seg2D)
            if(sl > segLen):
                segLen = sl
        segLens.append(segLen)

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes, \
        segLens = segLens, maxRes = maxRes)

def getLinesFromPts(pts):
    positions = []
    for i, pt in enumerate(pts):
        positions.append(pt)
        if(i > 0 and i < (len(pts)-1)):
            positions.append(pt)
    return positions

#see https://stackoverflow.com/questions/878862/drawing-part-of-a-b%c3%a9zier-curve-by-reusing-a-basic-b%c3%a9zier-curve-function/879213#879213
def getPartialSeg(seg, t0, t1):
    pts = [seg[0], seg[1], seg[2], seg[3]]

    if(t0 > t1):
        tt = t1
        t1 = t0
        t0 = tt

    u0 = 1.0 - t0
    u1 = 1.0 - t1

    qa = [pts[0][i]*u0*u0 + pts[1][i]*2*t0*u0 + pts[2][i]*t0*t0 for i in range(0, 3)]
    qb = [pts[0][i]*u1*u1 + pts[1][i]*2*t1*u1 + pts[2][i]*t1*t1 for i in range(0, 3)]
    qc = [pts[1][i]*u0*u0 + pts[2][i]*2*t0*u0 + pts[3][i]*t0*t0 for i in range(0, 3)]
    qd = [pts[1][i]*u1*u1 + pts[2][i]*2*t1*u1 + pts[3][i]*t1*t1 for i in range(0, 3)]

    pta = Vector([qa[i]*u0 + qc[i]*t0 for i in range(0, 3)])
    ptb = Vector([qa[i]*u1 + qc[i]*t1 for i in range(0, 3)])
    ptc = Vector([qb[i]*u0 + qd[i]*t0 for i in range(0, 3)])
    ptd = Vector([qb[i]*u1 + qd[i]*t1 for i in range(0, 3)])

    return [pta, ptb, ptc, ptd]

def getInterpolatedVertsCo(curvePts, numDivs):
    # Can be calculated only once
    curveLength = sum((curvePts[i] - curvePts[i-1]).length
        for i in range(1, len(curvePts)))

    if(floatCmpWithMargin(curveLength, 0)):
        return [curvePts[0]] * numDivs

    segLen = curveLength / numDivs
    vertCos = [curvePts[0]]

    actualLen = 0
    vertIdx = 0

    for i in range(1, numDivs):
        co = None
        targetLen = i * segLen

        while(not floatCmpWithMargin(actualLen, targetLen)
            and actualLen < targetLen):

            vertCo = curvePts[vertIdx]
            vertIdx += 1
            nextVertCo = curvePts[vertIdx]
            actualLen += (nextVertCo - vertCo).length

        if(floatCmpWithMargin(actualLen, targetLen)):
            co = curvePts[vertIdx]

        else:   #interpolate
            diff = actualLen - targetLen
            co = (nextVertCo - (nextVertCo - vertCo) * \
                (diff/(nextVertCo - vertCo).length))

            #Revert to last pt
            vertIdx -= 1
            actualLen -= (nextVertCo - vertCo).length
        vertCos.append(co)

    # ~ if(not vectCmpWithMargin(curvePts[0], curvePts[-1])):
    vertCos.append(curvePts[-1])

    return vertCos

#
# The following section is a Python conversion of the javascript
# a2c function at: https://github.com/fontello/svgpath
# (Copyright (C) 2013-2015 by Vitaly Puzrin)
#
# Note: Most of the comments are retained
######################## a2c start #######################

TAU = pi * 2

# eslint-disable space-infix-ops

# Calculate an angle between two unit vectors
#
# Since we measure angle between radii of circular arcs,
# we can use simplified math (without length normalization)
#
def unit_vector_angle(ux, uy, vx, vy):
    if(ux * vy - uy * vx < 0):
        sign = -1
    else:
        sign = 1

    dot  = ux * vx + uy * vy

    # Add this to work with arbitrary vectors:
    # dot /= sqrt(ux * ux + uy * uy) * sqrt(vx * vx + vy * vy)

    # rounding errors, e.g. -1.0000000000000002 can screw up this
    if (round(dot, 3) >=  1.0):
        dot =  1.0

    if (round(dot, 3) <= -1.0):
        dot = -1.0

    return sign * acos(dot)


# Convert from endpoint to center parameterization,
# see http:#www.w3.org/TR/SVG11/implnote.html#ArcImplementationNotes
#
# Return [cx, cy, theta1, delta_theta]
#
def get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi):
    # Step 1.
    #
    # Moving an ellipse so origin will be the middlepoint between our two
    # points. After that, rotate it to line up ellipse axes with coordinate
    # axes.
    #
    x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2
    y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2

    rx_sq  =  rx * rx
    ry_sq  =  ry * ry
    x1p_sq = x1p * x1p
    y1p_sq = y1p * y1p

    # Step 2.
    #
    # Compute coordinates of the centre of this ellipse (cx', cy')
    # in the new coordinate system.
    #
    radicant = (rx_sq * ry_sq) - (rx_sq * y1p_sq) - (ry_sq * x1p_sq)

    if (radicant < 0):
        # due to rounding errors it might be e.g. -1.3877787807814457e-17
        radicant = 0

    radicant /=   (rx_sq * y1p_sq) + (ry_sq * x1p_sq)
    factor = 1
    if(fa == fs):# Migration Note: note ===
        factor = -1
    radicant = sqrt(radicant) * factor #(fa === fs ? -1 : 1)

    cxp = radicant *  rx/ry * y1p
    cyp = radicant * -ry/rx * x1p

    # Step 3.
    #
    # Transform back to get centre coordinates (cx, cy) in the original
    # coordinate system.
    #
    cx = cos_phi*cxp - sin_phi*cyp + (x1+x2)/2
    cy = sin_phi*cxp + cos_phi*cyp + (y1+y2)/2

    # Step 4.
    #
    # Compute angles (theta1, delta_theta).
    #
    v1x =  (x1p - cxp) / rx
    v1y =  (y1p - cyp) / ry
    v2x = (-x1p - cxp) / rx
    v2y = (-y1p - cyp) / ry

    theta1 = unit_vector_angle(1, 0, v1x, v1y)
    delta_theta = unit_vector_angle(v1x, v1y, v2x, v2y)

    if (fs == 0 and delta_theta > 0):#Migration Note: note ===
        delta_theta -= TAU

    if (fs == 1 and delta_theta < 0):#Migration Note: note ===
        delta_theta += TAU

    return [ cx, cy, theta1, delta_theta ]

#
# Approximate one unit arc segment with bezier curves,
# see http:#math.stackexchange.com/questions/873224
#
def approximate_unit_arc(theta1, delta_theta):
    alpha = 4.0/3 * tan(delta_theta/4)

    x1 = cos(theta1)
    y1 = sin(theta1)
    x2 = cos(theta1 + delta_theta)
    y2 = sin(theta1 + delta_theta)

    return [ x1, y1, x1 - y1*alpha, y1 + x1*alpha, x2 + y2*alpha, y2 - x2*alpha, x2, y2 ]

def a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi, noSegs):
    sin_phi = sin(phi * TAU / 360)
    cos_phi = cos(phi * TAU / 360)

    # Make sure radii are valid
    #
    x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2
    y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2

    if (x1p == 0 and y1p == 0): # Migration Note: note ===
        # we're asked to draw line to itself
        return []

    if (rx == 0 or ry == 0): # Migration Note: note ===
        # one of the radii is zero
        return []

    # Compensate out-of-range radii
    #
    rx = abs(rx)
    ry = abs(ry)

    lmbd = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry)
    if (lmbd > 1):
        rx *= sqrt(lmbd)
        ry *= sqrt(lmbd)


    # Get center parameters (cx, cy, theta1, delta_theta)
    #
    cc = get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi)

    result = []
    theta1 = cc[2]
    delta_theta = cc[3]

    # Split an arc to multiple segments, so each segment
    # will be less than 90
    #
    segments = noSegs # int(max(ceil(abs(delta_theta) / (TAU / 4)), 1))
    delta_theta /= segments

    for i in range(0, segments):
        result.append(approximate_unit_arc(theta1, delta_theta))

        theta1 += delta_theta

    # We have a bezier approximation of a unit circle,
    # now need to transform back to the original ellipse
    #
    return getMappedList(result, rx, ry, sin_phi, cos_phi, cc)

def getMappedList(result, rx, ry, sin_phi, cos_phi, cc):
    mappedList = []
    for elem in result:
        curve = []
        for i in range(0, len(elem), 2):
            x = elem[i + 0]
            y = elem[i + 1]

            # scale
            x *= rx
            y *= ry

            # rotate
            xp = cos_phi*x - sin_phi*y
            yp = sin_phi*x + cos_phi*y

            # translate
            elem[i + 0] = xp + cc[0]
            elem[i + 1] = yp + cc[1]
            curve.append(complex(elem[i + 0], elem[i + 1]))
        mappedList.append(curve)
    return mappedList

######################## a2c end #######################

def get3DVector(cmplx, axisIdxs, z):
    vElems = [None] * 3
    vElems[axisIdxs[0]] = cmplx.real
    vElems[axisIdxs[1]] = cmplx.imag
    vElems[axisIdxs[2]] = z
    return Vector(vElems)

def getSegsForArc(start, radius, sweep, end, noSegs, axisIdxs, z):
    x1, y1 = start.real, start.imag
    x2, y2 = end.real, end.imag
    fa = 0
    fs = sweep
    rx, ry = radius.real, radius.imag
    phi = 0
    curvesPts = a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi, noSegs)
    newSegs = []
    for curvePts in curvesPts:
        newSegs.append([get3DVector(curvePts[0], axisIdxs, z), get3DVector(curvePts[1], axisIdxs, z), \
            get3DVector(curvePts[2], axisIdxs, z), get3DVector(curvePts[3], axisIdxs, z)])

    return newSegs

def getWSDataForSegs(segs):

    prevSeg = None
    wsData = []

    for j, seg in enumerate(segs):

        pt = seg[0]
        handleRight = seg[1]

        if(j == 0): handleLeft = pt
        else: handleLeft = prevSeg[2]

        ht = 'ALIGNED' if(vectCmpWithMargin(pt - handleLeft, handleRight - pt)) else 'FREE'
        wsData.append([handleLeft, pt, handleRight, ht, ht])
        prevSeg = seg

    if(prevSeg != None): wsData.append([seg[2], seg[3], seg[3], 'FREE', 'FREE'])
    else: return []

    return wsData


################### Common to Draw and Edit Flexi Bezier Ops ###################

# Some global constants

SEARCH_CURVE_RES = .5 # Per pixel seg divisions (.5 is one div per 2 pixel units)
DBL_CLK_DURN = 0.25
SNGL_CLK_DURN = 0.3

MAX_SEL_CURVE_RES = 1000
MAX_NONSEL_CURVE_RES = 100

EVT_NOT_CONS = 0
EVT_CONS = 1
EVT_META_OR_SNAP = 2

TOOL_TYPE_FLEXI_DRAW = 'Flexi Draw'
TOOL_TYPE_FLEXI_GREASE = 'Flexi Grease'
TOOL_TYPE_FLEXI_EDIT = 'Flexi Edit'
TOOL_TYPES_FLEXI_DRAW_COMMON = {TOOL_TYPE_FLEXI_DRAW, TOOL_TYPE_FLEXI_GREASE}
TOOL_TYPES_FLEXI_ALL = {TOOL_TYPE_FLEXI_DRAW, \
    TOOL_TYPE_FLEXI_GREASE, TOOL_TYPE_FLEXI_EDIT}

class BptDisplayInfo:
    # handleNos: 0: seg1-left, 1: seg1-right
    # tipColors: leftHdl, pt, rightHdl
    # Caller to make sure there are no tips without handle
    def __init__(self, pt, tipColors, handleNos = None):
        self.pt = pt
        if(len(tipColors) == 1):
            self.tipColors = [None, tipColors[0], None]
        elif(len(tipColors) == 3):
            self.tipColors = tipColors
        if(handleNos == None): self.handleNos = []
        else: self.handleNos = handleNos


class SegDisplayInfo:
    def __init__(self, segPts, segColor):
        self.segPts = segPts
        self.segColor = segColor

class RegionMouseXYInfo:

    def getRegionMouseXYInfo(event, exclInRgns):
        xyScreen = [event.mouse_x, event.mouse_y]
        idxs = getAreaRegionIdxs(xyScreen, exclInRgns)
        if(idxs == None):
            return None
        else:
            i, j = idxs
            area = bpy.context.screen.areas[i]
            region = area.regions[j]
            space3d = area.spaces[0]
            if(len(space3d.region_quadviews) > 0):
                qIdx = getWindowRegionIdx(area, j)
                rv3d = space3d.region_quadviews[qIdx]
            else:
                rv3d = space3d.region_3d
            xy = [xyScreen[0] - region.x, xyScreen[1] - region.y]
            return RegionMouseXYInfo(area, space3d, region, rv3d, xy, xyScreen)

    def __init__(self, area, space3d, region, rv3d, xy, xyScreen):
        self.area = area
        self.space3d = space3d
        self.region = region
        self.rv3d = rv3d
        self.xy = xy
        self.xyScreen = xyScreen

    def __eq__(self, other):
        if(other == None): return False
        return self.area == other.area and self.region == other.region and \
            self.rv3d == other.rv3d

def getLineShades(lineCos, baseColor, start, end, mid = True):
    if(len(lineCos) == 0 ): return [], []
    if(len(lineCos) == 1 ): return lineCos[0], [baseColor]
    if(mid): midPt = lineCos[0] + (lineCos[1] - lineCos[0]) / 2
    col1 = [start * c for c in baseColor]
    col2 = [end * c for c in baseColor]
    if(mid): return [lineCos[0], midPt, midPt, lineCos[1]], [col1, col2, col2, col1]
    else: return [lineCos[0], lineCos[1]], [col1, col2]


class BGLDrawInfo:
    def __init__(self, size, color, pts):
        self.size = size
        self.color = color
        self.pts = pts


class BGLDrawInfoLine(BGLDrawInfo):
    def __init__(self, size, color, pts, \
        gradientStart = None, gradientEnd = None, mid = True):
        super(BGLDrawInfoLine, self).__init__(size, color, pts)
        self.gradientStart = gradientStart
        self.gradientEnd = gradientEnd
        self.mid = mid


class BGLDrawMgr:
    def __init__(self, shader):
        self.lineInfoMap = {}
        self.ptInfoMap = {}
        self.shader = shader

    def addLineInfo(self, infoId, size, color, pts, \
        gradientStart = None, gradientEnd = None, mid = True):
        self.lineInfoMap[infoId] = BGLDrawInfoLine(size, color, pts, \
            gradientStart, gradientEnd, mid)

    def addPtInfo(self, infoId, size, color, pts):
        self.ptInfoMap[infoId] = BGLDrawInfo(size, color, pts)

    def redraw(self):
        lineInfos = sorted(self.lineInfoMap.values(), key = lambda x: (x.size))
        pos = []
        col = []
        batches = []
        for i, info in enumerate(lineInfos):
            if(i == 0 or info.size != lineInfos[i-1].size):
                if(i > 0):
                    bgl.glLineWidth(lineInfos[i-1].size)
                    batch = batch_for_shader(self.shader, \
                        'LINES', {"pos": pos, "color": col})
                    batch.draw(self.shader)
                pos = []
                col = []

            if(info.gradientEnd != None and info.gradientStart != None):
                if(len(info.pts) != 2 and len(info.color) != 1):
                    raise ValueError('Exactly two ' + \
                        'coordinates and one color for gradient line')
                linePos, lineCols = getLineShades(info.pts, info.color[0], \
                    info.gradientStart, info.gradientEnd, info.mid)
            else:
                if(len(info.pts) == 0): continue
                linePos = info.pts[:]
                lineCols = info.color[:]
                diff = len(linePos) - len(lineCols)
                if(diff >= 0):
                    for j in range(diff):
                        lineCols.append(info.color[-1])
                else:
                    for j in range(-diff):
                        lineCols.pop()
            pos += linePos
            col += lineCols

        if(len(pos) > 0):
            bgl.glLineWidth(lineInfos[-1].size)
            batch = batch_for_shader(self.shader, \
                'LINES', {"pos": pos, "color": col})
            batch.draw(self.shader)

        ptInfos = sorted(self.ptInfoMap.values(), key = lambda x: (x.size))
        for i, info in enumerate(ptInfos):
            if(i == 0 or info.size != ptInfos[i-1].size):
                if(i > 0):
                    bgl.glPointSize(ptInfos[i-1].size)
                    batch = batch_for_shader(self.shader, \
                        'POINTS', {"pos": pos, "color": col})
                    batch.draw(self.shader)
                pos = []
                col = []

            if(len(info.pts) == 0): continue
            ptCols = info.color[:]
            diff = len(info.pts) - len(info.color)
            if(diff >= 0):
                for j in range(diff):
                    ptCols.append(info.color[-1])
            else:
                for j in range(-diff):
                    ptCols.pop()

            pos += info.pts[:]
            col += ptCols

        if(len(pos) > 0):
            bgl.glPointSize(ptInfos[-1].size)
            batch = batch_for_shader(self.shader, \
                'POINTS', {"pos": pos, "color": col})
            batch.draw(self.shader)

    def resetLineInfo(self, infoId):
        drawInfo = self.lineInfoMap.get(infoId)
        if(drawInfo != None):
            drawInfo.pts = []

    def resetPtInfo(self, infoId):
        drawInfo = self.ptInfoMap.get(infoId)
        if(drawInfo != None):
            drawInfo.pts = []

    def reset(self):
        for key in list(self.ptInfoMap.keys()):
            self.ptInfoMap[key].pts = []
        for key in list(self.lineInfoMap.keys()):
            self.lineInfoMap[key].pts = []

# Return line batch for bezier line segments and handles and point batch for handle tips
def updateBezierBatches(bglDrawMgr, segDispInfos, bptDispInfos, areaRegionInfo, \
        defHdlType = 'ALIGNED'):

    lineCos = [] #segment is also made up of lines
    lineColors = []
    for i, info in enumerate(segDispInfos):
        segPts = info.segPts
        if(isStraightSeg(segPts)):
            lineCos += [segPts[0][1], segPts[1][1]]
            lineColors += [info.segColor, info.segColor]
        else:
            pts = getPtsAlongBezier2D(segPts, areaRegionInfo, \
                FTProps.dispCurveRes, maxRes = MAX_NONSEL_CURVE_RES)
            segLineCos = getLinesFromPts(pts)
            lineCos += segLineCos
            lineColors += [info.segColor for j in range(0, len(segLineCos))]

    tipCos = []
    tipColors = []
    for i, info in enumerate(bptDispInfos):
        pt = info.pt
        for hn in info.handleNos:
            lineCos += [pt[hn], pt[hn + 1]]

            if(len(pt) < 5): htype = defHdlType # For Draw
            else: htype = pt[3 + hn]

            lineColors += [ModalBaseFlexiOp.hdlColMap[htype], \
                    ModalBaseFlexiOp.hdlColMap[htype]]

        # Re-arrange tips so handles are on top of Bezier point
        tc = info.tipColors
        tc = [tc[1], tc[0], tc[2]]
        pt = [pt[1], pt[0], pt[2]]
        for j, tipColor in enumerate(tc):
            if(tipColor != None):
                tipCos.append(pt[j])
                tipColors.append(tipColor)

    bglDrawMgr.addLineInfo('bezLineBatch', FTProps.lineWidth, lineColors, lineCos)
    bglDrawMgr.addPtInfo('bezTipBatch', FTProps.drawPtSize, tipColors, tipCos)

def resetToolbarTool():
    win = bpy.context.window
    scr = win.screen
    areas3d  = [area for area in scr.areas if area.type == 'VIEW_3D']
    override = {'window':win,'screen':scr, 'scene' :bpy.context.scene}
    for a in areas3d:
        override['area'] = a
        regions   = [region for region in a.regions if region.type == 'WINDOW']
        for r in regions:
            override['region'] = r
            bpy.ops.wm.tool_set_by_index(override)

def updateMetaBtns(caller, event, keymap = None):
    if(keymap == None):
        keymap = {'LEFT_SHIFT': 'shift', 'RIGHT_SHIFT': 'shift',
            'LEFT_CTRL':'ctrl', 'RIGHT_CTRL':'ctrl',
            'LEFT_ALT': 'alt', 'RIGHT_ALT': 'alt'}

    var = keymap.get(event.type)

    if(var != None):
        expr = 'caller.' + var + ' = '
        if(event.value == 'PRESS'): exec(expr +'True')
        if(event.value == 'RELEASE'): exec(expr +'False')
        return True

    return False

unitMap = {'FEET': "'", 'METERS':'m'}

class FTHotKeyData:
    def __init__(self, id, key, label, description, isExclusive = False, \
        inclTools = None, exclTools = None):
        if(isExclusive and '+' in key):
            raise ValueError('Exclusive keys cannot be combined with meta keys')
        self.id = id
        self.key = key
        self.label = label
        self.description = description
        self.default = key
        # isExclusive means the hot key function will be invoked even if there are
        # meta keys (that are not part of the hot key combination) are
        # held down together with this key. This way user can have both
        # meta key related functionality (e.g. snapping) amd hot key functionality
        # (e. g. tweak position) simultaneously (Tweak from a snapped point)
        self.isExclusive = isExclusive
        self.exclTools = exclTools if(exclTools != None) else set()
        self.inclTools = inclTools if(inclTools != None) else TOOL_TYPES_FLEXI_ALL
        self.inclTools = self.inclTools - self.exclTools

class FTHotKeys:

    ##################### Key IDs #####################
    # Draw
    hkGrabRepos = 'hkGrabRepos'
    hkUndoLastSeg = 'hkUndoLastSeg'
    hkDissociateHdl = 'hkDissociateHdl'
    hkResetLastHdl = 'hkResetLastHdl'

    drawHotkeys = []

    drawHotkeys.append(FTHotKeyData(hkGrabRepos, 'G', 'Grab Bezier Point', \
            'Grab Bezier point while drawing', inclTools = TOOL_TYPES_FLEXI_DRAW_COMMON, \
                isExclusive = True))
    drawHotkeys.append(FTHotKeyData(hkDissociateHdl, 'V', 'Dissociate Draw Hendle', \
            'Dissociate right handle from left while drawing', \
                inclTools = TOOL_TYPES_FLEXI_DRAW_COMMON, isExclusive = True))
    drawHotkeys.append(FTHotKeyData(hkResetLastHdl, 'Shift+R', 'Reset Last Handle', \
            'Reset last handle so that new segment starts as straight line', \
                inclTools = TOOL_TYPES_FLEXI_DRAW_COMMON))
    drawHotkeys.append(FTHotKeyData(hkUndoLastSeg, 'BACK_SPACE', 'Undo Last Segment', \
            'Undo drawing of last segment', inclTools = TOOL_TYPES_FLEXI_DRAW_COMMON))

    # Edit
    hkUniSubdiv = 'hkUniSubdiv'
    hkBevelPt = 'hkBevelPt'
    hkAlignHdl = 'hkAlignHdl'
    hkDelPtSeg = 'hkDelPtSeg'
    hkToggleHdl = 'hkToggleHdl'
    hkSplitAtSel = 'hkSplitAtSel'
    hkMnHdlType = 'hkMnHdlType'
    hkMnSelect = 'hkMnSelect'
    hkMnDeselect = 'hkMnDeselect'

    editHotkeys = []
    editHotkeys.append(FTHotKeyData(hkUniSubdiv, 'W', 'Segment Uniform Subdivide', \
            'Hotkey to initiate Segment Uniform Subdiv op', \
                inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkBevelPt, 'Ctrl+B', 'Bevel Selected Points', \
            'Hotkey to initiate Bevel op', inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkAlignHdl, 'K', 'Align Handle', \
            'Hotkey to align one handle with the other', \
                inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkDelPtSeg, 'DEL', 'Delete Point / Seg', \
            'Delete selected Point / Segment, align selected handle with other point', \
                inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkToggleHdl, 'H', 'Hide / Unhide Handles', \
            'Toggle handle visibility', inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkSplitAtSel, 'L', 'Split At Selected Points', \
            'Split curve at selected Bezier points', inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkMnHdlType, 'S', 'Set Handle Type', \
            'Set type of selected handles', inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkMnSelect, 'A', 'Select', \
            'Select elements from existing spline selection', \
                inclTools = {TOOL_TYPE_FLEXI_EDIT}))
    editHotkeys.append(FTHotKeyData(hkMnDeselect, 'Alt+A', 'Deselect', \
            'Deselect elements from existing spline selection', \
                inclTools = {TOOL_TYPE_FLEXI_EDIT}))

    # Common
    hkToggleKeyMap = 'hkToggleKeyMap'
    hkSwitchOut = 'hkSwitchOut'
    hkTweakPos = 'hkTweakPos'
    hkToggleDrwEd = 'hkToggleDrwEd'
    hkReorient = 'hkReorient'

    commonHotkeys = []
    commonHotkeys.append(FTHotKeyData(hkToggleKeyMap, 'Ctrl+Shift+H', \
        'Hide / Unhide Keymap', \
            'Hide / Unhide Keymap Displayed When Flexi Tool Is Active'))
    commonHotkeys.append(FTHotKeyData(hkSwitchOut, 'F1', 'Exit Flexi Tool', \
            'Switch out of the Flexi Tool mode'))
    commonHotkeys.append(FTHotKeyData(hkTweakPos, 'P', 'Tweak Position', \
            'Tweak position or enter polar coordinates of the draw / edit point', \
                isExclusive = True))
    commonHotkeys.append(FTHotKeyData(hkToggleDrwEd, 'E', 'Toggle Draw / Edit', \
            'Toggle between Draw & Edit Flexi Tools', \
                exclTools = {TOOL_TYPE_FLEXI_GREASE}))
    commonHotkeys.append(FTHotKeyData(hkReorient, 'U', \
        'Origin to Active Face', \
            'Move origin / orientation to face under mouse cursor \" + \
                "if the origin / orientation is object face'))

    # Snapping
    hkSnapVert = 'hkSnapVert'
    hkSnapGrid = 'hkSnapGrid'
    hkSnapAngle = 'hkSnapAngle'

    # Snapping (Suffix important)
    hkSnapVertMeta = 'hkSnapVertMeta'
    hkSnapGridMeta = 'hkSnapGridMeta'
    hkSnapAngleMeta = 'hkSnapAngleMeta'

    snapHotkeys = [] # Order important
    snapHotkeys.append(FTHotKeyData(hkSnapVert, 'F5', 'Snap to Vert / Face', \
        'Key pressed with mouse click for snapping to Vertex or Face'))
    snapHotkeys.append(FTHotKeyData(hkSnapGrid, 'F6', 'Snap to Grid', \
        'Key pressed with mouse click for snapping to Grid'))
    snapHotkeys.append(FTHotKeyData(hkSnapAngle, 'F7', 'Snap to Angle', \
        'Key pressed with mouse click for snapping to Angle Increment'))

    snapHotkeysMeta = [] # Order should be same as snapHotkeys
    # These keys are not event.type (they can have an entry 'KEY')
    snapHotkeysMeta.append(FTHotKeyData(hkSnapVertMeta, 'ALT', 'Snap to Vert / Face', \
        'Key pressed with mouse click for snapping to Vertex or Face'))
    snapHotkeysMeta.append(FTHotKeyData(hkSnapGridMeta, 'CTRL', 'Snap to Grid', \
        'Key pressed with mouse click for snapping to Grid'))
    snapHotkeysMeta.append(FTHotKeyData(hkSnapAngleMeta, 'SHIFT', 'Snap to Angle', \
        'Key pressed with mouse click for snapping to Angle Increment'))

    # Metakeys not part of the map
    idDataMap = {h.id: h for h in \
        drawHotkeys + editHotkeys + commonHotkeys + snapHotkeys}

    # There can be multiple keydata with same keys (e.g. same key for draw and edit)
    # Not creating a separate map for now. May require later.
    # ~ keyDataMap = {h.key: h for h in \
        # ~ drawHotkeys + editHotkeys + commonHotkeys + snapHotkeys}

    exclKeys = {'RET', 'SPACE', 'ESC', 'X', 'Y', 'Z', \
        'ACCENT_GRAVE', 'COMMA', 'PERIOD', 'F2', 'F3'}

    metas = ['Alt', 'Ctrl', 'Shift']

    def getSnapHotKeys(kId):
        metaId = kId + 'Meta'
        keydatas = [k for k in FTHotKeys.snapHotkeysMeta if k.id == metaId]
        if(len(keydatas) > 0):
            keydata = keydatas[0]
            if(keydata.key != 'KEY'):
                return ['LEFT_' + keydata.key, 'RIGHT_' + keydata.key]
        keydata = FTHotKeys.idDataMap.get(kId)
        return [keydata.key]

    def getKey(key, metas):
        keyVal = ''
        for i, meta in enumerate(metas):
            if(meta): keyVal += FTHotKeys.metas[i] + '+'
        keyVal += key
        return keyVal

    def getHotKeyData(toolType, key, metas):
        # Map will be more efficient, but not so many keys... So ok for now
        kds = [kd for kd in FTHotKeys.idDataMap.values() \
            if(kd.key == FTHotKeys.getKey(key, metas) and toolType in kd.inclTools)]
        return kds[0] if len(kds) > 0 else None

    def isHotKey(id, key, metas):
        currKeyData = FTHotKeys.idDataMap[id]
        # Only compare part without meta for exclusive keys
        if(currKeyData.isExclusive):
            return currKeyData.key == key
        else:
            return currKeyData.key == FTHotKeys.getKey(key, metas)

    def haveCommonTool(keyData1, keyData2):
        return len(keyData1.inclTools.intersection(keyData2.inclTools)) > 0

    # The regular part of the snap keys is validated against assigned key without meta
    # So that if e.g. Ctrl+B is already assigned, B is not available as reg part
    def isAssignedWithoutMeta(kId, key):
        if(key.endswith(INVAL)): return True
        return any([kd.key.split('+')[-1] == key for kd in FTHotKeys.idDataMap.values() \
            if kd.id != kId and FTHotKeys.haveCommonTool(kd, FTHotKeys.idDataMap[kId])])

    def isAssigned(kId, key):
        if(key.endswith(INVAL)): return True
        currKeyData = FTHotKeys.idDataMap.get(kId)
        exclKeys = [kd for kd in FTHotKeys.idDataMap.values() \
            if(FTHotKeys.haveCommonTool(kd, currKeyData) and ((key == kd.key) or \
                (kd.isExclusive and (key.split('+')[-1] == kd.key)) or \
                    (currKeyData.isExclusive and (key == kd.key.split('+')[-1]))))]
        return len(exclKeys) > 0

    updating = False

    # Validation for single key text field (format: meta key + regular key)
    # UI Format: Key and 3 toggle buttons for meta
    # TODO: Separate update for each id?
    # TODO: Reset on any exception?
    def updateHotkeys(dummy, context):
        try:
            FTHotKeys.updateHKPropPrefs(context)
            # ~ FTHotKeys.keyDataMap = {h.key: h for h in \
                # ~ [k for k in (FTHotKeys.drawHotkeys + FTHotKeys.editHotkeys + \
                    # ~ FTHotKeys.commonHotkeys + FTHotKeys.snapHotkeys)]}
        except Exception as e:
            FTHotKeys.updateHKPropPrefs(context, reset = True)
            print("BezierUtils: Error fetching keymap preferences", e)

    # TODO: This method combines two rather different functionality based on reset flag
    # TODO: Metakey setting in updateSnapMetaKeys and reset here also
    def updateHKPropPrefs(context, reset = False):
        if(FTHotKeys.updating): return

        FTHotKeys.updating = True
        prefs = context.preferences.addons[__name__].preferences
        hmap = FTHotKeys.idDataMap
        # Checking entire map even if one key changed (TODO)
        for kId in hmap:
            if(reset): hmap[kId].key = hmap[kId].default
            # User pressed key
            prefKey = getattr(prefs, kId)
            combKey = ''
            for meta in FTHotKeys.metas:
                if(hasattr(prefs, kId + meta) and \
                    getattr(prefs, kId + meta) == True):
                    combKey += meta + '+'
            # key in map has format - meta1 + met2 ... + (event.type)
            prefKey = combKey + prefKey
            if(hmap[kId].key == prefKey): continue
            valid = (prefKey != INVAL and not reset)
            if(valid):
                if(kId in [k.id for k in FTHotKeys.snapHotkeys]):
                    valid = not FTHotKeys.isAssignedWithoutMeta(kId, prefKey)
                else:
                    valid = not FTHotKeys.isAssigned(kId, prefKey)
            # Revert to key from map if invalid
            if(not valid):
                comps = hmap[kId].key.split('+')
                setattr(prefs, kId, comps[-1])
                for meta in FTHotKeys.metas:
                    setattr(prefs, kId + meta, meta in comps[:-1])
            else:
                hmap[kId].key = prefKey

        if(reset):
            for i, metakeyData in enumerate(FTHotKeys.snapHotkeysMeta):
                metakeyData.key = metakeyData.default
                metaKeyId = metakeyData.id
                setattr(prefs, metaKeyId, metakeyData.key)

        FTHotKeys.updating = False
        ModalBaseFlexiOp.propsChanged()

    def getAvailKey(prefs):
        availKey = None
        oldKeys = set([m.upper() for m in FTHotKeys.metas])
        newKeys = set([getattr(prefs, kd.id) \
            for kd in FTHotKeys.snapHotkeysMeta])
        for availKey in (oldKeys - newKeys): break # Only one entry
        return availKey

    def initSnapMetaFromPref(context):
        prefs = context.preferences.addons[__name__].preferences
        for i, metakeyData in enumerate(FTHotKeys.snapHotkeysMeta):
            metaKeyId = metakeyData.id
            prefKeyMeta = getattr(prefs, metaKeyId)
            metakeyData.key = prefKeyMeta
            if(prefKeyMeta == 'KEY'):
                regKeyId = metaKeyId[0:metaKeyId.index('Meta')]
                regKeyData = FTHotKeys.idDataMap[regKeyId]
                prefKeyReg = getattr(prefs, regKeyId)
                regKeyData.key = prefKeyReg
        ModalBaseFlexiOp.propsChanged()

    # Validation for snap keys (format: EITHER meta key OR regular key)
    # (regular part validated by updateHotkeys)
    # UI Format: Drop-down with entries Ctrl, Alt, Shift, Keyboard and ...
    # a separate single key field activated with selection of 'Keyboard' in the dropdown
    def updateSnapMetaKeys(dummy, context):
        if(FTHotKeys.updating): return

        FTHotKeys.updating = True
        prefs = context.preferences.addons[__name__].preferences
        changedKeyIdx = -1
        prefKeyMeta = None

        for i, metakeyData in enumerate(FTHotKeys.snapHotkeysMeta):
            metaKeyId = metakeyData.id
            prefKeyMeta = getattr(prefs, metaKeyId)

            if(prefKeyMeta != metakeyData.key): # Changed entry
                changedKeyIdx = i
                break

        if(changedKeyIdx != -1):
            metaKeyId = metakeyData.id # loop broke at metakeyData
            metakeyData.key = prefKeyMeta
            if(prefKeyMeta == 'KEY'):
                regKeyId = metaKeyId[0:metaKeyId.index('Meta')]
                regKeyData = FTHotKeys.idDataMap[regKeyId]
                prefKeyReg = getattr(prefs, regKeyId)
                # Invalid regular key->assign available meta key to the dropdown entry
                if(prefKeyReg == INVAL or \
                    FTHotKeys.isAssignedWithoutMeta(regKeyId, prefKeyReg)):
                    setattr(prefs, regKeyId, regKeyData.key)
                else:
                    regKeyData.key = prefKeyReg
            else:
                # Change the other drop-down with the same key
                availKey = FTHotKeys.getAvailKey(prefs)
                for i, kd in enumerate(FTHotKeys.snapHotkeysMeta):
                    if(kd.key == prefKeyMeta and i != changedKeyIdx):
                        kd.key = availKey
                        setattr(prefs, kd.id, availKey)
                        break

        FTHotKeys.updating = False
        ModalBaseFlexiOp.propsChanged()

    def keyCodeMap():
        kcMap = {}
        digits = ['ZERO', 'ONE', 'TWO', 'THREE', 'FOUR', 'FIVE', \
                    'SIX', 'SEVEN', 'EIGHT', 'NINE']#, 'PERIOD']
        numpadDigits = ['NUMPAD_' + d for d in digits]
        # ~ 199: 'NUMPAD_PERIOD'

        for i in range(26):
            kcMap[i + 97] = chr(ord('A') + i)
        # Digits not available, used for keyboard input
        # ~ for i in range(10):
            # ~ kcMap[i + 48] = digits[i]
        for i in range(10):
            kcMap[i + 150] = numpadDigits[i]
        for i in range(12):
            kcMap[i + 300] = 'F' + str(i + 1)

        kcMap[161] = 'NUMPAD_SLASH'
        kcMap[160] = 'NUMPAD_ASTERIX'
        kcMap[162] = 'NUMPAD_MINUS'
        kcMap[163] = 'NUMPAD_ENTER'
        kcMap[164] = 'NUMPAD_PLUS'
        kcMap[165] = 'PAUSE'
        kcMap[166] = 'INSERT'
        kcMap[167] = 'HOME'
        kcMap[168] = 'PAGE_UP'
        kcMap[169] = 'PAGE_DOWN'
        kcMap[170] = 'END'
        kcMap[199] = 'NUMPAD_PERIOD'
        kcMap[219] = 'TAB'
        kcMap[220] = 'RET'
        kcMap[221] = 'SPACE'
        kcMap[223] = 'BACK_SPACE'
        kcMap[224] = 'DEL'
        kcMap[225] = 'SEMI_COLON'
        kcMap[226] = 'PERIOD'
        kcMap[227] = 'COMMA'
        kcMap[228] = "QUOTE"
        kcMap[229] = 'ACCENT_GRAVE'
        # ~ kcMap[230] = 'MINUS' # Reserved for keyboard input
        kcMap[232] = 'SLASH'
        kcMap[233] = 'BACK_SLASH'
        kcMap[234] = 'EQUAL'
        kcMap[235] = 'LEFT_BRACKET'
        kcMap[236] = 'RIGHT_BRACKET'

        return kcMap

    # Drop down item tuples with 400 entries
    # INVAL in all 3 fields for unavailable keycodes
    # TODO: Very big drop-down for every hotkey (because of default)
    def getKeyMapTupleStr():
        kcMap = FTHotKeys.keyCodeMap()
        tuples = []
        exclKeys = FTHotKeys.exclKeys

        for i in range(0, 400):
            if(kcMap.get(i) != None and kcMap[i] not in exclKeys):
                char = kcMap[i]
            else:
                char = INVAL
            tuples.append((char, char, char))
        return ''.join(["('" + t[0] + "','" + t[1] + "','" + t[2] +"')," \
            for t in tuples])

    def getHKFieldStr(keydata, addMeta):
        propName = keydata.id
        text = keydata.label
        description = keydata.description
        updateFn = 'FTHotKeys.updateHotkeys'
        default = keydata.default
        keys = default.split('+')
        key = keys[-1]
                # ~ "items = (('A', 'A','A'),('B', 'B', 'B')), " + \
        retVal = propName + ": EnumProperty(name = '" + text + "', " + \
                "items = (" + FTHotKeys.getKeyMapTupleStr() + "), " + \
                (("default = '"+ key + "', ") if default != INVAL else '') + \
                "update = " + updateFn + ", " + \
                "description = '"+ description +"') \n"
        if(addMeta):
            for meta in FTHotKeys.metas:
                retVal += propName + meta + ": BoolProperty(name='"+ meta +"', " + \
                        "description='"+ meta +" Key', " + \
                        "update = " + updateFn + ", " + \
                        "default = "+ ('True' if meta in keys[:-1] else 'False') +") \n"
        return retVal

    def getMetaHKFieldStr(keydataMeta):
        propName = keydataMeta.id
        text = keydataMeta.label
        description = keydataMeta.description
        updateFn = 'FTHotKeys.updateSnapMetaKeys'
        default = keydataMeta.default

        itemsStr = "('CTRL', 'Ctrl', 'Ctrl'), ('ALT', 'Alt', 'Alt'), \
                ('SHIFT', 'Shift', 'Shift'), ('KEY', 'Keyboard', 'Other keyboard key')"
        retVal = propName + ": EnumProperty(name = '" + text + "', " + \
                "items = (" + itemsStr + "), " + \
                "default = '"+ default + "', " + \
                "update = " + updateFn + ", " + \
                "description = '"+ description +"') \n"
        return retVal

    def getHKDispLines(toolType):
        hkData = [k for k in FTHotKeys.commonHotkeys if toolType not in k.exclTools]
        for i, hk in enumerate(FTHotKeys.snapHotkeysMeta):
            if(hk.key == 'KEY' and toolType not in hk.exclTools):
                hkData.append(FTHotKeys.snapHotkeys[i])
            else:
                hkData.append(hk)
        if(toolType in TOOL_TYPES_FLEXI_DRAW_COMMON):
            hkData += [k for k in FTHotKeys.drawHotkeys if toolType not in k.exclTools]
        if(toolType == TOOL_TYPE_FLEXI_EDIT):
            hkData += [k for k in FTHotKeys.editHotkeys if toolType not in k.exclTools]

        labels = []
        config = []
        keys = []
        stdKeylabels = []
        stdKeylabels.append(['Lock to YZ Plane', 'Shift+X'])
        stdKeylabels.append(['Lock to XZ Plane', 'Shift+Y'])
        stdKeylabels.append(['Lock to XY Plane', 'Shift+Z'])
        stdKeylabels.append(['Lock to Z-Axis', 'Z'])
        stdKeylabels.append(['Lock to Y-Axis', 'Y'])
        stdKeylabels.append(['Lock to X-Axis', 'X'])
        stdKeylabels.append(['Confirm Operation', 'Space / Enter'])
        for k in hkData:
            labels.append(k.label)
            config.append(True)
            keys.append(k.key)
        for k in stdKeylabels:
            labels.append(k[0])
            config.append(False)
            keys.append(k[1])

        return config, labels, keys


# Flexi Tool Constants
class FTProps:
    propUpdating = False

    def updateProps(dummy, context):
        FTProps.updatePropsPrefs(context)

    def updatePropsPrefs(context, resetPrefs = False):
        if(FTProps.propUpdating): return
        FTProps.propUpdating = True
        try:
            prefs = context.preferences.addons[__name__].preferences

            props = ['drawPtSize', 'lineWidth', 'axisLineWidth', 'editSubdivPtSize', \
            'greaseSubdivPtSize', 'colDrawSelSeg', 'colDrawNonHltSeg', 'colDrawHltSeg', \
            'colDrawMarker', 'colGreaseSelSeg', 'colGreaseNonHltSeg', 'colGreaseMarker', \
            'colHdlFree', 'colHdlVector', 'colHdlAligned', 'colHdlAuto', 'colSelTip', \
            'colHltTip', 'colBezPt', 'colHdlPtTip', 'colAdjBezTip', 'colEditSubdiv', \
            'colGreaseSubdiv', 'colGreaseBezPt', 'colKeymapText', 'colKeymapKey', 'snapDist', \
            'dispSnapInd', 'dispAxes', 'snapPtSize', 'liveUpdate', 'dispCurveRes', \
            'showKeyMap', 'keyMapFontSize', 'keyMapLocX', 'keyMapLocY', 'keyMapNextToTool', \
            'defBevelFact', 'maxBevelFact', 'minBevelFact', 'bevelIncr', 'numpadEntry', \
            'mathFnTxtFontSize', 'colMathFnTxt']

            if(resetPrefs):
                FTProps.initDefault()
                for prop in props:
                    exec('prefs.' + prop +' = FTProps.' + prop)
                    # ~ setattr(prefs, prop, getattr(FTProps, prop))
            else:
                for prop in props:
                    exec('FTProps.' + prop +' = prefs.' + prop)
                    # ~ setattr(FTProps, prop, getattr(prefs, prop))

        except Exception as e:
            print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
            FTProps.initDefault()

        ModalBaseFlexiOp.propsChanged()

        try: ModalBaseFlexiOp.opObj.refreshDisplaySelCurves()
        except: pass

        FTProps.propUpdating = False

    def initDefault():
        FTProps.drawPtSize = 5
        FTProps.lineWidth = 1.5
        FTProps.axisLineWidth = .25
        FTProps.editSubdivPtSize = 6
        FTProps.greaseSubdivPtSize = 4

        FTProps.colDrawSelSeg = (.6, .8, 1, 1)
        FTProps.colDrawNonHltSeg = (.1, .4, .6, 1)
        FTProps.colDrawHltSeg = (.2, .6, .9, 1)

        FTProps.colGreaseSelSeg = (0.2, .8, 0.2, 1)
        FTProps.colGreaseNonHltSeg = (0.2, .6, 0.2, 1)

        FTProps.colHdlFree = (.6, .05, .05, 1)
        FTProps.colHdlVector = (.4, .5, .2, 1)
        FTProps.colHdlAligned = (1, .3, .3, 1)
        FTProps.colHdlAuto = (.8, .5, .2, 1)

        FTProps.colDrawMarker = (.6, .8, 1, 1)
        FTProps.colGreaseMarker = (0.2, .8, 0.2, 1)

        FTProps.colSelTip = (.2, .7, .3, 1)
        FTProps.colHltTip = (.8, 1, .8, 1)
        FTProps.colBezPt = (1, 1, 0, 1)
        FTProps.colHdlPtTip = (.7, .7, 0, 1)
        FTProps.colAdjBezTip = (.1, .1, .1, 1)

        FTProps.colEditSubdiv = (.3, 0, 0, 1)

        FTProps.colGreaseSubdiv = (1, .3, 1, 1)
        FTProps.colGreaseBezPt = (1, .3, 1, 1)
        FTProps.colKeymapText = (1.0, 1.0, 1.0, 1.0)
        FTProps.colKeymapKey = (0.0, 1.0, 1.0, 1.0)

        FTProps.snapDist = 20
        FTProps.dispSnapInd = False
        FTProps.dispAxes = True
        FTProps.showKeyMap = False
        FTProps.keyMapFontSize = 10
        FTProps.keyMapLocX = 10
        FTProps.keyMapLocY = 10
        FTProps.keyMapNextToTool = True
        FTProps.snapPtSize = 3
        FTProps.liveUpdate = False
        FTProps.dispCurveRes = .4
        FTProps.defBevelFact = 4
        FTProps.maxBevelFact = 15
        FTProps.minBevelFact = -15
        FTProps.bevelIncr = .5
        FTProps.numpadEntry = False

        FTProps.mathFnTxtFontSize = 20
        FTProps.colMathFnTxt = (0.6, 1.0, 0.03, 1.0)

class FTMenuData:

    def __init__(self, hotkeyId, options, menuClassName, menuClassLabel, handler):
        self.hotkeyId = hotkeyId
        self.options = options
        self.menuClassName = menuClassName
        self.menuClassLabel = menuClassLabel
        self.handler = handler

# TODO: Move to end to avoid eval/getattr due to forward referencing
class FTMenu:

    propSuffix = 'ftMenuOpt'

    # Edit
    editMenus = []
    editMenus.append(FTMenuData(FTHotKeys.hkMnHdlType, \
        [['miHdlAuto', 'Auto', 'HANDLETYPE_AUTO_VEC'], \
         ['miHdlAligned', 'Aligned', 'HANDLETYPE_ALIGNED_VEC'], \
         ['miHdlFree', 'Free', 'HANDLETYPE_FREE_VEC'], \
         ['miHdlVector', 'Vector', 'HANDLETYPE_VECTOR_VEC']], \
            'VIEW3D_MT_FlexiEditHdlMenu', 'Set Handle Type', 'mnSetHdlType'))

    editMenus.append(FTMenuData(FTHotKeys.hkMnSelect, \
        [['miSelSegs', 'Segments', 'GP_SELECT_BETWEEN_STROKES'], \
         ['miSelAllSplines', 'All Splines', 'GP_MULTIFRAME_EDITING'], \
         ['miSelBezPts', 'Bezier Points', 'GP_ONLY_SELECTED'], \
         ['miSelHdls', 'Handles', 'CURVE_BEZCURVE'], \
         ['miSelObj', 'Curve Object', 'GP_SELECT_STROKES'], \
         ['miSelAll', 'Everything', 'SELECT_EXTEND']], \
            'VIEW3D_MT_FlexiEditSelMenu', 'Select', 'mnSelect'))

    editMenus.append(FTMenuData(FTHotKeys.hkMnDeselect, \
        [['miDeselSegs', 'Segments', 'GP_SELECT_BETWEEN_STROKES'], \
         ['miDeselBezPts', 'Bezier Points', 'GP_SELECT_POINTS'], \
         ['miDeselHdls', 'Handles', 'CURVE_BEZCURVE'], \
         ['miDeselObj', 'Curve Object', 'GP_SELECT_STROKES'], \
         ['miDeselInvert', 'Invert Selection', 'SELECT_SUBTRACT']], \
            'VIEW3D_MT_FlexiEditDeselMenu', 'Deselect', 'mnDeselect'))

    idDataMap = {m.hotkeyId: m for m in editMenus}
    toolClassMap = {'ModalFlexiEditBezierOp': set([m.hotkeyId for m in editMenus])}
    toolTypeMap = {'ModalFlexiEditBezierOp': TOOL_TYPE_FLEXI_EDIT}

    currMenuId = None
    abandoned = False

    def getMenuData(caller, hotkeyId):
        keyIds = FTMenu.toolClassMap.get(caller.__class__.__name__)
        if(keyIds != None and hotkeyId in keyIds):
            return FTMenu.idDataMap.get(hotkeyId)
        return None

    def procMenu(parent, context, event, outside):

        metakeys = parent.snapper.getMetakeys()
        evtType = event.type

        if(FTMenu.abandoned == True and evtType == 'ESC'):
            FTMenu.abandoned = False
            return True

        if(FTMenu.currMenuId != None):
            params = bpy.context.window_manager.bezierToolkitParams
            opt = FTMenu.getCurrMenuSel()
            if(opt != None or not evtType.startswith('TIMER')): # What's TIMER_REPORT?
                context.window_manager.event_timer_remove(parent.menuTimer)
                menuData = FTMenu.idDataMap.get(FTMenu.currMenuId)
                FTMenu.currMenuId = None
                if(evtType == 'TIMER'):
                    fn = getattr(parent, menuData.handler)
                    fn(opt)
                else:
                    FTMenu.abandoned = True
            return True

        if(outside): return False
        parentClassName = parent.__class__.__name__
        hkData = FTHotKeys.getHotKeyData(FTMenu.toolTypeMap.get(parentClassName), \
            evtType, metakeys)

        if(hkData == None): return False
        menuData = FTMenu.getMenuData(parent, hkData.id)
        if(menuData == None): return False
        if(event.value == 'RELEASE'):
            FTMenu.resetMenuOptions()
            FTMenu.currMenuId = menuData.hotkeyId
            parent.menuTimer = \
                context.window_manager.event_timer_add(time_step = 0.0001, \
                    window = context.window)
            ret = bpy.ops.wm.call_menu_pie(name = menuData.menuClassName)
        return True

    def getCurrMenuSel():
        menuData = FTMenu.idDataMap.get(FTMenu.currMenuId)
        if(menuData != None):
            params = bpy.context.window_manager.bezierToolkitParams
            for i, propName in enumerate(FTMenu.getAllOptPropNames()):
                if(getattr(params, propName)): return menuData.options[i]

        return None

    def resetMenuOptions():
        params = bpy.context.window_manager.bezierToolkitParams
        for propName in FTMenu.getAllOptPropNames():
            setattr(params, propName, False)

    def getAllOptPropNames():
        propNames = []
        maxMenuOpts = max(len(m.options) for m in FTMenu.idDataMap.values())
        for i in range(maxMenuOpts):
            propNames.append(FTMenu.propSuffix + str(i))
        return propNames

    def getMNClassDefStr(menuData):
        retStr = 'class ' + menuData.menuClassName + '(Menu):\n' + \
            '\tbl_label = "' + menuData.menuClassLabel + '"\n'+ \
            '\tdef draw(self, context): \n' + \
            '\t\tparams = bpy.context.window_manager.bezierToolkitParams\n' + \
            '\t\tlayout = self.layout\n' + \
            '\t\tpie = layout.menu_pie()\n'
        for i, opt in enumerate(menuData.options):
            retStr += '\t\top = pie.operator("object.ft_menu_options", text = "'+ \
                opt[1] + '"' + \
                    ((', icon = "'+ opt[2] +'"') if(opt[2] != '') else '') + ')\n'
            retStr += '\t\top.optIdx = ' + str(i) + '\n'
        return retStr

    def getMNPropDefStr(menuData):
        retStr = ''
        for propName in FTMenu.getAllOptPropNames():
            retStr += propName +': BoolProperty(default = False)\n'
        return retStr

class FTMenuOptionOp(Operator):
    bl_idname = "object.ft_menu_options"
    bl_label = "Set FT Menu Options"
    bl_description = "Set option"

    optIdx : IntProperty()

    def execute(self, context):
        FTMenu.resetMenuOptions()
        params = bpy.context.window_manager.bezierToolkitParams
        setattr(params, FTMenu.propSuffix + str(self.optIdx), True)
        return {'FINISHED'}


class SnapDigits:
    digitMap = {'ZERO':'0', 'ONE':'1', 'TWO':'2', 'THREE':'3', 'FOUR':'4', 'FIVE':'5', \
                'SIX':'6', 'SEVEN':'7', 'EIGHT':'8', 'NINE':'9', 'PERIOD':'.'}
    numpadDigitMap = {'NUMPAD_0':'0', 'NUMPAD_1':'1', 'NUMPAD_2':'2', 'NUMPAD_3':'3', \
        'NUMPAD_4':'4', 'NUMPAD_5':'5', 'NUMPAD_6':'6', 'NUMPAD_7':'7', 'NUMPAD_8':'8', \
            'NUMPAD_9':'9', 'NUMPAD_PERIOD': '.'}

    def getValidFloat(sign, digits):
        delta = 0
        valid = True
        for i in range(len(digits), 0, -1):
            td = digits[0:i]
            try:
                delta = float(sign + ''.join(td))
                return delta, valid
            except:
                valid = False
        return delta, valid

    def __init__(self, getFreeAxes, getEditCoPair):
        self.getFreeAxes = getFreeAxes
        self.getEditCoPair = getEditCoPair
        self.initialize()

    def initialize(self):
        self.axes = None # [0, 1, 2] etc.
        self.deltaVec = None # Always cartesian
        self.axisIdx = 0
        self.digitChars = []
        self.signChar = ''
        self.polar = False
        self.pDataIdx = 0 # Corresponding to axisIdx

    def hasVal(self):
        return self.axes != None

    # Theta in degrees
    def getPolarCos(self):
        axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
        val0, val1 = self.deltaVec[axis0], self.deltaVec[axis1]
        if(val0 == 0):
            if(val1 == 0): return 0, 0
            theta = (val1 / abs(val1)) * 90
        else:
            theta = degrees(atan(abs(val1) / abs(val0)))
            if(val0 < 0): theta = 180 - theta
            if(val1 < 0): theta = -theta
        r = sqrt(val1 * val1 + val0 * val0)
        return [r, theta]

    # Theta in degrees
    def getDeltaFromPolar(self, polCos):
        r, theta = polCos
        return r * cos(radians(theta)), r * sin(radians(theta))

    def addToPolar(self, delta):
        polCos = list(self.getPolarCos())
        polCos[self.pDataIdx] += delta
        val0, val1 = self.getDeltaFromPolar(polCos)
        return val0, val1

    def digitsToVec(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        if(not self.polar or self.pDataIdx == 0): delta /= getUnitScale()
        if(self.polar):
            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
            self.deltaVec[axis0], self.deltaVec[axis1] = self.addToPolar(delta)
        else:
            self.deltaVec[self.axisIdx] += delta

    def vecToDigits(self):
        if(self.polar):
            vals = self.getPolarCos()
            v = round(vals[self.pDataIdx], 4)
        else:
            v = round(self.deltaVec[self.axisIdx], 4)
        self.digitChars = list(str(abs(v))) if v != 0 else ''
        while(len(self.digitChars) > 0 and self.digitChars[-1] == '0'):
            self.digitChars.pop()
        if(len(self.digitChars) > 0 and self.digitChars[-1] == '.'):
            self.digitChars.pop()
        self.signChar = '-' if v < 0 else ''

    def updateDeltaFromEditCos(self):
        editCos = self.getEditCoPair()
        self.deltaVec = editCos[1] - editCos[0]

    def procEvent(self, context, event, metakeys):
        if(FTHotKeys.isHotKey(FTHotKeys.hkTweakPos, event.type, metakeys)):
            editCos = self.getEditCoPair()
            if(len(editCos) == 2):
                if(event.value == 'RELEASE'):
                    if(self.axes == None):
                        self.axes = self.getFreeAxes() # TODO: Always dynamic?
                        self.polar = False
                    # Only possible within a plane right now
                    elif(len(self.getFreeAxes()) == 2):
                        self.polar = not self.polar
                        self.pDataIdx = 0
                    self.updateDeltaFromEditCos()
                    self.digitChars = []
                    self.axisIdx = self.axes[0]
            return True

        dmap = self.digitMap.copy()
        if(FTProps.numpadEntry):
            dmap.update(self.numpadDigitMap)
        dval = dmap.get(event.type)
        if(dval != None):
            if(event.value == 'RELEASE'):
                self.digitChars.append(dval)
                if(self.axes == None):
                    self.axes = self.getFreeAxes() # TODO: Always dynamic?
                    self.deltaVec = Vector()
                    self.axisIdx = self.axes[0]
            return True

        if(not self.hasVal()): # No further processing if nothing to process
            return False

        retVal = True
        if(event.type == 'MINUS'):
            if(event.value == 'RELEASE'):
                self.signChar = '' if(self.signChar == '-') else '-'

        elif(event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                self.initialize()

        elif(event.type == 'BACK_SPACE'):
            if(event.value == 'RELEASE'):
                if(len(self.digitChars) > 0): self.digitChars.pop()
                else:
                    if(self.polar):
                        polarCos = self.getPolarCos()
                        if(not floatCmpWithMargin(polarCos[self.pDataIdx], 0)):
                            self.vecToDigits()

                            # TODO: Retain theta without this workaround
                            polarCos[self.pDataIdx] = 0.00001
                            delta = self.getDeltaFromPolar(polarCos)
                            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
                            self.deltaVec[axis0] = delta[0]
                            self.deltaVec[axis1] = delta[1]
                    else:
                        if(self.deltaVec[self.axisIdx] != 0):
                            self.vecToDigits()
                            self.deltaVec[self.axisIdx] = 0

        elif(event.type == 'TAB'):
            if(event.value == 'RELEASE'):
                self.digitsToVec()
                if(self.polar):
                    self.pDataIdx = (self.pDataIdx + 1) % 2
                else:
                    self.axisIdx = \
                        self.axes[(self.axes.index(self.axisIdx) + 1) % len(self.axes)]

                self.digitChars = []
                self.signChar = ''
        else:
            retVal = False

        return retVal

    def getCurrDelta(self):
        if(self.axes == None): return Vector()

        val = self.deltaVec.copy()
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        if(not self.polar or self.pDataIdx == 0): delta /= getUnitScale()
        if(self.polar):
            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
            val[axis0], val[axis1] = self.addToPolar(delta)
        else:
            val[self.axisIdx] += delta

        return val

    def getDeltaStrPolar(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        polCos = self.getPolarCos()
        polCos[0] *= getUnitScale()
        strs = ['['] * 2
        idx0 = self.pDataIdx
        idx1 = 1 - self.pDataIdx
        d = polCos[idx0]

        if(d != 0):
            strs[idx0] += str(round(d, 4)) + ('+' if (self.signChar == '') else '')

        strs[idx0] += self.signChar + ''.join(self.digitChars) + '] = ' \
            + str(round((delta + d), 4)) + ('' if valid else ' <Invalid>')

        strs[idx1] = '[' + str(round(polCos[1 - self.pDataIdx], 4)) + ']'

        return 'r: '+strs[0] + ' theta: '+ strs[1]

    def getCurrDeltaStr(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        retStr = '['
        d = self.deltaVec[self.axisIdx] * getUnitScale()
        if(d != 0):
            retStr += str(round(d, 4)) + ('+' if (self.signChar == '') else '')
        retStr += self.signChar + ''.join(self.digitChars) + '] = ' \
            + str(round((delta + d), 4)) + ('' if valid else ' <Invalid>')

        return retStr

class CustomAxis:

    def __init__(self):
        #TODO: What's better?
        if(bpy.data.scenes[0].get('btk_co1') == None):
            bpy.data.scenes[0]['btk_co1'] = LARGE_VECT
        if(bpy.data.scenes[0].get('btk_co2') == None):
            bpy.data.scenes[0]['btk_co2'] = LARGE_VECT
        self.axisPts = [Vector(bpy.data.scenes[0]['btk_co1']), \
            Vector(bpy.data.scenes[0]['btk_co2'])]
        if(bpy.data.scenes[0].get('btk_snapPtCnt') == None):
            bpy.data.scenes[0]['btk_snapPtCnt'] = 3
        self.snapCnt = bpy.data.scenes[0]['btk_snapPtCnt']
        self.inDrawAxis = False # User drawing the custom axis

    def length(self):
        pts = self.axisPts

        # Strange floating points!
        if(all(pt < (LARGE_NO - 1000) for pt in pts[0] + pts[1])):
            return (pts[1] - pts[0]).length
        else:
            return 0

    def set(self, idx, co):
        self.axisPts[idx] = co
        if(idx == 0): bpy.data.scenes[0]['btk_co1'] = [c for c in co]
        if(idx == 1): bpy.data.scenes[0]['btk_co2'] = [c for c in co]
        bpy.data.scenes[0]['btk_snapPtCnt'] = self.snapCnt

    def getSnapPts(self): # ptCnt excluding end points
        pts = self.axisPts
        if(self.length() == 0): return [pts[0], pts[1]]

        snapPts = [pts[0]]
        interval = self.snapCnt + 1
        incr = self.length() / interval
        diffV = pts[1] - pts[0]

        for i in range(1, interval):
           snapPts.append(pts[0] + (diffV * (incr * i)  / self.length()) )

        snapPts.append(pts[1])
        return snapPts

    def procDrawEvent(self, context, event, snapper, rmInfo):
        if(event.type == 'RIGHTMOUSE'):
            snapOrigin = bpy.context.window_manager.bezierToolkitParams.snapOrigin
            if(event.value == 'RELEASE' and snapOrigin == 'AXIS'):
                loc = snapper.get3dLocSnap(rmInfo, \
                    SnapParams(snapper, snapToAxisLine = False))
                if(not self.inDrawAxis): self.set(0, loc)
                else: self.set(1, loc)
                self.inDrawAxis = not self.inDrawAxis
            return True

        if(self.inDrawAxis):
            if(event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', 'NUMPAD_PLUS', \
                'NUMPAD_MINUS','PLUS', 'MINUS'}):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return True
                elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                    if(self.snapCnt < 20): self.snapCnt += 1
                elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                    if(self.snapCnt > 0): self.snapCnt -= 1

            if(event.type == 'MOUSEMOVE'):
                loc = snapper.get3dLocSnap(rmInfo, \
                    SnapParams(snapper, snapToAxisLine = False))
                self.set(1, loc)

            if(event.type == 'ESC'):
                self.set(0, LARGE_VECT)
                self.set(1, LARGE_VECT)
                self.inDrawAxis = False

            return True

        return False


# TODO: Make independent of snapper
class SnapParams:
    def __init__(self, \
        snapper, \
        xyDelta = [0, 0], \
        vec = None, \
        refreshStatus = True, \
        snapToAxisLine = True, \
        lastCo1Axis = False, \
        enableSnap = True, \
        vertSnap = None, \
        gridSnap = None, \
        angleSnap = None, \
        refLine = None, \
        refLineOrig = None, \
        selCo = None, \
        inEdit = None, \
        hasSel = None, \
        transType = None, \
        origType = None, \
        axisScale = None, \
        freeAxesN = None, \
        dispAxes = True,
        snapToPlane = None):

        self.xyDelta = xyDelta
        self.vec = vec
        self.refreshStatus = refreshStatus
        self.snapToAxisLine = snapToAxisLine
        self.lastCo1Axis = lastCo1Axis

        if(enableSnap):
            self.vertSnap = snapper.vertSnap if(vertSnap == None) else vertSnap
            self.gridSnap = snapper.gridSnap if(gridSnap == None) else gridSnap
            self.angleSnap = snapper.angleSnap if(angleSnap == None) else angleSnap
        else:
            self.vertSnap = False
            self.gridSnap = False
            self.angleSnap = False

        self.refLine = snapper.getRefLine() if(refLine == None) else refLine
        self.refLineOrig = snapper.getRefLineOrig() if(refLineOrig == None) \
            else refLineOrig
        self.selCo = snapper.getSelCo() if(selCo == None) else selCo

        self.inEdit = snapper.isEditing() if(inEdit == None) else inEdit
        self.hasSel = snapper.hasSelection() if(hasSel == None) else hasSel

        params = bpy.context.window_manager.bezierToolkitParams
        self.transType = params.snapOrient if(transType == None) else transType
        self.origType = params.snapOrigin if(origType == None) else origType
        self.axisScale = params.axisScale if(axisScale == None) else axisScale

        self.freeAxesN = snapper.getFreeAxesNormalized() \
            if(freeAxesN == None) else freeAxesN

        self.dispAxes = dispAxes

        if(snapToPlane != None):
            self.snapToPlane = snapToPlane
        else:
            if(showSnapToPlane(params)): # Only if Snap to Plane option is visible
                self.snapToPlane = params.snapToPlane
            else: self.snapToPlane = False


class Snapper:

    DEFAULT_ANGLE_SNAP_STEPS = 3
    MAX_SNAP_VERT_CNT = 1000
    MAX_SNAP_FACE_CNT = 1000

    def __init__(self, context, getSnapLocs, getRefLine, getRefLineOrig, getSelCo, \
        getCurrLine, hasSelection, isEditing):
        self.getSnapLocs = getSnapLocs
        self.getRefLine = getRefLine
        self.getCurrLine = getCurrLine
        self.getRefLineOrig = getRefLineOrig
        self.getSelCo = getSelCo
        self.hasSelection = hasSelection
        self.isEditing = isEditing
        self.angleSnapSteps = Snapper.DEFAULT_ANGLE_SNAP_STEPS
        self.customAxis = CustomAxis()
        self.snapDigits = SnapDigits(self.getSnapParamsFreeAxes, self.getEditCoPair)
        self.initialize()
        self.updateSnapKeyMap()
        # ~ bpy.context.space_data.overlay.show_axis_x = False
        # ~ bpy.context.space_data.overlay.show_axis_y = False

    def resetMetakeys(self):
        # metakeys for all functions
        self.shift = False
        self.ctrl = False
        self.alt = False

    def resetSnapKeys(self):
        self.angleSnap = False
        self.gridSnap = False
        self.vertSnap = False

    def initialize(self):
        self.resetMetakeys()
        self.resetSnapKeys()

        self.tm = None
        self.orig = None
        self.snapCo = None
        self.freezeOrient = False

        self.lastSnapTypes = set()

        self.resetSnap()

    def updateSnapKeyMap(self):
        self.snapKeyMap = {}
        kIds = [FTHotKeys.hkSnapVert, FTHotKeys.hkSnapAngle, FTHotKeys.hkSnapGrid]
        varMap = {FTHotKeys.hkSnapVert: 'vertSnap', \
            FTHotKeys.hkSnapAngle: 'angleSnap', FTHotKeys.hkSnapGrid: 'gridSnap'}
        for kId in kIds:
            keys = FTHotKeys.getSnapHotKeys(kId)
            for key in keys:
                self.snapKeyMap[key] = varMap[kId]

    def resetSnap(self): # Called even during isEditing

        self.freeAxes = [] # All axes free
        self.snapDigits.initialize()
        self.rmInfo = None
        self.snapParams = None

        # This variable lets caller know that return was pressed after digits were entered
        # Caller can reset snapper as per convenience
        self.digitsConfirmed = False
        self.lastSelCo = None

        # ~ self.snapStack = [] # TODO

    # Return user selection in header menu, [] for None
    # (Tightly coupled with strings in BezierToolkitParams)
    def getFreeAxesGlobal(self):
        constrAxes = bpy.context.window_manager.bezierToolkitParams.constrAxes
        if(constrAxes != 'NONE'):
            idx = constrAxes.find('-')
            axis = constrAxes[idx + 1]
            freeAxes = [ord(axis) - ord('X')]
            if(idx > 0): freeAxes = sorted(list({0, 1, 2} - set(freeAxes)))
            return freeAxes
        else: return []

    # Return actual free axes (menu + hotkey combined), [] for None
    def getFreeAxesCombined(self):
        constrAxes = self.getFreeAxesGlobal()
        if(len(constrAxes) == 0): return self.freeAxes
        if(len(self.freeAxes) == 0): return constrAxes
        freeAxes = sorted(list(set(constrAxes).intersection(set(self.freeAxes))))
        if(len(freeAxes) == 0): return constrAxes
        else: return freeAxes

    # Return actual free axes (menu + hotkey combined), [0, 1, 2] for None
    def getFreeAxesNormalized(self):
        if(len(self.getFreeAxesCombined()) == 0): return [0, 1, 2]
        else: return self.getFreeAxesCombined()

    def getSnapParamsFreeAxes(self):
        if(self.snapParams != None):
            return self.snapParams.freeAxesN
        else:
            return self.getFreeAxesNormalized()

    def getCurrOrig(self, rmInfo, obj, origType, refLineOrig, selCo):
        if(origType == 'AXIS'):
            if(self.customAxis.length() != 0): return self.customAxis.axisPts[0]
        elif(origType == 'REFERENCE'):
            if(refLineOrig != None): return refLineOrig
        elif(origType == 'CURR_POS'):
            if(selCo != None): return selCo
        elif(origType == 'OBJECT' and obj != None): return obj.location
        elif(origType == 'FACE' and rmInfo != None):
            selObj, location, normal, faceIdx = getSelFaceLoc(rmInfo.region, \
                rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT)
            if(faceIdx != None):
                return selObj.matrix_world @ selObj.data.polygons[faceIdx].center
        elif(origType == 'CURSOR'): return bpy.context.scene.cursor.location
        return Vector((0, 0, 0))

    def getTransMatsForOrient(self, rmInfo, obj, transType, axisScale):

        custAxis = self.customAxis
        if(abs(custAxis.length()) <= DEF_ERR_MARGIN): custAxis = None

        refLine = self.getRefLine()
        currLine = self.getCurrLine()
        if(refLine != None and len(refLine) < 2): refLine = None
        if(currLine != None and len(currLine) < 2): currLine = None

        tmScale = Matrix()
        if(custAxis != None and axisScale == 'AXIS'):
            unitD = custAxis.length() / (custAxis.snapCnt + 1)
            tmScale = Matrix.Scale(1 / unitD, 4)

        elif(refLine != None and axisScale == 'REFERENCE'):
            unitD = (refLine[1] - refLine[0]).length / 10
            if(unitD > DEF_ERR_MARGIN):
                tmScale = Matrix.Scale(1 / unitD, 4)

        tm = None
        if(transType == 'AXIS' and custAxis != None):
            tm, invTm = getLineTransMatrices(custAxis.axisPts[0], custAxis.axisPts[1])

        elif(transType == 'REFERENCE' and refLine != None):
            tm, invTm = getLineTransMatrices(refLine[0], refLine[1])

        elif(transType == 'CURR_POS' and currLine != None):
            tm, invTm = getLineTransMatrices(currLine[0], currLine[1])

        elif(transType == 'VIEW'):
            tm = rmInfo.rv3d.view_matrix

        elif(obj != None and transType == 'OBJECT'):
            tm = (obj.matrix_world).inverted_safe()

        elif(transType == 'FACE'):
            selObj, location, normal, faceIdx = getSelFaceLoc(rmInfo.region, \
                rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT)
            if(faceIdx != None):
                normal = selObj.data.polygons[faceIdx].normal
                quat = normal.to_track_quat('Z', 'X').to_matrix().to_4x4()
                tm = (selObj.matrix_world @ quat).inverted_safe()

        if(tm != None):
            trans, quat, scale = tm.decompose()
            tm = quat.to_matrix().to_4x4() @ tmScale
        else:
            tm = tmScale

        return tm, tm.inverted_safe()

    def isLocked(self):
        return len(self.freeAxes) > 0 or \
            (self.snapDigits.hasVal() and not self.digitsConfirmed)

    def getMetakeys(self):
        return [self.alt, self.ctrl, self.shift]

    # To be called in modal method of parent
    def procEvent(self, context, event):

        # update ctrl etc.
        retValMeta = updateMetaBtns(self, event)

        retValSnap = updateMetaBtns(self, event, self.snapKeyMap)

        if(retValMeta or retValSnap): return EVT_META_OR_SNAP

        metakeys = self.getMetakeys()
        refLineOrig = self.snapParams.refLineOrig if self.snapParams != None \
            else self.getRefLineOrig()

        if(refLineOrig != None):
            snapDProc = self.snapDigits.procEvent(context, event, metakeys)
            if(snapDProc):
                self.digitsConfirmed = False # Always reset if there was any digit entered
                return EVT_CONS

        if(FTHotKeys.isHotKey(FTHotKeys.hkReorient, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                self.freezeOrient = not self.freezeOrient
            return EVT_CONS

        if(not self.ctrl and event.type in {'X', 'Y', 'Z'}):
            self.digitsConfirmed = False # Always reset if there is any lock axis
            if(event.value == 'RELEASE'):
                self.freeAxes = [ord(event.type) - ord('X')]
                if(self.shift):
                    self.freeAxes = sorted(list({0, 1, 2} - set(self.freeAxes)))

                # if already part of global axes, don't store (no escape needed)
                if(self.getFreeAxesCombined() == self.getFreeAxesGlobal()):
                    self.freeAxes = []

            return EVT_CONS

        retVal = EVT_NOT_CONS
        # Consume escape or return / space only if there's something to process
        if(self.isLocked()):
            retVal = EVT_CONS
            if(event.type == 'RET' or event.type == 'SPACE'):
                if(event.value == 'RELEASE'):
                    # ~ self.resetSnap() # This is the responsibility of the caller
                    self.digitsConfirmed = True # Confirm first time
            elif(event.type == 'ESC'):
                if(event.value == 'RELEASE'): self.resetSnap()
            else:
                retVal = EVT_NOT_CONS

        return retVal

    def getStatusStr(self, unit, invTm, refPt, newPt):
        manualEntry = self.snapDigits.hasVal()

        if(manualEntry and self.snapDigits.polar):
            return self.snapDigits.getDeltaStrPolar()

        axes = self.getSnapParamsFreeAxes()
        diffV = self.snapDigits.getCurrDelta() \
            if manualEntry else (newPt - refPt)

        diffV *= getUnitScale()

        diffVActual = invTm @ diffV

        retStr = ''
        transformed = invTm != Matrix()

        axisDeltaFormat = 'D{axis}: {axisDelta}'
        axisDiffFormat = '{{{axisDiff}}}  '

        for i, d in enumerate(diffV):
            if(i not in axes): continue
            v1 = chr(ord('x') + i)

            if(manualEntry and i == self.snapDigits.axisIdx):
                v2 = self.snapDigits.getCurrDeltaStr()
            else:
                v2 = str(round(d, 4))

            retStr += axisDeltaFormat.format(axis = v1, axisDelta = v2)

            if(transformed):
                v3 = str(round(diffVActual[i], 4))
                retStr += axisDiffFormat.format(axisDiff = v3)

            retStr += '  '

        unitT = ''
        unitA = ''
        if(transformed): unitA = unit
        else: unitT = unit

        totalDeltaFormat = '({totalDelta}{unit})'
        diffVStr = str(round(diffV.length, 4))
        retStr += totalDeltaFormat.format(totalDelta = diffVStr, unit = unitT)

        totalDiffVFormat = '{{{totalDiffV}{unit}}}'
        if(transformed):
            diffVStr = str(round(diffVActual.length, 4))
            retStr += totalDiffVFormat.format(totalDiffV = diffVStr, unit = unitA)

        return retStr

    def getAllSnapLocs(self, snapToAxisLine):
        snapLocs = self.getSnapLocs()
        snapLocs.append(bpy.context.scene.cursor.location)
        snapLocs.append(Vector((0, 0, 0)))

        if(snapToAxisLine):
            snapLocs += self.customAxis.getSnapPts()

        vertCnt = 0
        aos = [bpy.context.object] if bpy.context.object != None else []
        objs = bpy.context.selected_objects + aos
        for obj in objs:
            snapLocs.append(obj.location)
            if(obj.type == 'MESH'):
                if(vertCnt + len(obj.data.vertices) < self.MAX_SNAP_VERT_CNT):
                    snapLocs += [obj.matrix_world @ v.co for v in obj.data.vertices]
                    vertCnt =+ len(obj.data.vertices)
                else:
                    break
        return snapLocs

    def getTMInfoAndOrig(self, rmInfo, transType, origType, freezeOrient, \
        axisScale, refLineOrig, selCo):
        obj = bpy.context.object

        if(self.tm != None and self.freezeOrient and transType == 'FACE'):
            tm, invTm = self.tm, self.tm.inverted_safe()
        else:
            tm, invTm = self.getTransMatsForOrient(rmInfo, obj, transType, axisScale)

        if(self.orig != None and freezeOrient and origType == 'FACE'):
            orig = self.orig
        else:
            orig = self.getCurrOrig(rmInfo, obj, origType, refLineOrig, selCo)

        return tm, invTm, orig

    def get3dLocSnap(self, rmInfo, snapParams = None):

        if(snapParams == None): snapParams = SnapParams(self)

        self.rmInfo = rmInfo
        self.snapParams = snapParams

        refreshStatus = snapParams.refreshStatus
        snapToAxisLine = snapParams.snapToAxisLine
        lastCo1Axis = snapParams.lastCo1Axis

        vertSnap = snapParams.vertSnap
        gridSnap = snapParams.gridSnap
        angleSnap = snapParams.angleSnap

        refLine = snapParams.refLine
        refLineOrig = snapParams.refLineOrig
        selCo = snapParams.selCo

        inEdit = snapParams.inEdit
        hasSel = snapParams.hasSel

        transType = snapParams.transType
        origType = snapParams.origType
        axisScale = snapParams.axisScale

        vec = snapParams.vec
        freeAxesN = snapParams.freeAxesN
        snapToPlane = snapParams.snapToPlane

        self.snapCo = None
        region = rmInfo.region
        rv3d = rmInfo.rv3d

        xy = [rmInfo.xy[0] - snapParams.xyDelta[0], \
            rmInfo.xy[1] - snapParams.xyDelta[1]]

        loc = None

        tm, invTm, orig = self.getTMInfoAndOrig(rmInfo, transType, \
            origType, self.freezeOrient, axisScale, refLineOrig, selCo)

        # Must be done after the call to getTMInfoAndOrig
        if(hasSel): self.freezeOrient = True

        vec = snapParams.vec if (snapParams.vec != None) else orig

        self.lastSnapTypes = set()

        unit = unitMap.get(getUnit())
        if(unit == None): unit = ''

        digitsValid = True
        # ~ freeAxesC = self.getFreeAxesCombined()
        # ~ freeAxesN = self.getFreeAxesNormalized()
        # ~ freeAxesG = self.getFreeAxesGlobal()

        if(FTProps.dispSnapInd or vertSnap):
            #TODO: Called very frequently (store the tree [without duplicating data])
            snapLocs = self.getAllSnapLocs((snapToAxisLine and \
                'AXIS' in {transType, origType, axisScale})) + [orig]

            kd = kdtree.KDTree(len(snapLocs))
            for i, l in enumerate(snapLocs):
                kd.insert(getCoordFromLoc(region, rv3d, l).to_3d(), i)
            kd.balance()

            coFind = Vector(xy).to_3d()
            searchResult = kd.find_range(coFind, FTProps.snapDist)

            if(len(searchResult) != 0):
                co, idx, dist = min(searchResult, key = lambda x: x[2])
                self.snapCo = snapLocs[idx]

        if(vertSnap):
            if(self.snapCo != None):
                loc = self.snapCo
            else:
                selObj, loc, normal, faceIdx = getSelFaceLoc(region, rv3d, xy, \
                    self.MAX_SNAP_FACE_CNT, checkEdge = True)

        if(loc != None):
            loc = tm @ loc
            self.lastSnapTypes.add('loc')
        else:
            loc = region_2d_to_location_3d(region, rv3d, xy, vec)
            loc = tm @ loc

            # TODO: Get gridSnap and angleSnap out of this if
            if((transType != 'GLOBAL' and inEdit) or \
                snapToPlane or gridSnap or \
                    self.snapDigits.hasVal() or \
                        (inEdit and (len(freeAxesN) < 3 or angleSnap))):

                # snapToPlane means global constrain axes selection is a plane
                # ~ if(snapToPlane or refLineOrig == None): refCo = orig
                # ~ else: refCo = refLineOrig

                refCo = tm @ orig

                if(self.snapDigits.hasVal()):
                    delta = self.snapDigits.getCurrDelta()
                    loc = tm @ orig + delta
                    self.lastSnapTypes.add('keyboard')
                else:
                    # Special condition for lock to single axis
                    # ~ if(len(freeAxesN) == 1 and refLineOrig != None):
                        # ~ refCo = tm @ refLineOrig
                    if(len(freeAxesN) == 2):
                        constrAxes = freeAxesN
                        loc = refCo.copy()
                        # Any other two points on the plane
                        ppt1 = loc.copy()
                        ppt1[constrAxes[0]] = loc[constrAxes[0]] + 10
                        ppt2 = loc.copy()
                        ppt2[constrAxes[1]] = loc[constrAxes[1]] + 10
                        ppt3 = ppt2.copy()
                        ppt2[constrAxes[0]] = loc[constrAxes[0]] + 10

                        # Raycast from 2d point onto the plane
                        pt = getPtProjOnPlane(region, rv3d, xy, invTm @ loc, \
                            invTm @ ppt1, invTm @ ppt2, invTm @ ppt3)

                        for axis in constrAxes:
                            # TODO: Better handling of boundary condition
                            if(pt == None or pt[axis] > 1000):
                                loc[axis] = refCo[axis]
                            else: loc[axis] = (tm @ pt)[axis]
                        self.lastSnapTypes.add('axis2')

                    if(len(freeAxesN) == 1):
                        if(lastCo1Axis): refCo = self.lastSelCo #TODO: More testing
                        axis = freeAxesN[0]
                        # Any one point on axis
                        ptOnAxis = refCo.copy()

                        # Convert everything to 2d
                        lastCo2d = getCoordFromLoc(region, rv3d, invTm @ refCo)

                        # Very small distance so that the point is on viewport
                        # TODO: This is not foolproof :(
                        ptOnAxis[axis] += .01
                        ptOnAxis2d = getCoordFromLoc(region, rv3d, invTm @ ptOnAxis)

                        # Find 2d projection (needed)
                        pt2d = geometry.intersect_point_line(xy, lastCo2d, ptOnAxis2d)[0]
                        # Any other 2 points on the plane, on which the axis lies
                        ppt1 = refCo.copy()
                        ppt1[axis] += 10
                        newAxis = [i for i in range(0, 3) if i != axis][0]
                        ppt2 = refCo.copy()
                        ppt2[newAxis] += 10

                        # Raycast from 2d point onto the plane
                        pt = getPtProjOnPlane(region, rv3d, pt2d, \
                            invTm @ refCo, invTm @ ppt1, invTm @ ppt2)
                        loc = refCo.copy()
                        if(pt == None or pt[axis] > 1000):
                            loc[axis] = refCo[axis]
                        else: loc[axis] = (tm @ pt)[axis]
                        self.lastSnapTypes.add('axis1')

                if(not self.snapDigits.hasVal() and gridSnap):
                    if(axisScale in {'AXIS' or 'REFERENCE'}):
                        # Independent of view distance
                        diffV = (loc - refCo)
                        loc = refCo + round(diffV.length) * (diffV / diffV.length)
                    else:
                        rounding = getViewDistRounding(rmInfo.space3d, rv3d)
                        loc = tm @ roundedVect(rmInfo.space3d, invTm @ loc, \
                            rounding, freeAxesN)
                    self.lastSnapTypes.add('grid')

                if(not self.snapDigits.hasVal() and angleSnap and len(refLine) > 0):
                    # ~ freeAxesC = [0, 1, 2] if len(freeAxesC) == 0 else freeAxesC
                    snapStart = tm @ orig
                    actualLoc = loc.copy()

                    #First decide the main movement axis
                    diff = [abs(v) for v in (actualLoc - snapStart)]
                    maxDiff = max(diff)
                    axis = diff.index(maxDiff)

                    loc = snapStart.copy()
                    loc[axis] = actualLoc[axis]

                    snapIncr = 45 / self.angleSnapSteps
                    snapAngles = [radians(snapIncr * a) \
                        for a in range(0, self.angleSnapSteps + 1)]

                    l1 =  actualLoc[axis] - snapStart[axis] #Main axis diff value

                    for i in range(0, 3):
                        if(i != axis and (i in freeAxesN)):
                            l2 =  (actualLoc[i] - snapStart[i]) #Minor axis value
                            angle = abs(atan(l2 / l1)) if l1 != 0 else 0
                            dirn = (l1 * l2) / abs(l1 * l2) if (l1 * l2) != 0 else 1
                            prevDiff = LARGE_NO
                            for j in range(0, len(snapAngles) + 1):
                                if(j == len(snapAngles)):
                                    loc[i] = snapStart[i] + \
                                        dirn * l1 * tan(snapAngles[-1])
                                    break
                                cmpAngle = snapAngles[j]
                                if(abs(angle - cmpAngle) > prevDiff):
                                    loc[i] = snapStart[i] + \
                                        dirn * l1 * tan(snapAngles[j-1])
                                    break
                                prevDiff = abs(angle - cmpAngle)

                    self.lastSnapTypes.add('angle')

        if(refreshStatus and inEdit):
            refPt = tm @ orig
            newPt = loc.copy()
            text = self.getStatusStr(unit, invTm, refPt, newPt)
        else:
            text = None

        self.setStatus(rmInfo.area, text)
        loc = invTm @ loc
        self.lastSelCo = loc
        self.tm = tm
        self.orig = orig

        return loc

    def getEditCoPair(self):
        refLineOrig = self.getRefLineOrig()
        if(self.lastSelCo == None or refLineOrig == None):
            return []
        if(self.snapParams == None):
            origType = bpy.context.window_manager.bezierToolkitParams.origType
            refLineOrig = self.getRefLineOrig()
            selCo = self.getSelCo()
        else:
            origType = self.snapParams.origType
            refLineOrig = self.snapParams.refLineOrig
            selCo = self.snapParams.selCo
        orig = self.getCurrOrig(self.rmInfo, bpy.context.object, origType, \
            refLineOrig, selCo)
        return (self.tm @ orig, self.tm @ self.lastSelCo)

    def setStatus(self, area, text): #TODO Global
        area.header_text_set(text)

    # Tightly coupled with get3dLocSnap
    def updateGuideBatches(self, bglDrawMgr):

        # Default values for resetting the bgl lines and points
        drawAxes = [0, 1, 2]
        axisLineCos = [[], [], []]
        axisLineCols = [[], [], []]

        snapLineCos = []
        snapLineCols = []

        custAxisLineCos = []
        custAxisLineCols = []
        custAxisGradStart = None
        custAxisGradEnd = None

        custAxisPtCos = []
        custAxisPtCols = []

        snapIndPtCos = []
        snapIndPtCols = []

        rmInfo = self.rmInfo

        if(rmInfo != None): # self.snapParams is also not None
            snapParams = self.snapParams
            refLine = snapParams.refLine
            refLineOrig = snapParams.refLineOrig
            selCo = snapParams.selCo
            freeAxesN = snapParams.freeAxesN

            transType = snapParams.transType
            origType = snapParams.origType
            axisScale = snapParams.axisScale
            dispAxes = snapParams.dispAxes

            tm, invTm, orig = self.getTMInfoAndOrig(rmInfo, transType, \
                origType, self.freezeOrient, axisScale, refLineOrig, selCo)

            if(dispAxes and FTProps.dispAxes and ((refLineOrig != None \
                or transType == 'VIEW' or len(freeAxesN) == 1) or (len(freeAxesN) > 0 \
                    and origType != 'REFERENCE'))):
                colors = [(.6, 0.2, 0.2, 1), (0.2, .6, 0.2, 1), (0.2, 0.4, .6, 1)]
                l = 2 * rmInfo.rv3d.view_distance

                if (self.lastSelCo != None and len(freeAxesN) == 1): orig = self.lastSelCo

                refCo = tm @ orig
                for axis in freeAxesN[:2]:
                    axisLineCols[axis] = [colors[axis]]
                    pt1 = refCo.copy()
                    pt2 = refCo.copy()
                    pt1[axis] = l + refCo[axis]
                    pt2[axis] = -l + refCo[axis]
                    axisLineCos[axis] = [invTm @ pt1, invTm @ pt2]

            if(refLineOrig != None and self.lastSelCo != None and \
                (self.angleSnap or ('keyboard' in self.lastSnapTypes \
                    and self.snapDigits.polar))):

                snapLineCos = [orig, self.lastSelCo]
                snapLineCols = [(.4, .4, .4, 1)]
                ptCol = (1, 1, 1, 1)

            if(self.customAxis.length() != 0 and (self.customAxis.inDrawAxis == True or \
                'AXIS' in {transType, origType, axisScale})):
                apts = self.customAxis.axisPts
                custAxisLineCos = [apts[0], apts[1]]
                custAxisLineCols = [(1, 1, 1, 1)]
                custAxisGradStart = .9
                custAxisGradEnd = .3

                custAxisPtCos = self.customAxis.getSnapPts()
                custAxisPtCols = [(1, .4, 0, 1)]

            if(FTProps.dispSnapInd and self.snapCo != None):
                snapIndPtCos = [self.snapCo]
                snapIndPtCols = [FTProps.colHltTip] #[(1, .4, 0, 1)]

        for i, axis in enumerate(drawAxes):
            axisGradStart = .2 if len(axisLineCos[i]) > 0 else None
            axisGradEnd = .9 if len(axisLineCos[i]) > 0 else None
            bglDrawMgr.addLineInfo('snapAxis'+str(axis), FTProps.axisLineWidth, \
                axisLineCols[i], axisLineCos[i], axisGradStart, axisGradEnd)

        bglDrawMgr.addLineInfo('SnapLine', FTProps.axisLineWidth, \
            snapLineCols, snapLineCos)

        bglDrawMgr.addLineInfo('CustAxisLine', FTProps.axisLineWidth, custAxisLineCols, \
            custAxisLineCos, custAxisGradStart, custAxisGradEnd, mid = False)

        bglDrawMgr.addPtInfo('CustAxisPt', FTProps.snapPtSize, custAxisPtCols, \
            custAxisPtCos)

        bglDrawMgr.addPtInfo('SnapPt', FTProps.snapPtSize, snapIndPtCols, snapIndPtCos)


################################## Flexi Tool Classes ##################################
#
#                                         Operator
#                                            |
#                                     ModalBaseFlexiOp
#                                            |
#                          -------------------------------------
#                          |                                   |
#                 ModalDrawBezierOp                   ModalFlexiEditBezierOp
#                          |
#          ---------------------------------
#          |                               |
#    ModalFlexiDrawBezierOp       ModalFlexiDrawGreaseOp

class ModalBaseFlexiOp(Operator):
    running = False
    drawHdlrRef = None
    drawTxtHdlrRef = None
    drawFunc = None
    shader = None
    bglDrawMgr = None
    opObj = None

    pointSize = 4 # For Draw (Marker is of diff size)

    def drawKeyMap():

        regions = [r for area in bpy.context.screen.areas if  area.type == 'VIEW_3D' \
            for r in area.regions if r.type == 'WINDOW']
        maxArea = max(r.width * r.height for r in regions)
        currRegion = [r for r in bpy.context.area.regions if r.type == 'WINDOW'][0]

        # Only display in window with max area
        if(currRegion.width * currRegion.height < maxArea): return
        toolRegion = [r for r in bpy.context.area.regions if r.type == 'TOOLS'][0]

        xOff1 = (10 + toolRegion.width) if(FTProps.keyMapNextToTool) else FTProps.keyMapLocX
        maxWidth = 0
        yStart = 10 if(FTProps.keyMapNextToTool) else FTProps.keyMapLocY
        yOff = yStart

        font_id = 0

        if(ModalBaseFlexiOp.opObj != None and FTProps.showKeyMap):

            descrCol = [0] + list(FTProps.colKeymapText)
            keyCol = [0] + list(FTProps.colKeymapKey)

            toolType = ModalBaseFlexiOp.opObj.getToolType()
            config, labels, keys = FTHotKeys.getHKDispLines(toolType)
            blf.size(font_id, FTProps.keyMapFontSize, 72)

            maxLabelWidth = max(blf.dimensions(font_id, l+'XXXXX')[0] for l in labels)
            xOff2 = xOff1 + maxLabelWidth

            lineHeight = 1.2 * max(blf.dimensions(font_id, l)[1] for l in labels)

            blf.position(font_id, xOff1, yOff, 0)
            blf.color(*descrCol)
            blf.draw(font_id, '*')
            dim = blf.dimensions(font_id, '*')
            blf.position(font_id, xOff1 + dim[0], yOff, 0)
            blf.draw(font_id, ' Indicates Configurable Hot Keys')
            yOff += 1.5 * lineHeight

            for i, label in enumerate(labels):
                blf.color(*descrCol)
                blf.position(font_id, xOff1, yOff, 0)
                blf.draw(font_id, label)
                if(config[i]):
                    dim = blf.dimensions(font_id, label)
                    blf.position(font_id, xOff1 + dim[0], yOff, 0)
                    blf.draw(font_id, '*')
                blf.color(*keyCol)
                blf.position(font_id, xOff2, yOff, 0)
                blf.draw(font_id, keys[i])
                yOff += lineHeight

            maxWidth = maxLabelWidth  + max(blf.dimensions(font_id, l)[0] for l in keys)
            header = toolType.title() + ' Keymap'
            headerX = xOff1 + (maxWidth - blf.dimensions(font_id, header)[0]) / 2
            blf.position(font_id, headerX, yOff + lineHeight * .5, 0)
            blf.color(*descrCol)
            blf.draw(font_id, header)

        mathFnTxts = MathFnDraw.getMathFnTxts()

        mathFnCol = [0] + list(FTProps.colMathFnTxt)
        blf.size(font_id, FTProps.mathFnTxtFontSize, 72)
        blf.color(*mathFnCol)
        lineHeight = blf.dimensions(font_id, 'yX')[1]

        yOff = lineHeight / 2
        if(mathFnTxts != None):
            for t in mathFnTxts:
                mathFnPos = (currRegion.width - blf.dimensions(font_id, t)[0]) / 2
                # ~ mathFnPos = xOff1 + maxWidth + 10
                blf.position(font_id, mathFnPos, yOff, 0)
                blf.draw(font_id, t)
                yOff += 1.5 * lineHeight

    def addDrawHandlers(context):
        hdlr = ModalBaseFlexiOp.opObj.__class__.drawHandler
        ModalBaseFlexiOp.drawHdlrRef = \
            bpy.types.SpaceView3D.draw_handler_add(hdlr, (), "WINDOW", "POST_VIEW")
        ModalBaseFlexiOp.drawTxtHdlrRef = \
            bpy.types.SpaceView3D.draw_handler_add(ModalBaseFlexiOp.drawKeyMap, \
                (), "WINDOW", "POST_PIXEL")

    def removeDrawHandlers():
        if(ModalBaseFlexiOp.drawHdlrRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(ModalBaseFlexiOp.drawHdlrRef, \
                "WINDOW")
            ModalBaseFlexiOp.drawHdlrRef = None
        if(ModalBaseFlexiOp.drawTxtHdlrRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(ModalBaseFlexiOp.drawTxtHdlrRef, \
                "WINDOW")
            ModalBaseFlexiOp.drawTxtHdlrRef = None

    def drawHandlerBase():
        if(ModalBaseFlexiOp.shader != None):
            ModalBaseFlexiOp.bglDrawMgr.redraw()

    def tagRedraw():
        areas = [a for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
        for a in areas:
            a.tag_redraw()

    # Called back after add-on preferences are changed
    def propsChanged():
        if(ModalBaseFlexiOp.opObj != None and \
            ModalBaseFlexiOp.opObj.snapper != None):
            ModalBaseFlexiOp.opObj.snapper.updateSnapKeyMap()
            ModalBaseFlexiOp.hdlColMap ={'FREE': FTProps.colHdlFree, \
                        'VECTOR': FTProps.colHdlVector,  \
                            'ALIGNED': FTProps.colHdlAligned, \
                                'AUTO': FTProps.colHdlAuto}

    def resetDisplayBase():
        ModalBaseFlexiOp.bglDrawMgr.reset()
        ModalBaseFlexiOp.tagRedraw()

    def refreshDisplayBase(segDispInfos, bptDispInfos, snapper):
        areaRegionInfo = getAllAreaRegions()

        updateBezierBatches(ModalDrawBezierOp.bglDrawMgr, segDispInfos, \
            bptDispInfos, areaRegionInfo)

        if(snapper != None):
            snapper.updateGuideBatches(ModalBaseFlexiOp.bglDrawMgr)

        ModalBaseFlexiOp.tagRedraw()

    @persistent
    def loadPostHandler(dummy):
        if(ModalBaseFlexiOp.shader != None):
            ModalBaseFlexiOp.resetDisplayBase()
        ModalBaseFlexiOp.running = False

    @persistent
    def loadPreHandler(dummy):
        ModalBaseFlexiOp.removeDrawHandlers()
        if(ModalBaseFlexiOp.drawFunc != None):
            bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc

    @classmethod
    def poll(cls, context):
        return not ModalBaseFlexiOp.running

    def getToolType(self):
        raise NotImplementedError('Call to abstract method.')

    def preInvoke(self, context, event):
        pass # placeholder

    def subInvoke(self, context, event):
        return {'RUNNING_MODAL'} # placeholder

    def invoke(self, context, event):
        ModalBaseFlexiOp.opObj = self
        ModalBaseFlexiOp.running = True
        self.preInvoke(context, event)
        ModalBaseFlexiOp.addDrawHandlers(context)

        ModalBaseFlexiOp.drawFunc = bpy.types.VIEW3D_HT_tool_header.draw
        bpy.types.VIEW3D_HT_tool_header.draw = drawSettingsFT
        context.space_data.show_region_tool_header = True

        self.snapper = Snapper(context, self.getSnapLocs, \
            self.getRefLine, self.getRefLineOrig, self.getSelCo, \
                self.getCurrLine, self.hasSelection, self.isEditing)

        self.rmInfo = None

        ModalBaseFlexiOp.shader = gpu.shader.from_builtin('3D_SMOOTH_COLOR')
        ModalBaseFlexiOp.bglDrawMgr = BGLDrawMgr(ModalBaseFlexiOp.shader)

        # ~ ModalBaseFlexiOp.shader.bind()
        context.window_manager.modal_handler_add(self)

        ModalBaseFlexiOp.ColGreaseHltSeg = (.3, .3, .3, 1) # Not used

        FTProps.updateProps(None, context)
        FTHotKeys.updateHotkeys(None, context)
        FTHotKeys.initSnapMetaFromPref(context)

        self.clickT, self.pressT = None, None
        self.click, self.doubleClick = False, False

        return self.subInvoke(context, event)

    def modal(self, context, event):

        snapper = self.snapper
        if(event.type == 'WINDOW_DEACTIVATE' and event.value == 'PRESS'):
            snapper.initialize()
            return {'PASS_THROUGH'}

        if(not self.isToolSelected(context)): # Subclass
            self.cancelOp(context)
            return {"CANCELLED"}

        self.click, self.doubleClick = False, False
        if(event.type == 'LEFTMOUSE'):
            if(event.value == 'PRESS'):
                self.pressT = time.time()
            elif(event.value == 'RELEASE'):
                t = time.time()
                if(self.clickT != None and (t - self.clickT) < DBL_CLK_DURN):
                    self.clickT = None
                    self.doubleClick = True
                elif(self.pressT != None and (t - self.pressT) < SNGL_CLK_DURN):
                    self.clickT = t
                    self.click = True
                self.pressT = None


        snapProc = snapper.procEvent(context, event)
        metakeys = snapper.getMetakeys()

        if(FTHotKeys.isHotKey(FTHotKeys.hkToggleKeyMap, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                try:
                    prefs = context.preferences.addons[__name__].preferences
                    prefs.showKeyMap = not prefs.showKeyMap
                except Exception as e:
                    print(e)
                    FTProps.showKeyMap = not FTProps.showKeyMap
            return {'RUNNING_MODAL'}

        # Special condition for case where user has configured a different snap key
        # In such case, pass through mouse clicks, if there's a meta key held down...
        # to allow e.g. alt + LMB to rotate view
        if(snapProc != EVT_CONS and any(metakeys) and \
            not any([snapper.angleSnap, snapper.gridSnap, snapper.vertSnap]) and \
                ((event.type == 'LEFTMOUSE' and event.value in {'PRESS', 'RELEASE'}) \
                    or event.type == 'MOUSEMOVE')):
            return {'PASS_THROUGH'}

        rmInfo = RegionMouseXYInfo.getRegionMouseXYInfo(event, self.exclToolRegion())

        ret = FTMenu.procMenu(self, context, event, rmInfo == None)
        if(ret):
            # Menu displayed on release, so retain metakeys till release
            if(event.value == 'RELEASE'):
                snapper.resetMetakeys()
                snapper.resetSnapKeys()
            return {'RUNNING_MODAL'}

        if(self.isEditing() and self.rmInfo != rmInfo):
            return {'RUNNING_MODAL'}
        if(rmInfo == None):
            return {'PASS_THROUGH'}

        self.rmInfo = rmInfo

        ret = snapper.customAxis.procDrawEvent(context, event, snapper, rmInfo)
        evtCons = (ret or snapProc == EVT_CONS)

        # Ignore all PRESS events if consumed, since action is taken only on RELEASE...
        # ...except 1) wheelup / down where there is no release & 2) snap / meta where...
        # ...refresh is needed even on press
        # TODO: Simplify the condition (Maybe return EVT values from all proc methods)
        if(evtCons and event.value == 'PRESS' and \
            event.type != 'MOUSEMOVE' and \
            not event.type.startswith('WHEEL') and (snapProc != EVT_META_OR_SNAP)):
            return {'RUNNING_MODAL'}

        if(FTHotKeys.isHotKey(FTHotKeys.hkSwitchOut, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                self.cancelOp(context)
                resetToolbarTool()
                return {"CANCELLED"}
            return {'RUNNING_MODAL'}

        evtCons = (evtCons or (snapProc == EVT_META_OR_SNAP))
        return self.subModal(context, event, evtCons)

    def cancelOpBase(self):
        for a in bpy.context.screen.areas:
            if (a.type == 'VIEW_3D'): a.header_text_set(None)

        ModalBaseFlexiOp.removeDrawHandlers()
        ModalBaseFlexiOp.running = False
        bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc
        self.snapper = None
        ModalBaseFlexiOp.opObj = None

    def getSnapLocs(self):
        return self.getSnapLocsImpl()

################################## Flexi Draw Bezier Curve ###############################
#
#                                         BaseDraw
#                                            |
#                                -------------------------------
#                                |                             |
#                           Primitive2DDraw                BezierDraw
#                                |
#                        -----------------------------
#                        |                           |
#               ClosedShapeDraw               MathFnDraw
#                        |
#        ---------------------------------
#        |               |               |
#  RectangleDraw    PolygonDraw     EllipseDraw

class BaseDraw:
    def __init__(self, parent):
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.curvePts = []

    def newPoint(self, loc, ltype, rtype):
        self.curvePts.append([loc, loc, loc, ltype, rtype])

    def setCurvePt(self, idx, pt):
        self.curvePts[idx] = pt

    def popCurvePt(self):
        self.curvePts.pop()

    def setCurvePts(self, curvePts):
        self.curvePts = curvePts

class Primitive2DDraw(BaseDraw):

    def __init__(self, parent):
        super(Primitive2DDraw, self).__init__(parent)
        self.shapeSegCnt = 4
        self.curveObjOrigin = None

    def initialize(self):
        super(Primitive2DDraw, self).initialize()
        self.bbStart = None
        self.bbEnd = None
        self.XYstart = None
        self.freeAxesN = None

    def set2DAxes(self):
        params = bpy.context.window_manager.bezierToolkitParams
        self.freeAxesN = self.parent.snapper.getFreeAxesNormalized()

        if(len(self.freeAxesN) > 2 and \
            params.snapOrient in {'GLOBAL', 'REFERENCE', 'CURR_POS'}):
            self.freeAxesN = getClosestPlaneToView(self.parent.rmInfo.rv3d)

    def getNumSegsLimits(self):
        return 2, 100

    def getShapePts(self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, \
        theta, axisIdxs, z):
        raise NotImplementedError('Call to abstract method.')

    # index: ((incr_key, decr_key), metakey, description)
    dynamicParams = {
        0:(['UP_ARROW', 'DOWN_ARROW'], 'Up (incr) / Down (decr)'),
        1:(['RIGHT_ARROW', 'LEFT_ARROW'], 'Left (incr) / Right (decr)'),
        2:(['PAGE_UP', 'PAGE_DOWN'], 'Page Up (incr) / Page Down (decr)'),
        3:(['W', 'S'], 'W (incr) / S (decr)'),
        4:(['A', 'D'], 'A (incr) / D (decr)'),
        5:(['LEFT_BRACKET', 'RIGHT_BRACKET'], '[ (incr) / ] (decr)'),
    }

    def getParamCnt():
        return len(Primitive2DDraw.dynamicParams)

    def getParamHotKeyDescriptions():
        return [Primitive2DDraw.dynamicParams[p][1] \
            for p in range(len(Primitive2DDraw.dynamicParams))]

    # default definition of all params
    # Format (can be overriden in subclass):
    #       def updateParam0(self, event, rmInfo, isIncr): # 0-5
    #           pass
    for i in range(len(dynamicParams)):
        exec('def updateParam' + str(i) + '(self, event, rmInfo, isIncr):\n\tpass')

    def afterShapeSegCnt(self): # Call back
        pass

    def getCurvePts(self, numSegs, axisIdxs, z = None):
        params = bpy.context.window_manager.bezierToolkitParams
        tm = self.parent.snapper.tm if self.parent.snapper.tm != None else Matrix()

        idx0, idx1, idx2 = axisIdxs
        startAngle = params.drawStartAngle
        sweep = params.drawAngleSweep

        bbEnd = tm @ self.bbEnd

        mode = params.drawObjMode

        if(mode == 'CENTER'):
            diffV = (self.bbEnd - self.bbStart)
            bbStart = tm @ (self.bbStart - diffV)
        else:
            bbStart = tm @ self.bbStart

        cX = (bbEnd[idx0] - bbStart[idx0]) / 2
        cY = (bbEnd[idx1] - bbStart[idx1]) / 2
        if(cX == 0 and cY == 0):
            return None

        if(z == None): z = bbStart[idx2]
        center2d = complex(cX, cY)

        snapOrigin = bpy.context.window_manager.bezierToolkitParams.snapOrigin
        orig = complex(bbStart[idx0], bbStart[idx1])
        self.curveObjOrigin = tm.inverted_safe() @ \
            get3DVector(orig + center2d, axisIdxs, z)

        curvePts = self.getShapePts(mode, numSegs, bbStart, bbEnd, center2d, \
            startAngle, sweep, axisIdxs, z)

        if(curvePts == None): return None

        curvePts = [[tm.inverted_safe() @ p if type(p) != str \
            else p for p in pts] for pts in curvePts]

        return curvePts

    def updateCurvePts(self):
        axisIdxs = self.freeAxesN + sorted(list({0, 1, 2} - set(self.freeAxesN)))

        curvePts = self.getCurvePts(axisIdxs = axisIdxs, \
            numSegs = self.shapeSegCnt)

        if(curvePts != None):
            self.setCurvePts(curvePts)
        else:
            self.setCurvePts([[self.bbStart, self.bbStart, self.bbStart], \
                [self.bbEnd, self.bbEnd, self.bbEnd]])

    def getPtLoc(self):
        parent = self.parent
        rmInfo = parent.rmInfo
        if(self.freeAxesN == None):
            self.set2DAxes()

        return self.parent.snapper.get3dLocSnap(rmInfo, \
            SnapParams(parent.snapper, lastCo1Axis = True, freeAxesN = self.freeAxesN, \
                snapToPlane = (len(self.freeAxesN) == 2)))

    def procDrawEvent(self, context, event, snapProc):
        parent = self.parent
        rmInfo = parent.rmInfo
        metakeys = parent.snapper.getMetakeys()

        if(len(self.curvePts) > 0 and not snapProc):
            if(event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', \
                    'NUMPAD_PLUS', 'NUMPAD_MINUS','PLUS', 'MINUS'}):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return {'RUNNING_MODAL'}
                self.updateSegCount(event, rmInfo, \
                    (event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')))
                return {'RUNNING_MODAL'}

            for i in range(len(self.dynamicParams)):
                hotkeys = self.dynamicParams[i][0]
                if(event.type in hotkeys):
                    if(event.value == 'RELEASE'):
                        exec('self.updateParam' + str(i) + \
                            '(event, rmInfo, (event.type == hotkeys[0]))')
                    return {'RUNNING_MODAL'}

            if(event.type == 'H' or event.type == 'h'):
                if(event.value == 'RELEASE'):
                    ModalDrawBezierOp.h = not ModalDrawBezierOp.h
                    self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
                return {"RUNNING_MODAL"}

        if(self.bbStart != None and (event.type == 'RET' or event.type == 'SPACE')):
            if(event.value == 'RELEASE'):
                self.updateCurvePts()
                self.parent.confirm(context, event, self.curveObjOrigin)
                self.parent.snapper.resetSnap()
                self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                self.parent.initialize() # should not access parent.initialize
                self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):
            if(self.bbStart == None):
                self.set2DAxes()
                loc = self.getPtLoc()
                self.XYstart = self.parent.rmInfo.xy
                self.bbStart = loc
                self.bbEnd = loc
                self.setCurvePts([[loc, loc, loc, 'FREE', 'FREE'], \
                    [loc, loc, loc, 'FREE', 'FREE']])
            else:
                loc = self.getPtLoc()
                self.bbEnd = loc
                self.updateCurvePts()
                parent.confirm(context, event, self.curveObjOrigin)
            return {"RUNNING_MODAL"}

        if (snapProc or event.type == 'MOUSEMOVE'):
            if(self.bbStart != None):
                self.bbEnd = self.getPtLoc()
                self.updateCurvePts()
            parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
            return {"RUNNING_MODAL"}

        return {"PASS_THROUGH"}

    #Reference point for restrict angle or lock axis
    def getRefLine(self):
        if(self.bbStart != None):
            if(self.bbEnd != None):
                return[self.bbStart, self.bbEnd]
            else:
                return [self.bbStart]

        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[0] if len(refLine) > 0 else None

class MathFnDraw(Primitive2DDraw):

    # Prefixes for param names generated dynamically for constants
    startPrefix = 'mathFnStart_'
    incrPrefix = 'mathFnIncr_'

    mathFnFileExt = 'mfn'

    # Default values
    defFnType = 'XY'
    defFNRes = 10

    defFNXYName = ''#'Sine Function'
    defFNXYDescr = ''#'Sine wave with A: Amplitude, B: Frequency, C: Phase shift'
    defFnXY = 'sin(x)'#'A * sin(B * x + C)'

    defFnParam1 = ''#'B * cos(t + C)'
    defFnParam2 = ''#'A * sin(t + D)'
    defTMapTo = 'HORIZONTAL'
    defTScale = 3
    defTStart = 0

    defXYMap = 'NORMAL_XY'
    defClipVal = 10

    defConstStart = 1
    defConstIncr = 0.1

    # XML tags / attributes
    xDocTag = "drawMathFn"
    xFnName = 'name'
    xFnDescr = 'descr'
    xFnCurveRes = 'curveRes'
    xFnType = 'type'
    xFns = 'functions'
    xEquation = 'equation'

    xXYFn = 'xyFn'
    xClipVal = 'clipValue'

    xParamFn1 = 'paramFn1'
    xParamFn2 = 'paramFn2'
    xTMapTo = 'tMapTo'
    xTScaleFact = 'tScaleFact'
    xTStart = 'tStartVal'

    xConstPrefix = 'constant_'
    xValue = 'value'
    xIncr = 'incr'

    mathFnDirty = True
    mathFnItems = None
    mathFnNoSelItem = ('NO_SEL', '', 'Function not saved')

    def __init__(self, parent, star = False):
        super(MathFnDraw, self).__init__(parent)
        params = bpy.context.window_manager.bezierToolkitParams
        self.shapeSegCnt = params.mathFnResolution

    def getNumSegsLimits(self):
        return 2, 99999

    def updateSegCount(self, event, rmInfo, isIncr):
        minSegs, maxSegs = self.getNumSegsLimits()
        params = bpy.context.window_manager.bezierToolkitParams
        if(isIncr and params.mathFnResolution < maxSegs): params.mathFnResolution += 1
        if(not isIncr and params.mathFnResolution > minSegs): params.mathFnResolution -= 1
        self.shapeSegCnt = params.mathFnResolution
        self.afterShapeSegCnt()
        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
        return True

    for i in range(len(Primitive2DDraw.dynamicParams)):
        fnStr = 'def updateParam' + str(i) + '(self, event, rmInfo, isIncr):\n\t' + \
            'params = bpy.context.window_manager.bezierToolkitParams\n\t' + \
            'incr = params.' + incrPrefix + str(i) + '\n\t' + \
            'params.'+ startPrefix + str(i) + ' += incr if(isIncr) else -incr\n\t' + \
            'self.updateCurvePts()\n\t' + \
            'self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)\n\t' + \
            'areas = [a for a in bpy.context.screen.areas]\n\t' + \
            'for a in areas:\n\t\t' + \
            'a.tag_redraw()\n\t' + \
            'return True\n\t'
        exec(fnStr)

    def afterShapeSegCnt(self):
        bpy.context.window_manager.bezierToolkitParams.mathFnResolution = self.shapeSegCnt
        areas = [a for a in bpy.context.screen.areas]
        for a in areas:
            a.tag_redraw()

    def testFn(expr, var): # var = 'x' or 'y'
        exec(var + ' = 1')
        # ~ exec(var.upper() + ' = 1') ??
        try:
            eval(expr)
            return True
        except Exception as e:
            return False

    # Returns None in case of invalid function
    def isInverted(expr): # var = 'x' or 'y'
        if(not MathFnDraw.testFn(expr, 'x')):
            if(not MathFnDraw.testFn(expr, 'y')): return None
            else: return True

        return False

    def getEvaluatedExpr(expr):
        params = bpy.context.window_manager.bezierToolkitParams
        for j in range(Primitive2DDraw.getParamCnt()):
            expr = expr.replace(str(chr(ord('A') + j)), \
                str(round(eval('params.'+ MathFnDraw.startPrefix + str(j)), 4)))
        return expr

    def getShapePts(self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, \
        theta, axisIdxs, z):

        params = bpy.context.window_manager.bezierToolkitParams
        curvePts = []
        idx0, idx1, idx2 = axisIdxs
        self.shapeSegCnt = params.mathFnResolution # Synch with params

        if(params.mathFnType == 'PARAMETRIC'):
            fn1 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric1)
            fn2 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric2)
            if(not MathFnDraw.testFn(fn1, 't')): return curvePts
            if(not MathFnDraw.testFn(fn2, 't')): return curvePts

            # Let the plot be always centered
            if(mode == 'CENTER'):
                bbStart[idx0] += center2d.real
                bbStart[idx1] += center2d.imag

            scaleFact = params.drawMathTScaleFact

            if(params.drawMathTMapTo in {'x', 'y', 'xy'}):
                xSpan = (bbEnd[idx0] - bbStart[idx0])
                ySpan = (bbEnd[idx1] - bbStart[idx1])
            else:
                xSpan = (self.parent.rmInfo.xy[0] - self.XYstart[0])
                ySpan = (self.parent.rmInfo.xy[1] - self.XYstart[1])
                scaleFact = scaleFact / 25 # Device dependent!

            if(params.drawMathTMapTo in {'X', 'HORIZONTAL'}):
                span = scaleFact * xSpan
            elif(params.drawMathTMapTo in {'Y', 'VERTICAL'}):
                span = scaleFact * ySpan
            else:
                span = scaleFact * sqrt(xSpan * xSpan + ySpan * ySpan)

            intervals = int(self.shapeSegCnt * abs(span))
            if(intervals == 0): intervals = self.shapeSegCnt
            incr = span / intervals

            t = params.drawMathTStart
            for step in range(intervals):
                try:
                    x = bbStart[idx0] + eval(fn1)
                    y = bbStart[idx1] + eval(fn2)
                    pt2d = complex(x, y)
                    pt = get3DVector(pt2d, axisIdxs, z)
                    curvePts.append([pt, pt, pt, 'VECTOR', 'VECTOR'])
                except Exception as e:
                    print(e, fn1, fn2)
                    pass

                t += incr

        else:
            expr = MathFnDraw.getEvaluatedExpr(params.drawMathFn)
            clip = params.mathFnclipVal

            inverted = MathFnDraw.isInverted(expr)

            if(inverted == None):
                return curvePts

            span = (bbEnd[idx1] - bbStart[idx1]) if(inverted) \
                else (bbEnd[idx0] - bbStart[idx0])

            intervals = int(self.shapeSegCnt * abs(span))
            if(intervals == 0): intervals = self.shapeSegCnt
            incr = span / intervals
            indep = bbStart[idx0] if(not inverted) else bbStart[idx1]

            for step in range(intervals):
                try:
                    if(inverted): y = indep
                    else: x = indep
                    dep = eval(expr)
                    if(abs(dep) > clip):
                        dep = clip * (dep / abs(dep))
                    if(inverted): x = dep + bbStart[idx0] + (center2d.real \
                        if(mode == 'CENTER') else 0)
                    else: y = dep + bbStart[idx1] + (center2d.imag \
                        if(mode == 'CENTER') else 0)
                    pt2d = complex(x, y)
                    pt = get3DVector(pt2d, axisIdxs, z)
                    curvePts.append([pt, pt, pt, 'VECTOR', 'VECTOR'])
                except Exception as e:
                    print(e, expr)
                    pass

                indep += incr

        return curvePts

    def getMathFnTxts():
        INVALID = '<Invalid Equation>'
        fnTxts = None
        params = bpy.context.window_manager.bezierToolkitParams
        # TODO: Should be somewhere else
        tool = bpy.context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)
        if(tool.idname == FlexiDrawBezierTool.bl_idname and params.drawObjType == 'MATH'):

            if(params.mathFnType == 'PARAMETRIC'):
                fn1 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric1)
                fn2 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric2)
                fnTxts = ['y = ' + (fn2 if MathFnDraw.testFn(fn2, 't') else INVALID), \
                    'x = '+ (fn1 if MathFnDraw.testFn(fn1, 't') else INVALID)]
            else:
                expr = MathFnDraw.getEvaluatedExpr(params.drawMathFn)
                inverted = MathFnDraw.isInverted(expr)
                if(inverted == None):
                    fnTxts = [INVALID]
                elif(inverted):
                    fnTxts = ['x = ' + expr]
                else:
                    fnTxts = ['y = ' + expr]
        return fnTxts

    def getMathFnFolder(create = True):
        userPath = bpy.utils.resource_path('USER')
        configPath = os.path.join(userPath, "config")
        mathFnFolder = configPath + '/mathFunctions'
        if(create and not os.path.isdir(mathFnFolder)):
            os.makedirs(mathFnFolder)
        return mathFnFolder

    def getMathFnList(dummy1 = None, dummy2 = None):
        if(MathFnDraw.mathFnItems == None or MathFnDraw.mathFnDirty):
            mathFnFolder = MathFnDraw.getMathFnFolder()
            fNames = sorted([fName for fName in os.listdir(mathFnFolder)
                if fName.endswith('.' + MathFnDraw.mathFnFileExt)], key=lambda s: s.lower())
            MathFnDraw.mathFnItems = [MathFnDraw.mathFnNoSelItem]
            for fName in fNames:
                with open(mathFnFolder + '/' + fName) as f:
                    doc = minidom.parse(f)
                fnName = doc.documentElement.getAttribute(MathFnDraw.xFnName)
                fnDescr = doc.documentElement.getAttribute(MathFnDraw.xFnDescr)
                MathFnDraw.mathFnItems.append((fnName, fnName, fnDescr))
            MathFnDraw.mathFnDirty = False

        return MathFnDraw.mathFnItems

    def refreshDefaultParams():
        params = bpy.context.window_manager.bezierToolkitParams
        # ~ params.mathFnDescr = MathFnDraw.defFNXYDescr
        # ~ params.mathFnResolution = MathFnDraw.defFNRes
        # ~ params.mathFnType = MathFnDraw.defFnType
        # ~ params.drawMathFn = MathFnDraw.defFnXY
        # ~ params.mathFnclipVal = MathFnDraw.defClipVal

        for i in range(Primitive2DDraw.getParamCnt()):
            char = chr(ord('A') + i)
            exec('params.' + MathFnDraw.startPrefix + str(i) + ' = ' + \
                str(MathFnDraw.defConstStart))
            exec('params.' + MathFnDraw.incrPrefix + str(i) + ' = ' + \
                str(MathFnDraw.defConstIncr))

    def refreshParamsFromFile(dummy1 = None, dummy2 = None):
        params = bpy.context.window_manager.bezierToolkitParams
        mathFnSel = params.mathFnList
        if(mathFnSel == MathFnDraw.mathFnNoSelItem[0]):
            refreshDefaultParams()
            return
        mathFnFolder = MathFnDraw.getMathFnFolder()
        filepath = mathFnFolder + '/' + mathFnSel + '.' + MathFnDraw.mathFnFileExt

        with open(filepath) as f:
            doc = minidom.parse(f)

        docElem = doc.documentElement

        fnName = docElem.getAttribute(MathFnDraw.xFnName)
        fnDescr = docElem.getAttribute(MathFnDraw.xFnDescr)
        fnType = docElem.getAttribute(MathFnDraw.xFnType)
        fnCurveRes = float(docElem.getAttribute(MathFnDraw.xFnCurveRes))

        params.mathFnName = fnName
        params.mathFnDescr = fnDescr
        params.mathFnType = fnType
        params.mathFnResolution = fnCurveRes

        fnElem = docElem.getElementsByTagName(MathFnDraw.xFns)[0]

        if(fnType == 'PARAMETRIC'):
            fnTMapTo = fnElem.getAttribute(MathFnDraw.xTMapTo)
            params.drawMathTMapTo = fnTMapTo

            fnTScaleFact = float(fnElem.getAttribute(MathFnDraw.xTScaleFact))
            params.drawMathTScaleFact = fnTScaleFact

            fnTStart = float(fnElem.getAttribute(MathFnDraw.xTStart))
            params.drawMathTStart = fnTStart

            elem = fnElem.getElementsByTagName(MathFnDraw.xParamFn1)[0]
            paramFn1 = elem.getAttribute(MathFnDraw.xEquation)
            params.drawMathFnParametric1 = paramFn1

            elem = fnElem.getElementsByTagName(MathFnDraw.xParamFn2)[0]
            paramFn2 = elem.getAttribute(MathFnDraw.xEquation)
            params.drawMathFnParametric2 = paramFn2
        else:
            fnYClip = float(fnElem.getAttribute(MathFnDraw.xClipVal))
            params.mathFnclipVal = fnYClip

            elem = fnElem.getElementsByTagName(MathFnDraw.xXYFn)[0]
            xyFn = elem.getAttribute(MathFnDraw.xEquation)
            params.drawMathFn = xyFn

        for i in range(Primitive2DDraw.getParamCnt()):
            char = chr(ord('A') + i)

            elem = fnElem.getElementsByTagName(MathFnDraw.xConstPrefix + char)[0]
            startVal = float(elem.getAttribute(MathFnDraw.xValue))
            incrVal = float(elem.getAttribute(MathFnDraw.xIncr))
            exec('params.' + MathFnDraw.startPrefix + str(i) + ' = ' + str(startVal))
            exec('params.' + MathFnDraw.incrPrefix + str(i) + ' = ' + str(incrVal))

        areas = [a for a in bpy.context.screen.areas]
        for a in areas:
            a.tag_redraw()


class ResetMathFn(Operator):
    bl_idname = "object.reset_math_fn"
    bl_label = "Reset"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        try:
            MathFnDraw.refreshParamsFromFile()
        except:
            MathFnDraw.refreshDefaultParams()

        return {'FINISHED'}

# TODO: Better validation and error handling
class SaveMathFn(Operator):
    bl_idname = "object.save_math_fn"
    bl_label = "Save"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        fnName = params.mathFnName
        fnDescr = params.mathFnDescr
        fnType = params.mathFnType
        fnCurveRes = params.mathFnResolution
        if(fnName.strip() == ''): # Error Condition #1
            return {'FINISHED'}

        doc = minidom.getDOMImplementation().createDocument(None, MathFnDraw.xDocTag, None)
        docElem = doc.documentElement

        docElem.setAttribute(MathFnDraw.xFnName, fnName)
        docElem.setAttribute(MathFnDraw.xFnDescr, fnDescr)
        docElem.setAttribute(MathFnDraw.xFnType, fnType)
        docElem.setAttribute(MathFnDraw.xFnCurveRes, str(fnCurveRes))

        fnElem = doc.createElement(MathFnDraw.xFns)
        docElem.appendChild(fnElem)
        if(fnType == 'PARAMETRIC'):
            paramFn1 = params.drawMathFnParametric1
            paramFn2 = params.drawMathFnParametric2

            if(paramFn1.strip() == '' or paramFn2.strip() == ''): # Error Condition #2
                return {'FINISHED'}

            fnTMapTo = params.drawMathTMapTo
            fnElem.setAttribute(MathFnDraw.xTMapTo, str(fnTMapTo))

            fnTScaleFact = params.drawMathTScaleFact
            fnElem.setAttribute(MathFnDraw.xTScaleFact, str(fnTScaleFact))

            fnTStart = params.drawMathTStart
            fnElem.setAttribute(MathFnDraw.xTStart, str(fnTStart))

            elem = doc.createElement(MathFnDraw.xParamFn1)
            elem.setAttribute(MathFnDraw.xEquation, paramFn1)
            fnElem.appendChild(elem)

            elem = doc.createElement(MathFnDraw.xParamFn2)
            elem.setAttribute(MathFnDraw.xEquation, paramFn2)
            fnElem.appendChild(elem)

        else:
            xyFn = params.drawMathFn
            if(xyFn.strip() == ''): # Error Condition #3
                return {'FINISHED'}

            elem = doc.createElement(MathFnDraw.xXYFn)
            elem.setAttribute(MathFnDraw.xEquation, xyFn)
            fnElem.appendChild(elem)

            fnYClip = params.mathFnclipVal
            fnElem.setAttribute(MathFnDraw.xClipVal, str(fnYClip))

        for i in range(Primitive2DDraw.getParamCnt()):
            char = chr(ord('A') + i)

            startVal = round(eval('params.' + MathFnDraw.startPrefix + str(i)), 4)
            incrVal = round(eval('params.' + MathFnDraw.incrPrefix + str(i)), 4)
            elem = doc.createElement(MathFnDraw.xConstPrefix + char)
            elem.setAttribute(MathFnDraw.xValue, str(startVal))
            elem.setAttribute(MathFnDraw.xIncr, str(incrVal))
            fnElem.appendChild(elem)

        mathFnFolder = MathFnDraw.getMathFnFolder()
        fnFile = mathFnFolder + '/' + fnName + '.' + MathFnDraw.mathFnFileExt

        try:
            with open(fnFile,"w") as f:
                doc.writexml(f)

            MathFnDraw.mathFnDirty = True
            params.mathFnList = fnName
        except: # Error Condition #4
            self.report({'ERROR'}, 'Error saving math function file')

        return {'FINISHED'}


# TODO: Better validation and error handling
class LoadMathFn(Operator):
    bl_idname = "object.load_math_fn"
    bl_label = "Load"
    bl_options = {'REGISTER', 'UNDO'}

    filter_glob : StringProperty(default = '*.' + MathFnDraw.mathFnFileExt, options={'HIDDEN'})
    filepath : StringProperty(subtype='FILE_PATH')

    def execute(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        try:
            with open(self.filepath) as f:
                doc = minidom.parse(f)
            fnName = doc.documentElement.getAttribute(MathFnDraw.xFnName)
            mathFnFolder = MathFnDraw.getMathFnFolder()
            destPath = mathFnFolder + '/' + fnName + '.' + MathFnDraw.mathFnFileExt
            copyfile(self.filepath, destPath)
            MathFnDraw.mathFnDirty = True
            params.mathFnList = fnName
        except:
            self.report({'ERROR'}, 'Error importing math function file')

        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class DeleteMathFn(bpy.types.Operator):
    bl_idname = "object.delete_math_fn"
    bl_label = "Remove function from list and delete it permanently?"
    bl_options = {'REGISTER', 'INTERNAL'}

    # ~ @classmethod
    # ~ def poll(cls, context):
        # ~ return True

    def execute(self, context):
        params = bpy.context.window_manager.bezierToolkitParams
        fnName = params.mathFnList
        if(fnName != MathFnDraw.mathFnNoSelItem[0]):
            mathFnFolder = MathFnDraw.getMathFnFolder()
            fnFile = mathFnFolder + '/' + fnName + '.' + MathFnDraw.mathFnFileExt
            try:
                os.remove(fnFile)
                MathFnDraw.mathFnDirty = True
                params.mathFnList = MathFnDraw.mathFnNoSelItem[0]
            except:
                self.report({'ERROR'}, 'Error deleting math function file')

        return {'FINISHED'}

    def invoke(self, context, event):
        params = bpy.context.window_manager.bezierToolkitParams
        fnName = params.mathFnList
        if(fnName != MathFnDraw.mathFnNoSelItem[0]):
            return context.window_manager.invoke_props_dialog(self)
        else:
            return {'FINISHED'}

class ClosedShapeDraw(Primitive2DDraw):
    def __init__(self, parent, star = False):
        super(ClosedShapeDraw, self).__init__(parent)

    def updateSegCount(self, event, rmInfo, isIncr):
        minSegs, maxSegs = self.getNumSegsLimits()
        if(isIncr and self.shapeSegCnt < maxSegs): self.shapeSegCnt += 1
        if(not isIncr and self.shapeSegCnt > minSegs): self.shapeSegCnt -= 1
        self.afterShapeSegCnt()
        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
        return True

    def updateParam1(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams
        theta = params.drawAngleSweep

        if(isIncr):
            if(theta >= -10 and theta <= 0): params.drawAngleSweep = 10
            elif(theta > 350): params.drawAngleSweep = -350
            else: params.drawAngleSweep = theta + 10
        else:
            if(theta <= 10 and theta >= 0): params.drawAngleSweep = -10
            elif(theta < -350): params.drawAngleSweep = 350
            else: params.drawAngleSweep = theta - 10
        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
        return True

class RectangleDraw(ClosedShapeDraw):
    def __init__(self, parent, star = False):
        super(RectangleDraw, self).__init__(parent)

    def getShapePts(self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, \
        theta, axisIdxs, z):
        idx0, idx1, idx2 = axisIdxs
        pt0 = complex(bbStart[idx0], bbStart[idx1])
        pt1 = complex(bbEnd[idx0], bbStart[idx1])
        pt2 = complex(bbEnd[idx0], bbEnd[idx1])
        pt3 = complex(bbStart[idx0], bbEnd[idx1])

        curvePts = []
        for pt2d in [pt0, pt1, pt2, pt3, pt0]:
            pt = get3DVector(pt2d, axisIdxs, z)
            curvePts.append([pt, pt, pt, 'VECTOR', 'VECTOR'])

        return curvePts

class PolygonDraw(ClosedShapeDraw):

    def updateParam0(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams
        offset = params.drawStarOffset

        params.drawStarOffset += 0.1 if(isIncr) else -0.1

        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
        return True

    def getNumSegsLimits(self):
        return 3, 100

    def __init__(self, parent, star = False):
        super(PolygonDraw, self).__init__(parent)
        self.star = star

    def getShapePts(self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, \
        theta, axisIdxs, z):
        params = bpy.context.window_manager.bezierToolkitParams
        params.drawSides = numSegs
        idx0, idx1, idx2 = axisIdxs
        cX, cY = center2d.real, center2d.imag
        radius = sqrt(cX * cX + cY * cY)
        offset = params.drawStarOffset
        orig = complex(bbStart[idx0], bbStart[idx1])

        if(cX == 0): startAngle = 90 * ((cY / abs(cY)) if cY != 0 else 1)
        else: startAngle = degrees(atan(cY / cX))
        if(cX < 0): startAngle += 180
        thetaIncr = theta / numSegs
        curvePts = []
        for segCnt in range(numSegs + 1):
            if(self.star):
                angle = startAngle + thetaIncr * segCnt - thetaIncr / 2
                shift = complex(offset * radius * cos(radians(angle)), \
                    offset * radius * sin(radians(angle)))
                pt = get3DVector(orig + center2d + shift, axisIdxs, z)
                curvePts.append([pt, pt, pt, 'VECTOR', 'VECTOR'])
            if(not self.star or segCnt < numSegs):
                angle = startAngle + thetaIncr * segCnt
                shift = complex(radius * cos(radians(angle)), radius * sin(radians(angle)))
                pt = get3DVector(orig + center2d + shift, axisIdxs, z)
                curvePts.append([pt, pt, pt, 'VECTOR', 'VECTOR'])

        return curvePts

class EllipseDraw(ClosedShapeDraw):

    # https://math.stackexchange.com/questions/22064/calculating-a-point-that-lies-on-an-ellipse-given-an-angle
    def getPtAtAngle(a, b, theta):
        if (theta < 0): theta += 2 * pi
        denom = sqrt(b * b + a * a * tan(theta) * tan(theta))
        num = (a * b)
        x = num / denom
        if(pi / 2 < theta <= 3 * pi / 2): x = -x
        y = x * tan(theta)
        return complex(x, y)

    def __init__(self, parent):
        super(EllipseDraw, self).__init__(parent)

    def updateParam0(self, event, rmInfo, isIncr):
        params = bpy.context.window_manager.bezierToolkitParams
        theta = params.drawStartAngle

        if(isIncr):
            if(theta > 350): params.drawStartAngle = 0
            else: params.drawStartAngle = theta + 10
        else:
            if(theta < -350): params.drawStartAngle = 0
            else: params.drawStartAngle = theta - 10

        self.updateCurvePts()
        self.parent.redrawBezier(rmInfo, hdlPtIdxs = {}, hltEndSeg = False)
        return True

    def getShapePts(self, mode, numSegs, bbStart, bbEnd, center2d, startAngle, \
        theta, axisIdxs, z):

        idx0, idx1, idx2 = axisIdxs
        cX, cY = center2d.real, center2d.imag

        if(cX == 0 or cY == 0):
            return None

        radius = complex(cX, cY) # Actually same as center2d
        orig = complex(bbStart[idx0], bbStart[idx1])

        large_arc = 0
        rotation = 0

        sweep = 1
        startIdx = 0
        endIdx = 1
        rvs = False

        if(theta < 0):
            sweep = 0
            startIdx = 1
            endIdx = 0
            rvs = True

        pt1 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(startAngle))
        a1 = startAngle + theta / 2
        if(a1 > 360): a1 = a1 - 360
        if(a1 < -360): a1 = a1 + 360
        a2 = startAngle + theta
        if(a2 > 360): a2 = a2 - 360
        if(a2 < -360): a2 = a2 + 360
        pt2 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(a1))
        pt3 = EllipseDraw.getPtAtAngle(abs(cX), abs(cY), radians(a2))

        pt1 = orig + pt1 + center2d
        pt2 = orig + pt2 + center2d
        pt3 = orig + pt3 + center2d

        endPts = [pt1, pt2]
        segs1 = getSegsForArc(endPts[startIdx], radius, 1, endPts[endIdx], \
            10, axisIdxs, z)

        endPts = [pt2, pt3]
        segs2 = getSegsForArc(endPts[startIdx], radius, 1, endPts[endIdx], \
            10, axisIdxs, z)

        segElems = [segs1, segs2]
        segs = segElems[startIdx] + segElems[endIdx]

        curvePts = getWSDataForSegs(segs)
        if(len(curvePts) < 2): return None

        pts = getInterpBezierPts(curvePts, subdivPerUnit = 100, segLens = None)
        if(len(pts) < 2): return None

        vertCos = getInterpolatedVertsCo(pts, numSegs)

        # TODO: A more efficient approach for dividing ellipse uniformly
        newSegs = []
        for i in range(1, len(vertCos)):
            segStart = complex(vertCos[i-1][idx0], vertCos[i-1][idx1])
            segEnd = complex(vertCos[i][idx0], vertCos[i][idx1])
            endPts = [segStart, segEnd]
            segs1 = getSegsForArc(endPts[startIdx], radius, sweep, \
                endPts[endIdx], 1, axisIdxs, z)

            newSegs += segs1

        if(rvs): newSegs = reversed(newSegs)
        curvePts = getWSDataForSegs(newSegs)

        if(len(curvePts) < 2): return None

        if(vectCmpWithMargin(curvePts[0][1], curvePts[-1][1])):
            ldiffV = curvePts[0][1] - curvePts[0][0]
            rdiffV = curvePts[-1][2] - curvePts[-1][1]
            if(vectCmpWithMargin(ldiffV, rdiffV)):
                curvePts[0][3] = 'ALIGNED'
                curvePts[0][4] = 'ALIGNED'
                curvePts[-1][3] = 'ALIGNED'
                curvePts[-1][4] = 'ALIGNED'

        return curvePts

class BezierDraw(BaseDraw):
    def __init__(self, parent):
        super(BezierDraw, self).__init__(parent)
        self.reset()

    def reset(self):
        self.capture = False
        self.grabRepos = False
        self.dissociateHdl = False

    def moveBezierPt(self, loc):
        if(len(self.curvePts) > 0):
            pt = self.curvePts[-1][:]
            self.setCurvePt(-1, [loc, loc, loc, pt[3], pt[4]])

    def movePointByDelta(self, delta):
        if(len(self.curvePts) > 0):
            pt = self.curvePts[-1][:]
            pt[1] += delta
            self.setCurvePt(-1, pt)

    def moveBptElem(self, handle, loc):
        idx = {'left':0, 'pt':1, 'right':2}[handle]
        if(len(self.curvePts) > 0):
            pt = self.curvePts[-1][:]
            pt[idx] = loc
            self.setCurvePt(-1, pt)

    def moveBptElemByDelta(self, handle, delta):
        idx = {'left':0, 'pt':1, 'right':2}[handle]
        if(len(self.curvePts) > 0):
            pt = self.curvePts[-1][:]
            pt[idx] += delta
            self.setCurvePt(-1, pt)

    def resetHandle(self, handle):
        if(len(self.curvePts) > 0):
            self.moveBptElem(handle, self.curvePts[-1][1])

    def isHandleSet(self):
        if(len(self.curvePts) == 0): return False
        co = self.curvePts[-1][1]
        lh = self.curvePts[-1][0]
        rh = self.curvePts[-1][2]
        if(not vectCmpWithMargin(co, lh) or not vectCmpWithMargin(co, rh)): return True
        return False

    def procDrawEvent(self, context, event, snapProc):
        rmInfo = self.parent.rmInfo
        snapper = self.parent.snapper
        metakeys = snapper.getMetakeys()

        if(self.capture and FTHotKeys.isHotKey(FTHotKeys.hkGrabRepos, \
            event.type, metakeys)):
            if(event.value == 'RELEASE'):
                self.dissociateHdl = False
                self.grabRepos = not self.grabRepos
            return {"RUNNING_MODAL"}

        if(self.capture and FTHotKeys.isHotKey(FTHotKeys.hkDissociateHdl, \
            event.type, metakeys)):
            if(event.value == 'RELEASE'):
                self.grabRepos = False
                self.dissociateHdl = not self.dissociateHdl
            return {"RUNNING_MODAL"}

        # This can happen only when space was entered and something was there
        # for Snapper to process
        if (snapProc and snapper.digitsConfirmed):
            snapper.resetSnap()

            # Because resetSnap sets this to False (TODO: Refactor resetSnap)
            snapper.digitsConfirmed = True

            # First space / enter is equivalent to mouse press without release
            if(not self.capture):
                self.capture = True
                snapper.setStatus(rmInfo.area, None)
                return {'RUNNING_MODAL'}
            else:
                # Second space / enter means it should be processed here,
                # set snapProc to False so this modal will process it
                snapProc = False

        if(not snapProc and event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                if(self.grabRepos):
                    self.grabRepos = False
                elif(self.dissociateHdl):
                    self.dissociateHdl = False
                elif(self.capture and self.isHandleSet()):
                    self.resetHandle('left')
                    self.resetHandle('right')

                    # Needed to indicate next space / entered to be processed here
                    snapper.digitsConfirmed = True
                    snapper.setStatus(rmInfo.area, None)
                else:
                    self.parent.initialize() # TODO: should not access parent.initialize
                self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and \
            FTHotKeys.isHotKey(FTHotKeys.hkResetLastHdl, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                if(len(self.curvePts) > 1):
                    if(len(self.curvePts) > 2):
                        idx = len(self.curvePts) - 2
                        pt = self.curvePts[idx][:]
                        pt[2] = pt[1].copy()
                        self.setCurvePt(idx, pt)
                    self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and \
            FTHotKeys.isHotKey(FTHotKeys.hkUndoLastSeg, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                snapper.resetSnap()
                if(not self.capture):
                    if(len(self.curvePts) > 0):
                        self.popCurvePt()

                #Because there is an extra point (the current one)
                if(len(self.curvePts) <= 1):
                    self.initialize()
                    self.reset()
                else:
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc)
                self.capture = False
                self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and (event.type == 'RET' or event.type == 'SPACE')):
            if(event.value == 'RELEASE'):
                if(snapper.digitsConfirmed):
                    self.reset()
                    snapper.digitsConfirmed = False
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.newPoint(loc, 'ALIGNED', 'ALIGNED')
                    self.parent.redrawBezier(rmInfo)
                else:
                    if(len(self.curvePts) > 0): self.popCurvePt()
                    self.parent.confirm(context, event)
                    snapper.resetSnap()
                    self.parent.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            if(len(self.curvePts) == 0):
                loc = snapper.get3dLocSnap(rmInfo)
                self.newPoint(loc, 'ALIGNED', 'ALIGNED')

            # Special condition for hot-key single axis lock (useful)
            if(len(snapper.freeAxes) == 1 and len(self.curvePts) > 1):
                snapper.resetSnap()

            if(self.capture): self.parent.pressT = None # Lock capture
            else: self.capture = True
            return {'RUNNING_MODAL'}

        if (not snapProc and event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):
            # ~ if(snapper.isLocked()):
                # ~ if(len(self.curvePts) == 1):
                    # ~ self.moveBptElem('right', \
                        # ~ snapper.get3dLocSnap(rmInfo))# changes only rt handle
                # ~ return {'RUNNING_MODAL'}

            # See the special condition above in event.value == 'PRESS'
            # ~ if(len(snapper.freeAxes) > 1):
            snapper.resetSnap()
            self.reset()

            # Rare condition: This happens e. g. when user clicks on header menu
            # like Object->Transform->Move. These ops consume press event but not release
            # So update the snap locations anyways if there was some transformation
            if(len(self.curvePts) == 0):
                self.parent.updateSnapLocs() # Subclass (TODO: have a relook)

            elif(self.parent.doubleClick):
                if(len(self.curvePts) > 0): self.popCurvePt()
                self.parent.confirm(context, event)
                self.parent.redrawBezier(rmInfo)

            else:
                if(self.parent.click):
                    loc = self.curvePts[-1][1]
                    self.moveBptElem('left', loc)
                    self.moveBptElem('right', loc)
                else:
                    loc = snapper.get3dLocSnap(rmInfo)

                # ~ if(len(self.curvePts) == 1):
                    # ~ self.moveBptElem('right', loc)# changes only rt handle

                self.newPoint(loc, 'ALIGNED', 'ALIGNED')
                self.parent.redrawBezier(rmInfo)

            return {'RUNNING_MODAL'}

        # Refresh also in case of snapper events
        # except when digitsConfirmed (to give user opportunity to draw a straight line)
        # ~ if ((snapProc and not snapper.digitsConfirmed) \
        if (snapProc or event.type == 'MOUSEMOVE'):

            bpy.context.window.cursor_set("DEFAULT")
            hdlPtIdxs = None
            if(len(self.curvePts) > 0):
                if(self.capture):
                    lastPt = self.curvePts[-1][:]
                    if(self.grabRepos):
                        rtHandle = lastPt[2].copy()
                        xy2 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, lastPt[1])
                        xy1 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, rtHandle)
                        loc = snapper.get3dLocSnap(rmInfo, SnapParams(snapper, \
                            xyDelta = [xy1[0] - xy2[0], xy1[1] - xy2[1]]))
                        delta = loc - lastPt[1]
                        self.moveBptElemByDelta('pt', delta)
                        self.moveBptElemByDelta('left', delta)
                        self.moveBptElemByDelta('right', delta)
                    else:
                        loc = snapper.get3dLocSnap(rmInfo)
                        delta = (loc - lastPt[1])
                        if(self.dissociateHdl):
                            lastPt[3] = 'FREE'
                            lastPt[4] = 'FREE'
                        else:
                            lastPt[0] = lastPt[1] - delta
                            lastPt[3] = 'ALIGNED'
                            lastPt[4] = 'ALIGNED'
                        lastPt[2] = lastPt[1] + delta
                        self.setCurvePt(-1, lastPt)
                    hdlPtIdxs = {len(self.curvePts) - 1}
                else:
                    loc = snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc)
                    hdlPtIdxs = {len(self.curvePts) - 2}

            self.parent.redrawBezier(rmInfo, hdlPtIdxs = hdlPtIdxs)
            return {'RUNNING_MODAL'}

        return {'PASS_THROUGH'} if not snapProc else {'RUNNING_MODAL'}

    def getRefLine(self):
        if(len(self.curvePts) > 0):
            idx = 0
            if(self.capture):
                if(self.grabRepos and len(self.curvePts) > 1): idx = -2
                else: idx = -1
            # There should always be min 2 pts if not capture, check anyway
            elif(len(self.curvePts) > 1):
                idx = -2
            if((len(self.curvePts) + (idx - 1)) >= 0):
                return[self.curvePts[idx-1][1], self.curvePts[idx][1]]
            else:
                return[self.curvePts[idx][1]]
        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[-1] if len(refLine) > 0 else None


class ModalDrawBezierOp(ModalBaseFlexiOp):

    # Static members shared by flexi draw and flexi grease
    markerSize = 8
    h = False

    drawObjMap = {}

    #static method
    def drawHandler():
        ModalBaseFlexiOp.drawHandlerBase()

    def updateDrawType(dummy, context):
        opObj = ModalDrawBezierOp.opObj
        if(opObj != None):
            opObj.setDrawObj()
            opObj.initialize()
            opObj.resetDisplay()
            ModalDrawBezierOp.updateDrawSides(dummy, context)

    def updateDrawSides(dummy, context):
        params = bpy.context.window_manager.bezierToolkitParams
        opObj = ModalDrawBezierOp.opObj
        if(opObj != None and opObj.drawType != 'BEZIER'):
            opObj.drawObj.shapeSegCnt = params.drawSides

    def getToolType(self):
        return TOOL_TYPE_FLEXI_DRAW

    def resetDisplay(self):
        ModalBaseFlexiOp.resetDisplayBase()

    def setDrawObj(self):
        params = bpy.context.window_manager.bezierToolkitParams
        self.drawObj = ModalDrawBezierOp.drawObjMap[params.drawObjType]
        self.drawType = params.drawObjType

    def preInvoke(self, context, event):
        self.bezierDrawObj = BezierDraw(self)
        self.rectangleDrawObj = RectangleDraw(self)
        self.ellipseDrawObj = EllipseDraw(self)
        self.polygonDrawObj = PolygonDraw(self)
        self.starDrawObj = PolygonDraw(self, star = True)
        self.mathDrawObj = MathFnDraw(self)

        ModalDrawBezierOp.drawObjMap = \
            {'BEZIER': self.bezierDrawObj, \
             'RECTANGLE': self.rectangleDrawObj, \
             'ELLIPSE': self.ellipseDrawObj, \
             'POLYGON': self.polygonDrawObj, \
             'STAR': self.starDrawObj, \
             'MATH': self.mathDrawObj, \
            }
        self.setDrawObj()

    #This will be called multiple times not just at the beginning
    def initialize(self):
        self.drawObj.initialize()
        self.snapper.initialize()

    def subInvoke(self, context, event):

        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)

        try:
            ModalDrawBezierOp.markerSize = \
                context.preferences.addons[__name__].preferences.markerSize
        except Exception as e:
            # ~ print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
            ModalDrawBezierOp.markerSize = 8

        self.updateDrawType(context)
        self.updateDrawSides(context)
        return {"RUNNING_MODAL"}

    def cancelOp(self, context):
        self.resetDisplay()
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)
        return self.cancelOpBase()

    def postUndoRedo(self, scene, dummy = None): # signature different in 2.8 and 2.81?
        self.updateSnapLocs() # subclass method

    def confirm(self, context, event, location = None):
        metakeys = self.snapper.getMetakeys()
        shift = self.snapper.angleSnap # Overloaded key op
        autoclose = (self.drawType == 'BEZIER' and shift and \
            (event.type == 'SPACE' or event.type == 'RET'))
        self.save(context, event, autoclose, location)
        self.resetDisplay()
        self.initialize()

    def exclToolRegion(self):
        return True

    def isEditing(self):
        return len(self.drawObj.curvePts) > 0

    def hasSelection(self):
        return self.isEditing()

    # Common subModal for Flexi Draw and Flexi Grease
    def baseSubModal(self, context, event, snapProc):
        return self.drawObj.procDrawEvent(context, event, snapProc)

    def refreshMarkerPos(self, rmInfo):
        colMap = self.getColorMap()
        colMarker = colMap['MARKER_COLOR']
        markerLoc = self.snapper.get3dLocSnap(rmInfo)

        self.resetDisplay()
        self.bglDrawMgr.addPtInfo('drawMarker', ModalDrawBezierOp.markerSize, \
            [colMarker], [markerLoc])

        ModalBaseFlexiOp.refreshDisplayBase(segDispInfos = [], bptDispInfos = [], \
            snapper = self.snapper)

    def redrawBezier(self, rmInfo, lastSegOnly = False, hdlPtIdxs = None, \
        hltEndSeg = True):
        curvePts = self.drawObj.curvePts

        ptCnt = len(curvePts)

        if(ptCnt == 0):
            self.refreshMarkerPos(rmInfo)
            return

        self.bglDrawMgr.resetPtInfo('drawMarker')

        colMap = self.getColorMap()
        colSelSeg = colMap['SEL_SEG_COLOR']
        colNonAdjSeg = colMap['NONADJ_SEG_COLOR']
        colTip = colMap['TIP_COLOR']
        colEndTip = colMap['ENDPT_TIP_COLOR']

        segColor = colSelSeg
        tipColors = [colTip, colEndTip, colTip] if (not ModalDrawBezierOp.h) \
            else [None, colEndTip, None]

        segDispInfos = []
        bptDispInfos = []

        if(hdlPtIdxs == None): hdlPtIdxs = {ptCnt - 2} # Default last but one
        elif(len(hdlPtIdxs) == 0): hdlPtIdxs = range(ptCnt)

        for hdlPtIdx in hdlPtIdxs:
            bptDispInfos.append(BptDisplayInfo(curvePts[hdlPtIdx], tipColors, \
                handleNos = [0, 1] if (not ModalDrawBezierOp.h) else []))

        startIdx = 0
        if(lastSegOnly and ptCnt > 1):
            startIdx = ptCnt - 2
        for i in range(startIdx, ptCnt - 1):
            if(not hltEndSeg or i == ptCnt - 2): segColor = colSelSeg
            else: segColor = colNonAdjSeg
            segDispInfos.append(SegDisplayInfo([curvePts[i], curvePts[i+1]], segColor))

        ModalBaseFlexiOp.refreshDisplayBase(segDispInfos, bptDispInfos, self.snapper)

    def getRefLine(self):
        return self.drawObj.getRefLine()

    def getRefLineOrig(self):
        return self.drawObj.getRefLineOrig()

    def getSelCo(self):
        return self.getRefLineOrig()

    def getCurrLine(self):
        return self.getRefLine()

class ModalFlexiDrawBezierOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves in object mode"
    bl_idname = "wm.flexi_draw_bezier_curves"
    bl_label = "Flexi Draw Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}

    def __init__(self):
        pass

    # For some curve-changing ops (like reset rotation); possible in draw
    def updateAfterGeomChange(self, scene = None, dummy = None): # 3 params in 2.81
        self.updateSnapLocs()

    def isToolSelected(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)

        if(tool == None or tool.idname != FlexiDrawBezierTool.bl_idname):
        # if(tool == None or tool.idname != 'flexi_bezier.draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {'SEL_SEG_COLOR': FTProps.colDrawSelSeg,
        'NONADJ_SEG_COLOR': FTProps.colDrawNonHltSeg,
        'TIP_COLOR': FTProps.colHdlPtTip,
        'ENDPT_TIP_COLOR': FTProps.colBezPt,
        'MARKER_COLOR': FTProps.colDrawMarker}

    def cancelOp(self, context):
        bpy.app.handlers.depsgraph_update_post.remove(self.updateAfterGeomChange)
        super(ModalFlexiDrawBezierOp, self).cancelOp(context)

    def preInvoke(self, context, event):
        super(ModalFlexiDrawBezierOp, self).preInvoke(context, event)
        # If the operator is invoked from context menu, enable the tool on toolbar
        if(not self.isToolSelected(context) and context.mode == 'OBJECT'):
            bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname)
            # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')

        # Object name -> [spline index, [pts]]
        # Not used right now (maybe in case of large no of curves)
        self.snapInfos = {}
        self.updateSnapLocs()
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()

        if(FTHotKeys.isHotKey(FTHotKeys.hkToggleDrwEd, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                bpy.ops.wm.tool_set_by_id(name = FlexiEditBezierTool.bl_idname)
                # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.edit_tool')
            return {"RUNNING_MODAL"}

        return self.baseSubModal(context, event, snapProc)

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if(len(self.drawObj.curvePts) > 0):
            locs += [pt[1] for pt in self.drawObj.curvePts[:-1]]
            # ~ locs += [self.curvePts[-1][0], self.curvePts[-1][1], self.curvePts[-1][2]]
        return locs

    def updateSnapLocs(self, objNames = None):
        updateCurveEndPtMap(self.snapInfos, addObjNames = objNames)

    def createCurveObj(self, context, startObj = None, \
        startSplineIdx = None, endObj = None, endSplineIdx = None, autoclose = False):
        # First create the new curve
        collection = context.collection
        if(collection == None):
            collection = context.scene.collection
        obj = createObjFromPts(self.drawObj.curvePts, '3D', collection, autoclose)

        # Undo stack in case the user does not want to join
        if(endObj != None or startObj != None):
            obj.select_set(True)
            # ~ bpy.context.view_layer.objects.active = obj
            bpy.ops.ed.undo_push()
        else:
            return obj

        endObjs = []

        # Connect the end curve (if exists) first
        splineIdx = endSplineIdx

        if(endObj != None and startObj != endObj):
            # first separate splines
            endObjs, changeCnt = splitCurve([endObj], split = 'spline', newColl = False)

            # then join the selected spline from end curve with new obj
            obj = joinSegs([endObjs[splineIdx], obj], optimized = True, \
                straight = False, srcCurve = endObjs[splineIdx])

            endObjs[splineIdx] = obj

            #Use this if there is no start curve
            objComps = endObjs

        if(startObj != None):
            # Repeat the above process for start curve
            startObjs, changeCnt = splitCurve([startObj], split = 'spline', \
                newColl = False)

            obj = joinSegs([startObjs[startSplineIdx], obj], \
                optimized = True, straight = False, srcCurve = startObjs[startSplineIdx])

            if(startObj == endObj and startSplineIdx != endSplineIdx):
                # If startSplineIdx == endSplineIdx the join call above would take care
                # but if they are different they need to be joined with a separate call
                obj = joinSegs([startObjs[endSplineIdx], obj], \
                    optimized = True, straight = False, \
                        srcCurve = startObjs[endSplineIdx])

                # can't replace the elem with new obj as in case of end curve
                # (see the seq below)
                startObjs.pop(endSplineIdx)
                if(endSplineIdx < startSplineIdx):
                    startSplineIdx -= 1

            # Won't break even if there were no endObjs
            objComps = startObjs[:startSplineIdx] + endObjs[:splineIdx] + \
                [obj] + endObjs[(splineIdx + 1):] + startObjs[(startSplineIdx + 1):]

        obj = joinCurves(objComps)

        if(any(p.co.z != 0 for s in obj.data.splines for p in s.bezier_points)):
            obj.data.dimensions = '3D'

        return obj

    #TODO: At least store in map instead of linear search
    def getSnapObjs(self, context, locs):
        retVals = [[None, 0, 0]] * len(locs)
        foundVals = 0
        for obj in bpy.data.objects:
            if(isBezier(obj)):
                mw = obj.matrix_world
                for i, s in enumerate(obj.data.splines):
                    if(s.use_cyclic_u or len(s.bezier_points) == 0): continue
                    for j, loc in enumerate(locs):
                        p = s.bezier_points[0]
                        if(vectCmpWithMargin(loc, mw @ p.co)):
                            retVals[j] = [obj, i, 0]
                            foundVals += 1
                        p = s.bezier_points[-1]
                        if(vectCmpWithMargin(loc, mw @ p.co)):
                            retVals[j] = [obj, i, -1]
                            foundVals += 1
                        if(foundVals == len(locs)):
                            return retVals
        return retVals

    def save(self, context, event, autoclose, location, align = True):
        curvePts = self.drawObj.curvePts
        if(len(curvePts) > 1):

            startObj, startSplineIdx, ptIdx2, endObj, endSplineIdx, ptIdx1 = \
                [x for y in self.getSnapObjs(context, [curvePts[0][1],
                    curvePts[-1][1]]) for x in y]

            ctrl = self.snapper.gridSnap # Overloaded key op

            # ctrl pressed and there IS a snapped end obj,
            # so user does not want connection

            # (no option to only connect to starting curve when end object exists)
            if(ctrl and endObj != None):
                obj = self.createCurveObj(context, autoclose = False)
            else:
                startObjName = startObj.name if(startObj != None) else ''
                endObjName = endObj.name if(endObj != None) else ''

                obj = self.createCurveObj(context, startObj, startSplineIdx, endObj, \
                    endSplineIdx, autoclose)

            if(align and startObj == None and endObj == None):
                alignToNormal(obj)
                bpy.context.evaluated_depsgraph_get().update()
                if(location == None): location = getObjBBoxCenter(obj)

            if(location != None):
                shiftOrigin(obj, location)
                obj.location = location
                bpy.context.evaluated_depsgraph_get().update()

            params = bpy.context.window_manager.bezierToolkitParams
            copyProperties(params.copyPropsObj, obj)

            #TODO: Why try?
            try:
                obj.select_set(True)
                # ~ bpy.context.view_layer.objects.active = obj
                self.updateSnapLocs([obj.name, startObjName, endObjName])
            except Exception as e:
                pass
        bpy.ops.ed.undo_push()

################### Flexi Draw Grease Bezier ###################

class ModalFlexiDrawGreaseOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves as grease pencil strokes"
    bl_idname = "wm.flexi_draw_grease_bezier_curves"
    bl_label = "Flexi Draw Grease Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}

    h = False

    def getToolType(self):
        return TOOL_TYPE_FLEXI_GREASE

    def __init__(self):
        # ~ curveDispRes = 200
        # ~ super(ModalFlexiDrawGreaseOp, self).__init__(curveDispRes)
        pass

    def isToolSelected(self, context):
        if(context.mode != 'PAINT_GPENCIL'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('PAINT_GPENCIL', \
            create = False)

        if(tool == None or tool.idname != FlexiGreaseBezierTool.bl_idname): 
        # if(tool == None or tool.idname != 'flexi_bezier.grease_draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {'SEL_SEG_COLOR': FTProps.colGreaseSelSeg,
        'NONADJ_SEG_COLOR': ModalBaseFlexiOp.ColGreaseHltSeg, #Not used
        'TIP_COLOR': FTProps.colHdlPtTip,
        'ENDPT_TIP_COLOR': FTProps.colGreaseBezPt,
        'MARKER_COLOR': FTProps.colGreaseMarker, }

    def preInvoke(self, context, event):
        super(ModalFlexiDrawGreaseOp, self).preInvoke(context, event)
        # If the operator is invoked from context menu, enable the tool on toolbar
        if(not self.isToolSelected(context) and context.mode == 'PAINT_GPENCIL'):
            bpy.ops.wm.tool_set_by_id(name = FlexiGreaseBezierTool.bl_idname)
            # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.grease_draw_tool')

        o = context.object
        if(o == None or o.type != 'GPENCIL'):
            d = bpy.data.grease_pencils.new('Grease Pencil Data')
            o = bpy.data.objects.new('Grease Pencil', d)
            context.scene.collection.objects.link(o)
        self.gpencil = o

        self.subdivCos = []
        self.interpPts = []
        self.subdivPerUnit = None

    # overridden
    def redrawBezier(self, rmInfo, hdlPtIdxs = None, hltEndSeg = True):
        curvePts = self.drawObj.curvePts
        ptCnt = len(curvePts)
        subdivCos = self.subdivCos if ptCnt > 1 else []
        self.bglDrawMgr.addLineInfo('gpSubdivLines', FTProps.lineWidth, \
            [FTProps.colGreaseNonHltSeg], getLinesFromPts(subdivCos))
        if(ModalFlexiDrawGreaseOp.h):
            self.bglDrawMgr.resetPtInfo('gpSubdivPts')
        else:
            self.bglDrawMgr.addPtInfo('gpSubdivPts', \
                FTProps.greaseSubdivPtSize, [FTProps.colGreaseSubdiv], subdivCos)

        if(self.drawType != 'BEZIER' and len(curvePts) > 0):
            ModalBaseFlexiOp.refreshDisplayBase(segDispInfos = [], bptDispInfos = [], \
                snapper = self.snapper)
        else:
            super(ModalFlexiDrawGreaseOp, self).redrawBezier(rmInfo, lastSegOnly = True, \
                hdlPtIdxs = {ptCnt-1}, hltEndSeg = hltEndSeg)

    def initialize(self):
        super(ModalFlexiDrawGreaseOp, self).initialize()
        self.subdivCos = []
        self.interpPts = []
        self.updateSnapLocs()

    def subModal(self, context, event, snapProc):
        curvePts = self.drawObj.curvePts
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()

        if(not metakeys[2]):
            cntIncr = 5 #if(self.isDrawShape) else 5

            if(self.drawType in {'BEZIER', 'ELLIPSE'} and \
                event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', 'NUMPAD_PLUS', \
                'NUMPAD_MINUS','PLUS', 'MINUS'} and len(curvePts) > 1):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return {'RUNNING_MODAL'}
                elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                    self.subdivAdd(cntIncr)
                elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                    self.subdivAdd(-cntIncr)

                self.redrawBezier(rmInfo)
                return {'RUNNING_MODAL'}

        if(event.type == 'H' or event.type == 'h'):
            if(event.value == 'RELEASE'):
                ModalFlexiDrawGreaseOp.h = not ModalFlexiDrawGreaseOp.h
                self.redrawBezier(self.rmInfo)
            return {"RUNNING_MODAL"}

        ptCnt = len(curvePts)

        retVal = self.baseSubModal(context, event, snapProc)

        newPtCnt = len(curvePts)
        # ~ if(newPtCnt - ptCnt != 0):
        if(len(curvePts) > 0):
            if(self.subdivPerUnit == None):
                viewDist = context.space_data.region_3d.view_distance
                self.initSubdivPerUnit = 5000.0 / viewDist # TODO: default configurable?
                self.subdivPerUnit = 0.02 * self.initSubdivPerUnit
                self.snapLocs.append(curvePts[0][1])
            if(len(curvePts) > 1):
                slens = self.getCurveSegLens()
                self.updateInterpPts(slens)
                self.updateSubdivCos(sum(slens))
                self.redrawBezier(rmInfo)

        return retVal

    def getCurveSegLens(self):
        clen = []
        curvePts = self.drawObj.curvePts

        for i in range(1, len(curvePts) - 1):
            clen.append(getSegLen([curvePts[i-1][1], curvePts[i-1][2], \
                curvePts[i][0], curvePts[i][1]]))
        return clen

    def updateSubdivCos(self, clen = None):
        if(self.drawType in {'POLYGON', 'STAR', 'RECTANGLE'}):
            self.subdivCos = [p[0] for p in self.drawObj.curvePts]
        elif(self.interpPts != []):
            if(clen == None): clen = sum(self.getCurveSegLens())
            cnt = round(self.subdivPerUnit * clen)
            if(cnt > 0):
                self.subdivCos = getInterpolatedVertsCo(self.interpPts, cnt)#[1:-1]
                return
            self.subdivCos = []

    def updateInterpPts(self, slens):
        curvePts = self.drawObj.curvePts[:] if(self.drawType != 'BEZIER') else \
            self.drawObj.curvePts[:-1]

        self.interpPts = getInterpBezierPts(curvePts, self.initSubdivPerUnit, slens)

        return self.interpPts

    def subdivAdd(self, addCnt):
        slens = self.getCurveSegLens()
        clen = sum(slens)
        if(clen == 0): return
        cnt = self.subdivPerUnit * clen + addCnt
        if(cnt < 1): cnt = 1

        self.subdivPerUnit = (cnt / clen)
        self.updateSubdivCos(clen)

    def getSnapLocsImpl(self):
        return self.snapLocs

    def updateSnapLocs(self):
        self.snapLocs = []
        gpencils = [o for o in bpy.data.objects if o.type == 'GPENCIL']
        for gpencil in gpencils:
            mw = gpencil.matrix_world
            for layer in gpencil.data.layers:
                for f in layer.frames:
                    for s in f.strokes:
                        if(len(s.points) > 0): # Shouldn't be needed, but anyway...
                            self.snapLocs += [mw @ s.points[0].co, mw @ s.points[-1].co]

    def save(self, context, event, autoclose, location):
        layer = self.gpencil.data.layers.active
        if(layer == None):
            layer = self.gpencil.data.layers.new('GP_Layer', set_active = True)
        if(len(layer.frames) == 0):
            layer.frames.new(0)
        frame = layer.frames[-1]

        invMw = self.gpencil.matrix_world.inverted_safe()
        if(len(self.subdivCos) > 0):
            brush = context.scene.tool_settings.gpencil_paint.brush
            lineWidth = brush.size
            strength = brush.gpencil_settings.pen_strength

            stroke = frame.strokes.new()
            stroke.display_mode = '3DSPACE'
            stroke.points.add(count = len(self.subdivCos))
            for i in range(0, len(self.subdivCos)):
                pt = self.subdivCos[i]
                stroke.points[i].co = self.gpencil.matrix_world.inverted_safe() @ pt
                stroke.points[i].strength = strength
            if(autoclose):
                stroke.points.add(count = 1)
                stroke.points[-1].co = stroke.points[0].co.copy()
                stroke.points[-1].strength = strength
            stroke.line_width = lineWidth
            self.snapLocs += [self.subdivCos[0][1], self.subdivCos[-1][1]]
        bpy.ops.ed.undo_push()

################### Flexi Edit Bezier Curve ###################

class EditSegDisplayInfo(SegDisplayInfo):

    def __init__(self, segPts, segColor, subdivCos):
        super(EditSegDisplayInfo, self).__init__(segPts, segColor)
        self.subdivCos = subdivCos

# fromMix True: points after shape key value / eval_time applied
def getBptData(obj, withShapeKey = True, shapeKeyIdx = None, fromMix = True, \
    updateDeps = False, local = False):
    # Less readable but more convenient than class
    # Format: [handle_left, co, handle_right, handle_left_type, handle_right_type]
    worldSpaceData = []
    mw = Matrix() if local else obj.matrix_world

    keydata = None
    dataIdx = 0
    shapeKey = obj.active_shape_key
    tmpsk = None
    if(withShapeKey and shapeKey != None):
        if(shapeKeyIdx == None):
            shapeKeyIdx = obj.active_shape_key_index
        if(fromMix):
            if(not obj.data.shape_keys.use_relative):
                val = obj.data.shape_keys.eval_time
            else:
                val = obj.data.shape_keys.key_blocks[obj.active_shape_key_index].value

            if(floatCmpWithMargin(val, 0)):
                keyBlock = obj.data.shape_keys.key_blocks[0]
            else:
                tmpsk = obj.shape_key_add(name = 'tmp', from_mix = True)
                keyBlock = obj.data.shape_keys.key_blocks[tmpsk.name]
        else:
            keyBlock = obj.data.shape_keys.key_blocks[shapeKeyIdx]
        keydata = keyBlock.data

    for spline in obj.data.splines:
        pts = []
        for pt in spline.bezier_points:
            lt, rt = pt.handle_left_type, pt.handle_right_type
            if(keydata != None):
                pt = keydata[dataIdx]
                dataIdx += 1

            pts.append([mw @ pt.handle_left, mw @ pt.co, mw @ pt.handle_right, lt, rt])
        worldSpaceData.append(pts)
    if(tmpsk != None):
        obj.shape_key_remove(tmpsk)
        obj.active_shape_key_index = shapeKeyIdx
        if(updateDeps): bpy.context.evaluated_depsgraph_get().update()
    return worldSpaceData

#TODO: splineIdx not needed if ptCnt given
def getAdjIdx(obj, splineIdx, startIdx, offset = 1, ptCnt = None):
    spline = obj.data.splines[splineIdx]
    if(ptCnt == None):
        ptCnt = len(spline.bezier_points)
    if(not spline.use_cyclic_u and
        ((startIdx + offset) >= ptCnt or (startIdx + offset) < 0)):
            return None
    return (ptCnt + startIdx + offset) % ptCnt # add ptCnt for negative offset

def getBezierDataForSeg(obj, splineIdx, segIdx, withShapeKey = True, shapeKeyIdx = None, \
    fromMix = True, updateDeps = False):
    wsData = getBptData(obj, withShapeKey, shapeKeyIdx, fromMix, updateDeps)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)
    if(segEndIdx == None):
        return []
    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0, pt1]

def getSegPtsInSpline(wsData, splineIdx, ptIdx, cyclic):
    splinePts = wsData[splineIdx]
    if(ptIdx < (len(splinePts) - 1) ): ptRange = [ptIdx, ptIdx + 1]
    elif(ptIdx == (len(splinePts) - 1)  and cyclic):
        ptRange = [-1, 0]
        if(splinePts[-1][4] == 'VECTOR'):
            splinePts[-1][2] = (splinePts[-1][1] + \
                1/3 * (splinePts[-1][1] - splinePts[0][1]))
        if(splinePts[0][3] == 'VECTOR'):
            splinePts[0][0] = (splinePts[0][1] + \
                1/3 * (splinePts[0][1] - splinePts[-1][1]))
    else: return []

    return [[splinePts[x][i] for i in range(5)] for x in ptRange]

def getInterpSegPts(wsData, splineIdx, ptIdx, cyclic, res, maxRes):
    segPts = getSegPtsInSpline(wsData, splineIdx, ptIdx, cyclic)
    areaRegionInfo = getAllAreaRegions() # TODO: To be passed from caller

    return getPtsAlongBezier2D(segPts, areaRegionInfo, res, maxRes)

# Wrapper for spatial search within segment
def getClosestPt2dWithinSeg(region, rv3d, coFind, selObj, selSplineIdx, selSegIdx, \
    withHandles, withBezPts):
    infos = {selObj: {selSplineIdx:[[selSegIdx],[]]}}

    # set selObj in objs for CurveBezPts
    return getClosestPt2d(region, rv3d, coFind, [selObj], infos, withHandles, \
        withBezPts, withObjs = False, maxSelObjRes = MAX_SEL_CURVE_RES)

def getClosestPt2d(region, rv3d, coFind, objs, selObjInfos, withHandles = True, \
    withBezPts = True, withObjs = True, maxSelObjRes = MAX_NONSEL_CURVE_RES, \
        withShapeKey = True):

    objLocMap = {}

    objLocList = [] # For mapping after search returns
    objInterpLocs = []
    objInterpCounts = []
    objBezPtCounts = []

    objSplineEndPts = []

    for obj in objs:
        #TODO: Check of shape key bounding box
        if(obj.active_shape_key == None and \
            not isPtIn2dBBox(obj, region, rv3d, coFind, FTProps.snapDist)):
            continue

        wsDataSK = None
        # Curve data with shape key value applied (if shape key exists)
        wsData = getBptData(obj, fromMix = True, updateDeps = True)
        if(withShapeKey and obj.active_shape_key != None):
            # active shape key data with value = 1
            wsDataSK = getBptData(obj, fromMix = False)

        for i, spline in enumerate(obj.data.splines):
            for j, pt in enumerate(spline.bezier_points):
                objLocList.append([obj, i, j])
                if(withObjs):
                    interpLocs = \
                        getInterpSegPts(wsData, i, j, spline.use_cyclic_u, \
                            res = SEARCH_CURVE_RES, maxRes = MAX_NONSEL_CURVE_RES)[1:-1]
                    if(wsDataSK != None):
                        interpLocs += \
                            getInterpSegPts(wsDataSK, i, j, spline.use_cyclic_u, \
                                res = SEARCH_CURVE_RES, \
                                    maxRes = MAX_NONSEL_CURVE_RES)[1:-1]

                    objInterpLocs += interpLocs
                    objInterpCounts.append(len(interpLocs))

                if(withBezPts):
                    cnt = 1
                    objSplineEndPts.append(wsData[i][j][1])#mw @ pt.co)
                    if(wsDataSK != None):
                        objSplineEndPts.append(wsDataSK[i][j][1])#mw @ pt.co)
                        cnt += 1
                    objBezPtCounts.append(cnt)

    selObjLocList = [] # For mapping after search returns
    selObjHdlList = [] # Better to create a new one, even if some redundancy

    segInterpLocs = []
    selObjInterpCounts = []
    selObjHdlCounts = []

    hdls = []

    for selObj in selObjInfos.keys():
        wsDataSK = None
        # Curve data with shape key value applied (if shape key exists)
        wsData = getBptData(selObj, fromMix = True, updateDeps = True)
        if(withShapeKey and selObj.active_shape_key != None):
            # active shape key data with value = 1
            wsDataSK = getBptData(selObj, fromMix = False)

        info = selObjInfos[selObj]
        for splineIdx in info.keys():
            cyclic = selObj.data.splines[splineIdx].use_cyclic_u
            segIdxs = info[splineIdx][0]
            for segIdx in segIdxs:
                selObjLocList.append([selObj, splineIdx, segIdx])
                interpLocs = getInterpSegPts(wsData, splineIdx, segIdx, cyclic, \
                        res = SEARCH_CURVE_RES * 5, maxRes = maxSelObjRes)[1:-1]
                if(wsDataSK != None):
                    interpLocs += \
                        getInterpSegPts(wsDataSK, splineIdx, segIdx, cyclic, \
                        res = SEARCH_CURVE_RES * 5, maxRes = maxSelObjRes)[1:-1]
                segInterpLocs += interpLocs
                selObjInterpCounts.append(len(interpLocs))

            if(withHandles):
                ptIdxs = info[splineIdx][1]
                for ptIdx in ptIdxs:
                    selObjHdlList.append([selObj, splineIdx, ptIdx])
                    hdlCnt = 2
                    if(wsDataSK != None):
                        pt = wsDataSK[splineIdx][ptIdx]
                        hdls += [pt[0], pt[2]]
                    else:
                        pt = wsData[splineIdx][ptIdx]
                        hdls += [pt[0], pt[2]]

                    selObjHdlCounts.append(hdlCnt)

    searchPtsList = [[], [], [], [], [], []]
    retStr = [[], [], [], [], [], []]

    #'SelHandles', 'SegLoc', 'CurveBezPt', 'CurveLoc'

    searchPtsList[0], retStr[0] = hdls, 'SelHandles'
    searchPtsList[1], retStr[1] = objSplineEndPts, 'CurveBezPt'
    searchPtsList[2], retStr[2] = segInterpLocs, 'SegLoc'
    searchPtsList[3], retStr[3] = objInterpLocs, 'CurveLoc'

    # TODO: Remove duplicates before sending for search?
    searchPtsList = [[getCoordFromLoc(region, rv3d, pt).to_3d() \
        for pt in pts] for pts in searchPtsList]

    srs = NestedListSearch(searchPtsList).findInLists(coFind, \
        searchRange = FTProps.snapDist)

    if(len(srs) == 0):
        return None

    sr = min(srs, key = lambda x: x[3])

    if(sr[0] > 1):
        # If seg loc then first priority to the nearby handle, end pt (even if farther)
        sr = min(srs, key = lambda x: (x[0], x[3]))

    idx = sr[1]
    retId = retStr[sr[0]]

    if(sr[0] == 0): # SelHandles
        listIdx = NestedListSearch.findListIdx(selObjHdlCounts, idx)
        obj, splineIdx, ptIdx = selObjHdlList[listIdx]
        # ~ obj, splineIdx, ptIdx = selObjHdlList[int(idx / 2)]
        return  retId, obj, splineIdx, ptIdx, 2 * (idx % 2)

    elif(sr[0]  == 1): # CurveBezPt
        listIdx = NestedListSearch.findListIdx(objBezPtCounts, idx)
        obj, splineIdx, ptIdx = objLocList[listIdx]
        # ~ obj, splineIdx, ptIdx = objLocList[int(idx / ptIdxCnt)]
        return retId, obj, splineIdx, ptIdx, 1 # otherInfo = segIdx

    elif(sr[0] == 2): # SegLoc
        listIdx = NestedListSearch.findListIdx(selObjInterpCounts, idx)
        obj, splineIdx, segIdx = selObjLocList[listIdx]
        return  retId, obj, splineIdx, segIdx, segInterpLocs[idx]

    else: # CurveLoc
        listIdx = NestedListSearch.findListIdx(objInterpCounts, idx)
        obj, splineIdx, segIdx = objLocList[listIdx]
        return retId, obj, splineIdx, segIdx, objInterpLocs[idx]

class NestedListSearch:
    # Find the list element containing the given idx from flattened list
    # return the index of the list element containing the idx
    def findListIdx(counts, idx):
        cumulCnt = 0
        cntIdx= 0
        while(idx >= cumulCnt):
            cumulCnt += counts[cntIdx] # cntIdx can never be >= len(counts)
            cntIdx += 1
        return cntIdx - 1

    def __init__(self, ptsList):
        self.ptsList = ptsList
        self.kd = kdtree.KDTree(sum(len(pts) for pts in ptsList))
        idx = 0
        self.counts = []
        for i, pts in enumerate(ptsList):
            self.counts.append(len(pts))
            for j, pt in enumerate(pts):
                self.kd.insert(pt, idx)
                idx += 1
        self.kd.balance()

    def findInLists(self, coFind, searchRange):
        if(searchRange == None):
            foundVals = [self.kd.find(coFind)]
        else:
            foundVals = self.kd.find_range(coFind, searchRange)
            foundVals = sorted(foundVals, key = lambda x: x[2])

        searchResults = []
        for co, idx, dist in foundVals:
            listIdx = NestedListSearch.findListIdx(self.counts, idx)
            ptIdxInList = idx - sum(len(self.ptsList[i]) for i in range(0, listIdx))
            searchResults.append([listIdx, ptIdxInList, co, dist])

        return searchResults


class SelectCurveInfo:
    def __init__(self, obj, splineIdx):
        self.obj = obj
        self.splineIdx = splineIdx
        self.updateWSData()

        # User Selection (mouse click); format ptIdx: set(sel)...
        # where sel: -1->seg, 0->left hdl, 1->bezier pt, 2->right hdl
        self.ptSels = {}

        # Highlighted point (mouse move)
        # 'ptIdx': ptIdx, 'hltIdx':hltIdx {-1, 0, 1, 2} (just as in sel above)
        self.hltInfo = {}

        # obj.name gives exception if obj is not in bpy.data.objects collection,
        # so keep a copy
        self.objName = obj.name
        self.interpPts = {}

        # Format 'ptIdx': segIdx, 'hdlIdx': hdlIdx, 'loc':loc, 't':t
        # hdlIdx - {-1, 0, 1, 2} similar to sel in ptSels
        self.clickInfo = {}

    def __hash__(self):
        return hash((self.objName, self.splineIdx))

    def updateWSData(self):
        self.hasShapeKey = (self.obj.active_shape_key != None)
        self.shapeKeyIdx = self.obj.active_shape_key_index if self.hasShapeKey else -1

        # for shape keys
        self.keyStartIdx = sum(len(self.obj.data.splines[i].bezier_points) \
            for i in range(self.splineIdx))

        # WS Data of the shape key (if exists)
        self.wsData = getBptData(self.obj, fromMix = False)[self.splineIdx]

    # For convenience
    def getAdjIdx(self, ptIdx, offset = 1):
        return getAdjIdx(self.obj, self.splineIdx, ptIdx, offset)

    def getBezierPt(self, ptIdx):
        return self.obj.data.splines[self.splineIdx].bezier_points[ptIdx]

    def getShapeKeyData(self, ptIdx, keyIdx = None):
        if(not self.hasShapeKey): return None
        if(keyIdx == None): keyIdx = self.obj.active_shape_key_index
        keydata = self.obj.data.shape_keys.key_blocks[keyIdx].data
        keyIdx = self.keyStartIdx + ptIdx
        return keydata[keyIdx] if(keyIdx < len(keydata)) else None

    def getAllShapeKeysData(self, ptIdx):
        if(not self.hasShapeKey): return None
        pts = []
        for keyIdx in range(len(self.obj.data.shape_keys.key_blocks)):
            pts.append(self.getShapeKeyData(ptIdx, keyIdx))
        return pts

    def getSegPtsInfo(self, ptIdx):
        nextIdx = self.getAdjIdx(ptIdx)
        pt0 = self.wsData[ptIdx]
        pt1 = self.wsData[nextIdx]
        return nextIdx, pt0, pt1

    def getSegPts(self, ptIdx):
        nextIdx, pt0, pt1 = self.getSegPtsInfo(ptIdx)
        return (pt0, pt1)

    # All selected points which have handles displayed (for example for snapping)
    def getAllPtsWithHdls(self):
        ptIdxs = set(self.ptSels.keys())
        nextIdxs = set(self.getSegPtsInfo(p)[0] for p in ptIdxs if -1 in self.ptSels[p])
        ptIdxs = ptIdxs.union(nextIdxs)
        return sorted(self.wsData[p] for p in ptIdxs)

    def setClickInfo(self, ptIdx, hdlIdx, clickLoc, lowerT = 0.001, higherT = .999):
        self.clickInfo = None
        if(clickLoc != None):
            nextIdx, pt0, pt1 = self.getSegPtsInfo(ptIdx)
            t = getTForPt([pt0[1], pt0[2], pt1[0], pt1[1]], clickLoc)

            if(t != None and (t < lowerT or t > higherT)):
                hdlIdx = 1
                if(t > higherT): ptIdx = nextIdx
            else:
                self.clickInfo = {'ptIdx': ptIdx, 'hdlIdx': hdlIdx, \
                    'loc':clickLoc, 't':t}

        if(self.clickInfo == None):
            self.clickInfo = {'ptIdx': ptIdx, 'hdlIdx': hdlIdx}

    def addSel(self, ptIdx, sel, toggle = False):
        self.addSels(ptIdx, set([sel]), toggle)

    def addSels(self, ptIdx, sels, toggle = False):
        # TODO: Check this condition at other places
        if( -1 in sels and self.getAdjIdx(ptIdx) == None): sels.remove(-1)
        if(len(sels) == 0): return

        currSels = self.ptSels.get(ptIdx)
        if(currSels == None): currSels = set()
        modSels = currSels.union(sels)
        if(toggle): modSels -= currSels.intersection(sels)
        self.ptSels[ptIdx] = modSels
        if(len(self.ptSels[ptIdx]) == 0 ): self.ptSels.pop(ptIdx)

    def removeSel(self, ptIdx, sel):
        self.removeSels(ptIdx, {sel})

    def removeSels(self, ptIdx, sels):
        for sel in sels:
            currSels = self.ptSels.get(ptIdx)
            if(currSels != None and sel in currSels):
                currSels.remove(sel)
                if(len(currSels) == 0): self.ptSels.pop(ptIdx)

    def resetClickInfo(self):
        self.clickInfo = {}

    def resetPtSel(self):
        self.ptSels = {}

    def resetHltInfo(self):
        self.hltInfo = {}

    def getHltInfo(self):
        return self.hltInfo

    def setHltInfo(self, ptIdx, hltIdx):
        self.hltInfo = {'ptIdx': ptIdx, 'hltIdx':hltIdx}

    def getClickLoc(self):
        return self.clickInfo.get('loc')

    def getSelCo(self):
        if(len(self.clickInfo) > 0):
            hdlIdx = self.clickInfo['hdlIdx']
            if(hdlIdx == -1):
                return self.clickInfo['loc']
            else:
                ptIdx = self.clickInfo['ptIdx']
                pt0 = self.wsData[ptIdx]
                return pt0[hdlIdx]

        return None

    def subdivSeg(self, subdivCnt):
        if(self.hasShapeKey): return False
        if(subdivCnt > 1):
            invMw = self.obj.matrix_world.inverted_safe()
            ts = []
            addCnt = 0
            for ptIdx in sorted(self.ptSels.keys()):
                if(-1 in self.ptSels[ptIdx]):
                    vertCos = getInterpolatedVertsCo(self.interpPts[ptIdx], \
                        subdivCnt)[1:-1]
                    changedIdx = ptIdx + addCnt
                    insertBezierPts(self.obj, self.splineIdx, changedIdx, \
                        [invMw @ v for v in vertCos], 'FREE')
                    addCnt += len(vertCos)
        return addCnt > 0

    def bevelPts(self, bevelCnt, deltaPos):
        if(self.hasShapeKey): return False
        pts, ptSels = self.getBevelPts(bevelCnt, self.wsData, deltaPos)
        spline = self.obj.data.splines[self.splineIdx]
        newPtCnt = len(pts) - len(self.wsData)
        spline.bezier_points.add(newPtCnt)
        for pt in spline.bezier_points:
            pt.handle_left_type = 'FREE'
            pt.handle_right_type = 'FREE'

        invMw = self.obj.matrix_world.inverted_safe()
        for i, pt in enumerate(spline.bezier_points):
            pt.handle_left = invMw @ pts[i][0]
            pt.co = invMw @ pts[i][1]
            pt.handle_right = invMw @ pts[i][2]
            pt.handle_left_type = pts[i][3]
            pt.handle_right_type = pts[i][4]

        self.updateWSData()
        self.ptSels = ptSels
        return True

    def initSubdivMode(self, rv3d):
        if(self.hasShapeKey): return False
        changed = False
        for ptIdx in self.ptSels.keys():
            if(-1 in self.ptSels[ptIdx]):
                self.interpPts[ptIdx] = getPtsAlongBezier3D(self.getSegPts(ptIdx), rv3d,
                    curveRes = 1000, minRes = 1000)
                changed = True
        return changed

    def isBevelabel(self, rv3d):
        if(self.hasShapeKey): return False
        changed = False
        for ptIdx in self.ptSels.keys():
            ptIdxs = [ptIdx]
            if(-1 in self.ptSels[ptIdx]): ptIdxs.append(self.getAdjIdx(ptIdx))
            elif(1 not in self.ptSels[ptIdx]): continue # only pt and seg selection
            for idx in ptIdxs:
                prevIdx = self.getAdjIdx(idx, -1)
                nextIdx = self.getAdjIdx(idx)
                if(nextIdx != None and prevIdx != None and \
                    not hasAlignedHandles(self.wsData[idx])):
                    changed = True
                    break
        return changed

    def getLastSegIdx(self):
        return getLastSegIdx(self.obj, self.splineIdx)

    # Remove all selected segments
    # Returns map with spline index and seg index change after every seg removal
    def removeSegs(self):
        changedSelMap = {}
        if(self.hasShapeKey): return changedSelMap
        segSels = [p for p in self.ptSels if -1 in self.ptSels[p]]
        cumulSegIdxIncr = 0
        changedSplineIdx = self.splineIdx
        segIdxIncr = 0

        for segIdx in sorted(segSels):
            changedSegIdx = segIdx + cumulSegIdxIncr
            splineIdxIncr, segIdxIncr = removeBezierSeg(self.obj, \
                changedSplineIdx, changedSegIdx)
            changedSplineIdx += splineIdxIncr
            cumulSegIdxIncr += segIdxIncr
            changedSelMap[segIdx] = [splineIdxIncr, segIdxIncr]
        return changedSelMap

    def straightenHandle(self, ptIdx, hdlIdx, allShapekeys = False):
        bpt = self.getBezierPt(ptIdx)
        prevIdx = self.getAdjIdx(ptIdx, -1)
        nextIdx = self.getAdjIdx(ptIdx)
        prevPts = None
        nextPts = None
        if(self.hasShapeKey):
            if(allShapekeys):
                pts = self.getAllShapeKeysData(ptIdx)
                if(prevIdx != None): prevPts = self.getAllShapeKeysData(prevIdx)
                if(nextIdx != None): nextPts = self.getAllShapeKeysData(nextIdx)
            else:
                pts = [self.getShapeKeyData(ptIdx)]
                if(prevIdx != None): prevPts = [self.getShapeKeyData(prevIdx)]
                if(nextIdx != None): nextPts = [self.getShapeKeyData(nextIdx)]

        else:
            pts = [bpt]
            if(prevIdx != None): prevPts = [self.getBezierPt(prevIdx)]
            if(nextIdx != None): nextPts = [self.getBezierPt(nextIdx)]

        if(hdlIdx == 0):
            if(bpt.handle_left_type != 'VECTOR'): bpt.handle_left_type = 'FREE'
            for i in range(len(pts)):
                pt = pts[i]
                if(prevPts != None): diffV = (pt.co - prevPts[i].co)
                else: diffV = (nextPts[i].co - pt.co)
                pt.handle_left = pt.co - diffV / 3
        elif(hdlIdx == 2):
            if(bpt.handle_right_type != 'VECTOR'):  bpt.handle_right_type = 'FREE'
            for i in range(len(pts)):
                pt = pts[i]
                if(nextPts != None): diffV = (nextPts[i].co - pt.co)
                else: diffV = (pt.co - prevPts[i].co)
                pt.handle_right = pt.co + diffV / 3

    def straightenSelHandles(self):
        changed = False
        for ptIdx in self.ptSels:
            for hdlIdx in self.ptSels[ptIdx]:
                self.straightenHandle(ptIdx, hdlIdx)
                changed = True
        return changed

    def alignHandle(self, ptIdx, hdlIdx, allShapekeys = False):
        if (hdlIdx == -1): return False
        oppIdx = 2 - hdlIdx
        if(self.hasShapeKey):
            if(allShapekeys): pts = self.getAllShapeKeysData(ptIdx)
            else: pts = [self.getShapeKeyData(ptIdx)]
        else: pts = [self.getBezierPt(ptIdx)]
        bpt = self.getBezierPt(ptIdx)
        if(hdlIdx == 0 and bpt.handle_left_type != 'ALIGNED'):
            bpt.handle_left_type = 'FREE'
        if(hdlIdx == 2 and bpt.handle_right_type != 'ALIGNED'):
            bpt.handle_right_type = 'FREE'

        for pt in pts:
            if(hdlIdx == 0): pt.handle_left = pt.co - \
                (pt.co - pt.handle_left).length * (pt.handle_right - pt.co).normalized()
            else: pt.handle_right = pt.co + \
                (pt.co - pt.handle_right).length * (pt.co - pt.handle_left).normalized()
        return True

    def alignSelHandles(self):
        changed = False
        invMw = self.obj.matrix_world.inverted_safe()
        for ptIdx in self.ptSels:
            sels = self.ptSels[ptIdx]
            for hdlIdx in sels:
                changed = self.alignHandle(ptIdx, hdlIdx) or changed
        return changed

    def insertNode(self, handleType, select = True):
        if(self.hasShapeKey): return False
        invMw = self.obj.matrix_world.inverted_safe()
        insertBezierPts(self.obj, self.splineIdx, \
            self.clickInfo['ptIdx'], [invMw @ self.clickInfo['loc']], handleType)
        return True

    def removeNode(self):
        if(self.hasShapeKey): return False
        changed = False

        toRemove = set() # Bezier points to remove from object
        toRemoveSel = set() # Selection entry to remove from ptSels

        nodeSels = [p for p in self.ptSels if 1 in self.ptSels[p]]

        for ptIdx in nodeSels:
            self.ptSels.pop(ptIdx)

        if(len(nodeSels) > 0):
            removeBezierPts(self.obj, self.splineIdx, nodeSels)
            changed = True

        if(changed):
            selIdxs = sorted(self.ptSels.keys())
            cnt = 0
            for ptIdx in nodeSels:
                cIdxs = [i for i in selIdxs if i >= (ptIdx - cnt)]
                for idx in cIdxs:
                    sels = self.ptSels.pop(idx - cnt)
                    self.ptSels[idx - cnt - 1] = sels
                cnt += 1

        return changed

    def getBevelPts(self, bevelCnt, pts, deltaPos):
        deltaLen = deltaPos.length
        if(floatCmpWithMargin(deltaLen, DEF_ERR_MARGIN)):
            return pts, self.ptSels

        # http://launchpadlibrarian.net/12692602/rcp.svg
        kFact = (bevelCnt/3) * (sqrt(2) - 1)
        maxT = .5

        # Deep copy
        pts = [[p if isinstance(p, str) else p.copy() for p in pt] for pt in pts]

        newPts = []
        newSelPtIdxs = []

        # Add both points of the selected segments in selection
        bevelPtIdxs = set()
        for ptIdx in self.ptSels.keys():
            if(1 in self.ptSels[ptIdx]): bevelPtIdxs.add(ptIdx)
            if(-1 in self.ptSels[ptIdx]):
                adjIdx = self.getAdjIdx(ptIdx)
                bevelPtIdxs.add(ptIdx)
                bevelPtIdxs.add(adjIdx)

        ptSels = {k:self.ptSels[k].copy() for k in self.ptSels.keys()}

        # Extra loop because next points need to be determined beforehand
        for ptIdx in bevelPtIdxs:
            if(ptSels.get(ptIdx) == None): ptSels[ptIdx] = {1}
            else: ptSels[ptIdx].add(1)
            nextIdx = self.getAdjIdx(ptIdx)
            prevIdx = self.getAdjIdx(ptIdx, -1)
            if(prevIdx != None and nextIdx != None and \
                not hasAlignedHandles(pts[ptIdx])):
                newSelPtIdxs.append(ptIdx)

        for ptIdx, pt in enumerate(pts):
            if(ptIdx in newSelPtIdxs):
                nextIdx = self.getAdjIdx(ptIdx)
                prevIdx = self.getAdjIdx(ptIdx, -1)
                prevPt = pts[prevIdx][:]
                diffV = (pt[1] - prevPt[1])
                segLen = diffV.length
                if(segLen < .001):
                    newPts.append(pt)
                else:
                    t = deltaLen / segLen
                    if(t > maxT and (prevIdx in newSelPtIdxs)):
                        t = maxT
                        k = kFact * (segLen / 2)
                    elif(t > 1):
                        t = 1
                        k = kFact * segLen
                    else:
                        k = kFact * deltaLen

                    seg = [pts[prevIdx][1], pts[prevIdx][2], pt[0], pt[1]]
                    partialSeg = getPartialSeg(seg, t0 = 0, t1 = 1 - t)
                    newPt = partialSeg[3]

                    tangent0 = getTangentAtT(pts[prevIdx][1], pts[prevIdx][2], \
                        pt[0], pt[1], 1 - t)
                    pt0_2 = newPt + k * (tangent0.normalized())

                    pt0 = [partialSeg[2], newPt, pt0_2, 'FREE', 'FREE']
                    newPts.append(pt0)

                    prevPt[2] = partialSeg[1]

                nextPt = pts[nextIdx][:]
                diffV = (nextPt[1] - pt[1])
                segLen = diffV.length
                if(floatCmpWithMargin(segLen, 0)):
                    newPts.append(pt)
                else:
                    t = deltaLen / segLen
                    if(t > maxT and (nextIdx in newSelPtIdxs)):
                        t = maxT
                        k = kFact * (segLen / 2)
                    elif(t > 1):
                        t = 1
                        k = kFact * segLen
                    else:
                        k = kFact * deltaLen

                    seg = [pt[1], pt[2], pts[nextIdx][0], pts[nextIdx][1]]
                    partialSeg = getPartialSeg(seg, t0 = t, t1 = 1)
                    newPt = partialSeg[0]

                    tangent1 = getTangentAtT(pt[1], pt[2], \
                        pts[nextIdx][0], pts[nextIdx][1], t)

                    pt1_0 = newPt - k * (tangent1.normalized())
                    pt1 = [pt1_0, newPt, partialSeg[1], 'FREE', 'FREE']

                    newPts.append(pt1)

                    nextPt[0] = partialSeg[2]
            else:
                newPts.append(pt)

        newPtSels = {}
        cnt = 0
        for ptIdx in sorted(ptSels.keys()):
            newPtSels[ptIdx + cnt] = ptSels[ptIdx].copy()
            if(ptIdx in newSelPtIdxs):
                newPtSels[ptIdx + cnt].add(-1)
                newPtSels[ptIdx + cnt + 1] = {1}
                cnt += 1

        return newPts, newPtSels

    def getDisplayInfos(self, hideHdls = False, subdivCnt = 0, \
        bevelCnt = 0, newPos = None, deltaPos = None):

        # Making long short
        cHltTip = FTProps.colHltTip
        cBezPt = FTProps.colBezPt
        cHdlPt = FTProps.colHdlPtTip
        cAdjBezTip = FTProps.colAdjBezTip
        cNonHltSeg = FTProps.colDrawNonHltSeg

        segDispInfos = []
        bptDispInfos = []

        pts = self.wsData[:]
        ptSels = self.ptSels
        if(newPos != None):
            # TODO: This method is in EditCurveInfo
            nPtIdxs, nPts = self.getOffsetSegPts(newPos)
            # Update list with new position (editing)
            for i, ptIdx in enumerate(nPtIdxs): pts[ptIdx] = nPts[i]
        elif(deltaPos != None):
            pts, ptSels = self.getBevelPts(bevelCnt, pts, deltaPos)

        # Default display of spline
        for i, pt in enumerate(pts):
            bptDispInfos.append(BptDisplayInfo(pt, [cAdjBezTip]))
            if(i > 0):
                segDispInfos.append(SegDisplayInfo([pts[i-1], pt], cNonHltSeg))
        lastIdx = self.getAdjIdx(len(self.wsData) - 1) # In case cyclic...
        if(lastIdx != None):
            segDispInfos.append(SegDisplayInfo([pts[-1], pts[0]], cNonHltSeg))

        hltInfo = self.getHltInfo()
        hltPtIdx = hltInfo.get('ptIdx')
        hltIdx = hltInfo.get('hltIdx')

        # Process highlighted segments before selected ones because...
        # selected segments take priority over highlighted
        if(hltIdx == -1):
            segDispInfos[hltPtIdx].segColor = FTProps.colDrawHltSeg
            bptDispInfos[hltPtIdx].tipColors[1] = cBezPt
            nextIdx = self.getAdjIdx(hltPtIdx)
            bptDispInfos[nextIdx].tipColors[1] = cBezPt

        # Process selections
        for ptIdx in sorted(ptSels.keys()):
            sels = ptSels[ptIdx]

            if(hideHdls):
                tipColors = [None, cBezPt, None]
                handleNos = []
            else:
                tipColors = [cHdlPt, cBezPt, cHdlPt]
                handleNos = [0, 1]

            bptDispInfos[ptIdx].tipColors = tipColors[:]
            bptDispInfos[ptIdx].handleNos = handleNos

            for hdlIdx in sorted(sels): # Start with seg selection i. e. -1
                if(hdlIdx == -1):
                    nextIdx = getAdjIdx(self.obj, self.splineIdx, ptIdx, ptCnt = len(pts))
                    segPts = [pts[ptIdx], pts[nextIdx]]

                    # process next only if there are no selection pts with that idx
                    if(nextIdx not in ptSels.keys()):
                        bptDispInfos[nextIdx].tipColors = tipColors[:]
                        bptDispInfos[nextIdx].handleNos = handleNos

                    vertCos = []
                    if(subdivCnt > 1):
                        vertCos = getInterpolatedVertsCo(self.interpPts[ptIdx], \
                            subdivCnt)[1:-1]

                    selSegDispInfo = EditSegDisplayInfo(segPts, \
                        FTProps.colDrawSelSeg, vertCos)
                    segDispInfos[ptIdx] = selSegDispInfo
                elif(hdlIdx == 1 or not hideHdls):
                    bptDispInfos[ptIdx].tipColors[hdlIdx] = FTProps.colSelTip

        # Process highlighted points after selected ones because...
        # highlighted points take priority over selected
        if(hltIdx in {0, 1, 2}):
            bptDispInfos[hltPtIdx].tipColors[hltIdx] = cHltTip

        return [segDispInfos, bptDispInfos]

class EditCurveInfo(SelectCurveInfo):
    def __init__(self, obj, splineIdx, ptSels = None):
        super(EditCurveInfo, self).__init__(obj, splineIdx)
        if(ptSels != None):
            self.ptSels = ptSels

    def syncAlignedHdl(self, pt, ctrlPLoc, hdlIdx):
        typeIdx = 3 if hdlIdx == 0 else 4
        if(pt[typeIdx] == 'ALIGNED'):
            oppTypeIdx = 4 if hdlIdx == 0 else 3
            if(pt[oppTypeIdx] in {'VECTOR', 'ALIGNED'}):
                oppHdlIdx = 2 if hdlIdx == 0 else 0
                oppHdlV = ctrlPLoc - pt[oppHdlIdx]
                if(oppHdlV.length != 0):
                    currL = (ctrlPLoc - pt[hdlIdx]).length
                    pt[hdlIdx] = ctrlPLoc + currL * oppHdlV / oppHdlV.length

    def setAlignedHdlsCo(self, pt, hdlIdx, ctrlPLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        if(pt[typeIdx] == 'ALIGNED'):
            oppTypeIdx = 4 if hdlIdx == 0 else 3
            if(pt[oppTypeIdx] != 'VECTOR'):
                pt[hdlIdx] += (ctrlPLoc - pt[1])
            else:
                self.syncAlignedHdl(pt, ctrlPLoc, hdlIdx)

    def setFreeHdlsCo(self, pt, hdlIdx, newLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        if(pt[typeIdx] == 'FREE'):
            pt[hdlIdx] += (newLoc - pt[1])

    def setVectHdlsCo(self, pt, newLoc, hdlIdx, prevPt, nextPt):
        typeIdx = 3 if hdlIdx == 0 else 4
        if(pt[typeIdx] == 'VECTOR'):
            typeIdx = 3 if hdlIdx == 0 else 4
            pts = [prevPt, nextPt] if hdlIdx == 0 else [nextPt, prevPt]
            diffV = None
            if(pts[0] != None): diffV = pts[0][1] - newLoc
            if(diffV == None and pts[1] != None):
                diffV = newLoc - pts[1][1]
            if(diffV == None): pt[hdlIdx] = newLoc
            else: pt[hdlIdx] = newLoc + diffV * 1 / 3

    # Calculate both handle and adjacent pt handles in case of Vector type
    # TODO: AUTO has a separate logic set to ALIGNED for now
    def syncCtrlPtHdls(self, ptIdx, newLoc):
        wsData = getBptData(self.obj, fromMix = False)
        pt = wsData[self.splineIdx][ptIdx]
        prevIdx = self.getAdjIdx(ptIdx, -1)
        prevPt = None if prevIdx == None else wsData[self.splineIdx][prevIdx]
        nextIdx = self.getAdjIdx(ptIdx)
        nextPt = None if nextIdx == None else wsData[self.splineIdx][nextIdx]

        ptIdxs = [ptIdx]
        pts = [pt]

        for typeIdx in [3, 4]:
            if(pt[typeIdx] == 'AUTO'): pt[typeIdx] = 'ALIGNED'
        for hdlIdx in [0, 2]:
            self.setVectHdlsCo(pt, newLoc, hdlIdx, prevPt, nextPt)
        for hdlIdx in [0, 2]:
            self.setFreeHdlsCo(pt, hdlIdx, newLoc)
        for hdlIdx in [0, 2]:
            self.setAlignedHdlsCo(pt, hdlIdx, newLoc)

        pt[1] = newLoc

        if(prevPt != None and prevPt[4] == 'VECTOR'):
            pPrevIdx = self.getAdjIdx(prevIdx, -1)
            pPrevPt = None if pPrevIdx == None else wsData[self.splineIdx][pPrevIdx]
            self.setVectHdlsCo(prevPt, prevPt[1], 2, pPrevPt, pt)
            self.setAlignedHdlsCo(prevPt, 0, prevPt[1])

            ptIdxs.append(prevIdx)
            pts.append(prevPt)

        if(nextPt != None and nextPt[3] == 'VECTOR'):
            nNextIdx = self.getAdjIdx(nextIdx)
            nNextPt = None if nNextIdx == None else wsData[self.splineIdx][nNextIdx]
            self.setVectHdlsCo(nextPt, nextPt[1], 0, pt, nNextPt)
            self.setAlignedHdlsCo(nextPt, 2, nextPt[1])

            ptIdxs.append(nextIdx)
            pts.append(nextPt)

        return ptIdxs, pts

    # Calculate the opposite handle values in case of ALIGNED and AUTO handles
    # Also set the type(s) of current (opposite) handle(s)
    def syncHdls(self, pt, hdlIdx, newLoc):
        typeIdx = 3 if hdlIdx == 0 else 4
        oppTypeIdx = 4 if hdlIdx == 0 else 3

        if(pt[typeIdx] == 'VECTOR'): pt[typeIdx] = 'FREE'
        if(pt[typeIdx] == 'AUTO'): pt[typeIdx] = 'ALIGNED'
        if(pt[oppTypeIdx] == 'AUTO' and pt[typeIdx] != 'FREE'): pt[oppTypeIdx] = 'ALIGNED'

        pt[hdlIdx] = newLoc

        self.syncAlignedHdl(pt, pt[1], 2 - hdlIdx) # First opposite

        self.syncAlignedHdl(pt, pt[1], hdlIdx)


    # Get seg points after change in position of handles or drag curve
    # The only function called on all 3 events: grab curve pt, grab handle, grab Bezier pt
    def getOffsetSegPts(self, newLoc):
        inf = self.clickInfo
        ptIdx = inf['ptIdx']
        hdlIdx = inf['hdlIdx']
        wsData = getBptData(self.obj, fromMix = False)
        pt = wsData[self.splineIdx][ptIdx]

        if(hdlIdx == -1): # Grab point on curve
            adjIdx = self.getAdjIdx(ptIdx)
            adjPt = wsData[self.splineIdx][adjIdx]

            ptIdxs = [ptIdx, adjIdx]
            pts = [pt, adjPt]

            delta = newLoc - inf['loc']
            if(delta == 0):
                return ptIdxs, pts
            t = inf['t']

            #****************************************************************
            # Magic Bezier Drag Equations (Courtesy: Inkscape)             #*
            #****************************************************************
                                                                           #*
            if (t <= 1.0 / 6.0):                                           #*
                weight = 0                                                 #*
            elif (t <= 0.5):                                               #*
                weight = (pow((6 * t - 1) / 2.0, 3)) / 2                   #*
            elif (t <= 5.0 / 6.0):                                         #*
                weight = (1 - pow((6 * (1-t) - 1) / 2.0, 3)) / 2 + 0.5     #*
            else:                                                          #*
                weight = 1                                                 #*
                                                                           #*
            offset0 = ((1 - weight) / (3 * t * (1 - t) * (1 - t))) * delta #*
            offset1 = (weight / (3 * t * t * (1 - t))) * delta             #*
                                                                           #*
            #****************************************************************

            # If the segment is edited, the 1st pt right handle...
            pts[0][2] += offset0

            if(pts[0][4] == 'VECTOR'): pts[0][4] = 'FREE'
            if(pts[0][4] == 'AUTO'): pts[0][4] = 'ALIGNED'
            # opposite handle must be changed if this is not FREE
            if(pts[0][3] == 'VECTOR' and pts[0][4] != 'FREE'): pts[0][3] = 'FREE'
            if(pts[0][3] == 'AUTO' and pts[0][4] != 'FREE'): pts[0][3] = 'ALIGNED'

            self.syncAlignedHdl(pts[0], pts[0][1], hdlIdx = 0)

            # ...and 2nd pt left handle impacted
            pts[1][0] += offset1

            if(pts[1][3] == 'VECTOR'): pts[1][3] = 'FREE'
            if(pts[1][3] == 'AUTO'): pts[1][3] = 'ALIGNED'
            # opposite handle must be changed if this is not FREE
            if(pts[1][4] == 'VECTOR' and  pts[1][3] != 'FREE'): pts[1][4] = 'FREE'
            if(pts[1][4] == 'AUTO' and  pts[1][3] != 'FREE'): pts[1][4] = 'ALIGNED'

            self.syncAlignedHdl(pts[1], pts[1][1], hdlIdx = 2)
            return ptIdxs, pts

        elif(hdlIdx in {0, 2}): # Grab one of the handles
            self.syncHdls(pt, hdlIdx, newLoc)
            return [ptIdx], [pt]

        else: # Grab the Bezier point
            return self.syncCtrlPtHdls(ptIdx, newLoc)

    def moveSeg(self, newPos):
        ptIdxs, pts = self.getOffsetSegPts(newPos)

        invMw = self.obj.matrix_world.inverted_safe()
        spline = self.obj.data.splines[self.splineIdx]
        bpts = [spline.bezier_points[idx] for idx in ptIdxs]

        for i, bpt in enumerate(bpts):
            bpt.handle_right_type = 'FREE'
            bpt.handle_left_type = 'FREE'

        if(self.hasShapeKey):
            for i, ptIdx in enumerate(ptIdxs):
                keydata = self.getShapeKeyData(ptIdx)
                keydata.handle_left = invMw @ pts[i][0]
                keydata.co = invMw @ pts[i][1]
                keydata.handle_right = invMw @ pts[i][2]
                if(pts[i][3] == 'AUTO'): pts[i][3] = 'ALIGNED'
                if(pts[i][4] == 'AUTO'): pts[i][4] = 'ALIGNED'
                impIdxs = [ptIdx, self.getAdjIdx(ptIdx, -1), self.getAdjIdx(ptIdx, 1)]
                for idx in impIdxs:
                    if(idx == None): continue
                    if(spline.bezier_points[idx].handle_left_type == 'AUTO'):
                        spline.bezier_points[idx].handle_left_type = 'ALIGNED'
                    if(spline.bezier_points[idx].handle_right_type == 'AUTO'):
                        spline.bezier_points[idx].handle_right_type = 'ALIGNED'
        else:
            for i, bpt in enumerate(bpts):
                bpt.handle_left = invMw @ pts[i][0]
                bpt.co = invMw @ pts[i][1]
                bpt.handle_right = invMw @ pts[i][2]

        for i, bpt in enumerate(bpts):
            bpt.handle_left_type = pts[i][3]
            bpt.handle_right_type = pts[i][4]

        self.updateWSData()

class ModalFlexiEditBezierOp(ModalBaseFlexiOp):
    bl_description = "Flexi editing of Bezier curves in object mode"
    bl_idname = "wm.modal_flexi_edit_bezier"
    bl_label = "Flexi Edit Curve"
    bl_options = {'REGISTER', 'UNDO'}

    h = False

    def drawHandler():
        ModalBaseFlexiOp.drawHandlerBase()

    def resetDisplay():
        ModalBaseFlexiOp.resetDisplayBase()

    # static method
    def refreshDisplay(segDispInfos, bptDispInfos, locOnCurve = None, snapper = None):

        ptCos = [co for d in segDispInfos if type(d) == EditSegDisplayInfo
            for co in d.subdivCos]

        # ~ if(locOnCurve != None): ptCos.append(locOnCurve) # For debugging

        ModalBaseFlexiOp.bglDrawMgr.addPtInfo('editSubdiv', FTProps.editSubdivPtSize, \
            [FTProps.colEditSubdiv], ptCos)

        ModalBaseFlexiOp.refreshDisplayBase(segDispInfos, bptDispInfos, snapper)

    def getToolType(self):
        return TOOL_TYPE_FLEXI_EDIT

    # Refresh display with existing curves (nonstatic)
    def refreshDisplaySelCurves(self, hltSegDispInfos = None, hltBptDispInfos = None, \
        locOnCurve = None, refreshPos = False):

        if(self.rmInfo == None): return # Possible in updateAfterGeomChange

        newPos = None
        if(FTProps.liveUpdate and self.editCurveInfo != None):
            newPos = self.getNewPos(refreshStatus = True)
            self.editCurveInfo.moveSeg(newPos)
            clickInfo = self.editCurveInfo.clickInfo
            if(clickInfo['hdlIdx'] == -1):
                self.editCurveInfo.setClickInfo(clickInfo['ptIdx'], \
                    clickInfo['hdlIdx'], newPos)
            self.xyPress = self.rmInfo.xy[:]

        segDispInfos = []
        bptDispInfos = []
        # ~ curveInfos = self.selectCurveInfos.copy()
        # ~ if(self.editCurveInfo != None):
            # ~ curveInfos.add(self.editCurveInfo)
        if(self.bevelMode):
            deltaPos = self.getNewDeltaPos(refreshStatus = True)
        else:
            deltaPos = None
        for c in self.selectCurveInfos:
            if(refreshPos and c == self.editCurveInfo and newPos == None):
                newPos = self.getNewPos(refreshStatus = True)
            else:
                newPos = None
            info1, info2 = c.getDisplayInfos(hideHdls = ModalFlexiEditBezierOp.h, \
                subdivCnt = self.subdivCnt, bevelCnt = self.bevelCnt, newPos = newPos, \
                    deltaPos = deltaPos)
            segDispInfos += info1
            bptDispInfos += info2

        # Highlighted at the top
        if(hltSegDispInfos != None): segDispInfos += hltSegDispInfos
        if(hltBptDispInfos != None): bptDispInfos += hltBptDispInfos

        ModalFlexiEditBezierOp.refreshDisplay(segDispInfos, bptDispInfos, \
            locOnCurve, self.snapper)

    def reset(self):
        self.editCurveInfo = None
        self.selectCurveInfos = set()
        #TODO: freezeOrient logic should be internal to Snapper
        if(self.snapper != None): self.snapper.freezeOrient = False
        ModalFlexiEditBezierOp.resetDisplay()

    def postUndoRedo(self, scene, dummy = None): # signature different in 2.8 and 2.81?
        # ~ self.snapper.customAxis.reload()
        self.updateAfterGeomChange()
        for ci in self.selectCurveInfos: ci.resetPtSel()

    def cancelOp(self, context):
        self.reset()
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)
        bpy.app.handlers.depsgraph_update_post.remove(self.updateAfterGeomChange)
        return self.cancelOpBase()

    def isToolSelected(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)
        if(tool == None or tool.idname != FlexiEditBezierTool.bl_idname):
        # if(tool == None or tool.idname != 'flexi_bezier.edit_tool'):
            return False
        return True

    # Will be called after the curve is changed (by the tool or externally)
    # So handle all possible conditions
    def updateAfterGeomChange(self, scene = None, dummy = None): # 3 params in 2.81
        ciRemoveList = []

        removeObjNames = set() # For snaplocs
        addObjNames = set()
        self.htlCurveInfo = None

        # TODO: check if self.editCurveInfo is to be set to None
        if(not FTProps.liveUpdate):
            self.editCurveInfo = None # Reset if editing (capture == True)

        for ci in self.selectCurveInfos:
            if(bpy.data.objects.get(ci.objName) != None):
                ci.obj = bpy.data.objects.get(ci.objName) #refresh anyway
                splines = ci.obj.data.splines
                if(ci.splineIdx >= len(ci.obj.data.splines)):
                    ciRemoveList.append(ci)
                    continue
                spline = splines[ci.splineIdx]
                bpts = spline.bezier_points
                bptsCnt = len(bpts)
                # Don't keep a point object / spline
                if(bptsCnt <= 1):
                    if(len(splines) == 1):
                        ciRemoveList.append(ci)
                    else:
                        splines.remove(spline)
                        if(ci.splineIdx >= (len(splines))):
                            ci.splineIdx = len(splines) - 1
                        ci.resetPtSel()
                else:
                    # If any of the current selections is / are...
                    # greater than last idx (pt) or last but one (seg)...
                    # move it / them to last idx (pt) or last but one idx (seg)
                    changeSegSels = set()
                    changePtSels = set()
                    lastIdx = bptsCnt - 1
                    lastSegIdx = ci.getLastSegIdx()
                    for ptIdx in ci.ptSels.keys():
                        sels = ci.ptSels[ptIdx]
                        if(-1 in sels and ptIdx > lastSegIdx):
                            changeSegSels.add(ptIdx)
                            sels.remove(-1)
                        if(ptIdx > lastIdx and len(sels) > 0):
                            changePtSels.add(ptIdx)
                    for ptIdx in changePtSels:
                        sels = ci.ptSels.pop(ptIdx)
                        ci.addSels(lastSegIdx, sels)
                    for pt in changeSegSels:
                        ci.addSels((lastIdx - 1), set([-1]))

                addObjNames.add(ci.objName)
                ci.updateWSData()
            else:
                ciRemoveList.append(ci)
                removeObjNames.add(ci.objName)

        if(len(ciRemoveList) > 0):
            for c in ciRemoveList:
                self.selectCurveInfos.remove(c)

        if(self.editCurveInfo == None): # exclude live update condition
            self.updateSnapLocs(addObjNames, removeObjNames)
            self.refreshDisplaySelCurves()

    def subInvoke(self, context, event):
        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

        self.editCurveInfo = None
        self.htlCurveInfo = None
        self.selectCurveInfos = set()
        self.subdivCnt = 0
        self.bevelCnt = 4
        self.bevelMode = False

        # For double click (TODO: remove; same as editCurveInfo == None?)
        self.capture = False
        self.xyPress = None # ...to avoid jerky movement at the beginning
        self.xyLoc = None # for bevel

        self.snapInfos = {}
        self.updateSnapLocs()

        return  {"RUNNING_MODAL"}

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if(not ModalFlexiEditBezierOp.h):
            for ci in self.selectCurveInfos:
                pts = ci.getAllPtsWithHdls()
                for pt in pts: # Already world space
                    locs.append(pt[0])
                    locs.append(pt[2])
        return locs

    def updateSnapLocs(self, addObjNames = None, removeObjNames = None):
        updateCurveEndPtMap(self.snapInfos, addObjNames, removeObjNames)

    def getRefLine(self):
        if(self.editCurveInfo != None):
            ei = self.editCurveInfo
            ptIdx = ei.clickInfo['ptIdx']
            hdlIdx = ei.clickInfo['hdlIdx']
            pt0 = ei.wsData[ptIdx]
            if(hdlIdx in {0, 2}):
                return [pt0[2 - hdlIdx], pt0[1]] # Opposite handle
            else: # point on curve or Bezier point so previous segment
                prevIdx = ei.getAdjIdx(ptIdx, -1)
                pPrevIdx = ei.getAdjIdx(ptIdx, -2)
                if(prevIdx != None and  pPrevIdx != None):
                    return [ei.wsData[pPrevIdx][1], ei.wsData[prevIdx][1]]
                else:
                    nextIdx = ei.getAdjIdx(ptIdx, 1)
                    nNextIdx = ei.getAdjIdx(ptIdx, 2)
                    if(nextIdx != None and nNextIdx != None):
                        return [ei.wsData[nNextIdx][1], ei.wsData[nextIdx][1]]
        return self.getCurrLine()

    def getCurrLine(self):
        ei = self.editCurveInfo
        if(ei != None):
            ptIdx = ei.clickInfo['ptIdx']
            hdlIdx = ei.clickInfo['hdlIdx']
            pt0 = ei.wsData[ptIdx]
            clickLoc = ei.getClickLoc()
            if(clickLoc != None): return [pt0[1], clickLoc]
            if(hdlIdx in {0, 2}):
                return [pt0[1], pt0[hdlIdx]] # Current handle
            elif(hdlIdx == 1):
                adjIdx = ei.getAdjIdx(ptIdx, -1)
                if(adjIdx == None):
                    adjIdx = ei.getAdjIdx(ptIdx, 1)
                if(adjIdx == None): return [pt0[1]]
                else: return [ei.wsData[adjIdx][1], pt0[1]]
        return []

    def getRefLineOrig(self):
        ei = self.editCurveInfo
        refLine = self.getRefLine()
        if(ei != None and len(refLine) > 0):
            return refLine[-1]
        return None

    def getSelCo(self):
        if(self.editCurveInfo != None):
            return self.editCurveInfo.getSelCo()
        return None

    def getEditableCurveObjs(self):
        return [b for b in bpy.data.objects if isBezier(b) and b.visible_get() \
                and not b.hide_select and len(b.data.splines[0].bezier_points) > 1]

    def getSearchQueryInfo(self): # TODO: Simplify if possible
        queryInfo = {}
        for ci in self.selectCurveInfos:
            info = queryInfo.get(ci.obj)
            if(info == None):
                info = {}
                queryInfo[ci.obj] = info

            segPtIdxs = info.get(ci.splineIdx)
            if(segPtIdxs == None):
                # First is for seg search, second for handles
                segPtIdxs = [[], []]
                info[ci.splineIdx] = segPtIdxs

            for ptIdx in ci.ptSels.keys():
                sels = ci.ptSels[ptIdx]
                if(-1 in sels):
                    # TODO: Could be duplicate
                    adjIdx = ci.getAdjIdx(ptIdx)
                    segPtIdxs[0].append(ptIdx)
                    segPtIdxs[1].append(ptIdx)
                    segPtIdxs[1].append(adjIdx)
                else:
                    segPtIdxs[1].append(ptIdx)

        return queryInfo

    def getSelInfoObj(self, obj, splineIdx):
        for ci in self.selectCurveInfos:
            if(ci.obj == obj and ci.splineIdx == splineIdx):
                return ci
        return None

    # Delete selected segments and synchronize remaining selections
    # TODO: Way too complicated, maybe there exists a much simpler way to do this
    def delSelSegs(self):
        changed = False
        curveInfoList = sorted(self.selectCurveInfos, \
            key = lambda x: (x.objName, x.splineIdx))

        # Process one spline at a time
        for cIdx, c in enumerate(curveInfoList):
            c.resetHltInfo()

            spline = c.obj.data.splines[c.splineIdx]
            wasCyclic = spline.use_cyclic_u
            oldPtCnt = len(spline.bezier_points)

            changedSelMap = c.removeSegs()

            if(len(changedSelMap) == 0): continue
            changed = True

            # Shift all the splineIdxs after the changed one by spline incr count
            totalSplineIdxIncr = sum(x[0] for x in changedSelMap.values())

            # Order doesn't matter (different curveInfo)
            for i in range(cIdx + 1, len(curveInfoList)):
                if(curveInfoList[i].objName != c.objName):
                    break
                curveInfoList[i].splineIdx += totalSplineIdxIncr

            # TODO: Remove the try after sufficient testing (or better replacement)
            # Exception means no selected points will be deleted
            try:
                # Copy old selections as they will change
                oIdxs = sorted(c.ptSels.keys())
                ptSelsCopy = c.ptSels.copy()
                newSplineIdx = c.splineIdx

                c.resetPtSel()
                currCurveInfo = c

                # Reflects new selection after every seg removal
                modifiedSegIdxs = {idx:idx for idx in oIdxs}

                # First get the segment selections out of the way
                for i, segIdx in enumerate(sorted(changedSelMap.keys())):
                    ptSelsCopy[segIdx].remove(-1)
                    if(len(ptSelsCopy[segIdx]) == 0):
                        ptSelsCopy.pop(segIdx)

                # Each of the 'if' blocks in the loop iterate over all the selections
                # and update them iteratively for each seg removal from spline
                # The updated selected seg idx for each iteration is in modifiedSegIdxs

                # segIdx and oIdx don't change throughout
                # they always refer to the selections that were there before removal
                for i, segIdx in enumerate(sorted(changedSelMap.keys())):
                    splineIdxIncr = changedSelMap[segIdx][0]
                    segIdxIncr = changedSelMap[segIdx][1]

                    # This will be executed only once (if at all), at first iteration
                    if(wasCyclic and i == 0):
                        for j, oIdx in enumerate(oIdxs):
                            # First iteration, so no need to refer to modifiedSegIdxs
                            newSegIdx = oIdx + segIdxIncr
                            if(newSegIdx < 0): newSegIdx += oldPtCnt
                            modifiedSegIdxs[oIdx] = newSegIdx
                            if(ptSelsCopy.get(oIdx) != None):
                                currCurveInfo.ptSels[newSegIdx] = \
                                    ptSelsCopy[oIdx].copy()

                    # 'removed' segment at one of the either ends
                    elif(splineIdxIncr == 0):
                        ptCnt = len(c.obj.data.splines[newSplineIdx].bezier_points)
                        for j, oIdx in enumerate(oIdxs):
                            prevIdx = modifiedSegIdxs[oIdx]
                            # segIdxIncr: only two values possible: 0, -1
                            newSegIdx = prevIdx + segIdxIncr
                            # ~ if(currCurveInfo.ptSels.get(prevIdx) != None):
                                # ~ currCurveInfo.ptSels.pop(prevIdx)
                            if(ptSelsCopy.get(oIdx) != None and \
                                newSegIdx >=0 and newSegIdx < ptCnt):
                                currCurveInfo.ptSels[newSegIdx] = ptSelsCopy[oIdx].copy()
                            modifiedSegIdxs[oIdx] = newSegIdx

                    # Most likely condition
                    elif(splineIdxIncr > 0):
                        splineCnt = len(c.obj.data.splines)
                        prevCurveInfo = currCurveInfo
                        newSplineIdx += 1
                        # No overwriting since the higher splineIdxs already moved above
                        # But it's possible this spline was removed in subsequent
                        # iterations by removeSegs, so check...
                        if(newSplineIdx < splineCnt):
                            currCurveInfo = SelectCurveInfo(c.obj, newSplineIdx)
                            self.selectCurveInfos.add(currCurveInfo)
                        # idxs and prevCurve have to be updated so continue
                        else:
                            currCurveInfo = None # Fail fast

                        for oIdx in oIdxs:
                            prevIdx = modifiedSegIdxs[oIdx]
                            # If prevIdx itself is negative, this is previous to previous
                            # So won't change
                            if(prevIdx < 0): continue
                            newSegIdx = prevIdx + segIdxIncr
                            # newSegIdx negative... first part of the split spline
                            if(newSegIdx < 0 and ptSelsCopy.get(oIdx) != None):
                                prevCurveInfo.ptSels[prevIdx] = ptSelsCopy[oIdx].copy()
                            # newSegIdx positive... second part of the split spline
                            elif(ptSelsCopy.get(oIdx) != None and newSegIdx >=0):
                                if(newSplineIdx < splineCnt):
                                    currCurveInfo.ptSels[newSegIdx] = \
                                        ptSelsCopy[oIdx].copy()
                                if(prevCurveInfo.ptSels.get(prevIdx) != None):
                                    prevCurveInfo.ptSels.pop(prevIdx)
                            modifiedSegIdxs[oIdx] = newSegIdx

                    elif(splineIdxIncr < 0):
                        # This is not the same as c
                        # (could be a new spline added in between)
                        toRemList = [x for x in self.selectCurveInfos \
                            if x.splineIdx == newSplineIdx]
                        if(len(toRemList) > 0):
                            self.selectCurveInfos.remove(toRemList[0])

            except Exception as e:
                c.resetPtSel()

        return changed

    def mnSelect(self, opt):
        if(opt[0] == 'miSelAllSplines'):
            curves = self.getEditableCurveObjs()
            allCurveInfos = [SelectCurveInfo(curve, i) for curve in curves \
                for i in range(len(curve.data.splines))]
            for c in allCurveInfos:
                if(not c in self.selectCurveInfos):
                    self.selectCurveInfos.add(c)
            # ~ for c in allCurveInfos:
                # ~ for ptIdx in range(len(c.wsData)):
                    # ~ c.addSel(ptIdx, 1)
        else:
            h = ModalFlexiEditBezierOp.h
            self.selHltInfo(makeActive = True)
            for i, c in enumerate(self.selectCurveInfos):
                if(opt[0] == 'miSelObj'):
                    c.obj.select_set(True)
                    if(self.htlCurveInfo == None and i == len(self.selectCurveInfos)-1):
                        bpy.context.view_layer.objects.active = c.obj
                else:
                    for ptIdx in range(len(c.wsData)):
                        if(opt[0] == 'miSelSegs'): c.addSel(ptIdx, -1)
                        if(opt[0] == 'miSelBezPts'): c.addSel(ptIdx, 1)
                        if(opt[0] == 'miSelHdls' and not h): c.addSels(ptIdx, {0, 2})
                        if(opt[0] == 'miSelAll'):
                            c.addSels(ptIdx, {-1, 1}.union({0, 2} if not h else set()))
        self.htlCurveInfo = None

    def mnDeselect(self, opt):
        h = ModalFlexiEditBezierOp.h
        if(opt[0] == 'miDeselObj'):
            if(self.htlCurveInfo != None):
                self.htlCurveInfo.obj.select_set(False)
                self.htlCurveInfo = None
        for c in self.selectCurveInfos:
            if(opt[0] == 'miDeselObj'):
                c.obj.select_set(False)
            else:
                for ptIdx in range(len(c.wsData)):
                    if(opt[0] == 'miDeselSegs'): c.removeSel(ptIdx, -1)
                    if(opt[0] == 'miDeselBezPts'): c.removeSel(ptIdx, 1)
                    if(opt[0] == 'miDeselHdls' and not h): c.removeSels(ptIdx, {0, 2})
                    if(opt[0] == 'miDeselInvert'):
                        c.addSels(ptIdx, {-1, 1}.union({0, 2} if not h else set()), \
                            toggle = True)

    def mnSetHdlType(self, opt):
        if(ModalFlexiEditBezierOp.h): return
        self.selHltInfo(hltIdxs = {0, 1, 2}, selHdls = True, selEndPts = True)

        hdlType = opt[1].upper()
        for c in self.selectCurveInfos:
            # TODO: Support for Auto handles
            if(hdlType == 'AUTO' and c.hasShapeKey): continue
            for ptIdx in c.ptSels:
                sels = c.ptSels[ptIdx]
                for sel in sels:
                    bpt = c.obj.data.splines[c.splineIdx].bezier_points[ptIdx]
                    if(sel == 0):
                        bpt.handle_left_type = hdlType
                        # Following manual alignment required for shape keys
                        if(hdlType == 'ALIGNED'):
                            c.alignHandle(ptIdx, 0, allShapekeys = True)
                        if(hdlType == 'VECTOR'):
                            c.straightenHandle(ptIdx, 0, allShapekeys = True)
                            if(bpt.handle_right_type == 'ALIGNED'):
                                c.alignHandle(ptIdx, 2, allShapekeys = True)
                    if(sel == 2):
                        bpt.handle_right_type = hdlType
                        # Following manual alignment required for shape keys
                        if(hdlType == 'ALIGNED'):
                            c.alignHandle(ptIdx, 2, allShapekeys = True)
                        if(hdlType == 'VECTOR'):
                            c.straightenHandle(ptIdx, 2, allShapekeys = True)
                            if(bpt.handle_right_type == 'ALIGNED'):
                                c.alignHandle(ptIdx, 0, allShapekeys = True)
        bpy.ops.ed.undo_push()

    def exclToolRegion(self):
        return False

    def isEditing(self):
        return self.editCurveInfo != None

    def hasSelection(self):
        return len(self.selectCurveInfos) > 0

    # SnapParams object for bevel indicator
    def getBevelIndSnapParam(self, orig):
        # TODO: Maybe a more efficient way to find two major axes
        locs = [region_2d_to_location_3d(self.rmInfo.region, self.rmInfo.rv3d, \
            [[0, 0], [1000,1000]][i], Vector()) for i in range(2)]

        axisIdxs = [x[0] for x in sorted([(i, -abs(y)) for i, y in \
            enumerate(locs[1] - locs[0])], key = lambda z: z[1])]

        return SnapParams(self.snapper, enableSnap = False, \
            freeAxesN = sorted(axisIdxs[:2]), refLineOrig = orig, inEdit = True, \
                transType = 'GLOBAL', origType = 'REFERENCE', dispAxes = False, \
                    vec = Vector(), snapToPlane = True)

    def getNewDeltaPos(self, refreshStatus):
        if(self.xyLoc != None):
            loc = self.snapper.get3dLocSnap(self.rmInfo, \
                self.getBevelIndSnapParam(self.xyLoc))
            return loc - self.xyLoc
        else:
            return Vector()

    def getNewPos(self, refreshStatus):
        selCo = self.editCurveInfo.getSelCo()
        xySel = getCoordFromLoc(self.rmInfo.region, self.rmInfo.rv3d, selCo)
        if(self.xyPress != None):
            return self.snapper.get3dLocSnap(self.rmInfo, \
                SnapParams(self.snapper, vec = selCo, refreshStatus = refreshStatus, \
                    xyDelta = [self.xyPress[0] - xySel[0], self.xyPress[1] - xySel[1]]))
        else:
            return self.snapper.get3dLocSnap(self.rmInfo, \
                SnapParams(self.snapper, vec = selCo, refreshStatus = refreshStatus))

    def confirmCurveOp(self):
        if(self.bevelMode or self.subdivCnt > 0):
            changed = False
            for c in self.selectCurveInfos:
                if(self.bevelMode):
                    changed = c.bevelPts(self.bevelCnt, self.getNewDeltaPos(False)) \
                        or changed
                else:
                    changed = c.subdivSeg(self.subdivCnt) or changed
                    c.resetPtSel()
            if(changed): bpy.ops.ed.undo_push()
            self.bevelMode = False
            self.subdivCnt = 0
            self.xyLoc = None
            self.bglDrawMgr.resetLineInfo('bevelLine')
            bpy.context.window.cursor_set("DEFAULT")
            self.snapper.resetSnap()
            self.refreshDisplaySelCurves()
            return True
        return False

    def getHltIdxFromRes(self, resType, otherInfo):
        # return 1:bez pt, -1:segloc, 0:lefthandle, 2:righthandle (like ptSels format)
        if(resType in {'SegLoc', 'CurveLoc'}): return -1
        else: return otherInfo

    # Select highlighted element for cases where op needs to be initiated without
    # mouse click (just by mouse hover)
    def selHltInfo(self, hltIdxs = None, makeActive = False, \
        selHdls = False, selEndPts = False):
        hltCurve = self.htlCurveInfo
        if(hltCurve != None and \
            all(sum(len(sel) for sel in c.ptSels.values() if hltIdxs == None or \
                len(sel.intersection(hltIdxs)) > 0) == 0 for c in self.selectCurveInfos)):
                hltIdx = hltCurve.hltInfo['hltIdx']
                if(hltIdxs == None or hltIdx in hltIdxs):
                    currIdx = hltCurve.hltInfo['ptIdx']
                    hltCurve.ptSels[currIdx] = {hltIdx}
                    ptIdxs = [currIdx]
                    if(selHdls and hltIdx == -1):
                        nextIdx = hltCurve.getAdjIdx(currIdx)
                        if(nextIdx != None):
                            hltCurve.ptSels[currIdx].add(1)
                            hltCurve.ptSels[nextIdx] = {1}
                            ptIdxs.append(nextIdx)
                        hltIdx = 1
                    for ptIdx in ptIdxs:
                        if(selHdls and hltIdx == 1):
                            hltCurve.ptSels[ptIdx] = hltCurve.ptSels[ptIdx].union({0, 2})
                        else:
                            hltCurve.ptSels[ptIdx] = {hltIdx}
                    self.selectCurveInfos.add(hltCurve)
                    if(makeActive):
                        bpy.context.view_layer.objects.active = self.htlCurveInfo.obj

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo
        metakeys = self.snapper.getMetakeys()
        alt = metakeys[0]
        ctrl = metakeys[1]
        shift = metakeys[2]
        opMode = self.bevelMode or self.subdivCnt > 0

        if(snapProc): retVal = {"RUNNING_MODAL"}
        else: retVal = {'PASS_THROUGH'}

        if(not snapProc and event.type == 'ESC'):
            # Escape processing sequence:
            # 1) Come out bevel mode
            # 2) Come out of snapper / snapdigits (not 1)
            # 3) Reset position if captured (double click) (not 2)
            # 4) Reset selection if captured and position already reset (not 3)
            if(event.value == 'RELEASE'):
                if(self.editCurveInfo == None):
                    if(self.subdivCnt > 0):
                        self.subdivCnt = 0
                        self.refreshDisplaySelCurves()
                    elif(self.bevelMode):
                        self.bevelMode = False
                        self.bglDrawMgr.resetLineInfo('bevelLine')
                        bpy.context.window.cursor_set("DEFAULT")
                        self.snapper.resetSnap()
                        self.refreshDisplaySelCurves()
                    else:
                        self.reset()
                        ModalFlexiEditBezierOp.resetDisplay()
                else:
                    if(self.capture and self.snapper.lastSelCo != None and
                        not vectCmpWithMargin(self.snapper.lastSelCo, \
                            self.editCurveInfo.getSelCo())):
                        self.snapper.lastSelCo = self.editCurveInfo.getSelCo()
                    else:
                        self.capture = False
                        self.editCurveInfo = None
                        self.snapper.resetSnap()
                    self.refreshDisplaySelCurves()
            return {"RUNNING_MODAL"}

        if(self.bevelMode or ((ctrl or alt) and (self.editCurveInfo == None or \
            (self.pressT != None and not self.click)))):
            bpy.context.window.cursor_set("CROSSHAIR")
        else:
            bpy.context.window.cursor_set("DEFAULT")

        if(not opMode and \
            FTHotKeys.isHotKey(FTHotKeys.hkSplitAtSel, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                selPtMap = {}
                # ~ self.selHltInfo(hltIdxs = {-1}, selEndPts = True)
                for c in self.selectCurveInfos:
                    if(selPtMap.get(c.obj) == None):
                        selPtMap[c.obj] = {}
                    ptIdxs = {p for p in c.ptSels.keys() if 1 in c.ptSels[p] \
                        or -1 in c.ptSels[p]}
                    if(len(ptIdxs) > 0):
                        selPtMap[c.obj][c.splineIdx] = ptIdxs
                        ptIdxs = [p for p in c.ptSels.keys() if -1 in c.ptSels[p]]
                        selPtMap[c.obj][c.splineIdx].update({c.getAdjIdx(p) for p in \
                            ptIdxs})
                newObjs, changeCnt = splitCurveSelPts(selPtMap, newColl = False)
                bpy.ops.ed.undo_push()
                self.reset()
                for o in newObjs:
                    for i in range(len(o.data.splines)):
                        self.selectCurveInfos.add(SelectCurveInfo(o, i))
            return {"RUNNING_MODAL"}

        if(not opMode and \
            FTHotKeys.isHotKey(FTHotKeys.hkToggleDrwEd, event.type, metakeys)):
            if(event.value == 'RELEASE'):
                self.reset()
                bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname)
                # bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')
            return {"RUNNING_MODAL"}

        if(not opMode and FTHotKeys.isHotKey(FTHotKeys.hkBevelPt, event.type, metakeys)):
            # Allow beveling seg / pt just with mouse hover
            self.selHltInfo(hltIdxs = {1, -1})
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    changed = False
                    for c in self.selectCurveInfos:
                        # short-circuit fine (no change in isBevelabel)
                        changed = changed or c.isBevelabel(rmInfo.rv3d)
                    self.bevelMode = changed
                    if(changed):
                        self.xyPress = rmInfo.xy[:]
                        self.xyLoc = self.snapper.get3dLocSnap(rmInfo, \
                            self.getBevelIndSnapParam(orig = None))
                        bpy.context.window.cursor_set("CROSSHAIR")
                        self.bevelCnt = FTProps.defBevelFact
                        self.refreshDisplaySelCurves()
                        self.htlCurveInfo = None
                return {"RUNNING_MODAL"}

        if(not opMode and \
            FTHotKeys.isHotKey(FTHotKeys.hkUniSubdiv, event.type, metakeys)):
            self.selHltInfo(hltIdxs = {-1})
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    changed = False
                    for c in self.selectCurveInfos:
                        changed = c.initSubdivMode(rmInfo.rv3d) or changed
                    if(changed):
                        self.subdivCnt = 2
                        self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        confirmed = False
        if(not snapProc and event.type in {'SPACE', 'RET'}):
            if(self.bevelMode or self.subdivCnt > 0):
                if(event.value == 'RELEASE'):
                    self.confirmCurveOp()
                return {"RUNNING_MODAL"}
            elif(self.editCurveInfo != None):
                confirmed = True

        elif(not snapProc and event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', \
            'NUMPAD_PLUS', 'NUMPAD_MINUS','PLUS', 'MINUS'}):
            if(len(self.selectCurveInfos) > 0 and \
                (self.subdivCnt > 0 or self.bevelMode)):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return {'RUNNING_MODAL'}
                elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                    if(self.bevelMode and self.bevelCnt > FTProps.minBevelFact):
                        self.bevelCnt -= FTProps.bevelIncr
                    elif(self.subdivCnt > 2):
                        self.subdivCnt -= 1
                elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                    if(self.bevelMode and self.bevelCnt < FTProps.maxBevelFact):
                        self.bevelCnt += FTProps.bevelIncr
                    elif(self.subdivCnt < 100 and self.subdivCnt > 0):
                        self.subdivCnt += 1
                self.refreshDisplaySelCurves()
                return {'RUNNING_MODAL'}

        if(FTHotKeys.isHotKey(FTHotKeys.hkToggleHdl, event.type, metakeys)):
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    ModalFlexiEditBezierOp.h = not ModalFlexiEditBezierOp.h
                    self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        if(not opMode and \
            FTHotKeys.isHotKey(FTHotKeys.hkDelPtSeg, event.type, metakeys)):
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    changed = self.delSelSegs()
                    for c in self.selectCurveInfos:
                        c.resetHltInfo()
                        changed = c.removeNode() or changed
                        changed = c.straightenSelHandles() or changed

                    if(changed):
                        # will be taken care by depsgraph?
                        self.updateAfterGeomChange()
                        bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if(not opMode and \
            FTHotKeys.isHotKey(FTHotKeys.hkAlignHdl, event.type, metakeys)):
            self.selHltInfo(hltIdxs = {0, 1, 2}, selHdls = True, selEndPts = True)
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    changed = False
                    for c in self.selectCurveInfos:
                        changed = c.alignSelHandles() or changed #selected node
                    if(changed): bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if(not snapProc and not self.capture \
            and event.type == 'LEFTMOUSE' and event.value == 'PRESS'):

            if(self.subdivCnt > 0 or self.bevelMode):
                return {'RUNNING_MODAL'}

            for ci in self.selectCurveInfos.copy():
                if(len(ci.ptSels) == 0): self.selectCurveInfos.remove(ci)

            self.xyPress = rmInfo.xy[:]
            coFind = Vector(rmInfo.xy).to_3d()

            objs = self.getEditableCurveObjs()

            selObjInfos = self.getSearchQueryInfo()

            #TODO: Move to Snapper?
            searchResult = getClosestPt2d(rmInfo.region, rmInfo.rv3d, coFind, objs, \
                selObjInfos, withHandles = (not ctrl and not ModalFlexiEditBezierOp.h))

            if(searchResult != None):
                resType, obj, splineIdx, segIdx, otherInfo = searchResult

                ci = self.getSelInfoObj(obj, splineIdx)

                if(ci == None):
                    ci = EditCurveInfo(obj, splineIdx)
                    self.selectCurveInfos.add(ci)
                elif(type(ci) != EditCurveInfo):
                    self.selectCurveInfos.remove(ci)
                    ci = EditCurveInfo(obj, splineIdx, ci.ptSels)
                    self.selectCurveInfos.add(ci)

                ptIdx = segIdx
                clickLoc = None
                if(resType == 'SelHandles'):
                    hdlIdx = otherInfo
                elif(resType == 'CurveBezPt'):
                    hdlIdx = 1
                else:#if(resType == 'SegLoc'):
                    hdlIdx = -1
                    searchResult = getClosestPt2dWithinSeg(rmInfo.region, rmInfo.rv3d, \
                        coFind, selObj = obj, selSplineIdx = splineIdx, \
                            selSegIdx = segIdx, withHandles = False, withBezPts = False)

                    # ~ if(searchResult != None): #Must never be None
                    resType, obj, splineIdx, segIdx, otherInfo = searchResult
                    clickLoc = otherInfo
                # ~ ci.addSel(ptIdx, hdlIdx)
                ci.setClickInfo(segIdx, hdlIdx, clickLoc)
                # ~ if(ci._t == None): ci = None

                self.editCurveInfo = ci
                ci.setHltInfo(ptIdx = segIdx, \
                    hltIdx = self.getHltIdxFromRes(resType, otherInfo))
                # ~ self.pressT = time.time()
                return {'RUNNING_MODAL'}

            if(not shift):
                self.reset()

            return retVal

        if(confirmed or self.snapper.digitsConfirmed or \
            (event.type == 'LEFTMOUSE' and event.value == 'RELEASE')):
            if(self.confirmCurveOp()):
                return {"RUNNING_MODAL"}

            if(self.editCurveInfo == None):
                return retVal

            ei = self.editCurveInfo
            tm = time.time()

            if(self.doubleClick):
                self.capture = True
            else:
                if(self.click and not self.capture):
                    ptIdx = ei.clickInfo['ptIdx']
                    hdlIdx = ei.clickInfo['hdlIdx']
                    pt = ei.wsData[ptIdx]
                    if(ctrl and ei.clickInfo['hdlIdx'] == -1):
                        if(shift): handleType = 'ALIGNED'
                        elif(alt): handleType = 'VECTOR'
                        else: handleType = 'FREE'

                        changed = ei.insertNode(handleType)
                        bpy.ops.ed.undo_push()
                        ModalFlexiEditBezierOp.resetDisplay()
                    elif(alt and (hdlIdx == -1 or (hdlIdx == 1 and hasAlignedHandles(pt)))):
                        if(hdlIdx == -1):
                            pts = ei.getSegPts(ei.clickInfo['ptIdx'])
                            seg = [pts[0][1], pts[0][2], pts[1][0], pts[1][1]]
                            t = ei.clickInfo['t']
                            tangent = getTangentAtT(*seg, t)
                            fact = tangent.normalized()
                            clickLoc = ei.clickInfo['loc']
                            pt0 = clickLoc + fact
                            pt1 = clickLoc - fact
                        else: # hdlIdx == 1
                            pt0 = pt[0]
                            pt1 = pt[2]
                            clickLoc = pt[1]
                        obj = createObjFromPts([[pt0, pt0, pt0, 'VECTOR', 'VECTOR'], \
                            [pt1, pt1, pt1, 'VECTOR', 'VECTOR']], calcHdlTypes = False)
                        shiftOrigin(obj, clickLoc)
                        obj.location = clickLoc
                        # ~ bpy.context.evaluated_depsgraph_get().update()
                        obj.select_set(True)
                        self.selectCurveInfos = {SelectCurveInfo(obj, 0)}
                        bpy.ops.ed.undo_push()
                    # Gib dem Benutzer Zeit zum Atmen!
                    else:
                        if(not shift or ctrl):
                            for ci in self.selectCurveInfos.copy():
                                if(ci != ei): self.selectCurveInfos.remove(ci)
                            ei.resetPtSel()
                        ei.addSel(ptIdx, hdlIdx, toggle = True)
                        if(hdlIdx == 1):
                            if(ptIdx in ei.ptSels and 1 in ei.ptSels[ptIdx]):
                                ei.addSel(ptIdx, 0, toggle = False)
                                ei.addSel(ptIdx, 2, toggle = False)
                            else:
                                ei.removeSel(ptIdx, 0)
                                ei.removeSel(ptIdx, 2)

                        self.selectCurveInfos.add(ei)
                        # ~ self.refreshDisplaySelCurves()
                else:
                    ei.moveSeg(self.getNewPos(refreshStatus = False))
                    # ~ self.updateAfterGeomChange() # TODO: Really needed?
                    bpy.ops.ed.undo_push()

                self.capture = False
                self.editCurveInfo = None
                self.refreshDisplaySelCurves()
                self.snapper.resetSnap()

            # ~ self.pressT = None
            return {"RUNNING_MODAL"}

        elif(snapProc or event.type == 'MOUSEMOVE'):
            segDispInfos = None
            bptDispInfos = None
            ei = self.editCurveInfo
            locOnCurve = None # For debug

            if(self.bevelMode):
                loc = self.snapper.get3dLocSnap(rmInfo, \
                    self.getBevelIndSnapParam(self.xyLoc))
                lineCol = (1, 1, 0, 1)
                self.bglDrawMgr.addLineInfo('bevelLine', 1, [lineCol], \
                    [self.xyLoc, loc])

            elif(self.subdivCnt > 0):
                pass

            # ei != None taken care by refreshDisplaySelCurves(refreshPos = True)
            elif(ei == None):
                self.htlCurveInfo = None
                # ~ coFind = Vector(rmInfo.xy).to_3d()
                coFind = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, \
                    self.snapper.get3dLocSnap(rmInfo, \
                        SnapParams(self.snapper, enableSnap = False))).to_3d()

                objs = self.getEditableCurveObjs()

                #Sel obj: low res (highlight only seg)
                selObjInfos = self.getSearchQueryInfo()

                #TODO: Move to Snapper
                searchResult = getClosestPt2d(rmInfo.region, rmInfo.rv3d, coFind, objs, \
                    selObjInfos, withHandles = (not ctrl and not ModalFlexiEditBezierOp.h))

                for c in self.selectCurveInfos: c.resetHltInfo()
                if(searchResult != None):
                    resType, obj, splineIdx, segIdx, otherInfo = searchResult
                    ci = self.getSelInfoObj(obj, splineIdx)

                    if(resType not in {'SelHandles', 'CurveBezPt'}):
                        locOnCurve = otherInfo
                    if(ci == None):
                        ci = SelectCurveInfo(obj, splineIdx)
                        ci.setHltInfo(ptIdx = segIdx, \
                            hltIdx = self.getHltIdxFromRes(resType, otherInfo))
                        segDispInfos, bptDispInfos = \
                            ci.getDisplayInfos(ModalFlexiEditBezierOp.h, \
                                subdivCnt = self.subdivCnt, bevelCnt = self.bevelCnt)
                    else:
                        ci.setHltInfo(ptIdx = segIdx, \
                            hltIdx = self.getHltIdxFromRes(resType, otherInfo))
                    self.htlCurveInfo = ci
            self.refreshDisplaySelCurves(segDispInfos, bptDispInfos, \
                locOnCurve, refreshPos = True)

            return retVal

        if(snapProc or opMode):
            self.refreshDisplaySelCurves(refreshPos = True)
            return {'RUNNING_MODAL'}
        else:
            return retVal

###################### Global Params ######################

def getConstrAxisTups(scene = None, context = None):
    axesMap = {\
       0: ('NONE', 'None', "Constrain only on hotkey event"), \
       1: ('-X', 'X', "Constrain to only X axis"), \
       2: ('-Y', 'Y', "Constrain to only Y axis"), \
       3: ('-Z', 'Z', "Constrain to only Z axis"), \
       4: ('shift-Z', 'XY', "Constrain to XY plane"), \
       5: ('shift-Y', 'XZ', "Constrain to XZ plane"), \
       6: ('shift-X', 'YZ', "Constrain to YZ plane"), \
      }
    transType = bpy.context.window_manager.bezierToolkitParams.snapOrient

    if(transType in {'AXIS', 'GLOBAL', 'OBJECT', 'FACE'}): keyset = range(0, 7)
    elif(transType in {'VIEW', 'REFERENCE', 'CURR_POS'}): keyset = [0] + [i for i in range(4, 7)]

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

    intersectNonactive: BoolProperty(name="Only Non-active", \
        description="Action is not performed on active curve but " + \
            "only other selected curves", \
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

    intersectExpanded: BoolProperty(name = "Intersect Curves", default = False)

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

    snapOrient: EnumProperty(name = 'Orientation',#"Align contrained axes and snap angle to",
        items = (('GLOBAL', 'Global Axes', "Orient to world space"), \
        ('REFERENCE', 'Reference Line', "Orient to preceding segment or opposite handle"),
        ('CURR_POS', 'Current Segment', "Orient to current segment or current handle"),
        ('AXIS', 'Custom Axes', "Orient to custom axis (if available)"), \
        ('VIEW', 'View', "Orient to window"), \
        ('OBJECT', 'Active Object', "Orient to local space of active object"),
        ('FACE', 'Selected Object Face', \
         "Orient to normal of face of selected object under mouse pointer ")),
        default = 'GLOBAL',
        description='Orientation for Draw / Edit')


    snapOrigin: EnumProperty(name = 'Origin',#"Align contrained axes and snap angle to",
        items = (('GLOBAL', 'Global Origin', \
          "Draw / Edit with reference to global origin"), \
        ('CURSOR', '3D Cursor Location', \
          "Draw / Edit with reference to 3D Cursor location"), \
        ('AXIS', 'Custom Axis Start', \
          "Draw / Edit with reference to starting point of custom axis"), \
        ('REFERENCE', 'Reference Line Point', \
          "Draw / Edit with reference to the appropriate reference line point"), \
        ('OBJECT', 'Active Object Location', \
          "Draw / Edit with reference to active object location"), \
        ('FACE', 'Selected Object Face', \
          "Draw / Edit with reference to the center of " + \
          "Selected object face under mouse pointer"),
        ('CURR_POS', 'Current Position', \
          "Edit with reference to the current mouse position")), \
        default = 'REFERENCE',
        description='Origin for Draw / Edit')

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
    bl_context_mode='PAINT_GPENCIL'

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

def showSnapToPlane(params):
    return (params.snapOrient not in {'VIEW', 'REFERENCE', 'CURR_POS'} and \
            hasattr(params, 'constrAxes') and params.constrAxes.startswith('shift'))

def drawSettingsFT(self, context):
    params = bpy.context.window_manager.bezierToolkitParams
    self.layout.use_property_split = True
    self.layout.row(align=True).template_header()
    from bl_ui.space_toolsystem_common import ToolSelectPanelHelper
    toolHeader = ToolSelectPanelHelper.draw_active_tool_header(
        context, self.layout,
        tool_key=('VIEW_3D', context.mode),
    )

    toolObj = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)
    toolGP = context.workspace.tools.from_space_view3d_mode('PAINT_GPENCIL', create = False)

    self.layout.use_property_decorate = True

    gpMode = (context.mode == 'PAINT_GPENCIL' and \
            toolGP.idname == FlexiGreaseBezierTool.bl_idname)
            # toolGP.idname == 'flexi_bezier.grease_draw_tool')
    drawMode = (context.mode == 'OBJECT' and \
                toolObj.idname  == FlexiDrawBezierTool.bl_idname)
                # toolObj.idname  == 'flexi_bezier.draw_tool')
    if(drawMode or gpMode):
        if(gpMode):
            brush = context.scene.tool_settings.gpencil_paint.brush
            self.layout.prop(brush, 'size', text ='')
            self.layout.prop(brush.gpencil_settings, 'pen_strength', text ='')
        self.layout.prop(params, "drawObjType", text = '')
        if(params.drawObjType != 'BEZIER'):
            if(params.drawObjType == 'MATH'):
                self.layout.prop(params, "mathFnType", text = '')
                if(params.mathFnType == 'PARAMETRIC'):
                    self.layout.prop(params, "drawMathFnParametric1", text = '')
                    self.layout.prop(params, "drawMathFnParametric2", text = '')
                else:
                    self.layout.prop(params, "drawMathFn", text = '')
            self.layout.prop(params, "drawObjMode", text = '')
            if(params.drawObjType == 'ELLIPSE'):
                self.layout.prop(params, "drawStartAngle", text = '')
            if(params.drawObjType in {'POLYGON', 'STAR'}):
                self.layout.prop(params, "drawSides", text = '')
            if(params.drawObjType == 'STAR'):
                self.layout.prop(params, "drawStarOffset", text = '')
            if(params.drawObjType != 'RECTANGLE' and params.drawObjType != 'MATH'):
                self.layout.prop(params, "drawAngleSweep", text = '')

    self.layout.prop(params, "snapOrient", text = '')
    self.layout.prop(params, "snapOrigin", text = '')
    self.layout.prop(params, "constrAxes", text = '')
    if(params.constrAxes not in [a[0] for a in getConstrAxisTups()]):
        params.constrAxes = 'NONE'

    # Only available for planes not axis
    if(showSnapToPlane(params)):
        self.layout.prop(params, "snapToPlane")

    self.layout.prop(params, "axisScale", text = '')

    # if((context.mode == 'OBJECT' and toolObj.idname  == 'flexi_bezier.draw_tool')):
    if((context.mode == 'OBJECT' and toolObj.idname  == FlexiDrawBezierTool.bl_idname)):
        self.layout.prop(params, "copyPropsObj", text = '')


# ****************** Configurations In User Preferences ******************

def updatePanel(self, context):
    try:
        panel = BezierUtilsPanel
        if "bl_rna" in panel.__dict__:
            bpy.utils.unregister_class(panel)
    except Exception as e:
        print("BezierUtils: Unregistering Panel has failed", e)
        return

    try:
        panel.bl_category = context.preferences.addons[__name__].preferences.category
        bpy.utils.register_class(panel)

    except Exception as e:
        print("BezierUtils: Updating Panel locations has failed", e)
        panel.bl_category = 'Tool'
        bpy.utils.register_class(panel)


class ResetDefaultPropsOp(bpy.types.Operator):
    bl_idname = "object.reset_default_props"
    bl_label = "Reset Preferences"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Reset all property values to default"

    def execute(self, context):
        context.preferences.addons[__name__].preferences.category = 'Tool'
        FTProps.updatePropsPrefs(context, resetPrefs = True)
        return {'FINISHED'}


class ResetDefaultHotkeys(bpy.types.Operator):
    bl_idname = "object.reset_default_hotkeys"
    bl_label = "Reset Keymap"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Reset all hotkey values to default"

    def execute(self, context):
        FTHotKeys.updateHKPropPrefs(context, reset = True)
        return {'FINISHED'}


class BezierUtilsPreferences(AddonPreferences):
    bl_idname = __name__

    category: StringProperty(
            name = "Tab Category",
            description = "Choose a name for the category of the panel",
            default = "Tool",
            update = updatePanel
    )

    lineWidth: FloatProperty(
            name = "Line Thickness",
            description = "Thickness of segment & handle Lines of Flexi Draw and Edit",
            default = 1.5,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    drawPtSize: FloatProperty(
            name = "Handle Point Size",
            description = "Size of Flexi Draw and Edit Bezier handle points",
            default = 5,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    editSubdivPtSize: FloatProperty(
            name = "Uniform Subdiv Point Size",
            description = "Size of point marking subdivisions",
            default = 6,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    greaseSubdivPtSize: FloatProperty(
            name = "Flexi Grease Res Point Size",
            description = "Size of point marking resoulution in Flexi Grease Tool",
            default = 4,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    markerSize: FloatProperty(
            name = "Marker Size",
            description = "Size of Flexi Draw and Mark Starting Vertices",
            default = 6,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    axisLineWidth: FloatProperty(
            name = "Axis Line Thickness",
            description = "Thickness of Axis Lines for snapping & locking",
            default = 0.25,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    snapPtSize: FloatProperty(
            name = "Snap Point Size",
            description = "Size of snap point indicator",
            default = 5,
            min = 0.1,
            max = 20,
            update = FTProps.updateProps
    )

    defBevelFact: FloatProperty(
            name = "Default Bevel Factor",
            description = "Initial factor used for beveling points",
            default = 4,
            update = FTProps.updateProps
    )

    maxBevelFact: FloatProperty(
            name = "Maximum Bevel Factor",
            description = "Maximum bevel factor",
            default = 15,
            update = FTProps.updateProps
    )

    minBevelFact: FloatProperty(
            name = "Minimum Bevel Factor",
            description = "Minimum bevel factor",
            default = -15,
            update = FTProps.updateProps
    )

    bevelIncr: FloatProperty(
            name = "Bevel Factor Step",
            description = "Increment / decrement step of bevel factor",
            default = .5,
            update = FTProps.updateProps
    )

    liveUpdate: BoolProperty(
            name="Flexi Edit Live Update", \
            description='Update underlying object with editing', \
            default = False,
            update = FTProps.updateProps
    )

    dispCurveRes: FloatProperty(
            name = "Display Curve Resolution",
            description = "Segment divisions per pixel",
            default = .4,
            min = 0.001,
            max = 1,
            update = FTProps.updateProps
    )

    snapDist: FloatProperty(
            name = "Snap Distance",
            description = "Snapping distance (range) in pixels",
            default = 20,
            min = 1,
            max = 200,
            update = FTProps.updateProps
    )

    dispSnapInd: BoolProperty(
            name="Snap Indicator", \
            description='Display indicator when pointer within snapping range', \
            default = True,
            update = FTProps.updateProps
    )

    dispAxes: BoolProperty(
            name="Orientation / Origin Axis", \
            description='Display axes for selected orientation / origin', \
            default = False,
            update = FTProps.updateProps
    )

    numpadEntry: BoolProperty(
            name="Allow Numpad Entry", \
            description='Allow numpad entries as keyboard input', \
            default = False,
            update = FTProps.updateProps
    )

    showKeyMap: BoolProperty(
            name="Display Keymap", \
            description='Display Keymap When Flexi Tool Is Active', \
            default = True,
            update = FTProps.updateProps
    )

    keyMapNextToTool: BoolProperty(
            name="Display Next to Toolbar", \
            description='Display Keymap Next to Toolbar', \
            default = True,
            update = FTProps.updateProps
    )

    keyMapFontSize: IntProperty(
            name = "Keymap Font Size",
            description = "Font size of keymap text",
            default = 10,
            min = 1,
            max = 2000,
            update = FTProps.updateProps
    )

    mathFnTxtFontSize: IntProperty(
            name = "Math Function Font Size",
            description = "Font size of math function equations (Flexi Draw - Math Fn mode)",
            default = 20,
            min = 1,
            max = 2000,
            update = FTProps.updateProps
    )

    keyMapLocX: IntProperty(
            name = "Keymap Location X",
            description = "Horizontal starting position of keymap display",
            default = 10,
            min = 1,
            max = 2000,
            update = FTProps.updateProps
    )

    keyMapLocY: IntProperty(
            name = "Keymap Font Size",
            description = "Vertical starting position of keymap display",
            default = 10,
            min = 1,
            max = 2000,
            update = FTProps.updateProps
    )

    colDrawSelSeg: bpy.props.FloatVectorProperty(
        name="Selected Draw / Edit  Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.6, .8, 1, 1), \
        description = 'Color of the segment being drawn / edited',
        update = FTProps.updateProps
    )

    colDrawNonHltSeg: bpy.props.FloatVectorProperty(
        name="Adjacent Draw / Edit Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.1, .4, .6, 1), \
        description = 'Color of the segment adjacent to the' + \
        'one being drawn / edited', update = FTProps.updateProps
    )

    colDrawHltSeg: bpy.props.FloatVectorProperty(
        name="Highlighted Edit Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.2, .6, .9, 1), \
        description = 'Color of the segment under mouse curser in Flexi Edit', \
        update = FTProps.updateProps
    )

    colDrawMarker: bpy.props.FloatVectorProperty(
        name="Marker", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.6, .8, 1, 1), \
        description = 'Color of the marker', update = FTProps.updateProps
    )

    colGreaseSelSeg: bpy.props.FloatVectorProperty(
        name="Selected Grease Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(0.2, .8, 0.2, 1), \
        description = 'Color of the segment being drawn', \
        update = FTProps.updateProps
    )

    colGreaseNonHltSeg: bpy.props.FloatVectorProperty(
        name="Adjacent Grease Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (0.2, .6, 0.2, 1), \
        description = 'Color of the segment adjacent to the one being drawn', \
        update = FTProps.updateProps
    )

    colGreaseMarker: bpy.props.FloatVectorProperty(
        name="Marker", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (0.2, .8, 0.2, 1), \
        description = 'Color of the marker', update = FTProps.updateProps
    )

    colHdlFree: bpy.props.FloatVectorProperty(
        name="Free Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.6, .05, .05, 1), \
        description = 'Free handle color in all Flexi Tools', \
        update = FTProps.updateProps
    )

    colHdlVector: bpy.props.FloatVectorProperty(
        name="Vector Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.4, .5, .2, 1), \
        description = 'Vector handle color in all Flexi Tools', \
        update = FTProps.updateProps
    )

    colHdlAligned: bpy.props.FloatVectorProperty(
        name="Aligned Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(1, .3, .3, 1), \
        description = 'Aligned handle color in all Flexi Tools', \
        update = FTProps.updateProps
    )

    colHdlAuto: bpy.props.FloatVectorProperty(
        name="Auto Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(.8, .5, .2, 1), \
        description = 'Auto handle color in all Flexi Tools', \
        update = FTProps.updateProps
    )

    colSelTip: bpy.props.FloatVectorProperty(
        name="Selected Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (.2, .7, .3, 1), \
        description = 'Color of the selected Bezier or handle point', \
        update = FTProps.updateProps
    )

    colHltTip: bpy.props.FloatVectorProperty(
        name="Highlighted Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (.8, 1, .8, 1), \
        description = 'Color of Bezier or handle point under mouse pointer', \
        update = FTProps.updateProps
    )

    colBezPt: bpy.props.FloatVectorProperty(
        name="Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (1, 1, 0, 1), \
        description = 'Color of nonselected Bezier point', \
        update = FTProps.updateProps
    )

    colHdlPtTip: bpy.props.FloatVectorProperty(
        name="Handle Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (.7, .7, 0, 1), \
        description = 'Color of nonselected handle point', \
        update = FTProps.updateProps
    )

    colAdjBezTip: bpy.props.FloatVectorProperty(
        name="Adjacent Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (.1, .1, .1, 1), \
        description = 'Color of Bezier points of adjacent segments', \
        update = FTProps.updateProps
    )

    colEditSubdiv: bpy.props.FloatVectorProperty(
        name="Uniform Subdiv Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (.3, 0, 0, 1), \
        description = 'Color of point marking subdivisions', \
        update = FTProps.updateProps
    )

    colGreaseSubdiv: bpy.props.FloatVectorProperty(
        name="Resolution Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (1, .3, 1, 1), \
        description = 'Color of point marking curve resolution', \
        update = FTProps.updateProps
    )

    colGreaseBezPt: bpy.props.FloatVectorProperty(
        name="Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (1, .3, 1, 1), \
        description = 'Color of Bezier point', \
        update = FTProps.updateProps
    )

    colKeymapText: bpy.props.FloatVectorProperty(
        name="Keymap Description Text Color", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(1.0, 1.0, 1.0, 1.0), \
        description = 'Color of keymap description text',
        update = FTProps.updateProps
    )

    colKeymapKey: bpy.props.FloatVectorProperty(
        name="Keymap Key Text Color", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default=(0.0, 1.0, 1.0, 1.0), \
        description = 'Color of keymap key text',
        update = FTProps.updateProps
    )

    colMathFnTxt: bpy.props.FloatVectorProperty(
        name="Math Function Text Color", subtype="COLOR", size=4, min=0.0, max=1.0,\
        default = (0.6, 1.0, 0.03, 1.0), \
        description = "Color of math function equations (Flexi Draw - Math Fn mode)",
        update = FTProps.updateProps
    )

############################ Hotkeys ###############################
    for i, keySet in enumerate([FTHotKeys.drawHotkeys, \
        FTHotKeys.editHotkeys, FTHotKeys.commonHotkeys]):
        for j, keydata in enumerate(keySet):
            exec(FTHotKeys.getHKFieldStr(keydata, addMeta = (not keydata.isExclusive)))
            expStr = keydata.id + 'Exp'
            exec(expStr + ': BoolProperty(name="' + expStr + '", default = False)')
        expStr = 'hotKeySet'+ str(i)+'Exp'
        exec(expStr + ': BoolProperty(name="' + expStr + '", default = False)')

    for i, keydata in enumerate(FTHotKeys.snapHotkeys):
        keydataMeta = FTHotKeys.snapHotkeysMeta[i]
        exec(FTHotKeys.getMetaHKFieldStr(keydataMeta))
        exec(FTHotKeys.getHKFieldStr(keydata, addMeta = False))
        expStr = keydata.id + 'Exp'
        exec(expStr + ': BoolProperty(name="' + expStr + '", default = True)')
    hkSnapExp: BoolProperty(name="Snap Hotkey", default = False)

    colSizeExp: BoolProperty(name="Color & Size", default = False)
    keymapExp: BoolProperty(name="Keymap", default = False)
    snapOptExp: BoolProperty(name="Snapping", default = False)
    bevelingExp: BoolProperty(name="Beveling", default = False)
    othPrefExp: BoolProperty(name="Other Options", default = False)

    elemDimsExp: BoolProperty(name="Draw Dimensions", default = False)
    drawColExp: BoolProperty(name="Draw Colors", default = False)
    greaseColExp: BoolProperty(name="Grease Pencil Colors", default = False)
    handleColExp: BoolProperty(name="Handle Colors", default = False)

    def draw(self, context):
        layout = self.layout
        col = layout.column().split()
        col.label(text="Tab Category:")
        col.prop(self, "category", text="")

        ####################### Color & Sizes #######################

        row = layout.row()
        row.prop(self, "colSizeExp", icon = "TRIA_DOWN" \
            if self.colSizeExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Color & Sizes:")
        if self.colSizeExp:
            box = layout.box()
            row = box.row()
            row.prop(self, "elemDimsExp", icon = "TRIA_DOWN" \
                if self.elemDimsExp else "TRIA_RIGHT",  \
                    icon_only = True, emboss = False)
            row.label(text = "Draw / Edit Element Sizes (Common):")

            if self.elemDimsExp:
                col = box.column().split()
                col.label(text='Segment / Line Thickness:')
                col.prop(self, "lineWidth", text = '')
                col = box.column().split()
                col.label(text='Handle Point Size:')
                col.prop(self, "drawPtSize", text = '')
                col = box.column().split()
                col.label(text='Uniform Subdiv Point Size:')
                col.prop(self, "editSubdivPtSize", text = '')
                col = box.column().split()
                col.label(text='Flexi Grease Res Point Size:')
                col.prop(self, "greaseSubdivPtSize", text = '')
                col = box.column().split()
                col.label(text='Marker Size:')
                col.prop(self, "markerSize", text = '')
                col = box.column().split()
                col.label(text='Axis Line Thickness:')
                col.prop(self, "axisLineWidth", text = '')

            row = box.row()
            row.prop(self, "drawColExp", icon = "TRIA_DOWN" \
                if self.drawColExp else "TRIA_RIGHT",  icon_only = True, \
                    emboss = False)
            row.label(text = "Flexi Draw / Edit Colors:")

            if self.drawColExp:
                col = box.column().split()
                col.label(text="Selected Draw / Edit  Segment:")
                col.prop(self, "colDrawSelSeg", text = '')
                col = box.column().split()
                col.label(text="Adjacent Draw / Edit Segment:")
                col.prop(self, "colDrawNonHltSeg", text = '')
                col = box.column().split()
                col.label(text="Highlighted Edit Segment:")
                col.prop(self, "colDrawHltSeg", text = '')
                col = box.column().split()
                col.label(text="Draw Marker:")
                col.prop(self, "colDrawMarker", text = '')
                col = box.column().split()
                col.label(text="Bezier Point:")
                col.prop(self, "colBezPt", text = '')
                col = box.column().split()
                col.label(text="Subdivision Marker:")
                col.prop(self, "colEditSubdiv", text = '')

            row = box.row()
            row.prop(self, "greaseColExp", icon = "TRIA_DOWN" \
                if self.greaseColExp else "TRIA_RIGHT",  icon_only = True, \
                    emboss = False)
            row.label(text = "Flexi Grease Colors:")

            if self.greaseColExp:
                col = box.column().split()
                col.label(text="Selected Grease Segment:")
                col.prop(self, "colGreaseSelSeg", text = '')
                col = box.column().split()
                col.label(text="Adjacent Grease Segment:")
                col.prop(self, "colGreaseNonHltSeg", text = '')
                col = box.column().split()
                col.label(text="Draw Marker:")
                col.prop(self, "colGreaseMarker", text = '')
                col = box.column().split()
                col.label(text="Bezier Point:")
                col.prop(self, "colGreaseBezPt", text = '')
                col = box.column().split()
                col.label(text="Curve Resolution Marker:")
                col.prop(self, "colGreaseSubdiv", text = '')

            row = box.row()
            row.prop(self, "handleColExp", icon = "TRIA_DOWN" \
                if self.handleColExp else "TRIA_RIGHT",  icon_only = True, \
                    emboss = False)
            row.label(text = "Handle Colors (Common):")

            if self.handleColExp:
                col = box.column().split()
                col.label(text="Free Handle:")
                col.prop(self, "colHdlFree", text = '')
                col = box.column().split()
                col.label(text="Vector Handle:")
                col.prop(self, "colHdlVector", text = '')
                col = box.column().split()
                col.label(text="Aligned Handle:")
                col.prop(self, "colHdlAligned", text = '')
                col = box.column().split()
                col.label(text="Auto Handle:")
                col.prop(self, "colHdlAuto", text = '')
                col = box.column().split()
                col.label(text="Selected Point:")
                col.prop(self, "colSelTip", text = '')
                col = box.column().split()
                col.label(text="Highlighted Point:")
                col.prop(self, "colHltTip", text = '')
                col = box.column().split()
                col.label(text="Handle Point:")
                col.prop(self, "colHdlPtTip", text = '')
                col = box.column().split()
                col.label(text="Adjacent Bezier Point:")
                col.prop(self, "colAdjBezTip", text = '')

        ####################### Snapping Options #######################

        row = layout.row()
        row.prop(self, "snapOptExp", icon = "TRIA_DOWN" \
            if self.snapOptExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Snapping:")

        if self.snapOptExp:
            box = layout.box()
            col = box.column().split()
            col.label(text='Snap Distance:')
            col.prop(self, "snapDist", text = '')
            col = box.column().split()
            col.label(text='Snap Indicator:')
            col.prop(self, "dispSnapInd", text = '')
            col = box.column().split()
            col.label(text='Snap Point Size:')
            col.prop(self, "snapPtSize", text = '')

        ####################### Beveling Options #######################

        row = layout.row()
        row.prop(self, "bevelingExp", icon = "TRIA_DOWN" \
            if self.bevelingExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Beveling:")

        if self.bevelingExp:
            box = layout.box()
            col = box.column().split()
            col.label(text='Default Bevel Factor:')
            col.prop(self, "defBevelFact", text = '')
            col = box.column().split()
            col.label(text='Maximum Bevel Factor:')
            col.prop(self, "maxBevelFact", text = '')
            col = box.column().split()
            col.label(text='Minimum Bevel Factor:')
            col.prop(self, "minBevelFact", text = '')
            col = box.column().split()
            col.label(text='Bevel Factor Step:')
            col.prop(self, "bevelIncr", text = '')

        ####################### Other Options #######################

        row = layout.row()
        row.prop(self, "othPrefExp", icon = "TRIA_DOWN" \
            if self.othPrefExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Other Options:")

        if self.othPrefExp:
            box = layout.box()
            col = box.column().split()
            col.label(text='Flexi Edit Live Update (Experimental):')
            col.prop(self, "liveUpdate", text = '')
            col = box.column().split()
            col.label(text='Display Curve Resolution:')
            col.prop(self, "dispCurveRes", text = '')
            col = box.column().split()
            col.label(text='Display Orientation / Origin Axes:')
            col.prop(self, "dispAxes", text = '')
            col = box.column().split()
            col.label(text='Allow Numpad Entry:')
            col.prop(self, "numpadEntry", text = '')

            col = box.column().split()
            col.label(text='Math Function Text Size:')
            col.prop(self, "mathFnTxtFontSize", text = '')
            col = box.column().split()
            col.label(text='Math Function Text Color:')
            col.prop(self, "colMathFnTxt", text = '')

            box = box.box()
            col = box.column().split()
            col.label(text='Display Keymap:')
            col.prop(self, "showKeyMap", text = '')
            if(self.showKeyMap):
                col = box.column().split()
                col.label(text='Keymap Description Text Color:')
                col.prop(self, "colKeymapText", text = '')
                col = box.column().split()
                col.label(text='Keymap Key Text Color:')
                col.prop(self, "colKeymapKey", text = '')
                col = box.column().split()
                col.label(text='Font Size:')
                col.prop(self, "keyMapFontSize", text = '')
                col = box.column().split()
                col.label(text='Display Next To Toolbar:')
                col.prop(self, "keyMapNextToTool", text = '')
                if(not self.keyMapNextToTool):
                    col = box.column().split()
                    col.label(text='Location X:')
                    col.prop(self, "keyMapLocX", text = '')
                    col = box.column().split()
                    col.label(text='Location Y:')
                    col.prop(self, "keyMapLocY", text = '')


        ####################### Keymap #######################

        row = layout.row()
        row.prop(self, "keymapExp", icon = "TRIA_DOWN" \
            if self.keymapExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Keymap:")

        if self.keymapExp:
            box = layout.box()
            ####################### Common / Draw / Edit Hotkeys #######################
            labels = ["Flexi Tools Common:", "Flexi Draw / Grease:", \
                "Flexi Edit:"]

            for i, keySet in enumerate([FTHotKeys.commonHotkeys, FTHotKeys.drawHotkeys, \
                FTHotKeys.editHotkeys]):
                row = box.row()
                expStr = 'hotKeySet'+ str(i)+'Exp'
                expanded = getattr(self, expStr)
                row.prop(self, expStr, icon = "TRIA_DOWN" \
                    if expanded else "TRIA_RIGHT",  icon_only = True, emboss = False)
                row.label(text = labels[i])

                if expanded:
                    colM = box.column()
                    boxIn = colM.grid_flow(row_major = True, \
                        even_columns = True, columns = 3)
                    j = 0
                    for j, keydata in enumerate(keySet):
                        expStr = keydata.id + "Exp"
                        col = boxIn.box().column(align=True)
                        col.prop(self, expStr, icon = "TRIA_DOWN" \
                            if eval('self.' + expStr) else "TRIA_RIGHT", \
                                emboss = False, text = keydata.label + ':')
                        if  getattr(self, expStr):
                            col.prop(self, keydata.id, text = '', event = True)
                            if(not keydata.isExclusive):
                                rowC = col.row()
                                rowC.prop(self, keydata.id + 'Alt', \
                                    text = 'Alt', toggle = True)
                                rowC.prop(self, keydata.id + 'Ctrl', \
                                    text = 'Ctrl', toggle = True)
                                rowC.prop(self, keydata.id + 'Shift', \
                                    text = 'Shift', toggle = True)
                    for idx in range(j % 3, 2):
                        col = boxIn.column(align = True)
            ####################### Snap Hotkeys #######################

            row = box.row()
            row.prop(self, "hkSnapExp", icon = "TRIA_DOWN" \
                if self.hkSnapExp else "TRIA_RIGHT",  icon_only = True, emboss = False)
            row.label(text = "Snapping")

            if self.hkSnapExp:
                colM = box.column()
                box = colM.grid_flow(row_major = True, even_columns = True, columns = 3)
                for i, keydata in enumerate(FTHotKeys.snapHotkeys):
                    keydataMeta = FTHotKeys.snapHotkeysMeta[i]
                    expStr = keydata.id + "Exp"
                    col = box.box().column(align=True)
                    col.prop(self, expStr, icon = "TRIA_DOWN" if eval('self.' + expStr) \
                        else "TRIA_RIGHT", emboss = False, text = keydataMeta.label + ':')
                    if  getattr(self, expStr):
                        if(getattr(self, keydataMeta.id) != 'KEY'):
                            col.prop(self, keydataMeta.id, text = '')
                        else:
                            col = col.box().column()
                            col.prop(self, keydataMeta.id, text = '')
                            col.prop(self, keydata.id, text = '', event = True)

        col = layout.column().split()
        col.operator('object.reset_default_props')
        col.operator('object.reset_default_hotkeys')

classes = [
    ModalMarkSegStartOp,
    SeparateSplinesObjsOp,
    SplitBezierObjsOp,
    splitBezierObjsPtsOp,
    JoinBezierSegsOp,
    CloseSplinesOp,
    CloseStraightOp,
    OpenSplinesOp,
    ExportSVGOp,
    RemoveDupliVertCurveOp,
    IntersectCurvesOp,
    convertToMeshOp,
    SetHandleTypesOp,
    SetCurveColorOp,
    PasteLengthOp,
    AlignToFaceOp,
    RemoveCurveColorOp,
    SelectInCollOp,
    InvertSelOp,
    BezierUtilsPanel,
    ModalFlexiDrawBezierOp,
    ModalFlexiDrawGreaseOp,
    ModalFlexiEditBezierOp,
    BezierUtilsPreferences,
    BezierToolkitParams,
    ResetDefaultPropsOp,
    ResetDefaultHotkeys,
    FTMenuOptionOp,
    SaveMathFn,
    LoadMathFn,
    ResetMathFn,
    DeleteMathFn,
]

for menuData in FTMenu.editMenus:
    exec(FTMenu.getMNClassDefStr(menuData))
    classes.append(eval(menuData.menuClassName))


def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.WindowManager.bezierToolkitParams = \
        bpy.props.PointerProperty(type = BezierToolkitParams)

    BezierUtilsPanel.colorCurves(add = True)
    bpy.app.handlers.depsgraph_update_post.append(BezierUtilsPanel.colorCurves)

    bpy.utils.register_tool(FlexiDrawBezierTool)
    bpy.utils.register_tool(FlexiEditBezierTool)
    bpy.utils.register_tool(FlexiGreaseBezierTool)
    
    updatePanel(None, bpy.context)

    bpy.app.handlers.load_post.append(ModalBaseFlexiOp.loadPostHandler)
    bpy.app.handlers.load_pre.append(ModalBaseFlexiOp.loadPreHandler)

def unregister():
    BezierUtilsPanel.colorCurves(remove = True)
    bpy.app.handlers.depsgraph_update_post.remove(BezierUtilsPanel.colorCurves)

    try: ModalBaseFlexiOp.opObj.cancelOp(bpy.context)
    except: pass # If not invoked or already unregistered

    bpy.app.handlers.load_post.remove(ModalBaseFlexiOp.loadPostHandler)
    bpy.app.handlers.load_pre.remove(ModalBaseFlexiOp.loadPreHandler)

    del bpy.types.WindowManager.bezierToolkitParams

    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

    bpy.utils.unregister_tool(FlexiDrawBezierTool)
    bpy.utils.unregister_tool(FlexiEditBezierTool)
    bpy.utils.unregister_tool(FlexiGreaseBezierTool)
