#
#
# Blender add-on with tools to draw and edit Bezier curves along with other utility ops
#
# Supported Blender Version: 2.8 Beta
#
# Copyright (C) 2019  Shrinivas Kulkarni

# License: GPL (https://github.com/Shriinivas/blenderbezierutils/blob/master/LICENSE)
#

import bpy, bmesh, bgl, gpu
from bpy.props import BoolProperty, IntProperty, EnumProperty, StringProperty
from bpy.types import Panel, Operator, WorkSpaceTool, AddonPreferences
from mathutils import Vector, Matrix, geometry, kdtree
from math import log, atan, tan, pi, radians
from bpy_extras.view3d_utils import region_2d_to_vector_3d, region_2d_to_location_3d
from bpy_extras.view3d_utils import location_3d_to_region_2d
from gpu_extras.batch import batch_for_shader
import time
from bpy.app.handlers import persistent

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "version": (0, 8, 25),
    "location": "Properties > Active Tool and Workspace Settings > Bezier Utilities",
    "description": "Collection of Bezier curve utility ops",
    "category": "Object",
    "wiki_url": "https://github.com/Shriinivas/blenderbezierutils/blob/master/README.md",
    "blender": (2, 80, 0),
}

DEF_ERR_MARGIN = 0.0001

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
def copyBezierPt(src, target, freeHandles, srcMw = Matrix(), invDestMW = Matrix()):
    target.handle_left_type = 'FREE'
    target.handle_right_type = 'FREE'

    target.co = invDestMW @ (srcMw @ src.co)
    target.handle_left = invDestMW @ (srcMw @ src.handle_left)
    target.handle_right = invDestMW @ (srcMw @ src.handle_right)

    if(not freeHandles):
        target.handle_left_type = src.handle_left_type
        target.handle_right_type =  src.handle_right_type

def createSplineForSeg(curveData, bezierPts):
    spline = curveData.splines.new('BEZIER')
    spline.bezier_points.add(len(bezierPts)-1)
    spline.use_cyclic_u = False

    for i in range(0, len(bezierPts)):
        copyBezierPt(bezierPts[i], spline.bezier_points[i], freeHandles = True)

def createSpline(curveData, srcSpline, forceNoncyclic, freeHandles, excludePtIdxs = {}):
    spline = curveData.splines.new('BEZIER')
    spline.bezier_points.add(len(srcSpline.bezier_points) - len(excludePtIdxs) - 1)

    if(forceNoncyclic):
        spline.use_cyclic_u = False
    else:
        spline.use_cyclic_u = srcSpline.use_cyclic_u

    ptIdx = 0
    for i in range(0, len(srcSpline.bezier_points)):
        if(i not in excludePtIdxs):
            copyBezierPt(srcSpline.bezier_points[i], spline.bezier_points[ptIdx], freeHandles)
            ptIdx += 1

    if(forceNoncyclic == True and srcSpline.use_cyclic_u == True):
        spline.bezier_points.add(1)
        copyBezierPt(srcSpline.bezier_points[0], spline.bezier_points[-1], freeHandles)

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
            keyData.append([[d.handle_left, d.co, d.handle_right] for d in key.data])
            keyNames.append(key.name)

    return keyNames, keyData

def updateShapeKeyData(obj, keyData, keyNames, startIdx, cnt):
    if(obj.data.shape_keys == None):
        return

    removeShapeKeys(obj)

    for i, name in enumerate(keyNames):
        key = obj.shape_key_add(name = name)
        for j in range(0, cnt):
            keyIdx = j + startIdx
            key.data[j].handle_left = keyData[i][keyIdx][0]
            key.data[j].co = keyData[i][keyIdx][1]
            key.data[j].handle_right = keyData[i][keyIdx][2]

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

def addLastSeg(spline):
    if(spline.use_cyclic_u):
        spline.use_cyclic_u = False
        spline.bezier_points.add(1)
        copyObjAttr(spline.bezier_points[0], spline.bezier_points[-1])

def joinCurves(curves):
    obj = curves[0]
    invMW = obj.matrix_world.inverted()
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

def removeBezierPts(obj, splineIdx, removePtIdxs):
    oldSpline = obj.data.splines[splineIdx]
    bpts = oldSpline.bezier_points
    if(min(removePtIdxs) >= len(bpts)):
        return

    if(len(bpts) == 1 and 0 in removePtIdxs) :
        obj.data.splines.remove(oldSpline)
        if(len(obj.data.splines) == 0):
            safeRemoveObj(obj)
        return

    createSpline(obj.data, oldSpline, False, False, removePtIdxs)
    obj.data.splines.remove(oldSpline)
    splineCnt = len(obj.data.splines)
    nextIdx = splineIdx
    for idx in range(nextIdx, splineCnt - 1):
        oldSpline = obj.data.splines[nextIdx]
        createSpline(obj.data, oldSpline, False, False)
        obj.data.splines.remove(oldSpline)

def insertBezierPts(obj, splineIdx, startIdx, ts, handleType):

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
    addCnt = len(ts)
    
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

        t = ts[i]
        seg = [prevPt.co, prevPt.handle_right, nextPt.handle_left, nextPt.co]
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

#Change position of bezier points according to new matrix_world
def changeMW(obj, newMW):
    invMW = newMW.inverted()
    for spline in obj.data.splines:
        for pt in spline.bezier_points:
            pt.co = invMW @ (obj.mw @ pt.co)
            pt.handle_left = invMW @ (obj.mw @ pt.handle_left)
            pt.handle_right = invMW @ (obj.mw @ pt.handle_right)
    obj.matrix_world = newMW

#Round to logarithmic scale .1, 0, 10, 100 etc.
#(47.538, -1) -> 47.5; (47.538, 0) -> 48.0; (47.538, 1) -> 50.0; (47.538, 2) -> 0,
def roundedVect(vect, rounding):
    rounding += 1
    fact = (10 ** rounding) / 10
    return Vector([round(l / fact) * fact for l in vect])

###################### Screen functions ######################

def  getViewDistRounding(context):
    viewDist = context.space_data.region_3d.view_distance
    return int(log(viewDist, 10)) - 1

def getCoordFromLoc(context, loc):
    region = context.region
    rv3d = context.space_data.region_3d
    coord = location_3d_to_region_2d(region, rv3d, loc)
    # return a unlocatable pt if None to avoid errors
    return coord if(coord != None) else Vector((9e+99, 9e+99))

# To be called only from 3d view
def getCurrAreaRegion(context):
    a, r = [(a, r) for a in bpy.context.screen.areas if a.type == 'VIEW_3D' for r in a.regions \
        if(r == context.region)][0]
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

###################### Op Specific functions ######################

def closeSplines(curve, htype = None):
    for spline in curve.data.splines:
        if(htype != None):
            spline.bezier_points[0].handle_left_type = htype
            spline.bezier_points[-1].handle_right_type = htype
        spline.use_cyclic_u = True

#split value is one of {'spline', 'seg', 'point'} (TODO: Enum)
def splitCurve(selObjs, split, newColl = True):
    changeCnt = 0
    splineCnt = 0
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
                            newSpline.bezier_points[0], newWM.inverted(), mw)

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
                createSpline(objCopy.data, spline, forceNoncyclic = False, \
                    freeHandles = False)
                currSegCnt = len(objCopy.data.splines[0].bezier_points)
                updateShapeKeyData(objCopy, keyData, keyNames, segCnt, currSegCnt)
                newObjs.append(objCopy)
                segCnt += currSegCnt

        for collection in collections:
            if(newColl):
                collection.children.link(objGrp)
            collection.objects.unlink(obj)

        bpy.data.curves.remove(obj.data)
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
            closestCurve2, dist2 = getClosestCurve(srcMW, ncStart.co, remainingCurves, dist)
            if(closestCurve2 != None):
                reverseCurve(nextCurve)
                closestCurve = closestCurve2

        orderedCurves.append(closestCurve.name)
        nextCurve = closestCurve
        remainingCurves.remove(closestCurve)
    return [idMap[cn] for cn in orderedCurves]

def joinSegs(curves, optimized, straight, srcCurve = None):
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
    invSrcMW = srcMW.inverted()
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
        if(vectCmpWithMargin(srcMW @ currBezierPt.co, mw @ nextBezierPt.co)):
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
                vectCmpWithMargin(mw @ nextSpline.bezier_points[i].co, srcMW @ currSpline.bezier_points[0].co)):
                    currSpline.bezier_points[0].handle_left_type = 'FREE'
                    currSpline.bezier_points[0].handle_left = invSrcMW @ (mw @ nextSpline.bezier_points[i].handle_left)
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

def removeDupliVert(curve):
    newCurveData = curve.data.copy()
    newCurveData.splines.clear()
    dupliFound = False
    for spline in curve.data.splines:
        newCurveData.splines.new('BEZIER')
        currSpline = newCurveData.splines[-1]
        copyObjAttr(spline, currSpline)
        prevPt = spline.bezier_points[0]
        mextPt = None
        copyObjAttr(spline.bezier_points[0], currSpline.bezier_points[0])

        if(len(spline.bezier_points) == 1):
            continue

        cmpPts = spline.bezier_points[:]
        pt = spline.bezier_points[0]
        while(vectCmpWithMargin(cmpPts[-1].co, pt.co) and
            len(cmpPts) > 1):
            prevPt = cmpPts.pop()
            pt.handle_left_type = 'FREE'
            pt.handle_left = prevPt.handle_left
            currSpline.use_cyclic_u = True

        for pt in cmpPts:
            if(vectCmpWithMargin(prevPt.co, pt.co)):
                prevPt.handle_right_type = 'FREE'
                prevPt.handle_right = pt.handle_right
                dupliFound = True
                continue
            currSpline.bezier_points.add(1)
            copyObjAttr(pt, currSpline.bezier_points[-1])
            prevPt = pt

    if(dupliFound):
        curve.data = newCurveData
    else:
        bpy.data.curves.remove(newCurveData)

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
        return mesh

    planeVert = Vector([round(c, 5) for c in meshObj.data.vertices[0].co])
    mod = meshObj.modifiers.new('mod', type='SOLIDIFY')
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
                ("s" if(changeCnt > 1) else "") + " into " +str(len(newObjs)) + " new ones")

        return {'FINISHED'}


class SplitBezierObjsOp(Operator):

    bl_idname = "object.separate_segments"
    bl_label = "Separate Segments"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate segments of selected Bezier curves as new objects"

    def execute(self, context):

        selObjs = bpy.context.selected_objects
        newObjs, changeCnt = splitCurve(selObjs, split = 'seg')

        if(changeCnt > 0):
            bpy.context.view_layer.objects.active = newObjs[-1]
            self.report({'INFO'}, "Split "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + " into " + str(len(newObjs)) + " new objects")

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
                ("s" if(changeCnt > 1) else "") + " into " + str(len(newObjs)) + " new objects")

        return {'FINISHED'}


class JoinBezierSegsOp(Operator):
    bl_idname = "object.join_curves"
    bl_label = "Join"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Join selected curves (with new segments if required)"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]

        straight = bpy.context.scene.straight
        optimized = bpy.context.scene.optimized

        newCurve = joinSegs(curves, optimized = optimized, straight = straight)
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
            collections = bpy.context.active_object.users_collection

            for obj in bpy.context.selected_objects:
                obj.select_set(False)

            selectIntrvl = bpy.context.scene.selectIntrvl

            for collection in collections:
                objs = [o for o in collection.objects]
                idx = objs.index(bpy.context.active_object)
                objs = objs[idx:] + objs[:idx]
                for i, o in enumerate(objs):
                    if( i % (selectIntrvl + 1) == 0):
                        o.select_set(True)
                    else:
                        o.select_set(False)

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
                spline.bezier_points[0].handle_left_type = 'ALIGNED'
                spline.bezier_points[-1].handle_right_type = 'ALIGNED'
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
        ht = bpy.context.scene.handleType
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
            removeDupliVert(curve)

        return {'FINISHED'}


class convertTo2DMeshOp(Operator):
    bl_idname = "object.convert_2d_mesh"
    bl_label = "Convert"
    bl_description = "Convert 2D curve to mesh with quad faces"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        curve = context.object
        if(curve != None and isBezier(curve)):
            for spline in curve.data.splines:
                spline.use_cyclic_u = True
            curve.data.dimensions = '2D'
            curve.data.fill_mode = 'BOTH'
            meshObj = convertToMesh(curve)

            remeshDepth = bpy.context.scene.remeshDepth
            unsubdivide = bpy.context.scene.unsubdivide
            applyMeshModifiers(meshObj, remeshDepth)

            if(unsubdivide):
                unsubdivideObj(meshObj)

            meshObj.matrix_world = curve.matrix_world.copy()

            safeRemoveObj(curve)

            bpy.context.view_layer.objects.active = meshObj

        return {'FINISHED'}

def markVertHandler(self, context):
    if(self.markVertex):
        bpy.ops.wm.mark_vertex()


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
        self.shader.bind()

        MarkerController.drawHandlerRef = \
            bpy.types.SpaceView3D.draw_handler_add(self.drawHandler, \
                (), "WINDOW", "POST_VIEW")

        self.createBatch(context)

    def saveStartVerts(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            splines = curve.data.splines
            spMap = self.smMap[curveName]

            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                if(markerInfo[1] != 0):
                    pts = splines[splineIdx].bezier_points
                    loc, idx = markerInfo[0], markerInfo[1]
                    cnt = len(pts)

                    ptCopy = [[p.co.copy(), p.handle_right.copy(), \
                        p.handle_left.copy(), p.handle_right_type, \
                            p.handle_left_type] for p in pts]

                    for i, pt in enumerate(pts):
                        srcIdx = (idx + i) % cnt
                        p = ptCopy[srcIdx]

                        #Must set the types first
                        pt.handle_right_type = p[3]
                        pt.handle_left_type = p[4]
                        pt.co = p[0]
                        pt.handle_right = p[1]
                        pt.handle_left = p[2]

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
            states.append(s.overlay.show_curve_handles)
            s.overlay.show_curve_handles = False
        return states

    def resetShowHandleState(context, handleStates):
        spaces = MarkerController.getSpaces3D(context)
        for i, s in enumerate(spaces):
            s.overlay.show_curve_handles = handleStates[i]


class ModalMarkSegStartOp(bpy.types.Operator):
    bl_description = "Mark Vertex"
    bl_idname = "wm.mark_vertex"
    bl_label = "Mark Start Vertex"

    def cleanup(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        self.markerState.removeMarkers(context)
        MarkerController.resetShowHandleState(context, self.handleStates)
        bpy.context.scene.markVertex = False

    def modal (self, context, event):

        if(context.mode  == 'OBJECT' or event.type == "ESC" or\
            not bpy.context.scene.markVertex):
            self.cleanup(context)
            return {'CANCELLED'}

        elif(event.type == "RET"):
            self.markerState.saveStartVerts()
            self.cleanup(context)
            return {'FINISHED'}

        if(event.type == 'TIMER'):
            self.markerState.updateSMMap()
            self.markerState.createBatch(context)

        elif(event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}):
            self.ctrl = (event.value == 'PRESS')

        elif(event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}):
            self.shift = (event.value == 'PRESS')

            if(event.type not in {"MIDDLEMOUSE", "TAB", "LEFTMOUSE", \
                "RIGHTMOUSE", 'WHEELDOWNMOUSE', 'WHEELUPMOUSE'} and \
                not event.type.startswith("NUMPAD_")):
                return {'RUNNING_MODAL'}

        return {"PASS_THROUGH"}

    def execute(self, context):
        #TODO: Why such small step?
        self._timer = context.window_manager.event_timer_add(time_step = 0.0001, \
            window = context.window)
        self.ctrl = False
        self.shift = False

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

    bpy.types.Scene.markVertex = BoolProperty(name="Mark Starting Vertices", \
        description='Mark first vertices in all closed splines of selected curves', \
            default = False, update = markVertHandler)

    bpy.types.Scene.selectIntrvl = IntProperty(name="Selection Interval", \
        description='Interval between selected objects', \
            default = 0, min = 0)

    bpy.types.Scene.handleType = EnumProperty(name="Handle Type", items = \
        [("AUTO", 'Automatic', "Automatic"), \
         ('VECTOR', 'Vector', 'Straight line'), \
         ('ALIGNED', 'Aligned', 'Left and right aligned'), \
         ('FREE', 'Free', 'Left and right independent')], \
        description = 'Handle type of the control points',
        default = 'ALIGNED')

    bpy.types.Scene.remeshDepth = IntProperty(name="Remesh Depth", \
        description='Remesh depth for converting to mesh', \
            default = 4, min = 1, max = 10)

    bpy.types.Scene.unsubdivide = BoolProperty(name="Unsubdivide", \
        description='Unsubdivide to reduce the number of polygons', \
            default = False)

    bpy.types.Scene.straight = BoolProperty(name="Join With Straight Segments", \
        description='Join curves with straight segments', \
            default = False)

    bpy.types.Scene.optimized = BoolProperty(name="Join Optimized", \
        description='Join the nearest curve (reverse direction if necessary)', \
            default = True)

    bpy.types.Scene.splitExpanded = BoolProperty(name="Split Expanded State",
            default = False)

    bpy.types.Scene.joinExpanded = BoolProperty(name="Join Expanded State",
            default = False)

    bpy.types.Scene.selectExpanded = BoolProperty(name="Select Expanded State",
            default = False)

    bpy.types.Scene.convertExpanded = BoolProperty(name="Convert Expanded State",
            default = False)

    bpy.types.Scene.handleTypesExpanded = BoolProperty(name="Set Handle Types State",
            default = False)

    bpy.types.Scene.otherExpanded = BoolProperty(name="Other Expanded State",
            default = False)

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_CURVE'}

    def draw(self, context):
        layout = self.layout
        # ~ layout.use_property_split = True
        layout.use_property_decorate = False

        if(context.mode  == 'OBJECT'):
            row = layout.row()
            row.prop(context.scene, "splitExpanded",
                icon="TRIA_DOWN" if context.scene.splitExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False)
            row.label(text="Split Bezier Curves", icon = 'UNLINKED')

            if context.scene.splitExpanded:
                col = layout.column()
                col.operator('object.separate_splines')#icon = 'UNLINKED')
                col = layout.column()
                col.operator('object.separate_segments')#, icon = 'ORPHAN_DATA')
                col = layout.column()
                col.operator('object.separate_points')#, icon = 'ORPHAN_DATA')

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "joinExpanded",
                icon="TRIA_DOWN" if context.scene.joinExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text="Join Bezier Curves", icon = 'LINKED')

            if context.scene.joinExpanded:
                col = layout.column()
                col.prop(context.scene, 'straight')
                col = layout.column()
                col.prop(context.scene, 'optimized')
                col = layout.column()
                col.operator('object.join_curves')

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "selectExpanded",
                icon="TRIA_DOWN" if context.scene.selectExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Select Objects In Collection', icon='RESTRICT_SELECT_OFF')
            if context.scene.selectExpanded:
                col = layout.column()
                row = col.row()
                row.prop(context.scene, 'selectIntrvl')
                row.operator('object.select_in_collection')
                col = layout.column()
                col.operator('object.invert_sel_in_collection')

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "convertExpanded",
                icon="TRIA_DOWN" if context.scene.convertExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Convert Curve to Mesh', icon='MESH_DATA')

            if context.scene.convertExpanded:
                col = layout.column()
                row = col.row()
                row.prop(context.scene, 'remeshDepth')
                row.prop(context.scene, 'unsubdivide')
                col = layout.column()
                col.operator('object.convert_2d_mesh')

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "handleTypesExpanded",
                icon="TRIA_DOWN" if context.scene.convertExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Set Handle Type', icon='MOD_CURVE')

            if context.scene.handleTypesExpanded:
                col = layout.column()
                row = col.row()
                col.prop(context.scene, 'handleType')
                col = layout.column()
                col.operator('object.set_handle_types')

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "otherExpanded",
                icon="TRIA_DOWN" if context.scene.otherExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Other Tools', icon='TOOL_SETTINGS')

            if context.scene.otherExpanded:
                col = layout.column()
                col.operator('object.close_splines')
                col = layout.column()
                col.operator('object.close_straight')
                col = layout.column()
                col.operator('object.open_splines')
                col = layout.column()
                col.operator('object.remove_dupli_vert_curve')

        else:
            col = layout.column()
            col.prop(context.scene, 'markVertex', toggle = True)


################### Common Bezier Functions & Classes ###################

def getPtFromT(p0, p1, p2, p3, t):
    c = (1 - t)
    pt = (c ** 3) * p0 + 3 * (c ** 2) * t * p1 + \
        3 * c * (t ** 2) * p2 + (t ** 3) * p3
    return pt

# iterative brute force, not optimized, some iterations maybe redundant
def getTsForPt(p0, p1, p2, p3, co, coIdx, tolerance = 0.0001, maxItr = 1000):
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
def getTForPt(curve, testPt):
    minLen = 9e+99
    retT = None
    for coIdx in range(0, 3):
        ts = getTsForPt(curve[0], curve[1], curve[2], curve[3], testPt[coIdx], coIdx)
        for t in ts:
            pt = getPtFromT(curve[0], curve[1], curve[2], curve[3], t)
            pLen = (testPt - pt).length
            if(pLen < minLen):
                minLen = pLen
                retT = t
            # ~ if(all(abs(testPt[i] - pt[i]) < .1 or abs(testPt[i] - pt[i]) < pt[i] * .1 \
                # ~ for i in range(0, len(pt)))):
                # ~ return t
    return retT

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

# Get pt coords along curve defined by the four control pts (curvePts)
# divPerUnitLength: No of subdivisions per unit length
# (which is the same as no of pts excluding the end pts)
def getPtsAlongBezier(curvePts, curveRes, context = None):
    if(len(curvePts) < 2):
        return []

    pts = []
    minRes = 20

    if(context != None and hasattr(context, 'space_data') and \
        context.space_data != None and hasattr(context.space_data, 'region_3d')):
        viewDist = context.space_data.region_3d.view_distance

        # (the smaller the view dist (higher zoom level),
        # the higher the num of subdivisions
        curveRes = curveRes / viewDist

    for i in range(0, len(curvePts) - 1):
        pt0 = curvePts[i]
        pt1 = curvePts[i+1]
        curve = [pt0[1], pt0[2], pt1[0], pt1[1]]

        #Add 2 for start and end
        res = int((curve[2] - curve[1]).length * curveRes)
        if(res < minRes):
            res = minRes

        locs = geometry.interpolate_bezier(pt0[1], pt0[2], pt1[0], pt1[1], res)
        pts += locs

    return pts


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

    #Let's make at least the line segments of predictable length :)
    if(pts[0] == pts[1] and pts[2] == pts[3]):
        pt0 = Vector([(1 - t0) * pts[0][i] + t0 * pts[2][i] for i in range(0, 3)])
        pt1 = Vector([(1 - t1) * pts[0][i] + t1 * pts[2][i] for i in range(0, 3)])
        return [pt0, pt0, pt1, pt1]

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

################### Common to Draw and Edit Flexi Bezier Ops ###################

# Some static constants
SEL_SEG_COLOR = (.6, .8, 1, 1)
ADJ_SEG_COLOR = (.1, .4, .6, 1)
NONADJ_SEG_COLOR = (.2, .5, .8, 1)
SEL_TIP_COLOR = (.2, .7, .3, 1)
HLT_TIP_COLOR = (.2, 1, .9, 1)
ENDPT_TIP_COLOR = (1, 1, 0, 1)
ADJ_ENDPT_TIP_COLOR = (.1, .1, .1, 1)
TIP_COLOR = (.7, .7, 0, 1)
MARKER_COLOR = (.7, 1, .7, 1)

TIP_COL_PRIORITY = {ADJ_ENDPT_TIP_COLOR: 0, TIP_COLOR: 1,
    ENDPT_TIP_COLOR: 2, SEL_TIP_COLOR : 3, HLT_TIP_COLOR : 4, MARKER_COLOR: 5}

HANDLE_COLOR_MAP ={'FREE': (.6, .05, .05, 1), 'VECTOR': (.4, .5, .2, 1), \
    'ALIGNED': (1, .3, .3, 1), 'AUTO': (.8, .5, .2, 1)}

DEF_CURVE_RES = 100 # Per unit seg len for viewdist = 10
SNAP_DIST_PIXEL = 20
STRT_SEG_HDL_LEN_COEFF = 0.25
DBL_CLK_DURN = 0.25
SNGL_CLK_DURN = 0.3

SEL_CURVE_SEARCH_RES = 1000
NONSEL_CURVE_SEARCH_RES = 100
ADD_PT_CURVE_SEARCH_RES = 5000

class SegDisplayInfo:

    # handleNos: 0: seg1-left, 1: seg1-right, 2: seg2-left, 3: seg2-right
    # tipColors: leftHdl, pt, rightHdl for both ends (total 6) (None: don't show tip)
    # Caller to make sure there are no tips without handle
    def __init__(self, segPts, segColor, handleNos, tipColors):
        self.segPts = segPts
        self.segColor = segColor
        self.handleNos = handleNos
        self.tipColors = tipColors

# Return line batch for bezier line segments and handles and point batch for handle tips
def getBezierBatches(shader, displayInfos, context = None, defHdlType = 'ALIGNED'):
    lineCos = [] #segment is also made up of lines
    lineColors = []
    tipCos = []
    tipColors = []

    for i, displayInfo in enumerate(displayInfos):
        segPts = displayInfo.segPts
        pts = getPtsAlongBezier(segPts, DEF_CURVE_RES, context)
        segLineCos = getLinesFromPts(pts)
        lineCos += segLineCos
        lineColors += [displayInfo.segColor for j in range(0, len(segLineCos))]

    for i, displayInfo in enumerate(displayInfos):
        segPts = displayInfo.segPts
        for handleNo in displayInfo.handleNos:
            # [4] array to [2][2] array
            ptIdx = int(handleNo / 2)
            hdlIdx = handleNo % 2
            lineCos += [segPts[ptIdx][hdlIdx], segPts[ptIdx][hdlIdx + 1]]

            # Making exception for Flexi Draw, it does not need to store handle types
            # (all are aligned by default), so it can have segtPts with only 3 elements
            if(len(segPts[ptIdx]) < 5):
                htype = defHdlType
            else:
                htype = segPts[ptIdx][3 + hdlIdx]

            lineColors += [HANDLE_COLOR_MAP[htype], \
                HANDLE_COLOR_MAP[htype]]

        tipColInfo = []
        for j, tipColor in enumerate(displayInfo.tipColors):
            if(tipColor != None):
                # [6] array to [2][3] array (includes end points of curve)
                ptIdx = int(j / 3)
                hdlIdx = j % 3
                tipColInfo.append([tipColor, segPts[ptIdx][hdlIdx]])
        tipColInfo = sorted(tipColInfo, key = lambda x: TIP_COL_PRIORITY[x[0]])
        for ti in tipColInfo:
            tipColors.append(ti[0])
            tipCos.append(ti[1])

    lineBatch = batch_for_shader(shader, "LINES", {"pos": lineCos, "color": lineColors})
    tipBatch = batch_for_shader(shader, "POINTS", {"pos": tipCos, "color": tipColors})

    return lineBatch, tipBatch


################### Flexi Draw Bezier Curve ###################


class ModalDrawBezierOp(Operator):

    drawHandlerRef = None
    segBatch = None
    tipBatch = None
    shader = None
    defLineWidth = 1.5
    defPointSize = 5
    running = False

    @classmethod
    def poll(cls, context):
        return not ModalDrawBezierOp.running

    def addDrawHandler():
        ModalDrawBezierOp.drawHandlerRef = \
            bpy.types.SpaceView3D.draw_handler_add(ModalDrawBezierOp.drawHandler, \
                (), "WINDOW", "POST_VIEW")

    def removeDrawHandler():
        if(ModalDrawBezierOp.drawHandlerRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(ModalDrawBezierOp.drawHandlerRef, \
                "WINDOW")
            ModalDrawBezierOp.drawHandlerRef = None

    def cancelOp(self, context):
        ModalDrawBezierOp.createBatch(context, [], [])
        ModalDrawBezierOp.removeDrawHandler()
        ModalDrawBezierOp.running = False

    #static method
    def drawHandler():
        if(ModalDrawBezierOp.shader != None):
            bgl.glLineWidth(ModalDrawBezierOp.defLineWidth)
            bgl.glPointSize(ModalDrawBezierOp.defPointSize)

            if(ModalDrawBezierOp.segBatch != None):
                ModalDrawBezierOp.segBatch.draw(ModalDrawBezierOp.shader)

            if(ModalDrawBezierOp.tipBatch != None):
                ModalDrawBezierOp.tipBatch.draw(ModalDrawBezierOp.shader)

    @persistent
    def loadPostHandler(dummy):
        ModalDrawBezierOp.addDrawHandler()
        if(ModalDrawBezierOp.shader != None):
            ModalDrawBezierOp.createBatch(bpy.context, [], [])
        ModalDrawBezierOp.running = False

    @persistent
    def loadPreHandler(dummy):
        ModalDrawBezierOp.removeDrawHandler()

    #static method
    def createBatch(context, curvePts, endSegHandleNos, segColor = None, tipColor = None):
        if(tipColor == None): tipColor = TIP_COLOR
        ptCnt = len(curvePts)
        displayInfos = []
        for i in range(0, ptCnt - 1):
            segPts = [curvePts[i], curvePts[i+1]]
            if(i == ptCnt - 2):
                handleNos = endSegHandleNos
                tm = {0:0, 2:1, 3:2, 5:3}
                tipColors = [tipColor if tm.get(i) in handleNos else None \
                    for i in range (0, 6) ]
                tipColors[5] = MARKER_COLOR

                sc = SEL_SEG_COLOR
            else:
                handleNos = []
                tipColors = []
                sc = NONADJ_SEG_COLOR

            if(segColor != None):#overrides earlier setting
                sc = segColor

            displayInfos.append(SegDisplayInfo(segPts, sc, handleNos, tipColors))

        ModalDrawBezierOp.segBatch, ModalDrawBezierOp.tipBatch = \
            getBezierBatches(ModalDrawBezierOp.shader, displayInfos, context)

        if context.area:
            context.area.tag_redraw()

    def __init__(self, curveDispRes = 20):
        #No of subdivisions in the displayed curve with view dist 1
        self.curveDispRes = curveDispRes
        self.defaultSnapSteps = 3

    def invoke(self, context, event):
        ModalDrawBezierOp.running = True

        self.initialize()

        ModalDrawBezierOp.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
        ModalDrawBezierOp.shader.bind()

        context.window_manager.modal_handler_add(self)
        ModalDrawBezierOp.addDrawHandler()

        return {"RUNNING_MODAL"}

    #This will be called multiple times not just at the beginning
    def initialize(self):
        self.curvePts = []
        self.clickT = None #For double click
        self.pressT = None #For single click
        self.capture = False
        self.ctrl = False
        self.shift = False
        self.alt = False
        self.lockAxes = []
        self.snapSteps = self.defaultSnapSteps

    def confirm(self, context, event):
        self.save(context, event)
        self.curvePts = []
        self.capture = False
        ModalDrawBezierOp.createBatch(context, [], [])
        self.initialize()

    def modal(self, context, event):

        if(context.space_data == None):
            return {'PASS_THROUGH'}

        if(event.type == 'WINDOW_DEACTIVATE' and event.value == 'PRESS'):
            self.ctrl = False
            self.shift = False
            self.alt = False
            return {'PASS_THROUGH'}

        if(isOutside(context, event)):
            if(len(self.curvePts) > 0):
                return {'RUNNING_MODAL'}
            return {'PASS_THROUGH'}

        if((event.type == 'E' or event.type == 'e') and event.value == 'RELEASE'):
            # ~ bpy.ops.wm.tool_set_by_id(name = FlexiEditBezierTool.bl_idname) (T60766)
            bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.edit_tool')
            return {"RUNNING_MODAL"}

        if(event.type == 'RET' or event.type == 'SPACE'):
            self.confirm(context, event)
            self.redrawBezier(context, event)
            return {'RUNNING_MODAL'}

        if(event.type == 'ESC'):
            self.initialize()
            self.redrawBezier(context, event)
            return {'RUNNING_MODAL'}

        if(event.type == 'BACK_SPACE' and event.value == 'RELEASE'):
            self.lockAxes = []
            if(len(self.curvePts) > 0):
                self.curvePts.pop()
            if(len(self.curvePts) <= 1): #Because there is an extra point (the current one)
                self.curvePts = []
                self.capture = False
            self.redrawBezier(context, event)
            return {'RUNNING_MODAL'}

        if(event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}):
            self.ctrl = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}):
            self.shift = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type in {'LEFT_ALT', 'RIGHT_ALT'}):
            self.alt = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type == 'WHEELDOWNMOUSE' and self.shift):
            if(self.snapSteps < 10):
                self.snapSteps += 1
            return {'PASS_THROUGH'}

        if(event.type == 'WHEELUPMOUSE' and self.shift):
            if(self.snapSteps > 1):
                self.snapSteps -= 1
            return {'PASS_THROUGH'}

        if(event.type == 'MIDDLEMOUSE' and self.shift):
            self.snapSteps = self.defaultSnapSteps
            return {'PASS_THROUGH'}

        if (event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            if(len(self.curvePts) == 0):
                loc = self.get3dLocWithSnap(context, event)
                self.curvePts.append([loc, loc, loc])
            self.capture = True
            self.lockAxes = []
            self.pressT = time.time()
            return {'RUNNING_MODAL'}

        if (event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):
            self.capture = False
            self.lockAxes = []

            # Rare condition: This happens e. g. when user clicks on header menu
            # like Object->Transform->Move. These ops consume press event but not release
            # So update the snap locations anyways if there was some transformation
            if(len(self.curvePts) == 0):
                self.updateSnapLocs() # Subclass (TODO: have a relook)
                return {'RUNNING_MODAL'}

            #Looks like no 'DOUBLE_CLICK' event?
            t = time.time()
            if(self.clickT !=  None and (t - self.clickT) < DBL_CLK_DURN):
                self.confirm(context, event)
                self.redrawBezier(context, event)
                self.clickT = None
                return {'RUNNING_MODAL'}

            self.clickT = t

            if((self.pressT != None) and (t - self.pressT) < 0.2):
                loc = self.curvePts[-1][1]
                self.curvePts[-1][0] = loc
                self.curvePts[-1][2] = loc
            else:
                loc = self.get3dLocWithSnap(context, event)

            if(len(self.curvePts) == 1):
                self.curvePts[0][2] = loc # changes only rt handle

            self.curvePts.append([loc, loc, loc])
            self.redrawBezier(context, event)
            return {'RUNNING_MODAL'}

        if (event.type == 'MOUSEMOVE'):
            bpy.context.window.cursor_set("DEFAULT")
            self.redrawBezier(context, event)
            return {'RUNNING_MODAL'}

        if(len(self.curvePts) > 0 and event.type in {'X', 'Y', 'Z', 'U'}):
            if(event.type == 'U'):
                self.lockAxes = []
            else:
                self.lockAxes = [ord(event.type) - ord('X')]
                if(self.shift):
                    self.lockAxes = list({0, 1, 2} - set(self.lockAxes))
            return {'RUNNING_MODAL'}

        return {'PASS_THROUGH'}

    def redrawBezier(self, context, event):
        loc = self.get3dLocWithSnap(context, event)
        rtHandlePt = [] # additional things to draw other than curves and handles
        segColor = None
        tipColor = None
        if(not self.capture):
            if(len(self.curvePts) == 0):
                # Marker (dot), if drawing not started
                # (drawn as seg with all six pts at the same location)
                rtHandlePt = [[loc, loc, loc], [loc, loc, loc]]
                tipColor = MARKER_COLOR
            else:
                self.curvePts[-1] = [loc, loc, loc]
            handleNos = [0, 1]
        else:
            if(len(self.curvePts) == 1):
                # First handle (straight line), if user drags first pt
                # (first point from CurvePts, second the current location)
                rtHandlePt = [[loc, loc, loc]]
                segColor = HANDLE_COLOR_MAP['ALIGNED']
                handleNos  = [0, 2]
            elif(len(self.curvePts) > 1):
                end = self.curvePts[-1][1]
                ltHandle = end - (loc - end)
                rtHandle = end + (loc - end)
                self.curvePts[-1][0] = ltHandle
                self.curvePts[-1][2] = rtHandle
                handleNos  = [2, 3]
        ModalDrawBezierOp.createBatch(context, self.curvePts + rtHandlePt, \
            endSegHandleNos = handleNos, segColor = segColor, tipColor = tipColor)

    def get3dLocWithSnap(self, context, event):
        return self.get3dLoc(context, event, snapToObj = self.alt,
            snapToGrid = self.ctrl, restrict = self.shift, vec = None, fromActiveObj = False)

    #Reference point for restrict angle or lock axis
    def getPrevRefCo(self):
        if(len(self.curvePts) > 0):
            if(self.capture):
                return self.curvePts[-1][1]
            # There should always be min 2 pts if not capture, check anyway
            elif(len(self.curvePts) > 1):
                return self.curvePts[-2][1]
        return None

    # TODO Maybe too tightly coupled with the data structure
    def get3dLoc(self, context, event, snapToObj,
        snapToGrid, restrict, vec = None, fromActiveObj = True):

        rounding = getViewDistRounding(context)
        region = context.region
        rv3d = context.space_data.region_3d
        xy = event.mouse_region_x, event.mouse_region_y

        if(snapToObj):
            kd = kdtree.KDTree(len(self.getSnapLocs()))
            for i, l in enumerate(self.getSnapLocs()):
                kd.insert(getCoordFromLoc(context, l).to_3d(), i)
            kd.balance()

            coFind = Vector(xy).to_3d()
            searchResult = kd.find_range(coFind, SNAP_DIST_PIXEL)
            if(len(searchResult) != 0):
                co, idx, dist = min(searchResult, key = lambda x: x[2])
                return self.getSnapLocs()[idx]

        if(vec == None):
            if(fromActiveObj and context.active_object != None):
                vec = context.active_object.location
            else:
                # ~ vec = region_2d_to_vector_3d(region, rv3d, xy)
                vec = context.scene.cursor.location

        loc = region_2d_to_location_3d(region, rv3d, xy, vec)

        lastCo = self.getPrevRefCo()
        if(len(self.lockAxes) > 0 and lastCo != None):
            actualLoc = loc[:]
            loc = lastCo.copy()
            for axis in self.lockAxes: loc[axis] = actualLoc[axis]

        if(snapToGrid):
            loc = roundedVect(loc, rounding)

        if(restrict and lastCo != None):
            actualLoc = loc.copy()

            #First decide the main movement axis
            diff = [abs(v) for v in (actualLoc - lastCo)]
            maxDiff = max(diff)
            axis = 0 if abs(diff[0]) == maxDiff \
                else (1 if abs(diff[1]) == maxDiff else 2)

            loc = lastCo.copy()
            loc[axis] = actualLoc[axis]

            snapIncr = 45 / self.snapSteps
            snapAngles = [radians(snapIncr * a) for a in range(0, self.snapSteps + 1)]
            l1 =  actualLoc[axis] - lastCo[axis] #Main axis value

            for i in range(0, 3):
                if(i != axis):
                    l2 =  (actualLoc[i] - lastCo[i]) #Minor axis value
                    angle = abs(atan(l2 / l1)) if l1 != 0 else 0
                    dirn = (l1 * l2) / abs(l1 * l2) if (l1 * l2) != 0 else 1
                    prevDiff = 9e+99
                    for j in range(0, len(snapAngles) + 1):
                        if(j == len(snapAngles)):
                            loc[i] = lastCo[i] + dirn * l1 * tan(snapAngles[-1])
                            break
                        cmpAngle = snapAngles[j]
                        if(abs(angle - cmpAngle) > prevDiff):
                            loc[i] = lastCo[i] + dirn * l1 * tan(snapAngles[j-1])
                            break
                        prevDiff = abs(angle - cmpAngle)

        return loc


class ModalFlexiDrawBezierOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves in object mode"
    bl_idname = "wm.flexi_draw_bezier_curves"
    bl_label = "Flexi Draw Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}

    def __init__(self):
        curveDispRes = 200
        super(ModalFlexiDrawBezierOp, self).__init__(curveDispRes)

    def isDrawToolSelected(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)

        # ~ if(tool == None or tool.idname != FlexiDrawBezierTool.bl_idname): (T60766)
        if(tool == None or tool.idname != 'flexi_bezier.draw_tool'):
            return False

        return True

    def cancelOp(self, context):
        super(ModalFlexiDrawBezierOp, self).cancelOp(context)
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)

    def postUndoRedo(self, scene):
        self.updateSnapLocs()

    def invoke(self, context, event):

        # If the operator is invoked from context menu, enable the tool on toolbar
        if(not self.isDrawToolSelected(context) and context.mode == 'OBJECT'):
            # ~ bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname) (T60766)
            bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')

        # Object name -> [spline index, (startpt, endPt)]
        # Not used right now (maybe in case of large no of curves)
        self.snapInfos = {}

        self.updateSnapLocs()
        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)

        return super(ModalFlexiDrawBezierOp, self).invoke(context, event)

    def modal(self, context, event):
        if(not self.isDrawToolSelected(context)):
            self.cancelOp(context)
            return {"CANCELLED"}

        return super(ModalFlexiDrawBezierOp, self).modal(context, event)

    def getSnapLocs(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if(len(self.curvePts) > 0):
            locs.append(self.curvePts[0][1])

        return locs

    def updateSnapLocs(self, objNames = None):
        if(objNames == None):
            objNames = [o.name for o in bpy.data.objects]
            invalOs = self.snapInfos.keys() - set(objNames) # In case of redo
            for o in invalOs:
                del self.snapInfos[o]

        for objName in objNames:
            obj = bpy.data.objects.get(objName)
            if(obj != None and isBezier(obj) and obj.visible_get()):
                self.snapInfos[objName] = []
                mw = obj.matrix_world
                for i, s in enumerate(obj.data.splines):
                    if(s.use_cyclic_u == True):
                        continue
                    ends = [mw @ s.bezier_points[0].co, mw @ s.bezier_points[-1].co]
                    self.snapInfos[objName].append([i, ends])
            elif(self.snapInfos.get(objName) != None):
                del self.snapInfos[objName]

    def addPtToSpline(invMW, spline, idx, pt, handleType):
        spline.bezier_points[i].handle_left = invMW @pt[0]
        spline.bezier_points[i].co = invMW @pt[1]
        spline.bezier_points[i].handle_right = invMW @ pt[2]
        spline.bezier_points[i].handle_left_type = handleType
        spline.bezier_points[i].handle_right_type = handleType

    def createObjFromPts(self, context):
        data = bpy.data.curves.new('BezierCurve', 'CURVE')
        data.dimensions = '3D'
        obj = bpy.data.objects.new('BezierCurve', data)
        collection = context.collection
        if(collection == None):
            collection = context.scene.collection
        collection.objects.link(obj)
        obj.location = context.scene.cursor.location

        depsgraph = context.evaluated_depsgraph_get()
        depsgraph.update()

        invM = obj.matrix_world.inverted()

        spline = data.splines.new('BEZIER')
        spline.use_cyclic_u = False

        if(vectCmpWithMargin(self.curvePts[0][1], self.curvePts[-1][0])):
            spline.use_cyclic_u = True
            self.curvePts.pop()

        spline.bezier_points.add(len(self.curvePts) - 1)
        prevPt = None
        for i, pt in enumerate(self.curvePts):
            currPt = spline.bezier_points[i]
            currPt.co = invM @ pt[1]
            currPt.handle_right = invM @ pt[2]
            if(prevPt != None and prevPt.handle_right == prevPt.co \
                and pt[0] == pt[1] and currPt.co != prevPt.co): # straight line
                    diffV = (currPt.co - prevPt.co)
                    if(prevPt.handle_left_type == 'ALIGNED'):
                        prevPt.handle_left_type = 'FREE'
                    prevPt.handle_right_type = 'FREE'
                    prevPt.handle_right = prevPt.co +  STRT_SEG_HDL_LEN_COEFF * diffV
                    currPt.handle_left = currPt.co -  STRT_SEG_HDL_LEN_COEFF * diffV
                    currPt.handle_left_type = 'FREE'
                    currPt.handle_right_type = 'FREE'
            else:
                currPt.handle_left = invM @ pt[0]
                if(i == 0 or i == len(self.curvePts) -1):
                    currPt.handle_left_type = 'FREE'
                    currPt.handle_right_type = 'FREE'
                else:
                    currPt.handle_left_type = 'ALIGNED'
                    currPt.handle_right_type = 'ALIGNED'
            prevPt = currPt

        diffV = (spline.bezier_points[-1].co - spline.bezier_points[0].co)
        pt0 = spline.bezier_points[0]
        pt1 = spline.bezier_points[-1]
        if(diffV.length > 0 and pt0.handle_left == pt0.co and pt1.handle_right == pt1.co):
            pt0.handle_left = pt0.co + STRT_SEG_HDL_LEN_COEFF * diffV
            pt0.handle_left_type = 'FREE'
            pt1.handle_right = pt1.co - STRT_SEG_HDL_LEN_COEFF * diffV
            pt1.handle_right_type = 'FREE'

        return obj

    def getSnapObj(self, context, loc):
        for obj in bpy.data.objects:
            if(isBezier(obj)):
                mw = obj.matrix_world
                for i, s in enumerate(obj.data.splines):
                    p = s.bezier_points[0]
                    if(vectCmpWithMargin(loc, mw @ p.co)):
                        return obj, i, 0
                    p = s.bezier_points[-1]
                    if(vectCmpWithMargin(loc, mw @ p.co)):
                        return obj, i, -1

        return None, 0, 0

    def createCurveObj(self, context, startObj = None, \
        startSplineIdx = None, endObj = None, endSplineIdx = None):
        # First create the new curve
        obj = self.createObjFromPts(context)

        # Undo stack in case the user does not want to join
        if(endObj != None or startObj != None):
            obj.select_set(True)
            bpy.context.view_layer.objects.active = obj
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
            startObjs, changeCnt = splitCurve([startObj], split = 'spline', newColl = False)

            obj = joinSegs([startObjs[startSplineIdx], obj], \
                optimized = True, straight = False, srcCurve = startObjs[startSplineIdx])

            if(startObj == endObj and startSplineIdx != endSplineIdx):
                # If startSplineIdx == endSplineIdx the join call above would take care
                # but if they are different they need to be joined with a separate call
                obj = joinSegs([startObjs[endSplineIdx], obj], \
                    optimized = True, straight = False, srcCurve = startObjs[endSplineIdx])

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

    def save(self, context, event):
        if(len(self.curvePts) > 0):
            self.curvePts.pop()

        if(len(self.curvePts) > 1):

            endObj, endSplineIdx, ptIdx1 = \
                self.getSnapObj(context, self.curvePts[-1][1])
            startObj, startSplineIdx, ptIdx2 = \
                self.getSnapObj(context, self.curvePts[0][1])

            # ctrl pressed and there IS a snapped end obj, so user does not want connection
            # (no option to only connect to starting curve when end object exists)
            if(self.ctrl and endObj != None):
                obj = self.createCurveObj(context)
                endObj = None
            else:
                startObjName = startObj.name if(startObj != None) else ''
                endObjName = endObj.name if(endObj != None) else ''

                obj = self.createCurveObj(context, \
                    startObj, startSplineIdx, endObj, endSplineIdx)

            if(endObj == None  and self.shift \
                and (event.type == 'SPACE' or event.type == 'RET')):
                obj.data.splines[-1].use_cyclic_u = True

            #TODO: Why try?
            try:
                obj.select_set(True)
                bpy.context.view_layer.objects.active = obj
                self.updateSnapLocs([obj.name, startObjName, endObjName])
            except Exception as e:
                pass
        bpy.ops.ed.undo_push()

#(T60766)
# ~ class FlexiDrawBezierTool(WorkSpaceTool):
    # ~ bl_space_type='VIEW_3D'
    # ~ bl_context_mode='OBJECT'

    # ~ bl_idname = "flexi_bezier.draw_tool"
    # ~ bl_label = "Flexi Draw Bezier"
    # ~ bl_description = ("Flexible drawing of Bezier curves in object mode")
    # ~ bl_icon = "ops.gpencil.extrude_move"
    # ~ bl_widget = None
    # ~ bl_operator = "wm.flexi_draw_bezier_curves"
    # ~ bl_keymap = (
        # ~ ("wm.flexi_draw_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
         # ~ {"properties": []}),
    # ~ )

################### Flexi Edit Bezier Curve ###################

def getWSData(obj):
    # Less readable but more convenient than class
    # Format: [handle_left, co, handle_right, handle_left_type, handle_right_type]
    worldSpaceData = []
    mw = obj.matrix_world
    for spline in obj.data.splines:
        pts = []
        for pt in spline.bezier_points:
            pts.append([mw @ pt.handle_left, mw @ pt.co, mw @ pt.handle_right, \
                pt.handle_left_type, pt.handle_right_type])
        worldSpaceData.append(pts)
    return worldSpaceData

def getAdjIdx(obj, splineIdx, startIdx, offset = 1):
    spline = obj.data.splines[splineIdx]
    ptCnt = len(spline.bezier_points)
    if(not spline.use_cyclic_u and
        ((startIdx + offset) >= ptCnt or (startIdx + offset) < 0)):
            return None
    return (ptCnt + startIdx + offset) % ptCnt # add ptCnt for negative offset

def getCtrlPtsForSeg(obj, splineIdx, segIdx):
    wsData = getWSData(obj)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)
    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0[1], pt0[2], pt1[0], pt1[1]]

def getBezierPtsForSeg(obj, splineIdx, segIdx):
    wsData = getWSData(obj)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)
    if(segEndIdx == None):
        return []
    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0, pt1]

# Wrapper for spatial search within segment
def getClosestPt2dWithinSeg(context, coFind, selObj, selSplineIdx, selSegIdx, \
    selObjRes, withHandles, withEndPts):
    return getClosestPt2d(context, coFind, objs = [], objRes = 0, selObj = selObj, \
        selSplineIdx = selSplineIdx, selSegIdx = selSegIdx, selObjRes = selObjRes, \
            withHandles = withHandles, withEndPts = withEndPts)

# Wrapper for spatial search across objects
def getClosestPt2dAcrossObjs(context, coFind, objs, objRes = NONSEL_CURVE_SEARCH_RES):
    return getClosestPt2d(context, coFind, objs,
        objRes = NONSEL_CURVE_SEARCH_RES, selObj = None, selSplineIdx = None, \
            selSegIdx = None, selObjRes = 0)

# Requirement is more generic than geometry.interpolate_bezier
def getInterpSegPts(mw, spline, ptIdx, res, startT, endT):
    bpts = spline.bezier_points
    j = ptIdx
    if(j < (len(bpts) - 1) ):
        seg = [mw @ bpts[j].co, mw @ bpts[j].handle_right, \
            mw @ bpts[j+1].handle_left, mw @ bpts[j+1].co]
    elif(j == (len(bpts) - 1)  and spline.use_cyclic_u):
        seg = [mw @ bpts[-1].co, mw @ bpts[-1].handle_right, \
            mw @ bpts[0].handle_left, mw @ bpts[0].co]
    else:
        return []

    resProp = 2 + int(res * getSegLen(seg, error = .01)) # very rough
    if(resProp > 3):
        interpLocs = []
        interpIncr = float(endT - startT) / (resProp - 1)
        for x in range(0, resProp):
            interPt = getPtFromT(seg[0], seg[1], seg[2], seg[3],
                startT + interpIncr * x)
            interpLocs.append(interPt)
    else:
        interpLocs = [seg[0]]

    return interpLocs

# Find the list element containing the given idx from flattened list
# return the index of the list element containing the idx
def findListIdx(counts, idx):
    cumulCnt = 0
    cntIdx= 0
    while(idx >= cumulCnt):
        cumulCnt += counts[cntIdx] # cntIdx can never be >= len(counts)
        cntIdx += 1
    return cntIdx - 1

def getClosestPt2d(context, coFind, objs, objRes, selObj, selSplineIdx, selSegIdx, \
    selObjRes, withHandles = True, withEndPts = True, normalized = True, \
        objStartT = 0, objEndT = 1, selObjStartT = 0.05, selObjEndT = .95):

    objLocMap = {}

    if(normalized):
        viewDist = context.space_data.region_3d.view_distance
        objRes /= viewDist # inversely proportional
        selObjRes /= viewDist

    counts = []
    objLocList = []
    objLocs = []
    for obj in objs:
        mw = obj.matrix_world
        for i, spline in enumerate(obj.data.splines):
            for j in range(0, len(spline.bezier_points)):
                interpLocs = getInterpSegPts(mw, spline, j, objRes, objStartT, objEndT)
                counts.append(len(interpLocs))
                objLocList.append([obj, i, j])
                objLocs += interpLocs

    hdls = []
    endPts = []
    segInterpLocs = []

    if(selObj != None):
        mw = selObj.matrix_world

        spline = selObj.data.splines[selSplineIdx]
        segInterpLocs = getInterpSegPts(mw, spline, selSegIdx, selObjRes, \
            selObjStartT, selObjEndT)

        if(withHandles):
            hdls = [pt for pts in getBezierPtsForSeg(selObj, selSplineIdx, selSegIdx) \
                for pt in pts[:3]]

        if(withEndPts):
            bpts = spline.bezier_points
            endPts = [mw @ bpts[j].co for j in range(0, len(bpts))]

    searchPtsList = [[], [], [], []]
    searchPtsList[0] = hdls
    searchPtsList[1] = segInterpLocs
    searchPtsList[2] = endPts
    searchPtsList[3] = objLocs
    searchPtsList = [[getCoordFromLoc(context, pt).to_3d() \
        for pt in pts] for pts in searchPtsList]
        
    searchResults = search2dFromPtsList(searchPtsList, coFind, searchRange = 20)
    
    if(len(searchResults) == 0):
        return None

    searchResult = min(searchResults, key = lambda x: x[0])
    if(searchResult[0] == 0):
        return 'Handles', searchResult[1], searchResult[3]
    elif(searchResult[0] == 1):
        return 'SegLoc', segInterpLocs[searchResult[1]], searchResult[3]
    elif(searchResult[0] == 2):
        return 'EndPts', searchResult[1], searchResult[3]
    elif(searchResult[0] == 3):
        ptInListIdx = searchResult[1]
        listIdx = findListIdx(counts, ptInListIdx)
        obj, splineIdx, segIdx = objLocList[listIdx]
        return 'ObjLoc', obj, splineIdx, segIdx, objLocs[ptInListIdx]

def search2dFromPtsList(ptsList, coFind, searchRange):
    kd = kdtree.KDTree(sum(len(pts) for pts in ptsList))
    idx = 0
    counts = []
    for i, pts in enumerate(ptsList):
        counts.append(len(pts))
        for j, pt in enumerate(pts):
            kd.insert(pt, idx)
            idx += 1
    kd.balance()
    foundVals = kd.find_range(coFind, searchRange)
    foundVals = sorted(foundVals, key = lambda x: x[2])

    searchResults = []
    for co, idx, dist in foundVals:
        listIdx = findListIdx(counts, idx)
        ptInListIdx = idx - sum(len(ptsList[i]) for i in range(0, listIdx))
        searchResults.append([listIdx, ptInListIdx, co, dist])

    return searchResults

class EditCurveInfo:
    def __init__(self, obj, splineIdx, ptIdx):
        self.obj = obj
        self.splineIdx = splineIdx
        self.segIdx = ptIdx
        self._ctrlIdx = None # Handle and end points: 0, 1, 2 and 3, 4, 5
        self._clickLoc = None 

        # obj.name gives exception if obj is not in bpy.data.objects collection,
        # so keep a copy
        self.objName = obj.name

        # Can be derived from clickLoc, but stored to avoid repeated computation
        self._t = None

    def getClickLoc(self):
        return self._clickLoc

    def getLastSegIdx(self):
        spline = self.obj.data.splines[self.splineIdx]
        ptCnt = len(spline.bezier_points)
        # ~ if(ptCnt <= 1): return None # Condition not handled
        return ptCnt - 1 if(spline.use_cyclic_u) else ptCnt - 2

    def setClickLocSafe(self, clickLoc, lowerT = 0.05, higherT = .95):
        self._clickLoc = clickLoc
        self._t = getTForPt(self.getCtrlPts(), clickLoc)
        if(self._t < lowerT or self._t > higherT):
            if(self._t < lowerT):
                self._ctrlIdx = 1
            else:
                self._ctrlIdx = 4
            self._t = None
            self._clickLoc = None

    def getCtrlPtCoIdx(self):
        idx0 = int(self._ctrlIdx / 3)
        idx1 = self._ctrlIdx % 3
        pts = self.getSegPts()
        return pts[idx0][idx1], idx0, idx1

    def setCtrlIdxSafe(self, ctrlIdx):
        self._ctrlIdx = ctrlIdx
        if(ctrlIdx != None):
            idx0 = int(ctrlIdx / 3)
            idx1 = ctrlIdx % 3
            pts = self.getSegPts()
            # If handle pt too close to seg pt, select seg pt
            if(vectCmpWithMargin(pts[idx0][idx1], pts[idx0][1])):
                self._ctrlIdx = idx0 * 3 + 1
            self._clickLoc = None
            self.t = None

    def getSelCo(self):
        if(self._ctrlIdx != None):
            return self.getCtrlPtCoIdx()[0]
        return self._clickLoc

    # Calculate the opposite handle values in case of ALIGNED and AUTO handles
    # oldPts must not be None if hdlIdx is 1 (the end point)
    def getPtsAfterCtrlPtChange(self, pts, ptIdx, hdlIdx, oldPts = None):
        if(pts[ptIdx][3] in {'ALIGNED', 'AUTO'} and pts[ptIdx][4] in {'ALIGNED', 'AUTO'}):

            if(hdlIdx == 1): # edited the point itself
                pts[ptIdx][2] += pts[ptIdx][1] - oldPts[ptIdx][1]
                hdlIdx = 2 # Not good

            diffV = (pts[ptIdx][hdlIdx] - pts[ptIdx][1]) # Changed by user
            impIdx = (2 - hdlIdx) # 2's opposite handle is 0 and vice versa
            diffL = diffV.length
            if(diffL > 0 ):
                oldL = (pts[ptIdx][1] - pts[ptIdx][impIdx]).length #Affected
                if(round(oldL, 4) == 0.): oldL = 1
                pts[ptIdx][impIdx] = pts[ptIdx][1] - (oldL *  diffV/diffL)
                pts[ptIdx][3] = 'ALIGNED'
                pts[ptIdx][4] = 'ALIGNED'
        else:
            # Also move handles with end points to avoid weird curve shapes
            if(hdlIdx == 1):
                delta = pts[ptIdx][1] - oldPts[ptIdx][1]
                pts[ptIdx][2] += delta
                pts[ptIdx][0] += delta
            pts[ptIdx][3] = 'FREE'
            pts[ptIdx][4] = 'FREE'

        return pts

    # Get seg points after change in position of handles or drag curve
    def getOffsetSegPts(self, newPos):
        pts = self.getSegPts()

        # If ctrilIdx and clickLoc both are valid clickLoc gets priority
        if(self.getClickLoc() != None):
            delta = newPos - self.getClickLoc()
            if(delta == 0):
                return pts
            t = self._t

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

            pts[0][2] += offset0
            pts[1][0] += offset1

            # If the segment is edited, the 1st pt right handle and 2nd pt
            # left handle impacted (if handle type is aligned or auto)
            pts = self.getPtsAfterCtrlPtChange(pts, ptIdx = 0, hdlIdx  = 2)
            pts = self.getPtsAfterCtrlPtChange(pts, ptIdx = 1, hdlIdx  = 0)

        elif(self._ctrlIdx != None):
            co, ptIdx, hdlIdx = self.getCtrlPtCoIdx()
            oldPts = [p.copy() for p in pts[:3]] + pts[3:]
            pts[ptIdx][hdlIdx] = newPos
            pts = self.getPtsAfterCtrlPtChange(pts, ptIdx, hdlIdx, oldPts)

        return pts

    def getSegBezierPts(self):
        spline = self.obj.data.splines[self.splineIdx]
        return [spline.bezier_points[self.segIdx], \
            spline.bezier_points[self.getSegAdjIdx()]]

    def getPrevSegBezierPts(self):
        prevSegIdx = self.getSegAdjIdx(-1)
        if(prevSegIdx != None):
            spline = self.obj.data.splines[self.splineIdx]
            return [spline.bezier_points[prevSegIdx], \
                spline.bezier_points[self.segIdx]]
        return []

    def getNextSegBezierPts(self):
        nextSegIdx = self.getSegAdjIdx(1)
        if(nextSegIdx != None):
            spline = self.obj.data.splines[self.splineIdx]
            return [spline.bezier_points[self.segIdx], \
                spline.bezier_points[nextSegIdx]]
        return []

    def moveSeg(self, newPos):
        pts = self.getOffsetSegPts(newPos)

        invMw = self.obj.matrix_world.inverted()
        bpts = self.getSegBezierPts()

        for i, bpt in enumerate(bpts):
            bpt.handle_right_type = 'FREE'
            bpt.handle_left_type = 'FREE'

        for i, bpt in enumerate(bpts):
            bpt.handle_left = invMw @ pts[i][0]
            bpt.co = invMw @ pts[i][1]
            bpt.handle_right = invMw @ pts[i][2]

        for i, bpt in enumerate(bpts):
            bpt.handle_left_type = pts[i][3]
            bpt.handle_right_type = pts[i][4]

    def insertNode(self, handleType, select = True):
        if(self._t == None):
            return
        invMw = self.obj.matrix_world.inverted()
        bpts = self.getSegBezierPts()
        insertBezierPts(self.obj, self.splineIdx, self.segIdx, [self._t], handleType)

        if(select):
            self.setCtrlIdxSafe(4)

    def removeNode(self):
        if(self._ctrlIdx == None):
            return
        co, ptIdx, hdlIdx = self.getCtrlPtCoIdx()
        if(hdlIdx == 1):
            removeBezierPts(self.obj, self.splineIdx, {self.getSegAdjIdx(ptIdx)})
        else:
            spline = self.obj.data.splines[self.splineIdx]
            bptIdx = self.getSegAdjIdx(ptIdx)
            pt = spline.bezier_points[bptIdx]
            pt.handle_right_type = 'FREE'
            pt.handle_left_type = 'FREE'
            if(hdlIdx == 0):
                prevIdx = getAdjIdx(self.obj, self.splineIdx, bptIdx, -1)
                if(prevIdx != None):
                    ppt = spline.bezier_points[prevIdx]
                    diffV = (pt.co - ppt.co)
                else:
                    diffV = (pt.handle_right - pt.co)
                if(diffV.length == 0):
                    pt.handle_left = pt.co
                else:
                    pt.handle_left = pt.co - .2 * diffV
            else:
                nextIdx = getAdjIdx(self.obj, self.splineIdx, bptIdx, 1)
                if(nextIdx != None):
                    npt = spline.bezier_points[nextIdx]
                    diffV = (npt.co - pt.co)
                else:
                    diffV = (pt.co - pt.handle_left)
                if(diffV.length == 0):
                    pt.handle_right = pt.co
                else:
                    pt.handle_right = pt.co + .2 * diffV

            # Alway select the main point after this (should be done by caller actually)
            self._ctrlIdx = ptIdx * 3 + 1 

    # For convenience
    def getCtrlPts(self):
        return getCtrlPtsForSeg(self.obj, self.splineIdx, self.segIdx)

    def getSegPts(self):
        return getBezierPtsForSeg(self.obj, self.splineIdx, self.segIdx)

    def getSegAdjIdx(self, offset = 1):
        return getAdjIdx(self.obj, self.splineIdx, self.segIdx, offset)

    def getPrevSegPts(self):
        idx = getAdjIdx(self.obj, self.splineIdx, self.segIdx, -1)
        return getBezierPtsForSeg(self.obj, self.splineIdx, idx) if idx != None else []

    def getNextSegPts(self):
        idx = getAdjIdx(self.obj, self.splineIdx, self.segIdx)
        return getBezierPtsForSeg(self.obj, self.splineIdx, idx) if idx != None else []

    def getDisplayInfos(self, newPos = None, hltHdlIdx = None, hltEndPtIdx = None):
        displayInfos = []

        if(newPos == None):
            pts = self.getSegPts()
        else:
            pts = self.getOffsetSegPts(newPos)
            
        # Display of selected segment...
        tipColors = [TIP_COLOR, ENDPT_TIP_COLOR, TIP_COLOR, TIP_COLOR, \
            ENDPT_TIP_COLOR, TIP_COLOR]

        if(self._ctrlIdx != None):
            tipColors[self._ctrlIdx] = SEL_TIP_COLOR

        if(hltHdlIdx != None):
            tipColors[hltHdlIdx] = HLT_TIP_COLOR

        selSegDisplayInfo = SegDisplayInfo(pts, SEL_SEG_COLOR, \
            [0, 1, 2, 3], tipColors)

        # Display of non-selected segments...
        endPtTipList = [None, ADJ_ENDPT_TIP_COLOR, None]

        # Adj segs change in case of aligned and auto handles
        prevPts = self.getPrevSegPts()
        nextPts = self.getNextSegPts()
        if((len(prevPts) > 1) and (prevPts == nextPts)):
            prevPts[1] = pts[0][:]
            prevPts[0] = pts[1][:]
            displayInfos.append(SegDisplayInfo(prevPts, ADJ_SEG_COLOR, \
                [],  endPtTipList + [None, None, None]))            
        else:            
            if(len(prevPts) > 1):
                prevPts[1] = pts[0][:]
                displayInfos.append(SegDisplayInfo(prevPts, ADJ_SEG_COLOR, \
                    [],  endPtTipList + [None, None, None]))

            if(len(nextPts) > 1):
                nextPts[0] = pts[1][:]
                displayInfos.append(SegDisplayInfo(nextPts, ADJ_SEG_COLOR, \
                    [], [None, None, None] + endPtTipList))

        i = self.splineIdx
        spline = self.obj.data.splines[i]
        for j, pt in enumerate(spline.bezier_points):
            if(j == self.segIdx or j == self.getSegAdjIdx() or \
                j == self.getSegAdjIdx(-1)):
                    continue
            segPts = getBezierPtsForSeg(self.obj, i, j)

            if(segPts != []):
                startTips = endPtTipList[:] + endPtTipList[:]
                if(hltEndPtIdx == j):
                    startTips[1] = HLT_TIP_COLOR
                if(hltEndPtIdx == j + 1):
                    startTips[4] = HLT_TIP_COLOR
                di = SegDisplayInfo(segPts, ADJ_SEG_COLOR, [], startTips)
                displayInfos.append(di)
                
        # Append at the end so it's displayed on top of everything else
        displayInfos.append(selSegDisplayInfo)
        return displayInfos

class ModalFlexiEditBezierOp(Operator):
    bl_description = "Flexi editing of Bezier curves in object mode"
    bl_idname = "wm.modal_flexi_edit_bezier"
    bl_label = "Flexi Edit Curve"
    bl_options = {'REGISTER', 'UNDO'}

    running = False
    shader = None
    drawHandlerRef = None
    lineBatch = None
    pointBatch = None

    @classmethod
    def poll(cls, context):
        return not ModalFlexiEditBezierOp.running

    def addDrawHandler():
        ModalFlexiEditBezierOp.drawHandlerRef = \
            bpy.types.SpaceView3D.draw_handler_add(ModalFlexiEditBezierOp.drawHandler, \
                (), "WINDOW", "POST_VIEW")

    def removeDrawHandler(context):
        if(ModalFlexiEditBezierOp.drawHandlerRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(ModalFlexiEditBezierOp.drawHandlerRef, \
                "WINDOW")
            ModalFlexiEditBezierOp.drawHandlerRef = None

    def reset(self, context):
        self.ctrl = False
        self.alt = False
        self.shift = False
        ModalFlexiEditBezierOp.refreshDisplay(context, [])

    def cancelOp(self, context):
        ModalFlexiEditBezierOp.refreshDisplay(context, [])
        ModalFlexiEditBezierOp.removeDrawHandler(context)
        bpy.app.handlers.undo_post.remove(self.updateAfterGeomChange)
        bpy.app.handlers.redo_post.remove(self.updateAfterGeomChange)
        bpy.app.handlers.depsgraph_update_post.remove(self.updateAfterGeomChange)
        ModalFlexiEditBezierOp.running = False

    def isEditToolSelected(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)
        # ~ if(tool == None or tool.idname != FlexiEditBezierTool.bl_idname): (T60766)
        if(tool == None or tool.idname != 'flexi_bezier.edit_tool'):
            return False
        return True

    def drawHandler():
        defLineWidth = 1.5
        defPointSize = 4
        bgl.glLineWidth(defLineWidth)
        bgl.glPointSize(defPointSize)
        if(ModalFlexiEditBezierOp.lineBatch != None):
            ModalFlexiEditBezierOp.lineBatch.draw(ModalFlexiEditBezierOp.shader)
        if(ModalFlexiEditBezierOp.pointBatch != None):
            ModalFlexiEditBezierOp.pointBatch.draw(ModalFlexiEditBezierOp.shader)

    def refreshDisplay(context, displayInfos):
        #Only 3d space
        ModalFlexiEditBezierOp.lineBatch, ModalFlexiEditBezierOp.pointBatch = \
            getBezierBatches(ModalFlexiEditBezierOp.shader, displayInfos, context)

        if(hasattr(context, 'area') and context.area != None):
            context.area.tag_redraw()

    @persistent
    def loadPostHandler(dummy):
        ModalFlexiEditBezierOp.addDrawHandler()
        if(ModalFlexiEditBezierOp.shader != None):
            ModalFlexiEditBezierOp.refreshDisplay(bpy.context, [])
        ModalFlexiEditBezierOp.running = False

    @persistent
    def loadPreHandler(dummy):
        ModalFlexiEditBezierOp.removeDrawHandler(bpy.context)

    # Will be called after the curve is changed (by the tool or externally)
    # So handle all possible conditions
    def updateAfterGeomChange(self, scene = None, context = bpy.context):
        ei = self.editCurveInfo
        displayInfos = []
        if(ei != None and bpy.data.objects.get(ei.objName) != None):
            ei.obj = bpy.data.objects.get(ei.objName) #refresh anyway
            splines = ei.obj.data.splines
            spline = splines[ei.splineIdx]
            bpts = spline.bezier_points
            # Don't keep a point object / spline
            if(len(bpts) == 1):
                if(len(splines) == 1):
                    # ~ safeRemoveObj(ei.obj) # Don't remove (could be reult of external op)
                    self.editCurveInfo = None
                else:
                    splines.remove(spline)
                    if(ei.splineIdx >= (len(splines))):                    
                        ei.splineIdx = len(splines) - 1
            elif(ei.segIdx >= len(bpts) - 1):
                ei.segIdx = ei.getLastSegIdx()
        else:
            self.editCurveInfo = None

        if(self.editCurveInfo != None):
            displayInfos = self.editCurveInfo.getDisplayInfos()
        else:
            displayInfos = []

        ModalFlexiEditBezierOp.refreshDisplay(context, displayInfos)

    def invoke(self, context, event):
        ModalFlexiEditBezierOp.running = True
        ModalFlexiEditBezierOp.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
        ModalFlexiEditBezierOp.shader.bind()
        ModalFlexiEditBezierOp.addDrawHandler()

        context.window_manager.modal_handler_add(self)
        bpy.app.handlers.undo_post.append(self.updateAfterGeomChange)
        bpy.app.handlers.redo_post.append(self.updateAfterGeomChange)
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

        ModalFlexiEditBezierOp.lineBatch = None
        ModalFlexiEditBezierOp.pointBatch = None

        self.editCurveInfo = None
        self.clickT = None
        self.pressT = None
        self.isEditing = False
        self.ctrl = False
        self.shift = False
        self.alt = False
        return  {"RUNNING_MODAL"}

    def getEditableCurveObjs(self):
        return [b for b in bpy.data.objects if isBezier(b) and b.visible_get() \
                and len(b.data.splines[0].bezier_points) > 1]

    def modal(self, context, event):
        if(not self.isEditToolSelected(context)):
            self.cancelOp(context)
            return {"CANCELLED"}

        if(not self.isEditing and event.type == 'ESC' and event.value == 'RELEASE'):
            self.editCurveInfo = None
            self.reset(context)
            return {"RUNNING_MODAL"}

        if(isOutside(context, event, exclInRgns = False)):
            if(self.isEditing):
                return {'RUNNING_MODAL'}
            return {'PASS_THROUGH'}

        if(event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}):
            if(event.value == 'PRESS'): self.ctrl = True
            elif(event.value == 'RELEASE'): self.ctrl = False
            return {'PASS_THROUGH'}

        if(event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}):
            if(event.value == 'PRESS'): self.shift = True
            elif(event.value == 'RELEASE'): self.shift = False
            return {'PASS_THROUGH'}

        if(event.type in {'LEFT_ALT', 'RIGHT_ALT'}):
            if(event.value == 'PRESS'): self.alt = True
            elif(event.value == 'RELEASE'): self.alt = False
            return {'PASS_THROUGH'}

        if(self.ctrl and (not self.isEditing or (self.pressT != None and \
            time.time() - self.pressT) < SNGL_CLK_DURN)):
            bpy.context.window.cursor_set("CROSSHAIR")
        else:
            bpy.context.window.cursor_set("DEFAULT")

        if((event.type == 'E' or event.type == 'e') and event.value == 'RELEASE'):
            # ~ bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname) (T60766)
            bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')
            return {"RUNNING_MODAL"}

        if(event.type == 'DEL' and event.value == 'PRESS'):
            return {"RUNNING_MODAL"} if self.editCurveInfo != None and \
                self.editCurveInfo._ctrlIdx != None  else {"PASS_THROUGH"}

        if(event.type == 'DEL' and event.value == 'RELEASE'):
            ec = self.editCurveInfo
            if(ec != None and self.editCurveInfo._ctrlIdx != None):
                ec.removeNode()
                # will be taken care by depsgraph
                # ~ self.updateAfterGeomChange(context = context)
                bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}
            else:
                return {"PASS_THROUGH"}

        if(event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            coFind = Vector((event.mouse_region_x, event.mouse_region_y)).to_3d()
            newPos = self.get3dLoc(context, event)
            objs = self.getEditableCurveObjs()
            searchResult = 'NA' # No search done

            ci = self.editCurveInfo
            if(ci != None):
                #Search within segment
                # More precise for adding point
                res = ADD_PT_CURVE_SEARCH_RES if self.ctrl else SEL_CURVE_SEARCH_RES
                searchResult = getClosestPt2d(context, coFind, objs,
                    objRes = NONSEL_CURVE_SEARCH_RES, selObj = ci.obj, \
                        selSplineIdx = ci.splineIdx, selSegIdx = ci.segIdx, selObjRes = res, 
                            withHandles = not self.ctrl, withEndPts = not self.ctrl)
                if(searchResult != None):
                    resType = searchResult[0]
                    if(resType in {'Handles', 'SegLoc', 'EndPts' }):
                        if(resType == 'Handles'):
                            ci.setCtrlIdxSafe(searchResult[1])
                        elif(resType == 'SegLoc'):
                            ci.setClickLocSafe(searchResult[1])
                        else:
                            ptIdx = searchResult[1]
                            if(getAdjIdx(ci.obj, ci.splineIdx, ptIdx) != None):
                                ci.segIdx = ptIdx
                                ci.setCtrlIdxSafe(1)
                            # If this is the last pt, select seg ending on this pt
                            else:
                                ci.segIdx = ptIdx - 1
                                ci.setCtrlIdxSafe(4)

                        displayInfos = ci.getDisplayInfos()
                        ModalFlexiEditBezierOp.refreshDisplay(context, displayInfos)
                        self.isEditing = True
                        self.pressT = time.time()
                        return {"RUNNING_MODAL"}

            if(searchResult == None):
                return {"PASS_THROUGH"}

            if(searchResult == 'NA'):
                searchResult = getClosestPt2dAcrossObjs(context, coFind, objs,
                    objRes = NONSEL_CURVE_SEARCH_RES)

                if(searchResult == None):
                    return {"PASS_THROUGH"}

            resType, obj, splineIdx, ptIdx, dist = searchResult

            # Search for more precise clickLoc
            # Even more precise for adding point
            res = ADD_PT_CURVE_SEARCH_RES if self.ctrl else SEL_CURVE_SEARCH_RES
            searchResult = getClosestPt2dWithinSeg(context, coFind, selObj = obj, \
                    selSplineIdx = splineIdx, selSegIdx = ptIdx, \
                        selObjRes = res, withHandles = False, withEndPts = False)

            # rare condition, selection of one of the end pts (TODO: Handle better)
            if(searchResult == None): 
                return {"PASS_THROUGH"}

            co = searchResult[1]

            self.editCurveInfo = EditCurveInfo(obj, splineIdx, ptIdx)
            self.editCurveInfo.setClickLocSafe(co, lowerT = .1, higherT = .9)
            self.isEditing = True

            displayInfos = self.editCurveInfo.getDisplayInfos()
            ModalFlexiEditBezierOp.refreshDisplay(context, displayInfos)
            self.pressT = time.time()

            return {"RUNNING_MODAL"}

        if(event.type == 'MOUSEMOVE'):
            obj = None
            splineIdx = None
            ptIdx = None
            displayInfos = []
            selObjDispInfos = []
            ci = self.editCurveInfo

            if(ci != None):
                obj = ci.obj
                splineIdx = ci.splineIdx
                ptIdx = ci.segIdx
                # User is editing curve or control points (left mouse pressed)
                if(self.isEditing):
                    selCo = ci.getSelCo()
                    newPos = self.get3dLoc(context, event, selCo)

                    displayInfos = ci.getDisplayInfos(newPos)
                    ModalFlexiEditBezierOp.refreshDisplay(context, displayInfos)
                    return {"RUNNING_MODAL"}
                else:
                    selObjDispInfos = ci.getDisplayInfos()
            coFind = Vector((event.mouse_region_x, event.mouse_region_y)).to_3d()

            objs = self.getEditableCurveObjs()

            searchResult = getClosestPt2d(context, coFind, objs,
                objRes = NONSEL_CURVE_SEARCH_RES, selObj = obj, \
                    selSplineIdx = splineIdx, selSegIdx = ptIdx, \
                        selObjRes = NONSEL_CURVE_SEARCH_RES) # Low res (highlight only seg)

            if(searchResult != None):
                resType = searchResult[0]
                if(resType in {'Handles', 'EndPts'}):
                    hltIdx = searchResult[1]
                    if(resType == 'Handles'):
                        selObjDispInfos = ci.getDisplayInfos(newPos = None, \
                            hltHdlIdx = hltIdx, hltEndPtIdx = None)
                    else:
                        selObjDispInfos = ci.getDisplayInfos(newPos = None, \
                            hltHdlIdx = None, hltEndPtIdx = hltIdx)
                elif(resType == 'ObjLoc'):
                    resType, selObj, splineIdx, ptIdx, dist = searchResult
                    if(ci == None or selObj != ci.obj or splineIdx != ci.splineIdx \
                        or ptIdx != ci.segIdx):
                        editCurveInfo = EditCurveInfo(selObj, splineIdx, ptIdx)
                        displayInfos = [SegDisplayInfo(editCurveInfo.getSegPts(), \
                            NONADJ_SEG_COLOR, [], [])]

            ModalFlexiEditBezierOp.refreshDisplay(context, selObjDispInfos + displayInfos)
            return {"PASS_THROUGH"}

        if(event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):
            if(not self.isEditing):
                return {"PASS_THROUGH"}

            ei = self.editCurveInfo
            self.isEditing = False
            tm = time.time()

            if(self.clickT != None and (tm - self.clickT) < DBL_CLK_DURN):
                pass # Handle double click here

            self.clickT = tm
            if(self.pressT != None and (tm - self.pressT) < SNGL_CLK_DURN):
                if(self.ctrl):
                    if(self.shift):
                        handleType = 'ALIGNED'
                    elif(self.alt):
                        handleType = 'VECTOR'
                    else:
                        handleType = 'FREE'
                    ei.insertNode(handleType)
                    bpy.ops.ed.undo_push()
                    self.clickT = None
                    ModalFlexiEditBezierOp.refreshDisplay(context, ei.getDisplayInfos())
                    return {"RUNNING_MODAL"}
                # Gib dem Benutzer Zeit zum Atmen!
                else:
                    ModalFlexiEditBezierOp.refreshDisplay(context, ei.getDisplayInfos())
                    return {"RUNNING_MODAL"}

            self.pressT = None

            newPos = self.get3dLoc(context, event, self.editCurveInfo.getSelCo())
            self.editCurveInfo.moveSeg(newPos)
            # will be taken care by depsgraph
            # ~ self.updateAfterGeomChange(context = context) 
            bpy.ops.ed.undo_push()

        return {"PASS_THROUGH"}

    # TODO: Combine with snapping get3DLoc
    def get3dLoc(self, context, event, vec = None):
        region = context.region
        rv3d = context.space_data.region_3d
        xy = event.mouse_region_x, event.mouse_region_y
        if(vec == None):
            vec = region_2d_to_vector_3d(region, rv3d, xy)
        return region_2d_to_location_3d(region, rv3d, xy, vec)



# ~ class FlexiEditBezierTool(WorkSpaceTool):
    # ~ bl_space_type='VIEW_3D'
    # ~ bl_context_mode='OBJECT'

    # ~ bl_idname = "flexi_bezier.edit_tool"
    # ~ bl_label = "Flexi Edit Bezier"
    # ~ bl_description = ("Flexible editing of Bezier curves in object mode")
    # ~ bl_icon = "ops.pose.breakdowner"
    # ~ bl_widget = None
    # ~ bl_operator = "wm.modal_flexi_edit_bezier"
    # ~ bl_keymap = (
        # ~ ("wm.modal_flexi_edit_bezier", {"type": 'MOUSEMOVE', "value": 'ANY'},
         # ~ {"properties": []}),
    # ~ )

# ****** Temporary Workaround for Tool Not working on restart (T60766) *******

from bpy.utils.toolsystem import ToolDef
kmToolFlexiDrawBezier = "3D View Tool: Object, Flexi Draw Bezier"
kmToolFlexiEditBezier = "3D View Tool: Object, Flexi Edit Bezier"

@ToolDef.from_fn
def toolFlexiDraw():
    return dict(idname = "flexi_bezier.draw_tool",
        label = "Flexi Draw Bezier",
        description = "Flexible drawing of Bezier curves in object mode",
        icon = "ops.gpencil.extrude_move",
        widget = None,
        keymap = kmToolFlexiDrawBezier,)

@ToolDef.from_fn
def toolFlexiEdit():
    return dict(idname = "flexi_bezier.edit_tool",
        label = "Flexi Edit Bezier",
        description = "Flexible editing of Bezier curves in object mode",
        icon = "ops.pose.breakdowner",
        widget = None,
        keymap = kmToolFlexiEditBezier,)

def getToolList(spaceType, contextMode):
    from bl_ui.space_toolsystem_common import ToolSelectPanelHelper
    cls = ToolSelectPanelHelper._tool_class_from_space_type(spaceType)
    return cls._tools[contextMode]
    
def registerFlexiBezierTools():
    tools = getToolList('VIEW_3D', 'OBJECT')
    tools += None, toolFlexiDraw
    tools += None, toolFlexiEdit
    del tools

def unregisterFlexiBezierTools():
    tools = getToolList('VIEW_3D', 'OBJECT')
    index = tools.index(toolFlexiDraw) - 1 #None
    tools.pop(index)
    tools.remove(toolFlexiDraw)
    del index

    index = tools.index(toolFlexiEdit) - 1 #None
    tools.pop(index)
    tools.remove(toolFlexiEdit)

    del tools
    del index

keymapDraw = (kmToolFlexiDrawBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": [
            ("wm.flexi_draw_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
             {"properties": []}),
        ]},)

emptyKeymapDraw = (kmToolFlexiDrawBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": []},)

keymapEdit = (kmToolFlexiEditBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": [
            ("wm.modal_flexi_edit_bezier", {"type": 'MOUSEMOVE', "value": 'ANY'},
             {"properties": []}),
        ]},)

emptyKeymapEdit = (kmToolFlexiEditBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": []},)

def registerFlexiBezierKeymaps():
    keyconfigs = bpy.context.window_manager.keyconfigs
    kc_defaultconf = keyconfigs.default
    kc_addonconf = keyconfigs.addon

    from bl_keymap_utils.io import keyconfig_init_from_data
    keyconfig_init_from_data(kc_defaultconf, [emptyKeymapDraw, emptyKeymapEdit])    
    keyconfig_init_from_data(kc_addonconf, [keymapDraw, keymapEdit])

def unregisterFlexiBezierKeymaps():
    keyconfigs = bpy.context.window_manager.keyconfigs
    defaultmap = keyconfigs.get("blender").keymaps
    addonmap   = keyconfigs.get("blender addon").keymaps

    km_name, km_args, km_content = keymapDraw
    keymap = addonmap.find(km_name, **km_args)
    keymap_items = keymap.keymap_items
    for item in km_content['items']:
        item_id = keymap_items.find(item[0])
        if item_id != -1:
            keymap_items.remove(keymap_items[item_id])
    addonmap.remove(keymap)
    
    defaultmap.remove(defaultmap.find(km_name, **km_args))

    km_name, km_args, km_content = keymapEdit
    keymap = addonmap.find(km_name, **km_args)
    keymap_items = keymap.keymap_items
    for item in km_content['items']:
        item_id = keymap_items.find(item[0])
        if item_id != -1:
            keymap_items.remove(keymap_items[item_id])
    addonmap.remove(keymap)
    
    defaultmap.remove(defaultmap.find(km_name, **km_args))

# ****************** Category Configuration In User Preferences ******************

def updatePanel(self, context):
    try:
        panel = BezierUtilsPanel
        if "bl_rna" in panel.__dict__:
            bpy.utils.unregister_class(panel)

        panel.bl_category = context.preferences.addons[__name__].preferences.category
        bpy.utils.register_class(panel)

    except:
        print("BezierUtils: Updating Panel locations has failed")
        
class BezierUtilsPreferences(AddonPreferences):
    bl_idname = __name__

    category: StringProperty(
            name = "Tab Category",
            description = "Choose a name for the category of the panel",
            default = "Tool",
            update = updatePanel
    )

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        col = row.column()
        col.label(text="Tab Category:")
        col.prop(self, "category", text="")

classes = (
    ModalMarkSegStartOp,
    SeparateSplinesObjsOp,
    SplitBezierObjsOp,
    splitBezierObjsPtsOp,
    JoinBezierSegsOp,
    CloseSplinesOp,
    CloseStraightOp,
    OpenSplinesOp,
    RemoveDupliVertCurveOp,
    convertTo2DMeshOp,
    SetHandleTypesOp,
    SelectInCollOp,
    InvertSelOp,
    BezierUtilsPanel,
    ModalFlexiDrawBezierOp,
    ModalFlexiEditBezierOp,
    BezierUtilsPreferences,
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    
    # ~ bpy.utils.register_tool(FlexiDrawBezierTool) (T60766)
    # ~ bpy.utils.register_tool(FlexiEditBezierTool) (T60766)
    registerFlexiBezierTools()
    registerFlexiBezierKeymaps()
   
    bpy.app.handlers.load_post.append(ModalDrawBezierOp.loadPostHandler)
    bpy.app.handlers.load_pre.append(ModalDrawBezierOp.loadPreHandler)
    bpy.app.handlers.load_post.append(ModalFlexiEditBezierOp.loadPostHandler)
    bpy.app.handlers.load_pre.append(ModalFlexiEditBezierOp.loadPreHandler)

def unregister():
    bpy.app.handlers.load_post.remove(ModalDrawBezierOp.loadPostHandler)
    bpy.app.handlers.load_pre.remove(ModalDrawBezierOp.loadPreHandler)
    bpy.app.handlers.load_post.remove(ModalFlexiEditBezierOp.loadPostHandler)
    bpy.app.handlers.load_pre.remove(ModalFlexiEditBezierOp.loadPreHandler)

    unregisterFlexiBezierKeymaps()
    unregisterFlexiBezierTools()

    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
        
    # ~ bpy.utils.unregister_tool(FlexiDrawBezierTool) (T60766)
    # ~ bpy.utils.unregister_tool(FlexiEditBezierTool) (T60766)
    
