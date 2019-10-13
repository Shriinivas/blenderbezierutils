#
#
# Blender add-on with tools to draw and edit Bezier curves along with other utility ops
#
# Supported Blender Version: 2.8
#
# Copyright (C) 2019  Shrinivas Kulkarni

# License: GPL (https://github.com/Shriinivas/blenderbezierutils/blob/master/LICENSE)
#

import bpy, bmesh, bgl, gpu
from bpy.props import BoolProperty, IntProperty, EnumProperty, \
FloatProperty, StringProperty, CollectionProperty, FloatVectorProperty
from bpy.types import Panel, Operator, WorkSpaceTool, AddonPreferences
from mathutils import Vector, Matrix, geometry, kdtree
from math import log, atan, tan, sin, cos, pi, radians, degrees, sqrt
from bpy_extras.view3d_utils import region_2d_to_vector_3d, region_2d_to_location_3d, \
region_2d_to_origin_3d
from bpy_extras.view3d_utils import location_3d_to_region_2d
from gpu_extras.batch import batch_for_shader
import time
from bpy.app.handlers import persistent
from gpu_extras.presets import draw_circle_2d

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "version": (0, 9, 61),
    "location": "Properties > Active Tool and Workspace Settings > Bezier Utilities",
    "description": "Collection of Bezier curve utility ops",
    "category": "Object",
    "wiki_url": "https://github.com/Shriinivas/blenderbezierutils/blob/master/README.md",
    "blender": (2, 80, 0),
}

DEF_ERR_MARGIN = 0.0001
LARGE_NO = 9e+9

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
            copyBezierPt(srcSpline.bezier_points[i], \
                spline.bezier_points[ptIdx], freeHandles)
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

def insertBezierPts(obj, splineIdx, startIdx, cos, handleType):

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
        t = getTForPt(seg, co)
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
def roundedVect(vect, rounding, axes):
    rounding += 1
    fact = (10 ** rounding) / 10
    retVect = vect.copy()
    # ~ Vector([round(vect[i] / fact) * fact for i in axes])
    for i in axes: retVect[i] = round(vect[i] / fact) * fact
    return retVect

###################### Screen functions ######################

def get3dLoc(context, event, vec = None):
    region = context.region
    rv3d = context.space_data.region_3d
    xy = event.mouse_region_x, event.mouse_region_y
    if(vec == None):
        vec = region_2d_to_vector_3d(region, rv3d, xy)
    return region_2d_to_location_3d(region, rv3d, xy, vec)

def  getViewDistRounding(rv3d):
    viewDist = rv3d.view_distance
    return int(log(viewDist, 10)) - 1

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

def is3DVireport(context):
    return (context != None and hasattr(context, 'space_data') and \
        context.space_data != None and hasattr(context.space_data, 'region_3d'))

def getPtProjOnPlane(region, rv3d, xy, p1, p2, p3, p4 = None):
    vec = region_2d_to_vector_3d(region, rv3d, xy)
    orig = region_2d_to_origin_3d(region, rv3d, xy)

    pt = geometry.intersect_ray_tri(p1, p2, p3, vec, orig, False)#p4 != None)
    # ~ if(not pt and p4):
        # ~ pt = geometry.intersect_ray_tri(p2, p4, p3, vec, orig, True)
    return pt

def getLineTransMatrices(pt0, pt1):
    diffV = (pt1 - pt0)
    invTm = diffV.to_track_quat('X', 'Z').to_matrix().to_4x4()
    tm = invTm.inverted()
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
        return None, None, -1
    viewVect = region_2d_to_vector_3d(region, rv3d, xy)
    rayOrig = region_2d_to_origin_3d(region, rv3d, xy)
    mw = obj.matrix_world
    invMw = mw.inverted()

    rayTarget = rayOrig + viewVect
    rayOrigObj = invMw @ rayOrig
    rayTargetObj = invMw @ rayTarget
    rayDirObj = rayTargetObj - rayOrigObj

    success, location, normal, faceIdx = obj.ray_cast(rayOrigObj, rayDirObj)
    return mw @ location, normal, faceIdx

def getSelFaceLoc(region, rv3d, xy, maxFaceCnt):
    objCnt = 0
    aos = [bpy.context.object] if bpy.context.object != None else []
    objs = bpy.context.selected_objects + aos
    for obj in objs:
        if(obj.type == 'MESH'):
            if(not isPtIn2dBBox(obj, region, rv3d, xy)): continue
            loc, normal, faceIdx = getFaceUnderMouse(obj, region, rv3d, xy, maxFaceCnt)
            if(faceIdx >= 0):
                return obj, loc, normal, faceIdx
            objCnt += 1
            if(objCnt > maxFaceCnt): return None, None, None, -1
        # For edge snapping
        # ~ if(faceIdx >=0):
            # ~ eLoc = getEdgeUnderMouse(region, rv3d, vec, obj, faceIdx, loc)
            # ~ if(eLoc != None): loc = eLoc
    return None, None, None, -1


# ~ def getEdgeUnderMouse(region, rv3d, vec, obj, faceIdx, loc):
    # ~ p1 = (0, 0)
    # ~ p2 = (0, SNAP_DIST_PIXEL)
    # ~ minEdgeDist =
        # ~ (region_2d_to_location_3d(region, rv3d, p1, vec) -
            # ~ region_2d_to_location_3d(region, rv3d, p2, vec)).length
    # ~ closestLoc = None
    # ~ for ek in obj.data.polygons[faceIdx].edge_keys:
        # ~ if(len(ek) == 0):
            # ~ p1 = obj.data.vertices[ek[0]]
            # ~ p2 = obj.data.vertices[ek[1]]
            # ~ iLoc = geometry.intersect_point_line(loc, p1, p2)[0]
            # ~ ptDist = (iLoc - iLoc).length
            # ~ if(ptDist < minEdgeDist):
                # ~ minEdgeDist = ptDist
                # ~ closestLoc = iLoc
    # ~ if(closestLoc != None): loc = closestLoc

def get2dBBox(obj, region, rv3d):
    mw = obj.matrix_world
    bbox = obj.bound_box
    co2ds = [getCoordFromLoc(region, rv3d, mw @ Vector(b)) for b in bbox]
    minX = min(c[0] for c in co2ds)
    maxX = max(c[0] for c in co2ds)
    minY = min(c[1] for c in co2ds)
    maxY = max(c[1] for c in co2ds)

    return minX, minY, maxX, maxY

def isPtIn2dBBox(obj, region, rv3d, xy, extendBy = 0):
    minX, minY, maxX, maxY = get2dBBox(obj, region, rv3d)
    if(xy[0] > (minX - extendBy) and xy[0] < (maxX + extendBy) \
        and xy[1] > (minY - extendBy) and xy[1] < (maxY + extendBy)):
        return True
    else: return False

def isLocIn2dBBox(obj, region, rv3d, loc):
    pt = getCoordFromLoc(region, rv3d, loc)
    return isPtIn2dBBox(pt)

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
            closestCurve2, dist2 = getClosestCurve(srcMW, ncStart.co, \
                remainingCurves, dist)
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
                vectCmpWithMargin(mw @ nextSpline.bezier_points[i].co, \
                    srcMW @ currSpline.bezier_points[0].co)):

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
        return meshObj

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

        selObjs = bpy.context.selected_objects
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


class SetCurveColorOp(bpy.types.Operator):
    bl_idname = "object.set_curve_color"
    bl_label = "Set Color"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Set color of selected curves"

    def execute(self, context):
        curves = bpy.context.selected_objects

        for curve in curves:
            curve.data['curveColor'] = bpy.context.scene.curveColorPick

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

    #TODO: Move to params class
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

    bpy.types.Scene.curveColorExpanded = BoolProperty(name="Set Curve Colors State",
            default = False)

    bpy.types.Scene.otherExpanded = BoolProperty(name="Other Expanded State",
            default = False)

    bpy.types.Scene.curveColorPick = bpy.props.FloatVectorProperty(
        name="Color",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(1.0, 1.0, 1.0, 1.0)
    )

    bpy.types.Scene.applyCurveColor = BoolProperty(name="Apply Curve Color", \
        description='Apply color to all non selected curves ', \
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
                icon="TRIA_DOWN" if context.scene.handleTypesExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Set Handle Type', icon='MOD_CURVE')

            if context.scene.handleTypesExpanded:
                col = layout.column()
                row = col.row()
                col.prop(context.scene, 'handleType')
                col = layout.column()
                col.operator('object.set_handle_types')

            ######## Curve Color #########

            col = layout.column()
            col.separator()

            row = layout.row()
            row.prop(context.scene, "curveColorExpanded",
                icon="TRIA_DOWN" if context.scene.curveColorExpanded else "TRIA_RIGHT",
                icon_only=True, emboss=False
            )
            row.label(text='Set Curve Colors', icon='MATERIAL')

            if context.scene.curveColorExpanded:
                col = layout.column()
                row = col.row()
                row.prop(context.scene, "curveColorPick", text = 'Curve Color')
                row.operator('object.set_curve_color')
                row.operator('object.remove_curve_color')
                row = col.row()
                row.prop(context.scene, 'applyCurveColor', toggle = True)
                col = layout.column()

            ######## Other Tools #########

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
                col.operator('object.paste_length')
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

    ################ Stand-alone handler for changing curve colors #################

    drawHandlerRef = None
    shader = None
    lineBatch = None
    lineWidth = 1.5

    @persistent
    def colorCurves(scene = None, add = False, remove = False):
        def ccDrawHandler():
            bgl.glLineWidth(BezierUtilsPanel.lineWidth)
            if(BezierUtilsPanel.lineBatch != None):
                BezierUtilsPanel.lineBatch.draw(BezierUtilsPanel.shader)

        if(add and BezierUtilsPanel.drawHandlerRef == None):
            try:
                BezierUtilsPanel.lineWidth = \
                    bpy.context.preferences.addons[__name__].preferences.drawLineWidth
            except Exception as e:
                # ~ print("BezierUtils: Error fetching line width in ColorCurves: ", e)
                BezierUtilsPanel.lineWidth = 1.5

            BezierUtilsPanel.drawHandlerRef = \
                bpy.types.SpaceView3D.draw_handler_add(ccDrawHandler, \
                    (), "WINDOW", "POST_VIEW")
            BezierUtilsPanel.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
            BezierUtilsPanel.shader.bind()
            return

        elif(remove):
            if(BezierUtilsPanel.drawHandlerRef != None):
                bpy.types.SpaceView3D.draw_handler_remove(BezierUtilsPanel.drawHandlerRef, \
                    "WINDOW")
                BezierUtilsPanel.drawHandlerRef = None
                return

        if(bpy.context.scene.applyCurveColor):
            objs = [o for o in bpy.context.scene.objects if(isBezier(o) and \
                o.visible_get() and len(o.modifiers) == 0 and not o.select_get())]

            lineCos = []
            lineColors = []
            for o in objs:
                colorVal = o.data.get('curveColor')
                if(colorVal != None):
                    for i, spline in enumerate(o.data.splines):
                        for j in range(0, len(spline.bezier_points)):
                            segPts = getBezierDataForSeg(o, i, j)
                            if(segPts == None):
                                continue
                            pts = getPtsAlongBezier2D(segPts, getAllAreaRegions(), \
                                DEF_CURVE_RES_2D)
                            linePts = getLinesFromPts(pts)
                            lineCos += linePts
                            lineColors += [colorVal for i in range(0, len(linePts))]
            BezierUtilsPanel.lineBatch = batch_for_shader(BezierUtilsPanel.shader, \
                "LINES", {"pos": lineCos, "color": lineColors})
        else:
            BezierUtilsPanel.lineBatch = batch_for_shader(BezierUtilsPanel.shader, \
                "LINES", {"pos": [], "color": []})

        areas = [a for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
        for a in areas:
            a.tag_redraw()


################### Common Bezier Functions & Classes ###################

def getPtFromT(p0, p1, p2, p3, t):
    c = (1 - t)
    pt = (c ** 3) * p0 + 3 * (c ** 2) * t * p1 + \
        3 * c * (t ** 2) * p2 + (t ** 3) * p3
    return pt

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
def getTForPt(curve, testPt):
    minLen = LARGE_NO
    retT = None
    for coIdx in range(0, 3):
        ts = getTsForPt(curve[0], curve[1], curve[2], curve[3], testPt[coIdx], coIdx)
        for t in ts:
            pt = getPtFromT(curve[0], curve[1], curve[2], curve[3], t)
            pLen = (testPt - pt).length
            if(pLen < minLen):
                minLen = pLen
                retT = t
    return retT

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

# Get pt coords along curve defined by the four control pts (segPts)
# subdivPerUnit: No of subdivisions per unit length
# (which is the same as no of pts excluding the end pts)
def getInterpBezierPts(segPts, subdivPerUnit, segLens = None):
    if(len(segPts) < 2):
        return []

    curvePts = []
    for i in range(1, len(segPts)):
        seg = [segPts[i-1][1], segPts[i-1][2], segPts[i][0], segPts[i][1]]
        if(segLens != None and len(segLens) > (i-1)):
            res = int(segLens[i-1] * subdivPerUnit)
        else:
            res = int(getSegLen(seg) * subdivPerUnit)
        if(res > 1):
            curvePts += geometry.interpolate_bezier(*seg, res)

    return curvePts

# Used in functions where actual locs of pts on curve matter (like subdiv Bezier)
# (... kind of expensive)
def getPtsAlongBezier3D(segPts, rv3d, curveRes, minRes = 200):

    viewDist = rv3d.view_distance

    # (the smaller the view dist (higher zoom level),
    # the higher the num of subdivisions
    curveRes = curveRes / viewDist

    if(curveRes < minRes): curveRes = minRes

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes)

# Used in functions where only visual resolution of curve matters (like draw Bezier)
# (... not so expensive)
def getPtsAlongBezier2D(segPts, areaRegionInfo, curveRes):
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

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes, segLens = segLens)

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

    if(not vectCmpWithMargin(curvePts[0], curvePts[-1])):
        vertCos.append(curvePts[-1])

    return vertCos

################### Common to Draw and Edit Flexi Bezier Ops ###################

# Some global constants

DEF_CURVE_RES_2D = .5 # Per pixel seg divisions (.5 is one div per 2 pixel units)
SNAP_DIST_PIXEL = 20
STRT_SEG_HDL_LEN_COEFF = 0.25
DBL_CLK_DURN = 0.25
SNGL_CLK_DURN = 0.3

SEL_CURVE_SEARCH_RES = 1000
NONSEL_CURVE_SEARCH_RES = 100
ADD_PT_CURVE_SEARCH_RES = 5000

class PtDisplayInfo:
    # handleNos: 0: seg1-left, 1: seg1-right
    # tipColors: leftHdl, pt, rightHdl
    # Caller to make sure there are no tips without handle
    def __init__(self, pt, handleNos, tipColors):
        self.pt = pt
        self.handleNos = handleNos
        self.tipColors = tipColors

class SegDisplayInfo:
    # ~ def __init__(self, segPts, segColor):
        # ~ self.segPts = segPts
        # ~ self.segColor = segColor
    def __init__(self, segPts, segColor, handleNos = [], tipColors = []):
        self.segPts = segPts
        self.segColor = segColor
        self.handleNos = handleNos
        self.tipColors = tipColors

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
            if(len(area.spaces[0].region_quadviews) > 0):
                qIdx = getWindowRegionIdx(area, j)
                rv3d = area.spaces[0].region_quadviews[qIdx]
            else:
                rv3d = area.spaces[0].region_3d
            xy = [xyScreen[0] - region.x, xyScreen[1] - region.y]
            return RegionMouseXYInfo(area, region, rv3d, xy, xyScreen)

    def __init__(self, area, region, rv3d, xy, xyScreen):
        self.area = area
        self.region = region
        self.rv3d = rv3d
        self.xy = xy
        self.xyScreen = xyScreen

    def __eq__(self, other):
        if(other == None): return False
        return self.area == other.area and self.region == other.region and \
            self.rv3d == other.rv3d

def getSubdivBatches(shader, subdivCos, showSubdivPts):
        ptSubDivCos = [] if(not showSubdivPts) else subdivCos
        ptBatch = batch_for_shader(ModalDrawBezierOp.shader, \
            "POINTS", {"pos": ptSubDivCos, "color": [ModalBaseFlexiOp.colGreaseSubdiv \
                for i in range(0, len(ptSubDivCos))]})

        lineCos = getLinesFromPts(subdivCos)
        lineBatch = batch_for_shader(ModalDrawBezierOp.shader, \
            "LINES", {"pos": lineCos, "color": [ModalBaseFlexiOp.colGreaseNonHltSeg \
                for i in range(0, len(lineCos))]})

        return ptBatch, lineBatch

# Return line batch for bezier line segments and handles and point batch for handle tips
def getBezierBatches(shader, displayInfos, areaRegionInfo, \
    ptDispInfos = None, defHdlType = 'ALIGNED'):

    lineCos = [] #segment is also made up of lines
    lineColors = []
    for i, displayInfo in enumerate(displayInfos):
        segPts = displayInfo.segPts
        pts = getPtsAlongBezier2D(segPts, areaRegionInfo, DEF_CURVE_RES_2D)
        segLineCos = getLinesFromPts(pts)
        lineCos += segLineCos
        lineColors += [displayInfo.segColor for j in range(0, len(segLineCos))]

    tipColInfo = []
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

            lineColors += [ModalBaseFlexiOp.hdlColMap[htype], \
                ModalBaseFlexiOp.hdlColMap[htype]]

        for j, tipColor in enumerate(displayInfo.tipColors):
            if(tipColor != None):
                # [6] array to [2][3] array (includes end points of curve)
                ptIdx = int(j / 3)
                hdlIdx = j % 3
                tipColInfo.append([tipColor, segPts[ptIdx][hdlIdx]])
                
    if(ptDispInfos):
        for i, ptInfo in enumerate(ptDispInfos):
            pt = ptInfo.pt
            for hn in ptInfo.handleNos:
                lineCos += [pt[hn], pt[hn + 1]]
                if(len(pt) < 5):
                    htype = defHdlType
                else:
                    htype = segPts[ptIdx][3 + hdlIdx]
                lineColors += [ModalBaseFlexiOp.hdlColMap[htype], \
                        ModalBaseFlexiOp.hdlColMap[htype]]        
            for j, tipColor in enumerate(ptInfo.tipColors):
                if(tipColor != None): tipColInfo.append([tipColor, pt[j]])
        
    tipCos = []
    tipColors = []

    tipColInfo = sorted(tipColInfo, \
        key = lambda x: ModalBaseFlexiOp.tipColPriority[x[0]])

    for ti in tipColInfo:
        tipColors.append(ti[0])
        tipCos.append(ti[1])

    lineBatch = batch_for_shader(shader, "LINES", {"pos": lineCos, "color": lineColors})
    tipBatch = batch_for_shader(shader, "POINTS", {"pos": tipCos, "color": tipColors})

    return lineBatch, tipBatch

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

    if(keymap.get(event.type) != None):
        expr = 'caller.' + keymap[event.type] + ' = '
        if(event.value == 'PRESS'): exec(expr +'True')
        if(event.value == 'RELEASE'): exec(expr +'False')
        return True

    return False

unitMap = {'FEET': "'", 'METERS':'m'}

class SnapDigits():
    digitMap = {'ONE':'1', 'TWO':'2', 'THREE':'3', 'FOUR':'4', 'FIVE':'5', \
                'SIX':'6', 'SEVEN':'7', 'EIGHT':'8', 'NINE':'9', 'ZERO':'0', 'PERIOD':'.'}

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

    def procEvent(self, context, event):
        if(event.type == 'P'):
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

        dval = SnapDigits.digitMap.get(event.type)
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

        if(self.polar):
            axis0, axis1 = self.getFreeAxes()[0], self.getFreeAxes()[1]
            val[axis0], val[axis1] = self.addToPolar(delta)
        else:
            val[self.axisIdx] += delta

        return val

    def getDeltaStrPolar(self):
        delta, valid = SnapDigits.getValidFloat(self.signChar, self.digitChars)
        polCos = self.getPolarCos()
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
        d = self.deltaVec[self.axisIdx]
        if(d != 0):
            retStr += str(round(d, 4)) + ('+' if (self.signChar == '') else '')
        retStr += self.signChar + ''.join(self.digitChars) + '] = ' \
            + str(round((delta + d), 4)) + ('' if valid else ' <Invalid>')

        return retStr

class CustomAxis:

    def __init__(self):
        #TODO: What's better?
        params = bpy.context.window_manager.bezierToolkitParams
        if(bpy.data.scenes[0].get('btk_co1') == None):
            bpy.data.scenes[0]['btk_co1'] = [LARGE_NO, LARGE_NO, LARGE_NO]
        if(bpy.data.scenes[0].get('btk_co2') == None):
            bpy.data.scenes[0]['btk_co2'] = [LARGE_NO, LARGE_NO, LARGE_NO]
        self.axisPts = [Vector(bpy.data.scenes[0]['btk_co1']), \
            Vector(bpy.data.scenes[0]['btk_co2'])]
        self.snapCnt = params.customAxisSnapCnt

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

class Snapper():

    DEFAULT_ANGLE_SNAP_STEPS = 3
    MAX_SNAP_VERT_CNT = 1000
    MAX_SNAP_FACE_CNT = 1000

    def __init__(self, context, getSnapLocs, getRefLine, getRefLineOrig, \
        hasSelection, isEditing):
        self.getSnapLocs = getSnapLocs
        self.getRefLine = getRefLine
        self.getRefLineOrig = getRefLineOrig
        self.hasSelection = hasSelection
        self.isEditing = isEditing
        self.angleSnapSteps = Snapper.DEFAULT_ANGLE_SNAP_STEPS
        self.customAxis = CustomAxis()
        self.snapDigits = SnapDigits(self.getFreeAxesNormalized, self.getEditCoPair)
        self.initialize()
        # ~ bpy.context.space_data.overlay.show_axis_x = False
        # ~ bpy.context.space_data.overlay.show_axis_y = False

    def initialize(self):
        self.shift = False # For locking to plane
        self.ctrl = False # For excluding ctrl Z etc.
        self.alt = False # For future use

        # Duplicate of self.ctrl etc. but better for readability
        self.angleSnap = False
        self.gridSnap = False
        self.objSnap = False

        self.tm = None
        self.orig = None

        self.lastSnapTypes = set()

        self.resetSnap()

    def resetSnap(self): # Called even during isEditing

        self.freeAxes = [] # All axes free
        self.snapDigits.initialize()
        self.rmInfo = None

        # This variable lets caller know that return was pressed after digits were entered
        # Caller can reset snapper as per convenience
        self.digitsConfirmed = False
        self.snapCo = None

        self.inDrawAxis = False # User drawing the custom axis
        # ~ self.snapStack = [] # TODO

    def getGlobalOrient(self):
        return bpy.context.window_manager.bezierToolkitParams.snapOrient

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

    def getCurrOrig(self, obj, rmInfo):
        snapOrigin = bpy.context.window_manager.bezierToolkitParams.snapOrigin
        if(snapOrigin == 'AXIS'):
            if(self.customAxis.length() != 0): return self.customAxis.axisPts[0]
        elif(snapOrigin == 'REFERENCE'):
            refLineOrig = self.getRefLineOrig()
            if(refLineOrig != None): return refLineOrig
        elif(snapOrigin == 'OBJECT' and obj != None): return obj.location
        elif(snapOrigin == 'FACE' and rmInfo != None):
            selObj, location, normal, faceIdx = getSelFaceLoc(rmInfo.region, \
                rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT)
            if(faceIdx >= 0):
                    return selObj.matrix_world @ selObj.data.polygons[faceIdx].center
        elif(snapOrigin == 'CURSOR'): return bpy.context.scene.cursor.location
        return Vector((0, 0, 0))

    def getTransMatsForOrient(self, rmInfo, obj = None):
        tm = Matrix()
        params = bpy.context.window_manager.bezierToolkitParams
        orientType = params.snapOrient

        if(orientType == 'AXIS'):
            ca = self.customAxis
            if(abs(ca.length()) > DEF_ERR_MARGIN):
                pts = ca.axisPts
                tm, invTm = getLineTransMatrices(pts[0], pts[1])
                if(params.axisScale):
                    unitD = ca.length() / 10
                    tm = Matrix.Scale(1 / unitD, 4) @ tm
                    invTm = tm.inverted()

        elif(orientType == 'REFERENCE'):
            refLine = self.getRefLine()
            if(refLine != None and len(refLine) == 2):
                tm, invTm = getLineTransMatrices(refLine[0], refLine[1])
                if(params.axisScale):
                    unitD = (refLine[1] - refLine[0]).length / 10
                    tm = Matrix.Scale(1 / unitD, 4) @ tm
                    invTm = tm.inverted()

        elif(orientType == 'VIEW'):
            tm = rmInfo.rv3d.view_matrix

        if(obj != None and orientType == 'OBJECT'):
                tm = obj.matrix_world.inverted()

        if(orientType == 'FACE'):
            selObj, location, normal, faceIdx = getSelFaceLoc(rmInfo.region, \
                rmInfo.rv3d, rmInfo.xy, self.MAX_SNAP_FACE_CNT)
            if(faceIdx >= 0):
                normal = selObj.data.polygons[faceIdx].normal
                tm = normal.to_track_quat('Z', 'X').to_matrix().to_4x4().inverted()

        return tm, tm.inverted()

    # To be called in modal method of parent
    def procEvent(self, context, event, rmInfo):

        # update ctrl etc.
        updateMetaBtns(self, event)

        keymap = {'LEFT_SHIFT': 'angleSnap', 'RIGHT_SHIFT': 'angleSnap',
            'LEFT_CTRL':'gridSnap', 'RIGHT_CTRL':'gridSnap',
            'LEFT_ALT': 'objSnap', 'RIGHT_ALT': 'objSnap'}

        # will do the same as above call, but update different vars for readability
        updateMetaBtns(self, event, keymap)

        if(event.type == 'RIGHTMOUSE'):
            snapOrigin = bpy.context.window_manager.bezierToolkitParams.snapOrigin
            if(event.value == 'RELEASE' and snapOrigin == 'AXIS'):
                loc = self.get3dLocSnap(rmInfo, snapToAxisLine = False)
                if(not self.inDrawAxis): self.customAxis.set(0, loc)
                else: self.customAxis.set(1, loc)
                self.inDrawAxis = not self.inDrawAxis
            return True

        if(self.inDrawAxis):
            if(event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', 'NUMPAD_PLUS', \
                'NUMPAD_MINUS','PLUS', 'MINUS'}):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return True
                elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                    if(self.customAxis.snapCnt < 20): self.customAxis.snapCnt += 1
                elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                    if(self.customAxis.snapCnt > 0): self.customAxis.snapCnt -= 1

            if(event.type == 'MOUSEMOVE'):
                loc = self.get3dLocSnap(rmInfo, snapToAxisLine = False)
                self.customAxis.set(1, loc)

            if(event.type == 'ESC'):
                self.customAxis.set(0, Vector((LARGE_NO, LARGE_NO, LARGE_NO)))
                self.customAxis.set(1, Vector((LARGE_NO, LARGE_NO, LARGE_NO)))
                self.inDrawAxis = False

            return True

        retVal = False # Whether the event was 'consumed'

        snapDProc = False
        if(len(self.getRefLine()) > 0):
            snapDProc = self.snapDigits.procEvent(context, event)
            if(snapDProc):
                self.digitsConfirmed = False # Always reset if there was any digit entered
                return True

        if(not self.ctrl and event.type == 'U'):
            if(event.value == 'RELEASE'):
                self.tm = None # Force reorientation
                self.orig = None # Force origin shift
            return True

        if(not self.ctrl and event.type in {'X', 'Y', 'Z'}):
            self.digitsConfirmed = False # Always reset if there is any lock axis
            if(event.value == 'RELEASE'):
                self.freeAxes = [ord(event.type) - ord('X')]
                if(self.shift):
                    self.freeAxes = sorted(list({0, 1, 2} - set(self.freeAxes)))

                # if already part of global axes, don't store (no escape needed)
                if(self.getFreeAxesCombined() == self.getFreeAxesGlobal()):
                    self.freeAxes = []

            return True

        # ~ if(event.type == 'Q' and self.getGlobalOrient() == 'AXIS'):
            # ~ if(event.value == 'RELEASE'):
                # ~ self.customAxis.initialize()
            # ~ return True

        # ~ if(event.type == 'WHEELDOWNMOUSE' and self.angleSnap):# Event not consumed
            # ~ if(self.snapSteps < 10):
                # ~ self.snapSteps += 1

        # ~ if(event.type == 'WHEELUPMOUSE' and self.angleSnap):# Event not consumed
            # ~ if(self.snapSteps > 1):
                # ~ self.snapSteps -= 1

        # ~ if(event.type == 'MIDDLEMOUSE' and self.angleSnap):# Event not consumed
            # ~ self.snapSteps = self.defaultSnapSteps

        # Consume escape or return / space only if there's something to process
        if(len(self.freeAxes) > 0 or \
            (self.snapDigits.hasVal() and not self.digitsConfirmed)):
            retVal = True
            if(event.type == 'RET' or event.type == 'SPACE'):
                if(event.value == 'RELEASE'):
                    # ~ self.resetSnap() # This is the responsibility of the caller
                    self.digitsConfirmed = True # Confirm first time
            elif(event.type == 'ESC'):
                if(event.value == 'RELEASE'): self.resetSnap()
            else:
                retVal = snapDProc

        return retVal

    def getStatusStr(self, unit, invTm, refPt, newPt):
        manualEntry = self.snapDigits.hasVal()

        if(manualEntry and self.snapDigits.polar):
            return self.snapDigits.getDeltaStrPolar()

        axes = self.getFreeAxesNormalized()
        diffV = self.snapDigits.getCurrDelta() if manualEntry else (newPt - refPt)
        diffVActual = invTm @ diffV

        retStr = ''
        transformed = invTm != Matrix()

        for i, d in enumerate(diffV):
            if(i not in axes): continue

            retStr += 'D' + chr(ord('x') + i) + ': '
            if(manualEntry and i == self.snapDigits.axisIdx):
                retStr += self.snapDigits.getCurrDeltaStr()
            else:
                retStr += str(round(d, 4))

            if(transformed): retStr += '{'+ str(round(diffVActual[i], 4)) +'}'

            retStr += '  '

        unitT = ''
        unitA = ''
        if(transformed): unitA = unit
        else: unitT = unit

        retStr += '(' + str(round(diffV.length, 4)) + unitT +')'
        if(transformed): retStr += '{'+ str(round(diffVActual.length, 4))+ unitA +'}'

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

    def get3dLocSnap(self, rmInfo, vec = None, refreshStatus = True, \
        snapToAxisLine = True, xyDelta = [0, 0]):

        self.rmInfo = rmInfo
        obj = bpy.context.object
        xy = [rmInfo.xy[0] - xyDelta[0], rmInfo.xy[1] - xyDelta[1]]

        region = rmInfo.region
        rv3d = rmInfo.rv3d

        refLine = self.getRefLine()
        refLineOrig = self.getRefLineOrig()

        inEdit = self.isEditing()
        hasSel = self.hasSelection()

        params = bpy.context.window_manager.bezierToolkitParams
        transType = params.snapOrient
        origType = params.snapOrigin

        loc = None

        if(self.tm != None and hasSel and transType != 'REFERENCE'):
            tm, invTm = self.tm, self.tm.inverted()
        else:
            tm, invTm = self.getTransMatsForOrient(rmInfo, obj)

        if(self.orig != None and hasSel and origType != 'REFERENCE'):
            orig = self.orig
        else:
            orig = self.getCurrOrig(obj, rmInfo)
        
        if(vec == None): vec = orig

        self.lastSnapTypes = set()

        unit = unitMap.get(bpy.context.scene.unit_settings.length_unit)
        if(unit == None): unit = ''

        digitsValid = True
        freeAxesC = self.getFreeAxesCombined()
        freeAxesN = self.getFreeAxesNormalized()
        freeAxesG = self.getFreeAxesGlobal()

        if(self.objSnap):
            #TODO: Called very frequently (store the tree [without duplicating data])
            snapLocs = self.getAllSnapLocs((snapToAxisLine and \
                'AXIS' in {transType, origType})) + [orig]

            kd = kdtree.KDTree(len(snapLocs))
            for i, l in enumerate(snapLocs):
                kd.insert(getCoordFromLoc(region, rv3d, l).to_3d(), i)
            kd.balance()

            coFind = Vector(xy).to_3d()
            searchResult = kd.find_range(coFind, SNAP_DIST_PIXEL)

            if(len(searchResult) != 0):
                co, idx, dist = min(searchResult, key = lambda x: x[2])
                loc = snapLocs[idx]
                self.lastSnapTypes.add('loc')
            else:
                selObj, loc, normal, faceIdx = \
                    getSelFaceLoc(region, rv3d, xy, self.MAX_SNAP_FACE_CNT)

        if(loc != None):
            loc = tm @ loc
        else:
            loc = region_2d_to_location_3d(region, rv3d, xy, vec)
            loc = tm @ loc

            params = bpy.context.window_manager.bezierToolkitParams
            if(showSnapToPlane(params)):
                snapToPlane = params.snapToPlane
            else: snapToPlane = False

            # TODO: Get gridSnap and angleSnap out of this if
            if((transType != 'GLOBAL' and inEdit) or \
                snapToPlane or self.gridSnap or \
                    self.snapDigits.hasVal() or \
                        (inEdit and (len(freeAxesN) < 3 or self.angleSnap))):

                # snapToPlane means global constrain axes selection is a plane
                if(snapToPlane or refLineOrig == None): refCo = orig
                else: refCo = refLineOrig

                refCo = tm @ refCo

                if(self.snapDigits.hasVal()):
                    delta = self.snapDigits.getCurrDelta()
                    loc = tm @ orig + delta
                    self.lastSnapTypes.add('keyboard')
                else:
                    # Special condition for lock to single axis
                    if(len(freeAxesC) == 1 and refLineOrig != None):
                        refCo = tm @ refLineOrig
                    if(len(freeAxesC) == 2 or (len(freeAxesG) == 2 and snapToPlane)):
                        constrAxes = freeAxesG if (len(freeAxesG) == 2) else freeAxesC
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

                    if(len(freeAxesC) == 1 and len(refLine) > 0):

                        axis = freeAxesC[0]
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

                if(not self.snapDigits.hasVal() and self.gridSnap):
                    if(params.axisScale and (len(refLine) > 0 and transType == 'AXIS') \
                        or (refLineOrig and transType == 'REFERENCE')):
                        # Independent of view distance
                        diffV = (loc - tm @ refLineOrig)
                        loc = tm @ refLineOrig + \
                            round(diffV.length) * (diffV / diffV.length)
                    else:
                        rounding = getViewDistRounding(rv3d)
                        loc = tm @ roundedVect(invTm @ loc, rounding, freeAxesN)
                    self.lastSnapTypes.add('grid')

                if(not self.snapDigits.hasVal() and self.angleSnap and len(refLine) > 0):
                    freeAxesC = [0, 1, 2] if len(freeAxesC) == 0 else freeAxesC
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
                        if(i != axis and (i in freeAxesC)):
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
        self.snapCo = loc
        self.tm = tm
        self.orig = orig

        return loc

    def getEditCoPair(self):
        refLineOrig = self.getRefLineOrig()
        if(self.snapCo == None or refLineOrig == None):
            return []
        orig = self.getCurrOrig(bpy.context.object, self.rmInfo)
        return (self.tm @ orig, self.tm @ self.snapCo)

    def setStatus(self, area, text): #TODO Global
        area.header_text_set(text)

    def getGuideBatches(self, shader):
        rmInfo = self.rmInfo
        if(rmInfo == None): return []
        obj = bpy.context.object
        refLine = self.getRefLine()
        refLineOrig = self.getRefLineOrig()
        batches = []
        freeAxesC = self.getFreeAxesCombined()
        freeAxesN = self.getFreeAxesNormalized()

        if(self.tm != None):
            tm, invTm = self.tm, self.tm.inverted()
        else:
            tm, invTm = self.getTransMatsForOrient(rmInfo, obj)

        if(self.orig != None):
            orig = self.orig
        else:
            orig = self.getCurrOrig(obj, rmInfo)

        params = bpy.context.window_manager.bezierToolkitParams
        transType = params.snapOrient
        snapOrigin = params.snapOrigin

        lineCo = []
        if((refLineOrig != None or transType == 'VIEW' or len(freeAxesC) == 1)
            or (len(freeAxesN) > 0 and snapOrigin != 'REFERENCE')):
            colors = [(.6, 0.2, 0.2, 1), (0.2, .6, 0.2, 1), (0.2, 0.4, .6, 1)]
            l = 10 * rmInfo.rv3d.view_distance

            if (self.snapCo != None and len(freeAxesC) == 1): orig = self.snapCo

            refCo = tm @ orig

            for axis in freeAxesN[:2]:
                col = colors[axis]
                pt1 = refCo.copy()
                pt2 = refCo.copy()
                pt1[axis] = l + refCo[axis]
                pt2[axis] = -l + refCo[axis]
                slineCo, slineCol = getLineShades([invTm @ pt1, invTm @ pt2], col, .2, .9)
                batches.append(batch_for_shader(shader, \
                    "LINES", {"pos": slineCo, "color": slineCol}))

        if(refLineOrig != None and self.snapCo != None and \
            (self.angleSnap or ('keyboard' in self.lastSnapTypes \
                and self.snapDigits.polar))):

            slineCo = [orig, self.snapCo]
            slineCol = [(.4, .4, .4, 1)] * 2
            batches.append(batch_for_shader(shader, \
                "LINES", {"pos": slineCo, "color": slineCol}))
            batches.append(batch_for_shader(shader, \
                "POINTS", {"pos": lineCo, \
                    "color": [(1, 1, 1, 1) for i in range(0, len(lineCo))]}))

        if(self.customAxis.length() != 0 and \
            (self.inDrawAxis == True or transType == 'AXIS' or snapOrigin == 'AXIS')):
            apts = self.customAxis.axisPts
            lineCo = [apts[0], apts[1]]
            ptCo = self.customAxis.getSnapPts()
        else: ptCo = []

        # ~ if(self.snapCo != None and 'loc' in self.lastSnapTypes):
            # ~ ptCo.append(self.snapCo)

        # Axis Line
        slineCo, slineCol = getLineShades(lineCo, (1, 1, 1, 1), .9, .3, mid = False)
        batches.append(batch_for_shader(shader, \
            "LINES", {"pos": slineCo, "color": slineCol}))

        batches.append(batch_for_shader(shader, \
            "POINTS", {"pos": ptCo, \
                "color": [(1, 1, 1, 1) for i in range(0, len(ptCo))]}))

        return batches

def getLineShades(lineCos, baseColor, start, end, mid = True):
    if(len(lineCos) == 0 ): return [], []
    if(len(lineCos) == 1 ): return lineCos[0], [baseColor]
    if(mid): midPt = lineCos[0] + (lineCos[1] - lineCos[0]) / 2
    col1 = [start * c for c in baseColor]
    col2 = [end * c for c in baseColor]
    if(mid): return [lineCos[0], midPt, midPt, lineCos[1]], [col1, col2, col2, col1]
    else: return [lineCos[0], lineCos[1]], [col1, col2]


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
    drawHandlerRef = None
    drawFunc = None
    segBatch = None
    tipBatch = None
    shader = None
    snapperBatches = []
    opObj = None

    pointSize = 4 # For Draw (Marker is of diff size)

    def addDrawHandler(drawHandler):
        ModalBaseFlexiOp.drawHandlerRef = \
            bpy.types.SpaceView3D.draw_handler_add(drawHandler, (), "WINDOW", "POST_VIEW")

    def removeDrawHandler():
        if(ModalBaseFlexiOp.drawHandlerRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(ModalBaseFlexiOp.drawHandlerRef, \
                "WINDOW")
            ModalBaseFlexiOp.drawHandlerRef = None

    def drawHandlerBase():
        if(ModalBaseFlexiOp.shader != None):

            bgl.glLineWidth(ModalBaseFlexiOp.axisLineWidth)
            bgl.glPointSize(ModalBaseFlexiOp.snapPtSize)
            for batch in ModalBaseFlexiOp.snapperBatches:
                batch.draw(ModalBaseFlexiOp.shader)

            bgl.glLineWidth(ModalBaseFlexiOp.lineWidth)
            if(ModalBaseFlexiOp.segBatch != None):
                ModalBaseFlexiOp.segBatch.draw(ModalBaseFlexiOp.shader)

            bgl.glPointSize(ModalBaseFlexiOp.drawPtSize)
            if(ModalDrawBezierOp.tipBatch != None):
                ModalDrawBezierOp.tipBatch.draw(ModalBaseFlexiOp.shader)


    def resetDisplayBase():
        ModalBaseFlexiOp.segBatch = getResetBatch(ModalBaseFlexiOp.shader, "LINES")
        ModalBaseFlexiOp.tipBatch = getResetBatch(ModalBaseFlexiOp.shader, "POINTS")
        ModalBaseFlexiOp.snapperBatches = \
            [getResetBatch(ModalBaseFlexiOp.shader, "LINES"), \
                getResetBatch(ModalBaseFlexiOp.shader, "POINTS")]

        areas = [a for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
        for a in areas:
            a.tag_redraw()

    def refreshDisplayBase(displayInfos, snapper = None, ptDispInfos = None):
        areaRegionInfo = getAllAreaRegions()

        ModalBaseFlexiOp.segBatch, ModalBaseFlexiOp.tipBatch = \
            getBezierBatches(ModalDrawBezierOp.shader, displayInfos, \
                areaRegionInfo, ptDispInfos)

        if(snapper != None):
            ModalBaseFlexiOp.snapperBatches = \
                snapper.getGuideBatches(ModalBaseFlexiOp.shader)
        else:
            ModalBaseFlexiOp.snapperBatches = []

        areas = [a for a in bpy.context.screen.areas if a.type == 'VIEW_3D']
        for a in areas:
            a.tag_redraw()

    @persistent
    def loadPostHandler(dummy):
        if(ModalBaseFlexiOp.shader != None):
            ModalBaseFlexiOp.resetDisplayBase()
        ModalBaseFlexiOp.running = False

    @persistent
    def loadPreHandler(dummy):
        ModalBaseFlexiOp.removeDrawHandler()
        if(ModalBaseFlexiOp.drawFunc != None):
            bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc

    @classmethod
    def poll(cls, context):
        return not ModalBaseFlexiOp.running

    def preInvoke(self, context, event):
        pass # place holder

    def subInvoke(self, context, event):
        return {'RUNNING_MODAL'} # place holder

    def invoke(self, context, event):
        ModalBaseFlexiOp.opObj = self
        ModalBaseFlexiOp.running = True
        self.preInvoke(context, event)
        ModalBaseFlexiOp.addDrawHandler(self.__class__.drawHandler)
        ModalBaseFlexiOp.drawFunc = bpy.types.VIEW3D_HT_tool_header.draw
        bpy.types.VIEW3D_HT_tool_header.draw = drawSettingsFT
        context.space_data.show_region_tool_header = True

        self.snapper = Snapper(context, self.getSnapLocs, \
            self.getRefLine, self.getRefLineOrig, self.hasSelection, self.isEditing)

        self.rmInfo = None

        ModalBaseFlexiOp.shader = gpu.shader.from_builtin('3D_SMOOTH_COLOR')
        ModalBaseFlexiOp.shader.bind()
        context.window_manager.modal_handler_add(self)

        ModalBaseFlexiOp.ColGreaseHltSeg = (.3, .3, .3, 1) # Not used

        updateProps(None, context)

        return self.subInvoke(context, event)

    def modal(self, context, event):
        if(context.space_data == None):
            return {'PASS_THROUGH'}

        if(not is3DVireport(context)):
            self.cancelOp(context)
            resetToolbarTool()
            return {'CANCELLED'}

        if(event.type == 'WINDOW_DEACTIVATE' and event.value == 'PRESS'):
            self.resetMetaBtns() # Subclass
            self.snapper.initialize() # TODO: Check if needed
            return {'PASS_THROUGH'}

        if(not self.isToolSelected(context)): # Subclass
            self.cancelOp(context)
            return {"CANCELLED"}

        rmInfo = RegionMouseXYInfo.getRegionMouseXYInfo(event, self.exclToolRegion())
        if(self.isEditing() and self.rmInfo != rmInfo):
            return {'RUNNING_MODAL'}
        elif(rmInfo == None and not self.isEditing()):
            return {'PASS_THROUGH'}

        self.rmInfo = rmInfo
        snapProc = self.snapper.procEvent(context, event, rmInfo)

        if(not snapProc and not self.hasSelection() and event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                self.cancelOp(context)
                resetToolbarTool()
                return {"CANCELLED"}
            return {'RUNNING_MODAL'}

        return self.subModal(context, event, snapProc)

    def cancelOpBase(self):
        ModalBaseFlexiOp.removeDrawHandler()
        ModalBaseFlexiOp.running = False
        bpy.types.VIEW3D_HT_tool_header.draw = ModalBaseFlexiOp.drawFunc
        self.snapper = None

    def getSnapLocs(self):
        return self.getSnapLocsImpl()

################### Flexi Draw Bezier Curve ###################

class ModalDrawBezierOp(ModalBaseFlexiOp):

    # Static members shared by flexi draw and flexi grease
    subdivPtBatch = None
    subdivLineBatch = None

    markerBatch = None
    markerSize = 8

    #static method
    def drawHandler():
        if(ModalBaseFlexiOp.shader != None):
            bgl.glPointSize(ModalDrawBezierOp.markerSize)
            if(ModalDrawBezierOp.markerBatch != None):
                ModalDrawBezierOp.markerBatch.draw(ModalBaseFlexiOp.shader)

            # TODO: Move this to grease draw
            bgl.glPointSize(ModalBaseFlexiOp.greaseSubdivPtSize)
            if(ModalDrawBezierOp.subdivPtBatch != None):
                ModalDrawBezierOp.subdivPtBatch.draw(ModalBaseFlexiOp.shader)

            if(ModalDrawBezierOp.subdivLineBatch != None):
                ModalDrawBezierOp.subdivLineBatch.draw(ModalBaseFlexiOp.shader)

            ModalBaseFlexiOp.drawHandlerBase()

    def resetDisplay():
        ModalDrawBezierOp.subdivPtBatch = \
            getResetBatch(ModalBaseFlexiOp.shader, "POINTS")
        ModalDrawBezierOp.subdivLineBatch = \
            getResetBatch(ModalBaseFlexiOp.shader, "LINES")
        ModalDrawBezierOp.markerBatch = \
            getResetBatch(ModalBaseFlexiOp.shader, "POINTS")

        ModalBaseFlexiOp.resetDisplayBase()

    def refreshDisplay(displayInfos, subdivCos = [], \
        showSubdivPts = True, markerLoc = [], colMarker = None, snapper = None, \
            ptDispInfos = None):

        ModalDrawBezierOp.subdivPtBatch, ModalDrawBezierOp.subdivLineBatch = \
            getSubdivBatches(ModalBaseFlexiOp.shader, subdivCos, showSubdivPts)

        ModalDrawBezierOp.markerBatch = batch_for_shader(ModalBaseFlexiOp.shader, \
            "POINTS", {"pos": markerLoc, \
                "color": [colMarker for i in range(len(markerLoc))]})

        ModalBaseFlexiOp.refreshDisplayBase(displayInfos, snapper, ptDispInfos)

    def __init__(self, curveDispRes):
        pass

    #This will be called multiple times not just at the beginning
    def initialize(self):
        self.curvePts = []
        self.clickT = None #For double click
        self.pressT = None #For single click
        self.capture = False
        self.grabRepos = False
        self.resetMetaBtns()
        self.snapper.initialize()

    def subInvoke(self, context, event):
        self.initialize()

        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)

        try:
            ModalDrawBezierOp.markerSize = \
                context.preferences.addons[__name__].preferences.markerSize
        except Exception as e:
            # ~ print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
            ModalDrawBezierOp.markerSize = 8

        return {"RUNNING_MODAL"}

    def resetMetaBtns(self):
        self.ctrl = False
        self.shift = False
        self.alt = False

    def cancelOp(self, context):
        ModalDrawBezierOp.resetDisplay()
        bpy.app.handlers.undo_post.remove(self.postUndoRedo)
        bpy.app.handlers.redo_post.remove(self.postUndoRedo)
        return self.cancelOpBase()

    def postUndoRedo(self, scene, dummy = None): # signature different in 2.8 and 2.81?
        self.updateSnapLocs() # subclass method

    def confirm(self, context, event):
        self.save(context, event)
        self.curvePts = []
        self.capture = False
        ModalDrawBezierOp.resetDisplay()
        self.initialize()

    def newPoint(self, loc):
        self.curvePts.append([loc, loc, loc])

    def moveBezierPt(self, loc):
        if(len(self.curvePts) > 0): self.curvePts[-1] = [loc, loc, loc]

    def movePointByDelta(self, delta):
        if(len(self.curvePts) > 0):
            self.curvePts[-1][1] = self.curvePts[-1][1]  + delta

    def moveBptElem(self, handle, loc):
        idx = {'left':0, 'pt':1, 'right':2}[handle]
        if(len(self.curvePts) > 0):
            self.curvePts[-1][idx] = loc

    def moveBptElemByDelta(self, handle, delta):
        idx = {'left':0, 'pt':1, 'right':2}[handle]
        if(len(self.curvePts) > 0):
            self.curvePts[-1][idx] = self.curvePts[-1][idx] + delta

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

    def exclToolRegion(self):
        return True

    def isEditing(self):
        return len(self.curvePts) > 0

    def hasSelection(self):
        return self.isEditing()

    # Common subModal for Flexi Draw and Flexi Grease
    def baseSubModal(self, context, event, snapProc):
        rmInfo = self.rmInfo

        if(event.type == 'G'):
            if(event.value == 'RELEASE'):
                self.grabRepos = not self.grabRepos
            return {"RUNNING_MODAL"}

        # This can happen only when space was entered and something was there
        # for Snapper to process
        if (snapProc and self.snapper.digitsConfirmed):
            self.snapper.resetSnap()

            # Because resetSnap sets this to False (TODO: Refactor resetSnap)
            self.snapper.digitsConfirmed = True

            # First space / enter is equivalent to mouse press without release
            if(not self.capture):
                self.capture = True
                self.snapper.setStatus(rmInfo.area, None)
                return {'RUNNING_MODAL'}
            else:
                # Second space / enter means it should be processed here,
                # set snapProc to False so this modal will process it
                snapProc = False

        if(updateMetaBtns(self, event)):
            return {'RUNNING_MODAL'}

        if(not snapProc and event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                if(self.capture and self.isHandleSet()):
                    self.resetHandle('left')
                    self.resetHandle('right')

                    # Needed to indicate next space / entered to be processed here
                    self.snapper.digitsConfirmed = True
                    self.snapper.setStatus(rmInfo.area, None)
                else:
                    self.initialize()
                self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and event.type == 'BACK_SPACE'):
            if(event.value == 'RELEASE'):
                self.snapper.resetSnap()
                if(not self.capture):
                    if(len(self.curvePts) > 0):
                        self.curvePts.pop()

                    #Because there is an extra point (the current one)
                    if(len(self.curvePts) <= 1):
                        self.curvePts = []
                        self.capture = False
                        self.grabRepos = False
                    else:
                        loc = self.snapper.get3dLocSnap(rmInfo)
                        self.curvePts[-1] = [loc, loc, loc]
                self.capture = False
                self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and (event.type == 'RET' or event.type == 'SPACE')):
            if(event.value == 'RELEASE'):
                if(self.snapper.digitsConfirmed):
                    self.capture = False
                    self.snapper.digitsConfirmed = False
                    loc = self.snapper.get3dLocSnap(rmInfo)
                    self.newPoint(loc)
                    self.redrawBezier(rmInfo)
                else:
                    self.confirm(context, event)
                    self.snapper.resetSnap()
                    self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(not snapProc and event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            if(len(self.curvePts) == 0):
                loc = self.snapper.get3dLocSnap(rmInfo)
                self.newPoint(loc)

            # Special condition for hot-key single axis lock (useful)
            if(len(self.snapper.freeAxes) == 1 and len(self.curvePts) > 1):
                self.snapper.resetSnap()

            if(not self.capture):
                # ~ self.snapper.resetSnap()
                self.pressT = time.time()
                self.capture = True
            return {'RUNNING_MODAL'}

        if (not snapProc and event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):
            self.capture = False
            self.grabRepos = False

            # See the special condition above in event.value == 'PRESS'
            if(len(self.snapper.freeAxes) > 1):
                self.snapper.resetSnap()

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
                self.redrawBezier(rmInfo)
                self.clickT = None
                return {'RUNNING_MODAL'}

            self.clickT = t

            co = None
            if((self.pressT != None) and (t - self.pressT) < 0.2):
                loc = self.curvePts[-1][1]
                self.moveBptElem('left', loc)
                self.moveBptElem('right', loc)
            else:
                loc = self.snapper.get3dLocSnap(rmInfo)

            if(len(self.curvePts) == 1):
                self.moveBptElem('right', loc)# changes only rt handle

            self.newPoint(loc)
            self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        # Refresh also in case of snapper events
        # except when digitsConfirmed (to give user opportunity to draw a straight line)
        # ~ if ((snapProc and not self.snapper.digitsConfirmed) \
        if (snapProc\
            or (event.type == 'MOUSEMOVE')):

            # Unlock axes in case of pure mousemove (TODO: Can be better)
            # ~ if(snapProc and not self.snapper.digitsConfirmed): pass
            # ~ else: self.snapper.resetSnap()
            # ~ if(self.snapper.digitsConfirmed): self.snapper.resetSnap()

            bpy.context.window.cursor_set("DEFAULT")
            if(len(self.curvePts) > 1):
                if(self.capture):
                    if(self.grabRepos):
                        pt = self.curvePts[-1][1].copy()
                        rtHandle = self.curvePts[-1][2].copy()
                        xy2 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, pt)
                        xy1 = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, rtHandle)
                        loc = self.snapper.get3dLocSnap(rmInfo, \
                            xyDelta = [xy1[0] - xy2[0], xy1[1] - xy2[1]])
                        delta = loc.copy() - pt.copy()
                        self.moveBptElemByDelta('pt', delta)
                        self.moveBptElemByDelta('left', delta)
                        self.moveBptElemByDelta('right', delta)
                    else:
                        loc = self.snapper.get3dLocSnap(rmInfo)
                        pt = self.curvePts[-1][1]
                        delta = (loc - pt)
                        self.moveBptElem('left', pt - delta)
                        self.moveBptElem('right', pt + delta)
                else:
                    loc = self.snapper.get3dLocSnap(rmInfo)
                    self.moveBezierPt(loc)
            self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        return {'PASS_THROUGH'} if not snapProc else {'RUNNING_MODAL'}

    def redrawBezier(self, rmInfo, subdivCos = [], segIdxRange = None, \
        showSubdivPts = True):

        colMap = self.getColorMap()
        colSelSeg = colMap['SEL_SEG_COLOR']
        colNonAdjSeg = colMap['NONADJ_SEG_COLOR']
        colTip = colMap['TIP_COLOR']
        colEndTip = colMap['ENDPT_TIP_COLOR']
        colMarker = colMap['MARKER_COLOR']

        markerLoc = []
        handleNos  = [0, 1]

        if(self.capture): hdlPtIdx = 1
        else: hdlPtIdx = 0

        curvePts = self.curvePts

        if(not self.capture or len(self.curvePts) <= 1):
            loc = self.snapper.get3dLocSnap(rmInfo)
            if(self.capture): #curvePts len must be 1
                # First handle (straight line), if user drags first pt
                curvePts = curvePts + [[curvePts[0][0], curvePts[0][1], loc]]
                handleNos = [1] # display only right handle
                
            # Marker (dot), if drawing not started
            elif(len(curvePts) == 0): markerLoc = [loc]

        tipColors = [colTip, colEndTip, colTip]
        displayInfos = []
        ptDisplayInfos = []

        ptCnt = len(curvePts)
        idxRange = range(1, ptCnt) if segIdxRange == None else segIdxRange

        for i in idxRange:
            segPts = [curvePts[i-1], curvePts[i]]            
            if(i == ptCnt - 1):
                segColor = colSelSeg
                ptDisplayInfos.append(PtDisplayInfo(segPts[hdlPtIdx], \
                    handleNos, tipColors))
            else: 
                segColor = colNonAdjSeg
                
            displayInfos.append(SegDisplayInfo(segPts, segColor))

        ModalDrawBezierOp.refreshDisplay(displayInfos, subdivCos, showSubdivPts, \
            markerLoc, colMarker, self.snapper, ptDisplayInfos)

    #Reference point for restrict angle or lock axis
    def getRefLine(self):
        if(len(self.curvePts) > 0):
            idx = 0
            if(self.capture):
                idx = -1
                # ~ return [self.curvePts[-1][1]]
            # There should always be min 2 pts if not capture, check anyway
            elif(len(self.curvePts) > 1):
                idx = -2
                # ~ return [self.curvePts[-2][1]]
            if((len(self.curvePts) + (idx - 1)) >= 0):
                return[self.curvePts[idx-1][1], self.curvePts[idx][1]]
            else:
                return[self.curvePts[idx][1]]
        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[-1] if len(refLine) > 0 else None

class ModalFlexiDrawBezierOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves in object mode"
    bl_idname = "wm.flexi_draw_bezier_curves"
    bl_label = "Flexi Draw Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}

    def __init__(self):
        # ~ curveDispRes = 200
        # ~ super(ModalFlexiDrawBezierOp, self).__init__(curveDispRes)
        pass

    def isToolSelected(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)

        # ~ if(tool == None or tool.idname != FlexiDrawBezierTool.bl_idname): (T60766)
        if(tool == None or tool.idname != 'flexi_bezier.draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {'SEL_SEG_COLOR': ModalBaseFlexiOp.colDrawSelSeg,
        'NONADJ_SEG_COLOR': ModalBaseFlexiOp.colDrawNonHltSeg,
        'TIP_COLOR': ModalBaseFlexiOp.colHdlPtTip,
        'ENDPT_TIP_COLOR': ModalBaseFlexiOp.colBezPt,
        'MARKER_COLOR': ModalBaseFlexiOp.colDrawMarker}

    def preInvoke(self, context, event):

        # If the operator is invoked from context menu, enable the tool on toolbar
        if(not self.isToolSelected(context) and context.mode == 'OBJECT'):
            # ~ bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname) (T60766)
            bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')

        # Object name -> [spline index, [pts]]
        # Not used right now (maybe in case of large no of curves)
        self.snapInfos = {}
        self.updateSnapLocs()

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo
        if(event.type == 'E'):
            if(event.value == 'RELEASE'):
                # ~ bpy.ops.wm.tool_set_by_id(name = FlexiEditBezierTool.bl_idname) (T60766)
                bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.edit_tool')
            return {"RUNNING_MODAL"}

        return self.baseSubModal(context, event, snapProc)

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if(len(self.curvePts) > 0):
            locs += [pt[1] for pt in self.curvePts[:-1]]

        return locs

    def updateSnapLocs(self, objNames = None):
        updateCurveEndPtMap(self.snapInfos, addObjNames = objNames)

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

    def createCurveObj(self, context, startObj = None, \
        startSplineIdx = None, endObj = None, endSplineIdx = None):
        # First create the new curve
        obj = self.createObjFromPts(context)

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
                    if(s.use_cyclic_u): continue
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

    def save(self, context, event):
        if(len(self.curvePts) > 0):
            self.curvePts.pop()

        if(len(self.curvePts) > 1):

            startObj, startSplineIdx, ptIdx2, endObj, endSplineIdx, ptIdx1 = \
                [x for y in self.getSnapObjs(context, [self.curvePts[0][1],
                    self.curvePts[-1][1]]) for x in y]

            # ctrl pressed and there IS a snapped end obj,
            # so user does not want connection

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
                # ~ bpy.context.view_layer.objects.active = obj
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

################### Flexi Draw Grease Bezier ###################

class ModalFlexiDrawGreaseOp(ModalDrawBezierOp):
    bl_description = "Flexible drawing of Bezier curves as grease pencil strokes"
    bl_idname = "wm.flexi_draw_grease_bezier_curves"
    bl_label = "Flexi Draw Grease Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}
    h = False

    def __init__(self):
        # ~ curveDispRes = 200
        # ~ super(ModalFlexiDrawGreaseOp, self).__init__(curveDispRes)
        pass

    def isToolSelected(self, context):
        if(context.mode != 'PAINT_GPENCIL'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('PAINT_GPENCIL', \
            create = False)

        # ~ if(tool == None or tool.idname != FlexiDrawBezierTool.bl_idname): (T60766)
        if(tool == None or tool.idname != 'flexi_bezier.grease_draw_tool'):
            return False

        return True

    def getColorMap(self):
        return {'SEL_SEG_COLOR': ModalBaseFlexiOp.colGreaseSelSeg,
        'NONADJ_SEG_COLOR': ModalBaseFlexiOp.ColGreaseHltSeg, #Not used
        'TIP_COLOR': ModalBaseFlexiOp.colHdlPtTip,
        'ENDPT_TIP_COLOR': ModalBaseFlexiOp.colGreaseBezPt,
        'MARKER_COLOR': ModalBaseFlexiOp.colGreaseMarker, }

    def preInvoke(self, context, event):
        # If the operator is invoked from context menu, enable the tool on toolbar
        if(not self.isToolSelected(context) and context.mode == 'PAINT_GPENCIL'):
            # ~ bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname) (T60766)
            bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.grease_draw_tool')

        o = context.object
        if(o == None or o.type != 'GPENCIL'):
            d = bpy.data.grease_pencils.new('Grease Pencil Data')
            o = bpy.data.objects.new('Grease Pencil', d)
            context.scene.collection.objects.link(o)
        self.gpencil = o

        self.subdivCos = []
        self.interpPts = []

    # overridden
    def redrawBezier(self, rmInfo, subdivCos = []):
        ptCnt = len(self.curvePts)
        segIdxRange = range(ptCnt - 1, ptCnt) if(ptCnt > 1) else None

        super(ModalFlexiDrawGreaseOp, self).redrawBezier(rmInfo, \
            self.subdivCos if ptCnt > 1 else [], \
                segIdxRange, showSubdivPts = not ModalFlexiDrawGreaseOp.h)

    def initialize(self):
        super(ModalFlexiDrawGreaseOp, self).initialize()
        self.subdivCos = []
        self.interpPts = []
        self.updateSnapLocs()

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo

        if(event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', 'NUMPAD_PLUS', \
            'NUMPAD_MINUS','PLUS', 'MINUS'} and len(self.curvePts) > 1):
            if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                and event.value == 'PRESS'):
                return {'RUNNING_MODAL'}
            elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                self.subdivAdd(5)
            elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                self.subdivAdd(-5)

            self.redrawBezier(rmInfo)
            return {'RUNNING_MODAL'}

        if(event.type == 'H' or event.type == 'h'):
            if(event.value == 'RELEASE'):
                ModalFlexiDrawGreaseOp.h = not ModalFlexiDrawGreaseOp.h
                self.redrawBezier(rmInfo)
            return {"RUNNING_MODAL"}

        ptCnt = len(self.curvePts)

        retVal = self.baseSubModal(context, event, snapProc)

        newPtCnt = len(self.curvePts)
        if(newPtCnt - ptCnt != 0):
            if(newPtCnt == 1):
                viewDist = context.space_data.region_3d.view_distance
                self.initSubdivPerUnit = 5000.0 / viewDist # TODO: default configurable?
                self.subdivPerUnit = 0.02 * self.initSubdivPerUnit
                self.snapLocs.append(self.curvePts[0][1])
            else:
                slens = self.getCurveSegLens()
                self.updateInterpPts(slens)
                self.updateSubdivCos(sum(slens))
                self.redrawBezier(rmInfo)

        return retVal

    def getCurveSegLens(self):
        clen = []
        for i in range(1, len(self.curvePts) - 1):
            clen.append(getSegLen([self.curvePts[i-1][1], self.curvePts[i-1][2], \
                self.curvePts[i][0], self.curvePts[i][1]]))
        return clen

    def updateSubdivCos(self, clen = None):
        if(self.interpPts != []):
            if(clen == None): clen = sum(self.getCurveSegLens())
            cnt = round(self.subdivPerUnit * clen)
            if(cnt > 0):
                self.subdivCos = getInterpolatedVertsCo(self.interpPts, cnt)#[1:-1]
                return
        self.subdivCos = []

    def updateInterpPts(self, slens):
        self.interpPts = getInterpBezierPts(self.curvePts[:-1],
            self.initSubdivPerUnit, slens)

        return self.interpPts

    def subdivAdd(self, addCnt):
        slens = self.getCurveSegLens()
        clen = sum(slens)
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

    def save(self, context, event):
        layer = self.gpencil.data.layers.active
        if(layer == None):
            layer = self.gpencil.data.layers.new('GP_Layer', set_active = True)
        if(len(layer.frames) == 0):
            layer.frames.new(0)
        frame = layer.frames[-1]

        invMw = self.gpencil.matrix_world.inverted()
        if(len(self.subdivCos) > 0):
            self.curvePts.pop()
            stroke = frame.strokes.new()
            stroke.display_mode = '3DSPACE'
            stroke.points.add(count = len(self.subdivCos))
            for i in range(0, len(self.subdivCos)):
                pt = self.subdivCos[i]
                stroke.points[i].co = self.gpencil.matrix_world.inverted() @ pt
            self.snapLocs += [self.subdivCos[0][1], self.subdivCos[-1][1]]
        bpy.ops.ed.undo_push()

################### Flexi Edit Bezier Curve ###################

class EditSegDisplayInfo(SegDisplayInfo):

    def __init__(self, segPts, segColor, handleNos, tipColors, subdivCos):
        super(EditSegDisplayInfo, self).__init__(segPts, segColor, handleNos, tipColors)
        self.subdivCos = subdivCos

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

# segIdx -> startPtIdx
def getCtrlPtsForSeg(obj, splineIdx, segIdx):
    wsData = getWSData(obj)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)

    if(segEndIdx == None):
        return None

    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0[1], pt0[2], pt1[0], pt1[1]]

def getBezierDataForSeg(obj, splineIdx, segIdx):
    wsData = getWSData(obj)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)
    if(segEndIdx == None):
        return []
    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0, pt1]

# Requirement is more generic than geometry.interpolate_bezier
# TODO: revisit (maybe merge with getPtsAlongBezier2D)
def getInterpSegPts(mw, spline, ptIdx, res, startT, endT, maxRes):
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

    resProp = int(res * getSegLen(seg))
    if(resProp > 1):
        # Otherwise too slow, when zoom level very high
        if(resProp > maxRes): resProp = maxRes
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

# Wrapper for spatial search within segment
def getClosestPt2dWithinSeg(region, rv3d, coFind, selObj, selSplineIdx, selSegIdx, \
    selObjRes, withHandles, withBezPts):
    infos = {selObj: {selSplineIdx:[selSegIdx]}}
    # set selObj in objs for CurveBezPts
    return getClosestPt2d(region, rv3d, coFind, [selObj], 0, \
        infos, selObjRes, withHandles, withBezPts, withObjs = False)

def getClosestPt2d(region, rv3d, coFind, objs, objRes, selObjInfos, selObjRes, \
    withHandles = True, withBezPts = True, withObjs = True, normalized = True, \
        objStartT = 0, objEndT = 1, selObjStartT = 0, selObjEndT = 1):

    objLocMap = {}

    if(normalized):
        #TODO: Should be pixel based (2d)
        viewDist = rv3d.view_distance
        objRes /= viewDist # inversely proportional
        selObjRes /= viewDist

    objInterpCounts = []
    objLocList = []

    objInterpLocs = []
    objSplineEndPts = []

    for obj in objs:
        mw = obj.matrix_world
        if(not isPtIn2dBBox(obj, region, rv3d, coFind, SNAP_DIST_PIXEL)):
            continue
            
        for i, spline in enumerate(obj.data.splines):
            for j, pt in enumerate(spline.bezier_points):
                objLocList.append([obj, i, j])
                if(withObjs):
                    interpLocs = getInterpSegPts(mw, spline, j, objRes, \
                        objStartT, objEndT, maxRes = NONSEL_CURVE_SEARCH_RES)[1:-1]

                    objInterpLocs += interpLocs
                    objInterpCounts.append(len(interpLocs))

                if(withBezPts):
                    objSplineEndPts.append(mw @ pt.co)

    selObjInterpCounts = []
    selObjLocList = []

    segInterpLocs = []
    hdls = []

    for selObj in selObjInfos.keys():
        mw = selObj.matrix_world
        info = selObjInfos[selObj]
        for splineIdx in info.keys():
            spline = selObj.data.splines[splineIdx]
            segIdxs = info[splineIdx]
            for segIdx in segIdxs:
                interpLocs = getInterpSegPts(mw, spline, segIdx, selObjRes, \
                    selObjStartT, selObjEndT, maxRes = ADD_PT_CURVE_SEARCH_RES * 5)[1:-1]
                segInterpLocs += interpLocs
                selObjInterpCounts.append(len(interpLocs))
                selObjLocList.append([selObj, splineIdx, segIdx])

                if(withHandles):
                    bpts = getBezierDataForSeg(selObj, splineIdx, segIdx)
                    hdls += [pt for pts in bpts for pt in [pts[0], pts[2]]]

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

    srs = search2dFromPtsList(searchPtsList, coFind, searchRange = SNAP_DIST_PIXEL)

    if(len(srs) == 0):
        return None

    sr = min(srs, key = lambda x: x[3])

    if(sr[0] > 1):
        # If seg loc then first priority to the nearby handle, end pt (even if farther)
        sr = min(srs, key = lambda x: (x[0], x[3]))

    idx = sr[1]
    retId = retStr[sr[0]]

    if(sr[0] in {0, 2}):
        cntList = selObjInterpCounts if sr[0] == 2 else \
            [4 for i in range(0, len(selObjInterpCounts))]

        listIdx = findListIdx(cntList, idx)

        obj, splineIdx, segIdx = selObjLocList[listIdx]
        return  retId, obj, splineIdx, segIdx, idx % 4 if sr[0] == 0 \
            else segInterpLocs[idx]

    elif(sr[0] in {1, 3}):
        cntList = objInterpCounts if sr[0] == 3 else \
            [1 for i in range(0, len(objInterpCounts))]

        listIdx = findListIdx(cntList, idx)

        obj, splineIdx, segIdx = objLocList[listIdx]
        if(sr[0] == 3):
            return retId, obj, splineIdx, segIdx, objInterpLocs[idx]
        else:
            otherInfo = segIdx
            # ~ if(obj in selObjInfos.keys()):
                # ~ selInfo = selObjInfos[obj]
                # ~ if splineIdx in selInfo.keys():
                    # ~ if(segIdx in selInfo[splineIdx]): retId = 'SelSegEndPt'
                    # ~ else: retId = 'SelSplineEndPt'
            return retId, obj, splineIdx, segIdx, otherInfo

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
        ptIdxInList = idx - sum(len(ptsList[i]) for i in range(0, listIdx))
        searchResults.append([listIdx, ptIdxInList, co, dist])

    return searchResults

def getCtrlIdxFromSearchInfo(info):
    if(info[0] == 'SelHandles'):
        return {0:0, 1:2, 2:3, 3:5}[info[1]]
    return None


class SelectCurveInfo:
    def __init__(self, obj, splineIdx, segIdx1):
        self.obj = obj
        self.splineIdx = splineIdx
        self.ptIdxs = set()
        
        if(segIdx1 != None): self.addSegSel(segIdx1)

        # obj.name gives exception if obj is not in bpy.data.objects collection,
        # so keep a copy
        self.objName = obj.name
        self.subdivCnt = 0
        self.interpPts = None
        self._ctrlIdx = None # Handle and end points: 0, 1, 2 and 3, 4, 5
        self._clickLoc = None

        # Can be derived from clickLoc, but stored to avoid repeated computation
        self._t = None

    def addSegSel(self, segIdx1):        
        segIdx2 = self.getSegAdjIdx(segIdx1)
        if(segIdx2 != None): 
            self.ptIdxs.add(segIdx1)
            self.ptIdxs.add(segIdx2)

    def removeSegSel(self, segEndPtIdx):        
        adjIdx = self.getSegAdjIdx(segEndPtIdx)
        if(adjIdx != None):
            self.ptIdxs.remove(segEndPtIdx)
            self.ptIdxs.remove(adjIdx)

    def addPtSel(self, ptIdx):
        if(ptIdx != None): self.ptIdxs.add(ptIdx)

    def removePtSel(self, ptIdx):
        if(ptIdx != None): self.ptIdxs.remove(ptIdx)

    def resetPts(self):
        self.ptIdxs = set()

    def getSegIdxs(self):
        return [p[0] for p in self.getSegIdxPairs()]

    def getSegIdxPairs(self):
        segIdxs = []
        for idx in self.ptIdxs:
            nextIdx = getAdjIdx(self.obj, self.splineIdx, idx)
            if(nextIdx in self.ptIdxs):
                segIdxs.append([idx, nextIdx])
        return segIdxs

    def resetSel(self):
        self._clickLoc = None
        self._t = None
        self._ctrlIdx = None

    def getClickLoc(self):
        return self._clickLoc

    def setClickLocSafe(self, clickLoc, lowerT = 0.001, higherT = .999):
        self._clickLoc = clickLoc
        self._t = getTForPt(self.getCtrlPts(), clickLoc)
        if(self._t != None and (self._t < lowerT or self._t > higherT)):
            if(self._t < lowerT):
                self._ctrlIdx = 1
            else:
                self._ctrlIdx = 4
            self._t = None
            self._clickLoc = None

    def getCtrlPtCoIdx(self):
        if(self._ctrlIdx == None):
            return None

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
            self._t = None

    def getSelCo(self):
        if(self._ctrlIdx != None):
            return self.getCtrlPtCoIdx()[0]
        return self._clickLoc

    # Callback after subdivseg (in case seg from same object being subdivided)
    def updateSegIdx(self, objName, splineIdx, oldSegIdx, addCnt):
        if(objName == self.objName and splineIdx == self.splineIdx
            and oldSegIdx in self.getSegIdxs()):
                # ~ self.ptIdxs[0] += addCnt
                self.removeSegSel(oldSegIdx)
                self.addSegSel(oldSegIdx + addCnt)

    def subdivSeg(self):
        if(self.subdivCnt > 1):
            invMw = self.obj.matrix_world.inverted()
            ts = []
            vertCos = getInterpolatedVertsCo(self.interpPts, self.subdivCnt)[1:-1]
            insertBezierPts(self.obj, self.splineIdx, self.getSegIdxs()[0], \
                [invMw @ v for v in vertCos], 'FREE')
            self.subdivCnt = 0


    def subdivMode(self, rv3d):
        self.subdivCnt = 2
        self.interpPts = getPtsAlongBezier3D(self.getSegPts(), rv3d,
            curveRes = 1000, minRes = 1000)

    def subdivDecr(self):
        if(self.subdivCnt > 2):
            self.subdivCnt -= 1

    def subdivIncr(self):
        if(self.subdivCnt < 100):
            self.subdivCnt += 1

    def getLastSegIdx(self):
        spline = self.obj.data.splines[self.splineIdx]
        ptCnt = len(spline.bezier_points)
        # ~ if(ptCnt <= 1): return None # Condition not handled
        return ptCnt - 1 if(spline.use_cyclic_u) else ptCnt - 2

    def getSegBezierPts(self):
        spline = self.obj.data.splines[self.splineIdx]
        if(len(self.getSegIdxs()) > 0):
            segIdx = self.getSegIdxs()[0]
            return [spline.bezier_points[segIdx], \
                spline.bezier_points[self.getSegAdjIdx(segIdx)]]
        else: return []

    def getPrevSegBezierPts(self):
        if(len(self.getSegIdxs()) > 0):
            segIdx = self.getSegIdxs()[0]
            prevSegIdx = getAdjIdx(self.obj, self.splineIdx, segIdx, -1)
            if(prevSegIdx != None):
                spline = self.obj.data.splines[self.splineIdx]
                return [spline.bezier_points[prevSegIdx], spline.bezier_points[segIdx]]
        return []

    def getNextSegBezierPts(self):
        if(len(self.getSegIdxs()) > 0):
            segIdx = self.getSegIdxs()[0]
            nextSegIdx = getAdjIdx(self.obj, self.splineIdx, segIdx)
            if(nextSegIdx != None):
                spline = self.obj.data.splines[self.splineIdx]
                return [spline.bezier_points[segIdx], spline.bezier_points[nextSegIdx]]
        return []

    # For convenience
    def getCtrlPts(self):
        return getCtrlPtsForSeg(self.obj, self.splineIdx, self.getSegIdxs()[0])

    def getSegPts(self):
        if(len(self.ptIdxs) == 0): return []
        return getBezierDataForSeg(self.obj, self.splineIdx, self.getSegIdxs()[0])

    def getSegAdjIdx(self, ptIdx):
        adjIdx = getAdjIdx(self.obj, self.splineIdx, ptIdx)
        if(adjIdx == None):
            adjIdx = getAdjIdx(self.obj, self.splineIdx, ptIdx, -1)                
        return adjIdx

    def getPrevSegPts(self):
        idx = getAdjIdx(self.obj, self.splineIdx, self.getSegIdxs()[0], -1)
        return getBezierDataForSeg(self.obj, self.splineIdx, idx) if idx != None else []

    def getNextSegPts(self):
        idx = getAdjIdx(self.obj, self.splineIdx, self.getSegIdxs()[0])
        return getBezierDataForSeg(self.obj, self.splineIdx, idx) if idx != None else []

    def insertNode(self, handleType, select = True):
        if(self._t == None):
            return
        invMw = self.obj.matrix_world.inverted()
        bpts = self.getSegBezierPts()
        insertBezierPts(self.obj, self.splineIdx, \
            self.getSegIdxs()[0], [invMw @ self._clickLoc], handleType)

        if(select):
            self.setCtrlIdxSafe(4)

    def alignHandle(self):
        if(self._ctrlIdx == None):
            return

        co, ptIdx, hdlIdx = self.getCtrlPtCoIdx()
        if(hdlIdx == 1):
            return

        invMw = self.obj.matrix_world.inverted()
        oppIdx = 2 - hdlIdx
        pt = self.getSegPts()[ptIdx]
        diffV = (invMw @ pt[1] - invMw @ pt[oppIdx])

        if(diffV.length):
            co = diffV * ((invMw @ pt[1] - invMw @ pt[hdlIdx])).length / diffV.length
            bpt = self.getSegBezierPts()[ptIdx]
            bpt.handle_right_type = 'FREE'
            bpt.handle_left_type = 'FREE'
            if(hdlIdx == 0):
                bpt.handle_left = bpt.co + co
            else:
                bpt.handle_right = bpt.co + co

    def removeNode(self):
        if(self._ctrlIdx == None):
            return
            
        # TODO: Rename ptIdx in co, ptIdx, hdlIdx.. confusing
        co, ptIdx, hdlIdx = self.getCtrlPtCoIdx()
        segIdx = self.getSegIdxs()[0]
        bptIdx = getAdjIdx(self.obj, self.splineIdx, segIdx, ptIdx)
        if(hdlIdx == 1):
            removeBezierPts(self.obj, self.splineIdx, {bptIdx})
            self.removeSegSel(bptIdx)
        else:
            spline = self.obj.data.splines[self.splineIdx]
            pt = spline.bezier_points[bptIdx]
            pt.handle_right_type = 'FREE'
            pt.handle_left_type = 'FREE'
            if(hdlIdx == 0):
                prevIdx = getAdjIdx(self.obj, \
                    self.splineIdx, bptIdx, -1)
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
                nextIdx = getAdjIdx(self.obj, \
                    self.splineIdx, bptIdx, 1)
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

    # TODO: Redundant data structures
    # TODO: Better: Single Obj with multiples splineIdxs
    def getDisplayInfos(self, segPts = None, hltInfo = None, hideHdls = False,
        selSegCol = None, includeAdj = True):

        # Making long short
        cHltTip = ModalBaseFlexiOp.colHltTip
        cBezPt = ModalBaseFlexiOp.colBezPt
        cHdlPt = ModalBaseFlexiOp.colHdlPtTip
        cAdjBezTip = ModalBaseFlexiOp.colAdjBezTip
        cNonHltTip = ModalBaseFlexiOp.colDrawNonHltSeg

        def getTipList(hltIdx, idx):
            # Display of non-selected segments...
            tipList = [None, cAdjBezTip, None, \
                None, cAdjBezTip, None]
            if(hltIdx == idx): tipList[1] = cHltTip
            if(hltIdx == idx + 1): tipList[4] = cHltTip
            return tipList

        nextIdx = None
        prevIdx = None
        hltHdlIdx = None
        displayInfos = []

        if(hltInfo != None):
            # TODO: Unnecessarily complex
            hltHdlIdx = getCtrlIdxFromSearchInfo(hltInfo)
            if(hltHdlIdx == None):
                endPtIdx = hltInfo[1]
                if(len(self.ptIdxs) > 0):
                    idxs = self.getSegIdxPairs()[0]
                    if(endPtIdx == idxs[0]): hltHdlIdx = 1
                    if(endPtIdx == idxs[1]): hltHdlIdx = 4

        if(selSegCol == None): selSegCol = ModalBaseFlexiOp.colDrawSelSeg

        hltEndPtIdx = hltInfo[1] if(hltInfo != None and hltHdlIdx == None) else None

        if(len(self.ptIdxs) > 0):
            if(segPts == None):
                segPts = self.getSegPts()

            # Display of selected segment...
            tipColors = [cHdlPt, cBezPt, cHdlPt, cHdlPt, cBezPt, cHdlPt]

            if(self._ctrlIdx != None): tipColors[self._ctrlIdx] = \
                ModalBaseFlexiOp.colSelTip

            if(hltHdlIdx != None): tipColors[hltHdlIdx] = cHltTip

            hdlIdxs = [0, 1, 2, 3]

            if(hideHdls):
                hdlIdxs = []
                for i in [0, 2, 3, 5]: tipColors[i] = None

            vertCos = []
            if(self.subdivCnt > 1):
                vertCos = getInterpolatedVertsCo(self.interpPts, self.subdivCnt)[1:-1]

            selSegDisplayInfo = EditSegDisplayInfo(segPts, \
                selSegCol, hdlIdxs, tipColors, vertCos)

            if(not includeAdj):
                return [selSegDisplayInfo]

            # Adj segs change in case of aligned and auto handles
            prevPts = self.getPrevSegPts()
            nextPts = self.getNextSegPts()
            segIdx = self.getSegIdxs()[0]
            nextIdx = getAdjIdx(self.obj, self.splineIdx, segIdx)
            prevIdx = getAdjIdx(self.obj, self.splineIdx, segIdx, -1)
            if((len(prevPts) > 1) and (prevPts == nextPts)):
                prevPts[1] = segPts[0][:]
                prevPts[0] = segPts[1][:]
                displayInfos.append(SegDisplayInfo(prevPts, cNonHltTip, [],  []))
            else:
                if(len(prevPts) > 1):
                    prevPts[1] = segPts[0][:]
                    tipList = getTipList(hltEndPtIdx, prevIdx)
                    displayInfos.append(SegDisplayInfo(prevPts, cNonHltTip, [],  tipList))
                if(len(nextPts) > 1):
                    nextPts[0] = segPts[1][:]
                    tipList = getTipList(hltEndPtIdx, nextIdx)
                    displayInfos.append(SegDisplayInfo(nextPts, cNonHltTip, [], tipList))

        spline = self.obj.data.splines[self.splineIdx]
        for j, pt in enumerate(spline.bezier_points):
            if(len(self.ptIdxs) > 0 and \
                (j == self.getSegIdxs()[0] or j == nextIdx or j == prevIdx)):
                continue
            segPts = getBezierDataForSeg(self.obj, self.splineIdx, j)
            if(segPts != []):
                tipList = [] if(j == (len(spline.bezier_points) - 1)) \
                    else getTipList(hltEndPtIdx, j)
                displayInfos.append(SegDisplayInfo(segPts, cNonHltTip, [], tipList))

        if(len(self.ptIdxs) > 0 ):
            # Append at the end so it's displayed on top of everything else
            displayInfos.append(selSegDisplayInfo)

        return displayInfos

class EditCurveInfo():
    def __init__(self, selCurveInfo):
        self.selCurveInfo = selCurveInfo

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
        pts = self.selCurveInfo.getSegPts()
        if(newPos == None): return pts

        # If ctrilIdx and clickLoc both are valid clickLoc gets priority
        if(self.selCurveInfo.getClickLoc() != None):
            delta = newPos - self.selCurveInfo.getClickLoc()
            if(delta == 0):
                return pts
            t = self.selCurveInfo._t

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

        elif(self.selCurveInfo._ctrlIdx != None):
            co, ptIdx, hdlIdx = self.selCurveInfo.getCtrlPtCoIdx()
            oldPts = [p.copy() for p in pts[:3]] + pts[3:]
            pts[ptIdx][hdlIdx] = newPos
            pts = self.getPtsAfterCtrlPtChange(pts, ptIdx, hdlIdx, oldPts)

        return pts

    def moveSeg(self, newPos):
        pts = self.getOffsetSegPts(newPos)

        invMw = self.selCurveInfo.obj.matrix_world.inverted()
        bpts = self.selCurveInfo.getSegBezierPts()

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

class ModalFlexiEditBezierOp(ModalBaseFlexiOp):
    bl_description = "Flexi editing of Bezier curves in object mode"
    bl_idname = "wm.modal_flexi_edit_bezier"
    bl_label = "Flexi Edit Curve"
    bl_options = {'REGISTER', 'UNDO'}

    ptBatch = None
    h = False

    def drawHandler():
        if(ModalBaseFlexiOp.shader != None):
            ModalBaseFlexiOp.drawHandlerBase()
            bgl.glPointSize(ModalBaseFlexiOp.editSubdivPtSize)
            if(ModalFlexiEditBezierOp.ptBatch != None):
                ModalFlexiEditBezierOp.ptBatch.draw(ModalBaseFlexiOp.shader)

    def resetDisplay():
        ModalFlexiEditBezierOp.ptBatch = getResetBatch(ModalBaseFlexiOp.shader, "POINTS")
        ModalBaseFlexiOp.resetDisplayBase()

    # static method
    def refreshDisplay(displayInfos, locOnCurve = None, snapper = None):

        ptCos = [co for d in displayInfos if type(d) == EditSegDisplayInfo
            for co in d.subdivCos]

        # ~ if(locOnCurve != None): ptCos.append(locOnCurve) # For debugging

        ModalFlexiEditBezierOp.ptBatch = batch_for_shader(ModalBaseFlexiOp.shader, \
            "POINTS", {"pos": ptCos, "color": [ModalBaseFlexiOp.colEditSubdiv \
                for i in range(0, len(ptCos))]})

        ModalBaseFlexiOp.refreshDisplayBase(displayInfos, snapper)

    # Refresh display with existing curves (nonstatic)
    def refreshDisplaySelCurves(self, displayInfosMap = {}, locOnCurve = None):
        if(self.rmInfo == None): return # Possible in updateAfterGeomChange
        displayInfos = list(v for vs in displayInfosMap.values() for v in vs)
        dispCurveInfoObjs = displayInfosMap.keys()
        dispInfoObjs = [c.obj for c in dispCurveInfoObjs] # bpy Curve objects
        for c in self.selectCurveInfos:
            if(c in dispCurveInfoObjs):
                continue
            includeAdj = True
            if(c.obj in dispInfoObjs):
                includeAdj = False
            displayInfos += c.getDisplayInfos(hideHdls = ModalFlexiEditBezierOp.h, \
                includeAdj = includeAdj)

        displayInfos = sorted(displayInfos, key = lambda \
            x:ModalBaseFlexiOp.segColPriority[x.segColor])

        ModalFlexiEditBezierOp.refreshDisplay(displayInfos, locOnCurve, self.snapper)

    def reset(self):
        self.editCurveInfo = None
        self.selectCurveInfos = set()
        ModalFlexiEditBezierOp.resetDisplay()

    def postUndoRedo(self, scene, dummy = None): # signature different in 2.8 and 2.81?
        # ~ self.snapper.customAxis.reload()
        self.updateAfterGeomChange()

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
        # ~ if(tool == None or tool.idname != FlexiEditBezierTool.bl_idname): (T60766)
        if(tool == None or tool.idname != 'flexi_bezier.edit_tool'):
            return False
        return True

    # Will be called after the curve is changed (by the tool or externally)
    # So handle all possible conditions
    def updateAfterGeomChange(self, scene = None, dummy = None): # 3 params in 2.81
        toRemove = []

        removeObjNames = set() # For snaplocs
        addObjNames = set()

        #can never be called during editing, so don't consider editInfo
        for ci in self.selectCurveInfos:
            displayInfos = []
            if(bpy.data.objects.get(ci.objName) != None):
                ci.obj = bpy.data.objects.get(ci.objName) #refresh anyway
                splines = ci.obj.data.splines
                spline = splines[ci.splineIdx]
                bpts = spline.bezier_points
                # Don't keep a point object / spline
                if(len(bpts) == 1):
                    if(len(splines) == 1):
                        toRemove.append(ci)
                    else:
                        splines.remove(spline)
                        if(ci.splineIdx >= (len(splines))):
                            ci.splineIdx = len(splines) - 1
                elif(len(ci.ptIdxs) > 0 and ci.getSegIdxs()[0] >= len(bpts) - 1):
                    ci.addSegSel(ci.getLastSegIdx())
                addObjNames.add(ci.objName)
            else:
                toRemove.append(ci)
                removeObjNames.add(ci.objName)

        if(len(toRemove) > 0):
            for c in toRemove:
                self.selectCurveInfos.remove(c)

        self.updateSnapLocs(addObjNames, removeObjNames)

        self.refreshDisplaySelCurves()

    def subInvoke(self, context, event):
        bpy.app.handlers.undo_post.append(self.postUndoRedo)
        bpy.app.handlers.redo_post.append(self.postUndoRedo)
        bpy.app.handlers.depsgraph_update_post.append(self.updateAfterGeomChange)

        self.editCurveInfo = None
        self.selectCurveInfos = set()
        self.clickT = None
        self.pressT = None
        self.ctrl = False
        self.shift = False
        self.alt = False
        self.subdivMode = False

        # For double click (TODO: remove; same as editCurveInfo == None?)
        self.capture = False
        self.xyPress = None # ...to avoid jerky movement at the beginning

        self.snapInfos = {}
        self.updateSnapLocs()

        return  {"RUNNING_MODAL"}

    def getSnapLocsImpl(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        cis = list(self.selectCurveInfos)
        for i, c in enumerate(cis):
            mw = c.obj.matrix_world
            bpts = c.getSegBezierPts()
            for pt in bpts:
                locs.append(mw @ pt.handle_left)
                locs.append(mw @ pt.handle_right)
        return locs

    def updateSnapLocs(self, addObjNames = None, removeObjNames = None):
        updateCurveEndPtMap(self.snapInfos, addObjNames, removeObjNames)

    def getRefLine(self):
        if(self.editCurveInfo != None):
            ci = self.editCurveInfo.selCurveInfo
            segPts = ci.getSegPts()
            if(ci.getClickLoc() != None):
                return [ci.getClickLoc()]
            coIdx = ci.getCtrlPtCoIdx()
            if(coIdx != None):
                co, ptIdx, hdlIdx = coIdx
                if(hdlIdx in {0, 2}):
                    return [segPts[ptIdx][1], segPts[ptIdx][hdlIdx]]
                else:
                    return [segPts[1-ptIdx][1], segPts[ptIdx][hdlIdx]]
        return []

    def getRefLineOrig(self):
        refLine = self.getRefLine()
        return refLine[0] if len(refLine) > 0 else None

    def getEditableCurveObjs(self):
        return [b for b in bpy.data.objects if isBezier(b) and b.visible_get() \
                and len(b.data.splines[0].bezier_points) > 1]

    def getSearchQueryInfo(self):
        queryInfo = {}
        for ci in self.selectCurveInfos:
            info = queryInfo.get(ci.obj)
            if(info == None):
                info = {}
                queryInfo[ci.obj] = info

            segIdxs = info.get(ci.splineIdx)
            if(segIdxs == None):
                segIdxs = []
                info[ci.splineIdx] = segIdxs

            if(len(ci.ptIdxs) > 0): segIdxs.append(ci.getSegIdxs()[0])
        return queryInfo

    def getSelInfoObj(self, searchResult):
        resType, obj, splineIdx, segIdx, otherInfo = searchResult
        for ci in self.selectCurveInfos:
            if(ci.obj == obj and ci.splineIdx == splineIdx):
                if(len(ci.ptIdxs) > 0 and (segIdx == ci.getSegIdxs()[0])):
                    return ci
                elif(resType == 'CurveBezPt' and (segIdx in ci.ptIdxs)):
                    return ci                    
        return None

    def getSelInfoObjSpline(self, obj, splineIdx):
        for ci in self.selectCurveInfos:
            if(ci.obj == obj and ci.splineIdx == splineIdx):
                return ci
        return None

    def resetMetaBtns(self):
        self.ctrl = False
        self.shift = False
        self.alt = False

    def exclToolRegion(self):
        return False

    def isEditing(self):
        return self.editCurveInfo != None

    def hasSelection(self):
        return len(self.selectCurveInfos) > 0

    def getNewPos(self, refreshStatus):
        selCo = self.editCurveInfo.selCurveInfo.getSelCo()
        xySel = getCoordFromLoc(self.rmInfo.region, self.rmInfo.rv3d, selCo)
        if(self.xyPress != None):
            return self.snapper.get3dLocSnap(self.rmInfo, \
                vec = selCo, refreshStatus = refreshStatus, \
                    xyDelta = [self.xyPress[0] - xySel[0], self.xyPress[1] - xySel[1]])
        else:
            return self.snapper.get3dLocSnap(self.rmInfo, \
                vec = selCo, refreshStatus = refreshStatus)

    def subModal(self, context, event, snapProc):
        rmInfo = self.rmInfo

        if(snapProc): retVal = {"RUNNING_MODAL"}
        else: retVal = {'PASS_THROUGH'}

        if(not snapProc and event.type == 'ESC'):
            if(event.value == 'RELEASE'):
                if(self.editCurveInfo == None): self.reset()
                else:
                    self.capture = False
                    self.editCurveInfo = None
                    self.snapper.resetSnap()
                ModalFlexiEditBezierOp.resetDisplay()
            return {"RUNNING_MODAL"}

        updateMetaBtns(self, event)

        if(self.ctrl and (self.editCurveInfo == None or (self.pressT != None and \
            time.time() - self.pressT) < SNGL_CLK_DURN)):
            bpy.context.window.cursor_set("CROSSHAIR")
        else:
            bpy.context.window.cursor_set("DEFAULT")

        if(event.type == 'E' or event.type == 'e'):
            if(event.value == 'RELEASE'):
                # ~ bpy.ops.wm.tool_set_by_id(name = FlexiDrawBezierTool.bl_idname) (T60766)
                self.reset()
                bpy.ops.wm.tool_set_by_id(name = 'flexi_bezier.draw_tool')
            return {"RUNNING_MODAL"}

        if(event.type in {'W', 'w'}):
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    self.subdivMode = True
                    for c in self.selectCurveInfos: c.subdivMode(rmInfo.rv3d)
                    self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        confirmed = False
        if(not snapProc and event.type in {'SPACE', 'RET'}):
            if(self.subdivMode):
                if(event.value == 'RELEASE'):
                    cis = list(self.selectCurveInfos)
                    for i, c in enumerate(cis):
                        addCnt = (c.subdivCnt - 1)
                        c.subdivSeg()
                        for c1 in cis[(i + 1):]:
                            # if same obj multiple times in selection!!
                            c1.updateSegIdx(c.objName, c.splineIdx, \
                                c.getSegIdxs()[0], addCnt)
                    for c in cis: c.resetPts()
                    bpy.ops.ed.undo_push()
                    self.subdivMode = False
                return {"RUNNING_MODAL"}
            elif(self.editCurveInfo != None):
                confirmed = True

        elif(not snapProc and event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE', \
            'NUMPAD_PLUS', 'NUMPAD_MINUS','PLUS', 'MINUS'}):
            if(len(self.selectCurveInfos) > 0 and self.subdivMode):
                if(event.type in {'NUMPAD_PLUS', 'NUMPAD_MINUS', 'PLUS', 'MINUS'} \
                    and event.value == 'PRESS'):
                    return {'RUNNING_MODAL'}
                elif(event.type =='WHEELDOWNMOUSE' or event.type.endswith('MINUS')):
                    for c in self.selectCurveInfos: c.subdivDecr()
                elif(event.type =='WHEELUPMOUSE' or event.type.endswith('PLUS')):
                    for c in self.selectCurveInfos: c.subdivIncr()
                self.refreshDisplaySelCurves()
                return {'RUNNING_MODAL'}

        if(event.type == 'H' or event.type == 'h'):
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    ModalFlexiEditBezierOp.h = not ModalFlexiEditBezierOp.h
                    self.refreshDisplaySelCurves()
                return {"RUNNING_MODAL"}

        if(event.type == 'DEL'):
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    for c in self.selectCurveInfos: c.removeNode() #selected node
                    # will be taken care by depsgraph?
                    self.updateAfterGeomChange()
                    bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if(event.type in {'K', 'k'}):
            #TODO: check _ctrlIdx != None for any
            if(len(self.selectCurveInfos) > 0):
                if(event.value == 'RELEASE'):
                    for c in self.selectCurveInfos: c.alignHandle() #selected node
                    bpy.ops.ed.undo_push()
                return {"RUNNING_MODAL"}

        if(not snapProc and not self.capture \
            and event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            self.xyPress = rmInfo.xy[:]
            coFind = Vector(rmInfo.xy).to_3d()

            objs = self.getEditableCurveObjs()

            selObjInfos = self.getSearchQueryInfo()

            #TODO: Move to Snapper?
            searchResult = getClosestPt2d(rmInfo.region, rmInfo.rv3d, coFind, objs, \
                NONSEL_CURVE_SEARCH_RES, selObjInfos, NONSEL_CURVE_SEARCH_RES, \
                    withHandles = (not self.ctrl and not ModalFlexiEditBezierOp.h))

            if(searchResult != None):
                resType, obj, splineIdx, segIdx, otherInfo = searchResult

                for ci in self.selectCurveInfos: ci.resetSel()

                ci = self.getSelInfoObj(searchResult)

                if(ci == None):
                    ci = SelectCurveInfo(obj, splineIdx, segIdx)
                    if(not self.shift or self.ctrl): self.selectCurveInfos = set()
                    self.selectCurveInfos.add(ci)

                if(resType  == 'SelHandles'):
                    ci.setCtrlIdxSafe(getCtrlIdxFromSearchInfo([resType, otherInfo]))
                elif(resType  == 'CurveBezPt'):
                    if(len(ci.ptIdxs) > 0 and segIdx == ci.getSegIdxs()[0]): 
                        ci.setCtrlIdxSafe(1)
                    else: ci.setCtrlIdxSafe(4)
                else:
                    # More precise for adding point
                    selRes = ADD_PT_CURVE_SEARCH_RES if self.ctrl \
                        else SEL_CURVE_SEARCH_RES

                    searchResult = getClosestPt2dWithinSeg(rmInfo.region, rmInfo.rv3d, \
                        coFind, selObj = obj, selSplineIdx = splineIdx, \
                            selSegIdx = segIdx, selObjRes = selRes, \
                                withHandles = False, withBezPts = False)

                    if(searchResult != None): #Must never be None
                        resType, obj, splineIdx, segIdx, otherInfo = searchResult
                        ci.setClickLocSafe(otherInfo)
                        if(ci._t == None): ci = None

                self.editCurveInfo = EditCurveInfo(ci)

                self.refreshDisplaySelCurves()
                self.pressT = time.time()
                return {'RUNNING_MODAL'}

            if(not self.shift):
                self.reset()

            return retVal

        if(confirmed or self.snapper.digitsConfirmed or \
            (event.type == 'LEFTMOUSE' and event.value == 'RELEASE')):
            if(self.editCurveInfo == None):
                return retVal

            ei = self.editCurveInfo
            tm = time.time()

            if(self.clickT != None and (tm - self.clickT) < DBL_CLK_DURN):
                self.capture = True
                self.clickT = None
            else:
                self.capture = False
                if(self.pressT != None and (tm - self.pressT) < SNGL_CLK_DURN):
                    if(self.ctrl):
                        if(self.shift): handleType = 'ALIGNED'
                        elif(self.alt): handleType = 'VECTOR'
                        else: handleType = 'FREE'
                        for c in self.selectCurveInfos: c.insertNode(handleType)
                        bpy.ops.ed.undo_push()
                        ModalFlexiEditBezierOp.resetDisplay()
                    # Gib dem Benutzer Zeit zum Atmen!
                    else:
                        ModalFlexiEditBezierOp.resetDisplay()
                else:
                    ei.moveSeg(self.getNewPos(refreshStatus = False))
                    bpy.ops.ed.undo_push()

                self.clickT = tm
                self.editCurveInfo = None
                self.snapper.resetSnap()
                self.updateAfterGeomChange()

            self.pressT = None
            return {"RUNNING_MODAL"}

        elif(snapProc or event.type == 'MOUSEMOVE'):
            displayInfosMap = {}
            ei = self.editCurveInfo
            locOnCurve = None # For debug

            if(ei != None):
                # User is editing curve or control points (left mouse pressed)
                ci = ei.selCurveInfo
                segPts = ei.getOffsetSegPts(self.getNewPos(refreshStatus = True))
                displayInfosMap = {ci: ci.getDisplayInfos(segPts, \
                    hideHdls = ModalFlexiEditBezierOp.h)}
            else:
                # ~ coFind = Vector(rmInfo.xy).to_3d()
                coFind = getCoordFromLoc(rmInfo.region, rmInfo.rv3d, \
                    self.snapper.get3dLocSnap(rmInfo)).to_3d()

                objs = self.getEditableCurveObjs()

                #Sel obj: low res (highlight only seg)
                selObjInfos = self.getSearchQueryInfo()

                #TODO: Move to Snapper
                searchResult = getClosestPt2d(rmInfo.region, rmInfo.rv3d, coFind, objs, \
                    NONSEL_CURVE_SEARCH_RES, selObjInfos, NONSEL_CURVE_SEARCH_RES, \
                        withHandles = (not self.ctrl and not ModalFlexiEditBezierOp.h))

                if(searchResult != None):
                    hltInfo = None

                    resType, obj, splineIdx, segIdx, otherInfo = searchResult
                    ci = self.getSelInfoObj(searchResult)

                    if(resType in {'SelHandles', 'CurveBezPt'}):
                        if(resType == 'CurveBezPt'):
                            ci = self.getSelInfoObjSpline(obj, splineIdx)
                            if(ci == None): segIdx = None
                        hltInfo = [resType, otherInfo]
                    else:
                        locOnCurve = otherInfo

                    if(ci == None):
                        ci = SelectCurveInfo(obj, splineIdx, segIdx)
                        segColor = ModalBaseFlexiOp.colDrawHltSeg
                        hideHdl = True
                    else:
                        segColor = ModalBaseFlexiOp.colDrawSelSeg
                        hideHdl = ModalFlexiEditBezierOp.h

                    displayInfosMap = {ci: ci.getDisplayInfos(ci.getSegPts(), \
                        hltInfo, hideHdl, segColor)}

            self.refreshDisplaySelCurves(displayInfosMap, locOnCurve)

            return retVal

        return retVal

###################### Global Params ######################

def getSnapOrientTups(scene, context):
    orients = [\
       ('GLOBAL', 'Global Axes', "Orient to world space"), \
       ('AXIS', 'Custom Axes', "Orient to custom axis (if available)"), \
       ('VIEW', 'View', "Orient to window"), \
      ]

    tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)

    # ~ if(tool == None or tool.idname != FlexiDrawBezierTool.bl_idname): (T60766)
    if(tool != None and (tool.idname == 'flexi_bezier.draw_tool' \
        or tool.idname == 'flexi_bezier.grease_draw_tool')):
        orients.insert(1, ('REFERENCE', 'Reference Line', "Orient to preceding segment"))

    if(tool != None and tool.idname == 'flexi_bezier.edit_tool'):
        orients.insert(1, ('REFERENCE', 'Reference Line', \
            "Orient to preceding segment line or current handle"))

    # ~ if(context.active_object != None):
    orients.insert(3, ('OBJECT', 'Active Object', \
        "Orient to local space of active object"))
    orients.insert(4, ('FACE', 'Selected Object Face', \
        "Orient to normal of face of selected object under mouse pointer "))
    return orients

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
    elif(transType in {'VIEW', 'REFERENCE'}): keyset = [0] + [i for i in range(4, 7)]

    return [axesMap[key] for key in keyset]

def getSnapOriginTups(scene = None, context = None):
    return [\
       ('GLOBAL', 'Global Origin', \
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
            "Selected object face under mouse pointer"), \
      ]

class BezierToolkitParams(bpy.types.PropertyGroup):
    snapOrient: EnumProperty(name = 'Orientation',#"Align contrained axes and snap angle to",
            items = getSnapOrientTups,
            description='Orientation for Draw / Edit')

    snapOrigin: EnumProperty(name = 'Origin',#"Align contrained axes and snap angle to",
            items = getSnapOriginTups,
            description='Origin for Draw / Edit')

    constrAxes: EnumProperty(name = 'Constrain Axis', #"Constrain axis for draw and edit ops",
            items = getConstrAxisTups,
            description='Constrain Draw / Edit Axes')

    snapToPlane: BoolProperty(name="Snap to Plane",
        description='During draw / edit snap the point to the selected plane', \
                    default = False)

    axisScale: BoolProperty(name="Axis Scale", \
        description='Use custom axis scale for grid snap and transform values entered', \
                    default = False)

    customAxisCo1: FloatVectorProperty(default = (LARGE_NO, LARGE_NO, LARGE_NO))
    customAxisCo2: FloatVectorProperty(default = (LARGE_NO, LARGE_NO, LARGE_NO))
    customAxisSnapCnt: IntProperty(default = 3, min = 0)

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
kmToolFlexiGreaseDrawBezier = "3D View Tool: Object, Flexi Grease Draw Bezier"

def showSnapToPlane(params):
    return (params.snapOrient != 'VIEW' and params.snapOrient != 'REFERENCE' and\
        hasattr(params, 'constrAxes') and params.constrAxes.startswith('shift'))

def drawSettingsFT(self, context):
    params = bpy.context.window_manager.bezierToolkitParams
    self.layout.use_property_split = True
    self.layout.row(align=True).template_header()
    from bl_ui.space_toolsystem_common import ToolSelectPanelHelper
    tool = ToolSelectPanelHelper.draw_active_tool_header(
        context, self.layout,
        tool_key=('VIEW_3D', context.mode),
    )

    self.layout.use_property_decorate = True
    self.layout.prop(params, "snapOrient", text = '')
    self.layout.prop(params, "snapOrigin", text = '')
    self.layout.prop(params, "constrAxes", text = '')
    if(params.constrAxes not in [a[0] for a in getConstrAxisTups()]):
        params.constrAxes = 'NONE'

    # Only available for planes not axis
    if(showSnapToPlane(params)):
        self.layout.prop(params, "snapToPlane")

    if(params.snapOrient in {'AXIS', 'REFERENCE'}):
        self.layout.prop(params, "axisScale")

@ToolDef.from_fn
def toolFlexiDraw():

    return dict(idname = "flexi_bezier.draw_tool",
        label = "Flexi Draw Bezier",
        description = "Flexible drawing of Bezier curves in object mode",
        icon = "ops.gpencil.extrude_move",
        widget = None,
        keymap = kmToolFlexiDrawBezier,
        # ~ draw_settings = drawSettingsFT,
        )

@ToolDef.from_fn
def toolFlexiGreaseDraw():

    return dict(idname = "flexi_bezier.grease_draw_tool",
        label = "Flexi Grease Bezier",
        description = "Flexible drawing of Bezier curves as grease pencil strokes",
        icon = "ops.gpencil.extrude_move",
        widget = None,
        keymap = kmToolFlexiGreaseDrawBezier,
        # ~ draw_settings = drawSettingsFT,
        )

@ToolDef.from_fn
def toolFlexiEdit():

    return dict(idname = "flexi_bezier.edit_tool",
        label = "Flexi Edit Bezier",
        description = "Flexible editing of Bezier curves in object mode",
        icon = "ops.pose.breakdowner",
        widget = None,
        keymap = kmToolFlexiEditBezier,
        # ~ draw_settings = drawSettingsFT,
    )

def getToolList(spaceType, contextMode):
    from bl_ui.space_toolsystem_common import ToolSelectPanelHelper
    cls = ToolSelectPanelHelper._tool_class_from_space_type(spaceType)
    return cls._tools[contextMode]

def registerFlexiBezierTools():
    tools = getToolList('VIEW_3D', 'OBJECT')
    tools += None, toolFlexiDraw, toolFlexiEdit
    # ~ tools += None, toolFlexiEdit
    del tools

    tools = getToolList('VIEW_3D', 'PAINT_GPENCIL')
    tools += None, toolFlexiGreaseDraw
    del tools

def unregisterFlexiBezierTools():
    tools = getToolList('VIEW_3D', 'OBJECT')

    index = tools.index(toolFlexiDraw) - 1 #None
    tools.pop(index)
    tools.remove(toolFlexiDraw)

    index = tools.index(toolFlexiEdit) - 1 #None
    tools.pop(index)
    tools.remove(toolFlexiEdit)
    del tools

    tools = getToolList('VIEW_3D', 'PAINT_GPENCIL')
    index = tools.index(toolFlexiGreaseDraw) - 1 #None
    tools.pop(index)
    tools.remove(toolFlexiGreaseDraw)
    del tools


keymapDraw = (kmToolFlexiDrawBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": [
            ("wm.flexi_draw_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
             {"properties": []}),
        ]},)

emptyKeymapDraw = (kmToolFlexiDrawBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": []},)

keymapGreaseDraw = (kmToolFlexiGreaseDrawBezier,
        {"space_type": 'VIEW_3D', "region_type": 'WINDOW'},
        {"items": [
            ("wm.flexi_draw_grease_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
             {"properties": []}),
        ]},)

emptyKeymapGreaseDraw = (kmToolFlexiGreaseDrawBezier,
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
    keyconfig_init_from_data(kc_defaultconf, [emptyKeymapDraw, \
        emptyKeymapGreaseDraw, emptyKeymapEdit])

    keyconfig_init_from_data(kc_addonconf, [keymapDraw, keymapGreaseDraw, keymapEdit])

def unregisterFlexiBezierKeymaps():
    keyconfigs = bpy.context.window_manager.keyconfigs
    defaultmap = keyconfigs.get("blender").keymaps
    addonmap   = keyconfigs.get("blender addon").keymaps

    for km_name, km_args, km_content in [keymapDraw, keymapGreaseDraw, keymapEdit]:
        keymap = addonmap.find(km_name, **km_args)
        keymap_items = keymap.keymap_items
        for item in km_content['items']:
            item_id = keymap_items.find(item[0])
            if item_id != -1:
                keymap_items.remove(keymap_items[item_id])
        addonmap.remove(keymap)
        defaultmap.remove(defaultmap.find(km_name, **km_args))

# ****************** Configurations In User Preferences ******************

def updatePanel(self, context):
    try:
        panel = BezierUtilsPanel
        if "bl_rna" in panel.__dict__:
            bpy.utils.unregister_class(panel)

        panel.bl_category = context.preferences.addons[__name__].preferences.category
        bpy.utils.register_class(panel)

    except Exception as e:
        print("BezierUtils: Updating Panel locations has failed", e)

def updateProps(self, context):
    try:
        prefs = context.preferences.addons[__name__].preferences
        BezierUtilsPanel.bl_category = prefs.category
        ModalBaseFlexiOp.drawPtSize = prefs.drawPtSize
        ModalBaseFlexiOp.lineWidth = prefs.drawLineWidth

        ModalBaseFlexiOp.axisLineWidth = prefs.axisLineWidth
        ModalBaseFlexiOp.snapPtSize = prefs.snapPtSize
        ModalBaseFlexiOp.editSubdivPtSize = prefs.editSubdivPtSize
        ModalBaseFlexiOp.greaseSubdivPtSize = prefs.greaseSubdivPtSize

        ModalBaseFlexiOp.colDrawSelSeg = prefs.colDrawSelSeg
        ModalBaseFlexiOp.colDrawNonHltSeg = prefs.colDrawNonHltSeg
        ModalBaseFlexiOp.colDrawHltSeg = prefs.colDrawHltSeg
        ModalBaseFlexiOp.colDrawMarker = prefs.colDrawMarker

        ModalBaseFlexiOp.colGreaseSelSeg = prefs.colGreaseSelSeg
        ModalBaseFlexiOp.colGreaseNonHltSeg = prefs.colGreaseNonHltSeg
        ModalBaseFlexiOp.colGreaseMarker = prefs.colGreaseMarker

        ModalBaseFlexiOp.colHdlFree = prefs.colHdlFree
        ModalBaseFlexiOp.colHdlVector = prefs.colHdlVector
        ModalBaseFlexiOp.colHdlAligned = prefs.colHdlAligned
        ModalBaseFlexiOp.colHdlAuto = prefs.colHdlAuto

        ModalBaseFlexiOp.colSelTip = prefs.colSelTip
        ModalBaseFlexiOp.colHltTip = prefs.colHltTip
        ModalBaseFlexiOp.colBezPt = prefs.colBezPt
        ModalBaseFlexiOp.colHdlPtTip = prefs.colHdlPtTip
        ModalBaseFlexiOp.colAdjBezTip = prefs.colAdjBezTip

        ModalBaseFlexiOp.colEditSubdiv = prefs.colEditSubdiv
        ModalBaseFlexiOp.colGreaseSubdiv = prefs.colGreaseSubdiv
        ModalBaseFlexiOp.colGreaseBezPt = prefs.colGreaseBezPt

    except Exception as e:
        print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
        ModalBaseFlexiOp.drawPtSize = 4
        ModalBaseFlexiOp.lineWidth = 1.5
        ModalBaseFlexiOp.axisLineWidth = .25
        ModalBaseFlexiOp.snapPtSize = 2
        ModalBaseFlexiOp.editSubdivPtSize = 6
        ModalBaseFlexiOp.greaseSubdivPtSize = 4

        ModalBaseFlexiOp.colDrawSelSeg = (.6, .8, 1, 1) # DRAW_SEL_SEG_COLOR
        ModalBaseFlexiOp.colDrawNonHltSeg = (.1, .4, .6, 1) # DRAW_ADJ_SEG_COLOR
        ModalBaseFlexiOp.colDrawHltSeg = (.2, .6, .9, 1) # DRAW_NONADJ_SEG_COLOR

        ModalBaseFlexiOp.colGreaseSelSeg = (0.2, .8, 0.2, 1) # GREASE_SEL_SEG_COLOR
        ModalBaseFlexiOp.colGreaseNonHltSeg = (0.2, .6, 0.2, 1) # GREASE_ADJ_SEG_COLOR

        ModalBaseFlexiOp.colHdlFree = (.6, .05, .05, 1) # HANDLE_COLOR_FREE
        ModalBaseFlexiOp.colHdlVector = (.4, .5, .2, 1) # HANDLE_COLOR_VECTOR
        ModalBaseFlexiOp.colHdlAligned = (1, .3, .3, 1) # HANDLE_COLOR_ALIGNED
        ModalBaseFlexiOp.colHdlAuto = (.8, .5, .2, 1) # HANDLE_COLOR_AUTO

        ModalBaseFlexiOp.colDrawMarker = (.6, .8, 1, 1) # DRAW_MARKER_COLOR
        ModalBaseFlexiOp.colGreaseMarker = (0.2, .8, 0.2, 1) # GREASE_MARKER_COLOR

        ModalBaseFlexiOp.colSelTip = (.2, .7, .3, 1) # SEL_TIP_COLOR
        ModalBaseFlexiOp.colHltTip = (.2, 1, .9, 1) # HLT_TIP_COLOR
        ModalBaseFlexiOp.colBezPt = (1, 1, 0, 1) # ENDPT_TIP_COLOR
        ModalBaseFlexiOp.colHdlPtTip = (.7, .7, 0, 1) # TIP_COLOR
        ModalBaseFlexiOp.colAdjBezTip = (.1, .1, .1, 1) # ADJ_ENDPT_TIP_COLOR

        ModalBaseFlexiOp.colEditSubdiv = (.3, 0, 0, 1) # EDIT_SUBDIV_PT_COLOR

        ModalBaseFlexiOp.colGreaseSubdiv = (1, .3, 1, 1) # GREASE_SUBDIV_PT_COLOR
        ModalBaseFlexiOp.colGreaseBezPt = (1, .3, 1, 1) # GREASE_ENDPT_TIP_COLOR

    ModalBaseFlexiOp.segColPriority = {ModalBaseFlexiOp.colDrawNonHltSeg: 0, \
        ModalBaseFlexiOp.colDrawHltSeg: 1, ModalBaseFlexiOp.colDrawSelSeg: 2}

    ModalBaseFlexiOp.hdlColMap ={'FREE': ModalBaseFlexiOp.colHdlFree, \
        'VECTOR': ModalBaseFlexiOp.colHdlVector,  \
            'ALIGNED': ModalBaseFlexiOp.colHdlAligned, \
                'AUTO': ModalBaseFlexiOp.colHdlAuto}

    ModalBaseFlexiOp.tipColPriority = {ModalBaseFlexiOp.colAdjBezTip: 0, \
        ModalBaseFlexiOp.colHdlPtTip: 1, ModalBaseFlexiOp.colBezPt: 2, \
            ModalBaseFlexiOp.colGreaseBezPt: 3, ModalBaseFlexiOp.colSelTip : 4, \
                ModalBaseFlexiOp.colHltTip : 5, ModalBaseFlexiOp.colDrawMarker: 6, \
                    ModalBaseFlexiOp.colGreaseMarker: 7}

    if(updateProps.refreshDisp):
        try: ModalBaseFlexiOp.opObj.refreshDisplaySelCurves()
        except: pass

updateProps.refreshDisp = True


class ResetDefaultPropsOp(bpy.types.Operator):
    bl_idname = "object.reset_default_props"
    bl_label = "Reset"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Reset all values to default"

    def execute(self, context):
        panel = BezierUtilsPanel
        prefs = context.preferences.addons[__name__].preferences
        props = prefs.bl_rna.properties
        updateProps.refreshDisp = False
        for prop in props:
            try:
                if(not hasattr(prop, "default")): continue
                if(prop.type == 'STRING' and prop.identifier == 'category'):
                    execStr = 'prefs.' + prop.identifier + ' =  "' + \
                        str(prop.default) + '"'
                    exec(execStr)
                elif hasattr(prop, "default_array") and getattr(prop, "is_array", True):
                    for i, a in enumerate(prop.default_array):
                        execStr = 'prefs.' + prop.identifier + '[' + str(i) + '] = ' + str(a)
                        exec(execStr)
                else:
                    execStr = 'prefs.' + prop.identifier + ' =  ' + str(prop.default)
                    exec(execStr)
            except Exception as e:
                print(e)
        updateProps.refreshDisp = True
        updateProps(None, context)
        return {'FINISHED'}


class BezierUtilsPreferences(AddonPreferences):
    bl_idname = __name__

    category: StringProperty(
            name = "Tab Category",
            description = "Choose a name for the category of the panel",
            default = "Tool",
            update = updatePanel
    )

    drawLineWidth: FloatProperty(
            name = "Line Thickness",
            description = "Thickness of segment & handle Lines of Flexi Draw and Edit",
            default = 1.5,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    drawPtSize: FloatProperty(
            name = "Handle Point Size",
            description = "Size of Flexi Draw and Edit Bezier handle points",
            default = 4,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    editSubdivPtSize: FloatProperty(
            name = "Uniform Subdiv Point Size",
            description = "Size of point marking subdivisions",
            default = 6,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    greaseSubdivPtSize: FloatProperty(
            name = "Flexi Grease Res Point Size",
            description = "Size of point marking resoulution in Flexi Grease Tool",
            default = 4,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    markerSize: FloatProperty(
            name = "Marker Size",
            description = "Size of Flexi Draw and Mark Starting Vertices",
            default = 6,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    axisLineWidth: FloatProperty(
            name = "Axis Line Thickness",
            description = "Thickness of Axis Lines for snapping & locking",
            default = 0.25,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    snapPtSize: FloatProperty(
            name = "Snap Point Size",
            description = "Size of snap point indicator",
            default = 2,
            min = 0.1,
            max = 20,
            update = updateProps
    )

    colDrawSelSeg: bpy.props.FloatVectorProperty(
        name="Selected Draw / Edit  Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.6, .8, 1, 1), \
                description = 'Color of the segment being drawn / edited',
                    update = updateProps
    )

    colDrawNonHltSeg: bpy.props.FloatVectorProperty(
        name="Adjacent Draw / Edit Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.1, .4, .6, 1), \
                description = 'Color of the segment adjacent to the' + \
                    'one being drawn / edited', update = updateProps
    )

    colDrawHltSeg: bpy.props.FloatVectorProperty(
        name="Highlighted Edit Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.2, .6, .9, 1), \
                description = 'Color of the segment under mouse curser in Flexi Edit', \
                    update = updateProps
    )

    colDrawMarker: bpy.props.FloatVectorProperty(
        name="Marker", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.6, .8, 1, 1), \
                description = 'Color of the marker', update = updateProps
    )

    colGreaseSelSeg: bpy.props.FloatVectorProperty(
        name="Selected Grease Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(0.2, .8, 0.2, 1), \
                description = 'Color of the segment being drawn', \
                    update = updateProps
    )

    colGreaseNonHltSeg: bpy.props.FloatVectorProperty(
        name="Adjacent Grease Segment", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (0.2, .6, 0.2, 1), \
                description = 'Color of the segment adjacent to the one being drawn', \
                    update = updateProps
    )

    colGreaseMarker: bpy.props.FloatVectorProperty(
        name="Marker", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (0.2, .8, 0.2, 1), \
                description = 'Color of the marker', update = updateProps
    )

    colHdlFree: bpy.props.FloatVectorProperty(
        name="Free Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.6, .05, .05, 1), \
                description = 'Free handle color in all Flexi Tools', \
                        update = updateProps
    )

    colHdlVector: bpy.props.FloatVectorProperty(
        name="Vector Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.4, .5, .2, 1), \
                description = 'Vector handle color in all Flexi Tools', \
                    update = updateProps
    )

    colHdlAligned: bpy.props.FloatVectorProperty(
        name="Aligned Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(1, .3, .3, 1), \
                description = 'Aligned handle color in all Flexi Tools', \
                    update = updateProps
    )

    colHdlAuto: bpy.props.FloatVectorProperty(
        name="Auto Handle", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default=(.8, .5, .2, 1), \
                description = 'Auto handle color in all Flexi Tools', \
                    update = updateProps
    )

    colSelTip: bpy.props.FloatVectorProperty(
        name="Selected Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (.2, .7, .3, 1), \
                description = 'Color of the selected Bezier or handle point', \
                    update = updateProps
    )

    colHltTip: bpy.props.FloatVectorProperty(
        name="Highlighted Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (.2, 1, .9, 1), \
                description = 'Color of Bezier or handle point under mouse pointer', \
                    update = updateProps
    )

    colBezPt: bpy.props.FloatVectorProperty(
        name="Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (1, 1, 0, 1), \
                description = 'Color of nonselected Bezier point', \
                    update = updateProps
    )

    colHdlPtTip: bpy.props.FloatVectorProperty(
        name="Handle Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (.7, .7, 0, 1), \
                description = 'Color of nonselected handle point', \
                    update = updateProps
    )

    colAdjBezTip: bpy.props.FloatVectorProperty(
        name="Adjacent Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (.1, .1, .1, 1), \
                description = 'Color of Bezier points of adjacent segments', \
                    update = updateProps
    )

    colEditSubdiv: bpy.props.FloatVectorProperty(
        name="Uniform Subdiv Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (.3, 0, 0, 1), \
                description = 'Color of point marking subdivisions', \
                    update = updateProps
    )

    colGreaseSubdiv: bpy.props.FloatVectorProperty(
        name="Resolution Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (1, .3, 1, 1), \
                description = 'Color of point marking curve resolution', \
                    update = updateProps
    )

    colGreaseBezPt: bpy.props.FloatVectorProperty(
        name="Bezier Point", subtype="COLOR", size=4, min=0.0, max=1.0,\
            default = (1, .3, 1, 1), \
                description = 'Color of Bezier point', \
                    update = updateProps
    )

    elemDimsExpanded: BoolProperty(name="Draw Elem Expanded State", default = False)
    drawColExpanded: BoolProperty(name="Draw Col Expanded State", default = False)
    greaseColExpanded: BoolProperty(name="Grease Col Expanded State", default = False)
    handleColExpanded: BoolProperty(name="Handle Col Expanded State", default = False)

    def draw(self, context):
        layout = self.layout
        col = layout.column().split()
        col.label(text="Tab Category:")
        col.prop(self, "category", text="")

        layout.row().label(text="UI Preferences:")

        row = layout.row()
        row.prop(self, "elemDimsExpanded", icon = "TRIA_DOWN" \
            if self.elemDimsExpanded else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Draw / Edit Element Dimensions (Common)")#, icon = 'UNLINKED')

        if self.elemDimsExpanded:
            col = layout.column().split()
            col.label(text='Segment / Line Thickness:')
            col.prop(self, "drawLineWidth", text = '')
            col = layout.column().split()
            col.label(text='Handle Point Size:')
            col.prop(self, "drawPtSize", text = '')
            col = layout.column().split()
            col.label(text='Uniform Subdiv Point Size:')
            col.prop(self, "editSubdivPtSize", text = '')
            col = layout.column().split()
            col.label(text='Flexi Grease Res Point Size:')
            col.prop(self, "greaseSubdivPtSize", text = '')
            col = layout.column().split()
            col.label(text='Marker Size:')
            col.prop(self, "markerSize", text = '')
            col = layout.column().split()
            col.label(text='Axis Line Thickness:')
            col.prop(self, "axisLineWidth", text = '')
            col = layout.column().split()
            col.label(text='Snap Point Size:')
            col.prop(self, "snapPtSize", text = '')

        layout.column().separator()

        row = layout.row()
        row.prop(self, "drawColExpanded", icon = "TRIA_DOWN" \
            if self.drawColExpanded else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Flexi Draw / Edit Colors")#, icon = 'UNLINKED')

        if self.drawColExpanded:
            col = layout.column().split()
            col.label(text="Selected Draw / Edit  Segment:")
            col.prop(self, "colDrawSelSeg", text = '')
            col = layout.column().split()
            col.label(text="Adjacent Draw / Edit Segment:")
            col.prop(self, "colDrawNonHltSeg", text = '')
            col = layout.column().split()
            col.label(text="Highlighted Edit Segment:")
            col.prop(self, "colDrawHltSeg", text = '')
            col = layout.column().split()
            col.label(text="Draw Marker:")
            col.prop(self, "colDrawMarker", text = '')
            col = layout.column().split()
            col.label(text="Bezier Point:")
            col.prop(self, "colBezPt", text = '')
            col = layout.column().split()
            col.label(text="Subdivision Marker:")
            col.prop(self, "colEditSubdiv", text = '')

        layout.column().separator()

        row = layout.row()
        row.prop(self, "greaseColExpanded", icon = "TRIA_DOWN" \
            if self.greaseColExpanded else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Flexi Grease Colors")#, icon = 'UNLINKED')

        if self.greaseColExpanded:
            col = layout.column().split()
            col.label(text="Selected Grease Segment:")
            col.prop(self, "colGreaseSelSeg", text = '')
            col = layout.column().split()
            col.label(text="Adjacent Grease Segment:")
            col.prop(self, "colGreaseNonHltSeg", text = '')
            col = layout.column().split()
            col.label(text="Draw Marker:")
            col.prop(self, "colGreaseMarker", text = '')
            col = layout.column().split()
            col.label(text="Bezier Point:")
            col.prop(self, "colGreaseBezPt", text = '')
            col = layout.column().split()
            col.label(text="Curve Resolution Marker:")
            col.prop(self, "colGreaseSubdiv", text = '')

        layout.column().separator()

        row = layout.row()
        row.prop(self, "handleColExpanded", icon = "TRIA_DOWN" \
            if self.handleColExpanded else "TRIA_RIGHT",  icon_only = True, emboss = False)
        row.label(text = "Handle Colors (Common)")#, icon = 'UNLINKED')

        if self.handleColExpanded:
            col = layout.column().split()
            col.label(text="Free Handle:")
            col.prop(self, "colHdlFree", text = '')
            col = layout.column().split()
            col.label(text="Vector Handle:")
            col.prop(self, "colHdlVector", text = '')
            col = layout.column().split()
            col.label(text="Aligned Handle:")
            col.prop(self, "colHdlAligned", text = '')
            col = layout.column().split()
            col.label(text="Auto Handle:")
            col.prop(self, "colHdlAuto", text = '')
            col = layout.column().split()
            col.label(text="Selected Point:")
            col.prop(self, "colSelTip", text = '')
            col = layout.column().split()
            col.label(text="Highlighted Point:")
            col.prop(self, "colHltTip", text = '')
            col = layout.column().split()
            col.label(text="Handle Point:")
            col.prop(self, "colHdlPtTip", text = '')
            col = layout.column().split()
            col.label(text="Adjacent Bezier Point:")
            col.prop(self, "colAdjBezTip", text = '')
        col = layout.column()
        col.operator('object.reset_default_props')

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
    SetCurveColorOp,
    PasteLengthOp,
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
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.WindowManager.bezierToolkitParams = \
        bpy.props.PointerProperty(type = BezierToolkitParams)

    BezierUtilsPanel.colorCurves(add = True)
    bpy.app.handlers.depsgraph_update_post.append(BezierUtilsPanel.colorCurves)

    # ~ bpy.utils.register_tool(FlexiDrawBezierTool) (T60766)
    # ~ bpy.utils.register_tool(FlexiEditBezierTool) (T60766)
    registerFlexiBezierTools()
    registerFlexiBezierKeymaps()
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

    unregisterFlexiBezierKeymaps()
    unregisterFlexiBezierTools()

    del bpy.types.WindowManager.bezierToolkitParams

    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

    # ~ bpy.utils.unregister_tool(FlexiDrawBezierTool) (T60766)
    # ~ bpy.utils.unregister_tool(FlexiEditBezierTool) (T60766)
