#
#
# Blender add-on with Bezier Curve utility ops
#
# Supported Blender Version: 2.8 Beta
#
# Copyright (C) 2019  Shrinivas Kulkarni

# License: MIT (https://github.com/Shriinivas/blenderbezierutils/blob/master/LICENSE)
#
# Not yet pep8 compliant

import bpy, bmesh, bgl, gpu
from bpy.props import BoolProperty, IntProperty
from bpy.types import Panel, Operator, WorkSpaceTool
from mathutils import Vector, Matrix, geometry
from math import log, atan, tan, pi, radians
from bpy_extras.view3d_utils import region_2d_to_vector_3d, region_2d_to_location_3d
from bpy_extras.view3d_utils import location_3d_to_region_2d
from gpu_extras.batch import batch_for_shader
import time

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
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
        and bObj.data.splines[0].type == 'BEZIER'

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

def createSpline(curveData, srcSpline, forceNoncyclic, freeHandles):
    spline = curveData.splines.new('BEZIER')
    spline.bezier_points.add(len(srcSpline.bezier_points)-1)

    if(forceNoncyclic):
        spline.use_cyclic_u = False
    else:
        spline.use_cyclic_u = srcSpline.use_cyclic_u

    for i in range(0, len(srcSpline.bezier_points)):
        copyBezierPt(srcSpline.bezier_points[i], spline.bezier_points[i], freeHandles)

    if(forceNoncyclic == True and srcSpline.use_cyclic_u == True):
        spline.bezier_points.add(1)
        copyBezierPt(srcSpline.bezier_points[0], spline.bezier_points[-1], freeHandles)

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

#Change position of bezier points according to new matrix_world
def changeMW(obj, newMW):
    invMW = newMW.inverted()
    for spline in obj.data.splines:
        for pt in spline.bezier_points:
            pt.co = invMW @ (obj.mw @ pt.co)
            pt.handle_left = invMW @ (obj.mw @ pt.handle_left)
            pt.handle_right = invMW @ (obj.mw @ pt.handle_right)
    obj.matrix_world = newMW

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

        bpy.context.view_layer.objects.active = curve

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

        bpy.context.view_layer.objects.active = curve

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

        bpy.context.view_layer.objects.active = curve

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

        bpy.context.view_layer.objects.active = curve

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
        context.window_manager.AssignShapeKeyParams.markVertex = False

    def modal (self, context, event):
        params = context.window_manager.AssignShapeKeyParams

        if(context.mode  == 'OBJECT' or event.type == "ESC" or\
            not context.window_manager.AssignShapeKeyParams.markVertex):
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

########## Draw Flexi Bezier Curve - Start ###############

def  getViewDistRounding(context):
    viewDist = context.space_data.region_3d.view_distance
    return int(log(viewDist, 10)) - 1

def getCoordFromLoc(context, loc):
    region = context.region
    rv3d = context.space_data.region_3d
    return location_3d_to_region_2d(region, rv3d, loc)

#Round to logarithmic scale .1, 0, 10, 100 etc.
#(47.538, -1) -> 47.5; (47.538, 0) -> 48.0; (47.538, 1) -> 50.0; (47.538, 2) -> 0,
def roundedVect(vect, rounding):
    rounding += 1
    fact = (10 ** rounding) / 10
    return Vector([round(l / fact) * fact for l in vect])

def isOutside(context, event):
    x = event.mouse_region_x
    y = event.mouse_region_y
    region = context.region
    return (x < 0 or x > region.width or y < 0 or y > region.height)

# Get pt coords along curve defined by the four control pts (curvePts)
# divPerUnitLength: No of subdivisions per unit length
# (which is the same as no of pts excluding the end pts)
def getPtsAlongBezier(context, curvePts, curveRes, normalized):
    if(len(curvePts) < 2):
        return []

    pts = []

    if(normalized):
        viewDist = context.space_data.region_3d.view_distance

        # (the smaller the view dist (higher zoom level),
        # the higher the num of subdivisions
        curveRes = curveRes / viewDist

    for i in range(0, len(curvePts) - 1):
        pt0 = curvePts[i]
        pt1 = curvePts[i+1]
        curve = [pt0[1], pt0[2], pt1[0], pt1[1]]

        #Add 2 for start and end
        res = 2 + int((curve[2] - curve[1]).length * curveRes)
        incr = 1 / (res - 1)

        for j in range(0, res):
            t = j * incr
            c = (1 - t)
            loc = (c ** 3) * curve[0] + 3 * (c ** 2) * t * curve[1] + \
                3 * c * (t ** 2) * curve[2] + (t ** 3) * curve[3]
            pts.append(loc)

    return pts


class ModalDrawBezierOp(Operator):

    def __init__(self, curveDispRes = 20, defLineWidth = 1.5, defPointSize = 7):
        self.defLineWidth = defLineWidth
        self.defPointSize = defPointSize

        #No of subdivisions in the displayed curve with view dist 1
        self.curveDispRes = curveDispRes
        self.drawHandlerRef = None
        self.addHandle = True
        self.defaultSnapSteps = 3

    def invoke(self, context, event):
        self.cleanup(context)
        self.initialize()

        self.shader = gpu.shader.from_builtin('3D_FLAT_COLOR')
        self.shader.bind()
        
        context.window_manager.modal_handler_add(self)

        return {"RUNNING_MODAL"}

    #This will be called multiple times not just at the beginning
    def initialize(self):
        self.curvePts = []
        self.clickT = None #For double click
        self.capture = False
        self.ctrl = False
        self.shift = False
        self.alt = False
        self.axis = None
        self.snapSteps = self.defaultSnapSteps
        self.drawHandlerRef = bpy.types.SpaceView3D.draw_handler_add(self.drawHandler, \
            (), "WINDOW", "POST_VIEW")
        self.batch = None
        self.batch2 = None

    def cleanup(self, context):
        if(self.drawHandlerRef != None):
            bpy.types.SpaceView3D.draw_handler_remove(self.drawHandlerRef, "WINDOW")
            if(context.area and hasattr(context.space_data, 'region_3d')):
                context.area.tag_redraw()
            self.drawHandlerRef = None

    def confirm(self, context):
        self.save(context)
        self.curvePts = []
        self.capture = False
        self.cleanup(context)
        self.initialize()

    def modal(self, context, event):

        if(context.space_data == None):
            return {'PASS_THROUGH'}
            
        if(event.type == 'WINDOW_DEACTIVATE' and event.value == 'PRESS'):
            self.ctrl = False
            self.shift = False
            self.alt = False
            return {'PASS_THROUGH'}

        if(event.type == 'RET'):
            self.confirm(context)   #subclass
            return {'RUNNING_MODAL'}

        if(event.type == 'ESC'):
            self.cleanup(context)
            self.initialize()
            return {'PASS_THROUGH'}

        if(isOutside(context, event)):
            if(event.value == 'PRESS'):
                self.curvePts = []
                self.refreshAfterMove(context, event)
            return {'PASS_THROUGH'}

        if(event.type == 'BACK_SPACE' and event.value == 'RELEASE'):
            if(len(self.curvePts) > 1):
                self.curvePts.pop()
                self.refreshAfterMove(context, event)
            else:
                self.curvePts = []
                self.capture = False

        if (event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            if(len(self.curvePts) == 0):
                loc = self.get3dLocWithSnap(context, event)
                self.curvePts.append([loc, loc, loc])
            self.capture = True
            return {'PASS_THROUGH'}

        if(event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}):
            self.ctrl = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}):
            self.shift = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type in {'LEFT_ALT', 'RIGHT_ALT'}):
            self.alt = (event.value == 'PRESS')
            return {'RUNNING_MODAL'}

        if(event.type == 'WHEELDOWNMOUSE' and self.alt):
            if(self.snapSteps < 10):
                self.snapSteps += 1
            return {'RUNNING_MODAL'}

        if(event.type == 'WHEELUPMOUSE' and self.alt):
            if(self.snapSteps > 1):
                self.snapSteps -= 1
            return {'RUNNING_MODAL'}

        if(event.type == 'MIDDLEMOUSE' and self.alt):
            self.snapSteps = self.defaultSnapSteps
            return {'RUNNING_MODAL'}

        if (event.type == 'LEFTMOUSE' and event.value == 'RELEASE'):

            #Looks like no 'DOUBLE_CLICK' event?
            t = time.time()
            if(self.clickT !=  None):
                if((t - self.clickT) < 0.25):
                    self.confirm(context)
                    return {'RUNNING_MODAL'}
            self.clickT = t

            if(len(self.curvePts) == 1): # len can't be zero
                loc = self.get3dLocWithSnap(context, event)
                self.curvePts[0][2] = loc # changes only rt handle
            else:
                loc = self.get3dLocWithSnap(context, event)

            self.curvePts.append([loc, loc, loc])
            self.createBatch(context, event, handlePtIdx = -2, \
                addHandle = self.addHandle)
            self.capture = False
            return {'RUNNING_MODAL'}

        if (event.type == 'MOUSEMOVE'):
            self.refreshAfterMove(context, event)
            if(self.capture):
                return {'RUNNING_MODAL'}
            else:
                return {'PASS_THROUGH'}

        return {'PASS_THROUGH'}

    def refreshAfterMove(self, context, event):
        loc = self.get3dLocWithSnap(context, event)
        if(not self.capture and len(self.curvePts) > 0):
            self.curvePts[-1] = [loc, loc, loc]
            handlePtIdx = -2
        else:
            if(len(self.curvePts) > 1):
                end = self.curvePts[-1][1]
                ltHandle = end - (loc - end)
                rtHandle = end + (loc - end)
                self.curvePts[-1][0] = ltHandle
                self.curvePts[-1][2] = rtHandle
            handlePtIdx  = -1
        self.createBatch(context, event, handlePtIdx = handlePtIdx, \
            addHandle = self.addHandle)

    def getLinesFromPts(self, pts):
        positions = []
        for i, pt in enumerate(pts):
            positions.append(pt)
            if(i > 0 and i < (len(pts)-1)):
                positions.append(pt)
        return positions


    def getPositions(self, context, event, handlePtIdx, addHandle):
        if(len(self.curvePts) <= 1):
            loc = self.get3dLocWithSnap(context, event)
            if(len(self.curvePts) == 1):
                #Draw straight line (no change in pts)
                positions = [self.curvePts[0][1], loc]
            else:
                positions = [loc, loc] #Will draw the dot twice.. ok for now.
        else:
            pts = getPtsAlongBezier(context, self.curvePts,
                self.curveDispRes, normalized = True)

            positions = self.getLinesFromPts(pts)
            if(addHandle == True):
                positions += [self.curvePts[handlePtIdx][0],
                    self.curvePts[handlePtIdx][2]]

        return positions

    def createBatch(self, context, event, handlePtIdx, addHandle):
        positions = self.getPositions(context, event, handlePtIdx, addHandle)
        finishedColor = (1, 0, 0, 1)
        handleColor = (0, 1, 1, 1)
        tipColor = (.8, 1, 0, 1)

        colors = [finishedColor for i in range(0, len(positions) - 1)]
        if(addHandle):
            colors.append(handleColor)
        else:
            colors.append(finishedColor)

        self.batch = batch_for_shader(self.shader, \
            "LINES", {"pos": positions, "color": colors})

        if(addHandle and len(positions) >= 2):
            pos2 = [positions[-1], positions[-2]]
            self.batch2 = batch_for_shader(self.shader, \
                "POINTS", {"pos": pos2, "color": [tipColor, tipColor]})

        if context.area:
            context.area.tag_redraw()

    def get3dLocWithSnap(self, context, event):
        return self.get3dLoc(context, event, snapToObj = self.alt,
            snapToGrid = self.ctrl, restrict = self.shift, vec = None, fromActiveObj = True)

    def get3dLoc(self, context, event, snapToObj,
        snapToGrid, restrict, vec = None, fromActiveObj = True):

        rounding = getViewDistRounding(context)
        region = context.region
        rv3d = context.space_data.region_3d
        xy = event.mouse_region_x, event.mouse_region_y

        self.axis = None

        if(snapToObj):
            minDist = 9e+99
            matchSnapLoc = None
            for snapLoc in self.getSnapLocs():
                loc = region_2d_to_location_3d(region, rv3d, xy, snapLoc)
                dist = (loc - snapLoc).length
                if(dist < (10 ** rounding)):
                    if(dist < minDist):
                        minDist = dist
                        matchSnapLoc = snapLoc
            if(matchSnapLoc != None):
                return matchSnapLoc

        if(vec == None):
            if(fromActiveObj and context.active_object != None):
                vec = context.active_object.location
            else:
                vec = region_2d_to_vector_3d(region, rv3d, xy)

        loc = region_2d_to_location_3d(region, rv3d, xy, vec)

        if(snapToGrid):
            loc = roundedVect(loc, rounding)

        if(restrict and len(self.curvePts) > 0):
            if(len(self.curvePts) == 1 and not self.capture):
                return loc

            actualLoc = loc.copy()
            lastCo = self.curvePts[-1][1] if(self.capture) else self.curvePts[-2][1]

            #First decide the main movement axis
            if(self.axis == None):
                diff = [abs(v) for v in (actualLoc - lastCo)]
                maxDiff = max(diff)
                self.axis = 0 if abs(diff[0]) == maxDiff \
                    else (1 if abs(diff[1]) == maxDiff else 2)

            loc = lastCo.copy()
            loc[self.axis] = actualLoc[self.axis]

            snapIncr = 45 / self.snapSteps
            snapAngles = [radians(snapIncr * a) for a in range(0, self.snapSteps + 1)]
            l1 =  actualLoc[self.axis] - lastCo[self.axis] #Main axis value

            for i in range(0, 3):
                if(i != self.axis):
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

    def drawHandler(self):
        if(not hasattr(self, 'batch') or self.batch == None):
            return

        bgl.glLineWidth(self.defLineWidth)
        bgl.glPointSize(self.defPointSize)

        self.batch.draw(self.shader)

        if(hasattr(self, 'batch2') and self.batch2 != None):
            self.batch2.draw(self.shader)


class ModalFlexiBezierOp(ModalDrawBezierOp):
    bl_description = "Draw Bezier curves by manipulating end point handles"
    bl_idname = "wm.draw_flexi_bezier_curves"
    bl_label = "Draw Bezier Curves"
    bl_options = {'REGISTER', 'UNDO'}

    running = False

    def __init__(self):
        curveDispRes = 100
        defLineWidth = 1.5
        defPointSize = 7
        super(ModalFlexiBezierOp, self).__init__(curveDispRes, defLineWidth, defPointSize)

    @classmethod
    def poll(cls, context):
        return not ModalFlexiBezierOp.running

    def isValidContext(self, context):
        if(context.mode != 'OBJECT'):
            return False

        tool = context.workspace.tools.from_space_view3d_mode('OBJECT', create = False)

        if(tool == None or tool.idname != DrawFlexiBezierTool.bl_idname):
            return False
            
        return True
    
    def cancelOp(self, context):
        self.cleanup(context)
        ModalFlexiBezierOp.running = False
        return {"CANCELLED"}

    def invoke(self, context, event):
        if(ModalFlexiBezierOp.running):
            return {"CANCELLED"}

        if(not self.isValidContext(context)):
            return self.cancelOp(context)

        ModalFlexiBezierOp.running = True

        #Object name -> [spline index, (startpt, endPt)]
        self.snapInfos = {}

        self.updateSnapPts([o.name for o in bpy.data.objects])
        return super(ModalFlexiBezierOp, self).invoke(context, event)

    def modal(self, context, event):
        if(not self.isValidContext(context)):
            return self.cancelOp(context)

        return super(ModalFlexiBezierOp, self).modal(context, event)

    def getSnapLocs(self):
        locs = []
        infos = [info for values in self.snapInfos.values() for info in values]
        for info in infos:
            locs += info[1]

        if(len(self.curvePts) > 0):
            locs.append(self.curvePts[0][1])

        return locs

    def updateSnapPts(self, objNames):
        for objName in objNames:
            obj = bpy.data.objects.get(objName)
            if(obj == None or not isBezier(obj) or not obj.visible_get()):
                if(self.snapInfos.get(objName) != None):
                    del self.snapInfos[objName] #Just in case
                continue

            self.snapInfos[objName] = []
            mw = obj.matrix_world
            for i, s in enumerate(obj.data.splines):
                if(s.use_cyclic_u == True):
                    continue
                ends = [mw @ s.bezier_points[0].co, mw @ s.bezier_points[-1].co]
                self.snapInfos[objName].append([i, ends])

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
        spline = data.splines.new('BEZIER')
        spline.use_cyclic_u = False

        if(vectCmpWithMargin(self.curvePts[0][1], self.curvePts[-1][0])):
            spline.use_cyclic_u = True
            self.curvePts.pop()

        spline.bezier_points.add(len(self.curvePts) - 1)
        for i, pt in enumerate(self.curvePts):
            spline.bezier_points[i].handle_left = pt[0]
            spline.bezier_points[i].co = pt[1]
            spline.bezier_points[i].handle_right = pt[2]
            spline.bezier_points[i].handle_left_type = 'ALIGNED'
            spline.bezier_points[i].handle_right_type = 'ALIGNED'

        depsgraph = context.evaluated_depsgraph_get()
        depsgraph.update()
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

    def createCurveObj(self, context, startObj, startSplineIdx, endObj, endSplineIdx):
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
                # but it they are different they need to be joined with a separate call
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

    def save(self, context):
        if(len(self.curvePts) > 0):
            self.curvePts.pop()

        if(len(self.curvePts) > 1):
            endObj, endSplineIdx, ptIdx1 = \
                self.getSnapObj(context, self.curvePts[-1][1])
            startObj, startSplineIdx, ptIdx2 = \
                self.getSnapObj(context, self.curvePts[0][1])

            startObjName = startObj.name if(startObj != None) else ''
            endObjName = endObj.name if(endObj != None) else ''

            obj = self.createCurveObj(context, \
                startObj, startSplineIdx, endObj, endSplineIdx)

            #TODO: Why try?
            try:
                obj.select_set(True)
                bpy.context.view_layer.objects.active = obj
                self.updateSnapPts([obj.name, startObjName, endObjName])
            except Exception as e:
                pass
        bpy.ops.ed.undo_push()


class DrawFlexiBezierTool(WorkSpaceTool):
    bl_space_type='VIEW_3D'
    bl_context_mode='OBJECT'

    bl_idname = "flexi_bezier.draw_tool"
    bl_label = "Flexi Bezier"
    bl_description = ("Draw flexi Bezier curve")
    bl_icon = "ops.gpencil.extrude_move"
    bl_widget = None
    bl_operator = "wm.draw_flexi_bezier_curves"
    bl_keymap = (
        ("wm.draw_flexi_bezier_curves", {"type": 'MOUSEMOVE', "value": 'ANY'},
         {"properties": []}),
    )

########## Draw Flexi Bezier Curve - End ###############


def register():
    bpy.utils.register_class(ModalMarkSegStartOp)
    bpy.utils.register_class(SeparateSplinesObjsOp)
    bpy.utils.register_class(SplitBezierObjsOp)
    bpy.utils.register_class(splitBezierObjsPtsOp)
    bpy.utils.register_class(JoinBezierSegsOp)
    bpy.utils.register_class(CloseSplinesOp)
    bpy.utils.register_class(CloseStraightOp)
    bpy.utils.register_class(OpenSplinesOp)
    bpy.utils.register_class(RemoveDupliVertCurveOp)
    bpy.utils.register_class(convertTo2DMeshOp)
    bpy.utils.register_class(SelectInCollOp)
    bpy.utils.register_class(InvertSelOp)
    bpy.utils.register_class(BezierUtilsPanel)
    bpy.utils.register_class(ModalFlexiBezierOp)
    bpy.utils.register_tool(DrawFlexiBezierTool)

def unregister():
    bpy.utils.unregister_class(ModalMarkSegStartOp)
    bpy.utils.unregister_class(SeparateSplinesObjsOp)
    bpy.utils.unregister_class(SplitBezierObjsOp)
    bpy.utils.unregister_class(splitBezierObjsPtsOp)
    bpy.utils.unregister_class(JoinBezierSegsOp)
    bpy.utils.unregister_class(CloseSplinesOp)
    bpy.utils.unregister_class(CloseStraightOp)
    bpy.utils.unregister_class(OpenSplinesOp)
    bpy.utils.unregister_class(RemoveDupliVertCurveOp)
    bpy.utils.unregister_class(convertTo2DMeshOp)
    bpy.utils.unregister_class(SelectInCollOp)
    bpy.utils.unregister_class(InvertSelOp)
    bpy.utils.unregister_class(BezierUtilsPanel)
    bpy.utils.unregister_class(ModalFlexiBezierOp)
    bpy.utils.unregister_tool(DrawFlexiBezierTool)
