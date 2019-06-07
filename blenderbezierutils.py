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

import bpy, bmesh
from bpy.props import BoolProperty, IntProperty
from bpy.types import Panel
from mathutils import Vector, Matrix, geometry

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "location": "Properties > Active Tool and Workspace Settings > Bezier Utilities",
    "category": "Object",
    "blender": (2, 80, 0),    
}

DEF_ERR_MARGIN = 0.0001 #rather high

def floatCmpWithMargin(float1, float2, margin = DEF_ERR_MARGIN):
    return abs(float1 - float2) < margin 

def vectCmpWithMargin(v1, v2, margin = DEF_ERR_MARGIN):
    return all(floatCmpWithMargin(v1[i], v2[i]) for i in range(0, len(v1)))

class MarkerController:
    scale = Vector([1, 1, 1])
    
    def createMarkers(self, context):
        objs = context.selected_objects
        smMap = {}
        for curve in objs:
            if(not isBezier(curve)):
                continue
                
            smMap[curve.name] = {}
            for splineIdx, spline in enumerate(curve.data.splines):
                if(not spline.use_cyclic_u):
                    continue
                    
                #TODO: Maybe share the data
                bm = bmesh.new()
                bmesh.ops.create_uvsphere(bm, u_segments=10, v_segments=10, diameter = .05)
                d = bpy.data.meshes.new('___tmp')
                bm.to_mesh(d)
                marker = bpy.data.objects.new('___tmp', d)
                bpy.context.scene.collection.objects.link(marker)
                marker.show_in_front = True
                marker.scale = MarkerController.scale
                smMap[curve.name][splineIdx] = [marker, 0]#marker and curr start vert idx
                
                for pt in spline.bezier_points:
                    pt.select_control_point = False
                    
            if(len(smMap[curve.name]) == 0):
                del smMap[curve.name]
                
        return smMap
        
    def removeMarkers(self):
        for spMap in self.smMap.values():
            for markerInfo in spMap.values():
                marker, idx = markerInfo[0], markerInfo[1]
                safeRemoveObj(marker)
        
    def __init__(self, context):
        self.smMap = self.createMarkers(context)           
        self.moveMarkersToVerts()
    
    def saveStartVerts(self):        
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            splines = curve.data.splines
            spMap = self.smMap[curveName]
            
            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                if(markerInfo[1] != 0):
                    pts = splines[splineIdx].bezier_points
                    marker, idx = markerInfo[0], markerInfo[1]
                    cnt = len(pts)
                    
                    ptCopy = [[p.co.copy(), p.handle_right.copy(), \
                        p.handle_left.copy(), p.handle_right_type, p.handle_left_type] \
                            for p in pts]

                    for i, pt in enumerate(pts):
                        srcIdx = (idx + i) % cnt
                        p = ptCopy[srcIdx]
                        #Must set the types first
                        pt.handle_right_type = p[3]
                        pt.handle_left_type = p[4]
                        pt.co = p[0]
                        pt.handle_right = p[1]
                        pt.handle_left = p[2]
                            
    def moveMarkersToVerts(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            spMap = self.smMap[curveName]
            mw = curve.matrix_world
            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                marker, idx = markerInfo[0], markerInfo[1]
                pts = curve.data.splines[splineIdx].bezier_points
                selIdxs = [x for x in range(0, len(pts)) \
                    if pts[x].select_control_point == True]        
                if(len(selIdxs) > 0 ):
                    selIdx = selIdxs[0]
                    markerInfo[1] = selIdx
                    pt = pts[selIdx] 
                else:
                    pt = pts[idx]
                co = mw @ pt.co
                if(marker.location != co):
                    marker.location = co
                    
    def scaleMarker(self):
        markers = [m[0] for spMap in self.smMap.values() for m in spMap.values()]
        if(len(markers) == 0):
            return
        for m in markers:
            m.scale = MarkerController.scale
        
    def getSpaces3D(context):
        areas3d  = [area for area in context.window.screen.areas if area.type == 'VIEW_3D']
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


def isBezier(bObj):
    return bObj.type == 'CURVE' and len(bObj.data.splines) > 0 \
        and bObj.data.splines[0].type == 'BEZIER'

def safeRemoveObj(obj):
    try:
        collections = obj.users_collection

        for c in collections:
            c.objects.unlink(obj)
            
        if(obj.data.users == 1):
            bpy.data.meshes.remove(obj.data) #This also removes object?        
        else:
            bpy.data.objects.remove(obj)
    except:
        pass

def markVertHandler(self, context):
    if(self.markVertex):
        bpy.ops.wm.mark_vertex()

def copyBezierPt(src, target, freeHandles):
    target.co = src.co
    target.handle_left = src.handle_left
    target.handle_right = src.handle_right
    if(freeHandles):
        target.handle_left_type = 'FREE'
        target.handle_right_type = 'FREE'
    else:
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
            
def createSkeletalCurve(objGrp, obj):
    objCopy = obj.copy()
    objCopy.name = obj.name 
    dataCopy = obj.data.copy()
    dataCopy.splines.clear()
    objCopy.data = dataCopy
    objGrp.objects.link(objCopy)
    return objCopy

def removeShapeKeys(obj):
    keyblocks = reversed(obj.data.shape_keys.key_blocks)
    for sk in keyblocks:
        obj.shape_key_remove(sk)
    
def getKeyInfo(obj):
    keyData = []
    keyNames = []

    if(obj.data.shape_keys != None):
        keyblocks = obj.data.shape_keys.key_blocks
        for key in keyblocks:
            keyData.append([[d.co, d.handle_left, d.handle_right] for d in key.data])
            keyNames.append(key.name)
    return keyNames, keyData
        
def updateKeyData(obj, keyNames, keyData, retainIdxStart, retainCnt):
    if(obj.data.shape_keys == None):
        return 
        
    removeShapeKeys(obj)
    
    for i, name in enumerate(keyNames):
        key = obj.shape_key_add(name = name)
        for j in range(0, retainCnt):
            key.data[j].co = keyData[i][j + retainIdxStart][0]
            key.data[j].handle_left = keyData[i][j + retainIdxStart][1]
            key.data[j].handle_right = keyData[i][j + retainIdxStart][2]

def splitCurve(context, splitSegs):
    selObjs = bpy.context.selected_objects
    changeCnt = 0
    splineCnt = 0
    newObjs = []

    if(len(selObjs) == 0):
        self.report({'WARNING'}, "No Curve Objects Selected")
        return newObjs, changeCnt
    
    for obj in selObjs:

        if(not isBezier(obj) or (len(obj.data.splines) <= 1 and not splitSegs) or \
            (len(obj.data.splines) == 1 and len(obj.data.splines[0].bezier_points) == 2)):
            continue

        keyNames, keyData = getKeyInfo(obj)
        collections = obj.users_collection
        objGrp = bpy.data.collections.new(obj.name)
        segCnt = 0
        
        for i, spline in enumerate(obj.data.splines):
            if(splitSegs):
                for j in range(0, len(spline.bezier_points) - 1):
                    objCopy = createSkeletalCurve(objGrp, obj)
                    createSplineForSeg(objCopy.data, \
                        spline.bezier_points[j:j+2])
                    updateKeyData(objCopy, keyNames, keyData, len(newObjs), 2)
                    newObjs.append(objCopy)
                if(spline.use_cyclic_u):
                    objCopy = createSkeletalCurve(objGrp, obj)
                    createSplineForSeg(objCopy.data, \
                        [spline.bezier_points[-1], spline.bezier_points[0]])
                    updateKeyData(objCopy, keyNames, keyData, -1, 2)
                    newObjs.append(objCopy)
            else:
                objCopy = createSkeletalCurve(objGrp, obj)
                createSpline(objCopy.data, spline, forceNoncyclic = False, \
                    freeHandles = False)
                currSegCnt = len(objCopy.data.splines[0].bezier_points)
                updateKeyData(objCopy, keyNames, keyData, segCnt, currSegCnt)
                newObjs.append(objCopy)
                segCnt += currSegCnt

        for collection in collections:
            collection.children.link(objGrp)
            collection.objects.unlink(obj)
            
        bpy.data.curves.remove(obj.data)
        changeCnt += 1

    return newObjs, changeCnt

#TODO: Fix this hack if possible
def copyObjAttr(src, dest, invSrcMW = Matrix(), mw = Matrix()):
    for att in dir(src):
        try:
            if(att not in ['co', 'handle_left', 'handle_right', 'handle_left_type', 'handle_right_type']):
                setattr(dest, att, getattr(src, att))
        except Exception as e:            
            pass
    try:
        lt = src.handle_left_type
        rt = src.handle_right_type
        dest.handle_left_type = 'FREE'
        dest.handle_right_type = 'FREE'
        dest.co = invSrcMW @ (mw @ src.co)
        dest.handle_left = invSrcMW @ (mw @ src.handle_left)
        dest.handle_right = invSrcMW @ (mw @ src.handle_right)
        dest.handle_left_type = lt
        dest.handle_right_type = rt
        pass
    except Exception as e:
        pass
            
def reverseCurve(curve):
    cp = curve.data.copy()
    curve.data.splines.clear()
    for s in reversed(cp.splines):
        ns = curve.data.splines.new('BEZIER')
        copyObjAttr(s, ns)
        ns.bezier_points.add(len(s.bezier_points) - 1)
        for i, p in enumerate(reversed(s.bezier_points)):
            copyObjAttr(p, ns.bezier_points[i])
            ns.bezier_points[i].handle_left = p.handle_right
            ns.bezier_points[i].handle_left_type = p.handle_right_type
            ns.bezier_points[i].handle_right = p.handle_left
            ns.bezier_points[i].handle_right_type = p.handle_left_type
    bpy.data.curves.remove(cp)
    
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

def getArrangedCurves(curves):
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
    

def addLastSeg(spline):
    if(spline.use_cyclic_u):
        spline.bezier_points[0].handle_left_type = 'FREE'
        spline.bezier_points[0].handle_right_type = 'FREE'        
        spline.use_cyclic_u = False
        spline.bezier_points.add(1)
        copyObjAttr(spline.bezier_points[0], spline.bezier_points[-1])
        
def closeSplines(curve, force = False):
    for spline in curve.data.splines:
        if(not spline.use_cyclic_u or force):
            spline.bezier_points[0].handle_left_type = 'FREE'
            spline.bezier_points[0].handle_left = spline.bezier_points[0].co
            spline.bezier_points[-1].handle_right_type = 'FREE'
            spline.bezier_points[-1].handle_right = spline.bezier_points[-1].co
            spline.use_cyclic_u = True

def cmpPts(mw1, pt1, mw2, pt2, cmpMargin = 0.1):
    return (vectCmpWithMargin(mw1 @ pt1.co, mw2 @ pt2.co, cmpMargin) and
        vectCmpWithMargin(mw1 @ pt1.handle_right, mw2 @ pt2.handle_right, cmpMargin) and 
            vectCmpWithMargin(mw1 @ pt1.handle_left, mw2 @ pt2.handle_left, cmpMargin))
    
def joinSegs(curves, optimized, straight):

    if(len(curves) == 0):
        return None
    if(len(curves) == 1):    
        return curves[0]
    
    firstCurve = curves[0]
    
    if(optimized):
        curves = getArrangedCurves(curves)
    
    newCurve = firstCurve.copy()
    newCurveData = firstCurve.data.copy()
    newCurve.data = newCurveData
    
    collections = firstCurve.users_collection
    for collection in collections:
        collection.objects.link(newCurve)
        
    srcMW = firstCurve.matrix_world
    invSrcMW = srcMW.inverted()

    for curve in curves[1:]:
        mw = curve.matrix_world
        
        currSpline = newCurveData.splines[-1]
        nextSpline = curve.data.splines[0]

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
            currBezierPt.handle_right_type = 'FREE'
            currBezierPt.handle_right = currBezierPt.co
        
        for i in range(ptIdx, len(nextSpline.bezier_points)):
            if((i == len(nextSpline.bezier_points) - 1) and 
                vectCmpWithMargin(mw @ nextSpline.bezier_points[i].co, srcMW @ currSpline.bezier_points[0].co)):
                    currSpline.bezier_points[0].handle_right_type = 'FREE'
                    currSpline.bezier_points[0].handle_right = invSrcMW @ (mw @ nextSpline.bezier_points[i].handle_right)
                    currSpline.use_cyclic_u = True
                    break
            currSpline.bezier_points.add(1)
            currBezierPt = currSpline.bezier_points[-1]
            copyObjAttr(nextSpline.bezier_points[i], currBezierPt, invSrcMW, mw)
            if(straight and i ==  0):
                currBezierPt.handle_left_type = 'FREE'
                currBezierPt.handle_left = currBezierPt.co

        #Simply add the remaining splines
        for spline in curve.data.splines[1:]:
            newSpline = newCurveData.splines.new('BEZIER')
            copyObjAttr(spline, newSpline)
            for i, pt in enumerate(spline.bezier_points):
                if(i > 0):
                    newSpline.bezier_points.add(1)
                copyObjAttr(pt, newSpline.bezier_points[-1], invSrcMW, mw)
    
        safeRemoveObj(curve)
    safeRemoveObj(firstCurve)
    return newCurve

# ~ def scaledBezierData(curveData, scale):
    # ~ newData = curveData.copy()    
    # ~ for spline in newData.splines:
        # ~ for pt in spline.bezier_points:
            # ~ pt.co = [scale[i] * pt.co[i] for i in range(0, 3)]
            # ~ pt.handle_left = [scale[i] * pt.handle_left[i] for i in range(0, 3)]
            # ~ pt.handle_right = [scale[i] * pt.handle_right[i] for i in range(0, 3)]
    # ~ return newData
        
# ~ def transformedBezierData(curveData, mat):
    # ~ newData = curveData.copy()    
    # ~ for spline in newData.splines:
        # ~ for pt in spline.bezier_points:
            # ~ pt.co = mat @ pt.co
            # ~ pt.handle_left = mat @ pt.handle_left
            # ~ pt.handle_right = mat @ pt.handle_right
    # ~ return newData
        
# ~ def copyLength(fromCurve, toCurves):
    # ~ srcMW = fromCurve.matrix_world
    # ~ srcLen = sum(s.calc_length() for s in transformedBezierData(fromCurve.data, srcMW).splines)
    # ~ for curve in toCurves:        
        # ~ mw = curve.matrix_world
        # ~ destLen = sum(s.calc_length() for s in transformedBezierData(curve.data, mw).splines)
        # ~ fact = srcLen / destLen
        # ~ curve.scale *= fact
        
def removeDupli(curve):
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
    
class SeparateSplinesObjsOp(bpy.types.Operator):

    bl_idname = "object.separate_splines"
    bl_label = "Separate Bezier Splines"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate splines of selected Bezier curves as new objects"

    def execute(self, context):
        
        
        newObjs, changeCnt = splitCurve(context, splitSegs = False)

        if(changeCnt > 0):
            self.report({'INFO'}, "Separated "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + " into " +str(len(newObjs)) + " new ones")
        
        return {'FINISHED'}
        
class SplitBezierObjsOp(bpy.types.Operator):

    bl_idname = "object.separate_segments"
    bl_label = "Split Bezier Segments"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Separate segments of selected Bezier curves as new objects"

    def execute(self, context):
        newObjs, changeCnt = splitCurve(context, splitSegs = True)

        bpy.context.view_layer.objects.active = newObjs[-1]
        self.report({'INFO'}, "Split "+ str(changeCnt) + " curve object" + \
            ("s" if(changeCnt > 1) else "") + " into " + str(len(newObjs)) + " new objects")
        
        return {'FINISHED'}

class JoinBezierSegsOp(bpy.types.Operator):
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

class InvertSelOp(bpy.types.Operator):
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

class SelectInCollOp(bpy.types.Operator):
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

class CloseStraightOp(bpy.types.Operator):
    bl_idname = "object.close_straight"
    bl_label = "Close With Straight Segment"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Close selected curves with straight segmennt"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]
        
        for curve in curves:
            closeSplines(curve, force = True)
        
        bpy.context.view_layer.objects.active = curve
            
        return {'FINISHED'}


class RemoveDupliVertCurveOp(bpy.types.Operator):
    bl_idname = "object.remove_dupli_vert_curve"
    bl_label = "Remove Duplicate Curve Vertices"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Remove duplicate vertices and mark splines as cyclic if applicable"

    def execute(self, context):
        curves = [o for o in bpy.data.objects \
            if o in bpy.context.selected_objects and isBezier(o)]
        
        for curve in curves:
            removeDupli(curve)
        
        bpy.context.view_layer.objects.active = curve
            
        return {'FINISHED'}

class convertTo2DMeshOp(bpy.types.Operator):
    bl_idname = "object.convert_2d_mesh"
    bl_label = "Convert"
    bl_description = "Convert 2D curve to mesh with quad faces"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        curve = context.object
        if(curve != None and isBezier(curve) and curve.data.dimensions == '2D'):
            closeSplines(curve, force = False)            
            curve.data.fill_mode = 'BOTH'    
            meshObj = convertToMesh(curve)
            remeshDepth = bpy.context.scene.remeshDepth            
            applyMeshModifiers(meshObj, remeshDepth)
            meshObj.matrix_world = curve.matrix_world
            safeRemoveObj(curve)
            bpy.context.view_layer.objects.active = meshObj

        return {'FINISHED'}

class ModalMarkSegStartOp(bpy.types.Operator):
    bl_description = "Mark starting vertex of all the closed splines in the selection"
    bl_idname = "wm.mark_vertex"
    bl_label = "Mark Start Vertex"

    def cleanup(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        self.markerState.removeMarkers()
        MarkerController.resetShowHandleState(context, self.handleStates)
        context.scene.markVertex = False

    def modal (self, context, event):
        
        if(context.mode  == 'OBJECT' or event.type == "ESC" or\
            not context.scene.markVertex):
            self.cleanup(context)
            return {'CANCELLED'}
        
        elif(event.type == "RET"):
            self.markerState.saveStartVerts()
            self.cleanup(context)
            return {'FINISHED'}
            
        if(event.type == 'TIMER'):
            self.markerState.moveMarkersToVerts()
        elif(event.type in {'LEFT_CTRL', 'RIGHT_CTRL'}):
            self.ctrl = (event.value == 'PRESS')
        elif(event.type in {'LEFT_SHIFT', 'RIGHT_SHIFT'}):
            self.shift = (event.value == 'PRESS')
        else:
            if(event.type not in {"MIDDLEMOUSE", "TAB", "LEFTMOUSE", \
                "RIGHTMOUSE", 'WHEELDOWNMOUSE', 'WHEELUPMOUSE'} and \
                not event.type.startswith("NUMPAD_")):
                return {'RUNNING_MODAL'}
            elif(event.type in {'WHEELDOWNMOUSE', 'WHEELUPMOUSE'} \
                and self.shift and self.ctrl):
                up = event.type == 'WHEELUPMOUSE'                    
                MarkerController.scale *= 1.1 if(up) else .9                    
                self.markerState.scaleMarker()
                return {'RUNNING_MODAL'}
            elif(event.type == 'MIDDLEMOUSE' and self.shift and self.ctrl):
                MarkerController.scale = Vector([1,1,1])
                self.markerState.scaleMarker()                 
                return {'RUNNING_MODAL'}
            
        return {"PASS_THROUGH"}

    def execute(self, context):
        self._timer = context.window_manager.event_timer_add(time_step = 0.01, \
            window = context.window)
        self.ctrl = False
        self.shift = False
        context.window_manager.modal_handler_add(self)
        self.markerState = MarkerController(context)
        self.handleStates = MarkerController.hideHandles(context)

        return {"RUNNING_MODAL"}

class BezierUtilsPanel(Panel):    
    bl_label = "Bezier Utilities"
    bl_idname = "CURVE_PT_bezierutils"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Tool'    
    
    #Mark Vert
    bpy.types.Scene.markVertex = BoolProperty(name="Mark Starting Vertices", \
        description='Mark first vertices in all closed splines of selected curves', \
            default = False, update = markVertHandler)
            
    #Split Seg
    bpy.types.Scene.selectIntrvl = IntProperty(name="Selection Interval", \
        description='Selection interval of the split segments', \
            default = 0, min = 0)
                
    #Convert to mesh
    bpy.types.Scene.remeshDepth = IntProperty(name="Remesh Depth", \
        description='Remesh depth for converting to mesh', \
            default = 4, min = 1, max = 10)
                
    ###Join segs
    bpy.types.Scene.straight = BoolProperty(name="Join With Straight Segments", \
        description='Join curves with straight segments', \
            default = True)
            
    bpy.types.Scene.optimized = BoolProperty(name="Join Optimized", \
        description='Join the nearest curve (reverse direction if necessary)', \
            default = True)
            
    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_CURVE'} 

    def draw(self, context):
        layout = self.layout
        # ~ layout.use_property_split = True
        
        if(context.mode  == 'OBJECT'):
            col = layout.column()
            col.operator('object.separate_splines', icon = 'UNLINKED')
            col = layout.column()
            col.operator('object.close_straight', icon = 'IPO_LINEAR')
            col = layout.column()
            col.operator('object.remove_dupli_vert_curve', icon = 'X' )
            col = layout.column()
            col.operator('object.separate_segments', icon = 'ORPHAN_DATA')
            
            col.separator()
            col.separator()
            col.label(text='Select Objects In Collection', icon='HAIR')            
            col = layout.column()
            row = col.row()
            row.prop(context.scene, 'selectIntrvl')
            row.operator('object.select_in_collection')
            col = layout.column()
            col.operator('object.invert_sel_in_collection')
            
            col.separator()
            col.separator()
            col.label(text='Convert to Mesh', icon='MESH_DATA')            
            row = col.row()
            row.prop(context.scene, 'remeshDepth')
            row.operator('object.convert_2d_mesh')
            
            col.separator()
            col.separator()
            col.label(text='Join Bezier Curves', icon='LINKED')            
            col = layout.column()
            col.prop(context.scene, 'straight')
            col = layout.column()
            col.prop(context.scene, 'optimized')
            col = layout.column()
            col.operator('object.join_curves')
        else:
            col = layout.column()
            col.prop(context.scene, 'markVertex', toggle = True)
            
def register():
    bpy.utils.register_class(ModalMarkSegStartOp)
    bpy.utils.register_class(SeparateSplinesObjsOp)
    bpy.utils.register_class(SplitBezierObjsOp)
    bpy.utils.register_class(JoinBezierSegsOp)
    bpy.utils.register_class(CloseStraightOp)
    bpy.utils.register_class(RemoveDupliVertCurveOp)
    bpy.utils.register_class(convertTo2DMeshOp)
    bpy.utils.register_class(SelectInCollOp)
    bpy.utils.register_class(InvertSelOp)
    bpy.utils.register_class(BezierUtilsPanel)
    
def unregister():
    bpy.utils.unregister_class(ModalMarkSegStartOp)
    bpy.utils.unregister_class(SeparateSplinesObjsOp)
    bpy.utils.unregister_class(SplitBezierObjsOp)
    bpy.utils.unregister_class(JoinBezierSegsOp)
    bpy.utils.unregister_class(CloseStraightOp)
    bpy.utils.unregister_class(RemoveDupliVertCurveOp)
    bpy.utils.unregister_class(convertTo2DMeshOp)
    bpy.utils.unregister_class(SelectInCollOp)
    bpy.utils.unregister_class(InvertSelOp)
    bpy.utils.unregister_class(BezierUtilsPanel)
    
