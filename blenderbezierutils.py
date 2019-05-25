#
#
# This Blender add-on with Bezier Curve utility operations 
#
# Supported Blender Version: 2.8 Beta
#
# Copyright (C) 2019  Shrinivas Kulkarni

# License: MIT (https://github.com/Shriinivas/blenderbezierutils/blob/master/LICENSE)
#
# Not yet pep8 compliant 

import bpy, bmesh
from bpy.props import BoolProperty
from bpy.types import Panel
from mathutils import Vector

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "location": "Properties > Active Tool and Workspace Settings > Bezier Utilities",
    "category": "Object",
    "blender": (2, 80, 0),    
}

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
            
def createSkeletalCurve(objGrp, obj):
    objCopy = obj.copy()
    objCopy.name = obj.name 
    dataCopy = obj.data.copy()
    dataCopy.splines.clear()
    objCopy.data = dataCopy
    objGrp.objects.link(objCopy)
    return objCopy

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
        
    keyblocks = reversed(obj.data.shape_keys.key_blocks)
    for sk in keyblocks:
        obj.shape_key_remove(sk)
    
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
    segCnt = 0
    
    if(len(selObjs) == 0):
        self.report({'WARNING'}, "No Curve Objects Selected")
        return None, None, None
    
    for obj in selObjs:

        if(not isBezier(obj) or (len(obj.data.splines) <= 1 and not splitSegs) or \
            (len(obj.data.splines) == 1 and len(obj.data.splines[0].bezier_points) == 2)):
                
            continue

        keyNames, keyData = getKeyInfo(obj)
        collections = obj.users_collection
        objGrp = bpy.data.collections.new(obj.name)
        
        for i, spline in enumerate(obj.data.splines):
            if(splitSegs):
                for j in range(0, len(spline.bezier_points) - 1):
                    objCopy = createSkeletalCurve(objGrp, obj)
                    createSplineForSeg(objCopy.data, \
                        spline.bezier_points[j:j+2])
                    updateKeyData(objCopy, keyNames, keyData, segCnt, 2)
                    segCnt += 1
                if(spline.use_cyclic_u):
                    objCopy = createSkeletalCurve(objGrp, obj)
                    createSplineForSeg(objCopy.data, \
                        [spline.bezier_points[-1], spline.bezier_points[0]])
                    updateKeyData(objCopy, keyNames, keyData, -1, 2)
                    segCnt += 1
            else:
                objCopy = createSkeletalCurve(objGrp, obj)
                createSpline(objCopy.data, spline, forceNoncyclic = False, \
                    freeHandles = False)
                currSegCnt = len(objCopy.data.splines[0].bezier_points)
                updateKeyData(objCopy, keyNames, keyData, segCnt, currSegCnt)
                segCnt += currSegCnt 

        for collection in collections:
            collection.children.link(objGrp)
            collection.objects.unlink(obj)
            
        bpy.data.curves.remove(obj.data)
        changeCnt += 1
        splineCnt += (i + 1)
        
    return changeCnt, splineCnt, segCnt

class SeparateSplinesObjsOp(bpy.types.Operator):

    bl_idname = "object.separate_splines"
    bl_label = "Separate Bezier Splines"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        
        changeCnt, splineCnt, segCnt = splitCurve(context, splitSegs = False)

        if(changeCnt != None):
            self.report({'INFO'}, "Separated "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + " into " +str(splineCnt) + " new ones")
        
        return {'FINISHED'}
        
class SplitBezierObjsOp(bpy.types.Operator):

    bl_idname = "object.separate_segments"
    bl_label = "Split Bezier Segments"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        
        changeCnt, splineCnt, segCnt = splitCurve(context, splitSegs = True)
        
        if(changeCnt != None):
            self.report({'INFO'}, "Split "+ str(changeCnt) + " curve object" + \
                ("s" if(changeCnt > 1) else "") + " into " + str(segCnt) + " new objects")
        
        return {'FINISHED'}

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

def markVertHandler(self, context):
    if(self.markVertex):
        bpy.ops.wm.mark_vertex()

class ModalMarkSegStartOp(bpy.types.Operator):
    bl_description = "Mark Vertex"
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
    
    bpy.types.Scene.markVertex = BoolProperty(name="Mark Starting Vertices", \
        description='Mark first vertices in all closed splines of selected curves', \
            default = False, update = markVertHandler)
            
    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_CURVE'} 

    def draw(self, context):
        layout = self.layout
        if(context.mode  == 'OBJECT'):
            col = layout.column()
            col.operator("object.separate_splines")
            col = layout.column()
            col.operator("object.separate_segments")
        else:
            col = layout.column()
            col.prop(context.scene, "markVertex", toggle = True)
            
def register():
    bpy.utils.register_class(ModalMarkSegStartOp)
    bpy.utils.register_class(SeparateSplinesObjsOp)
    bpy.utils.register_class(SplitBezierObjsOp)
    bpy.utils.register_class(BezierUtilsPanel)
    
def unregister():
    bpy.utils.unregister_class(ModalMarkSegStartOp)
    bpy.utils.unregister_class(SeparateSplinesObjsOp)
    bpy.utils.unregister_class(SplitBezierObjsOp)
    bpy.utils.unregister_class(BezierUtilsPanel)
    
