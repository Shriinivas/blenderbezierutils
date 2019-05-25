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

import bpy
from bpy.types import Panel

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


class BezierUtilsPanel(Panel):    
    bl_label = "Bezier Utilities"
    bl_idname = "CURVE_PT_bezierutils"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Tool'    
    bl_context = '.objectmode'

    def draw(self, context):
        layout = self.layout
        col = layout.column()
        col.operator("object.separate_splines")
        col = layout.column()
        col.operator("object.separate_segments")

def register():
    bpy.utils.register_class(SeparateSplinesObjsOp)
    bpy.utils.register_class(SplitBezierObjsOp)
    bpy.utils.register_class(BezierUtilsPanel)
    
def unregister():
    bpy.utils.unregister_class(SeparateSplinesObjsOp)
    bpy.utils.unregister_class(SplitBezierObjsOp)
    bpy.utils.unregister_class(BezierUtilsPanel)
    
