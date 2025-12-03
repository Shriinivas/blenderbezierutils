# bezier_utils/drawing/math_fn.py

import bpy
import os
from xml.dom import minidom
from shutil import copyfile
from bpy.types import Operator
from bpy.props import StringProperty
from math import sqrt
from .primitives import Primitive2DDraw
from ..utils.bezier_math import get3DVector
from ..constants import GP_CONTEXT_MODE

class MathFnDraw(Primitive2DDraw):
    
    # Prefixes for param names generated dynamically for constants
    startPrefix = 'mathFnStart_'
    incrPrefix = 'mathFnIncr_'
    
    mathFnFileExt = 'mfn'
    
    # Default values
    defFnType = 'XY'
    defFNRes = 10
    
    defFNXYName = ''
    defFNXYDescr = ''
    defFnXY = 'sin(x)'
    
    defFnParam1 = ''
    defFnParam2 = ''
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
        fnStr = '''def updateParam{0}(self, event, rmInfo, isIncr):
\tparams = bpy.context.window_manager.bezierToolkitParams
\tincr = params.{1}{0}
\tparams.{2}{0} += incr if(isIncr) else -incr
\tself.updateCurvePts()
\tself.parent.redrawBezier(rmInfo, hdlPtIdxs = {{}}, hltEndSeg = False)
\tareas = [a for a in bpy.context.screen.areas]
\tfor a in areas:
\t\ta.tag_redraw()
\treturn True'''.format(i, incrPrefix, startPrefix)
        exec(fnStr)
    
    def afterShapeSegCnt(self):
        bpy.context.window_manager.bezierToolkitParams.mathFnResolution = self.shapeSegCnt
        areas = [a for a in bpy.context.screen.areas]
        for a in areas:
            a.tag_redraw()
    
    def testFn(expr, var):
        # Use a shared namespace for exec and eval
        # Import math functions into the namespace for function evaluation
        from math import sin, cos, tan, asin, acos, atan, sinh, cosh, tanh
        from math import sqrt, pow, exp, log, log10, pi, e, radians, degrees
        
        namespace = {
            'sin': sin, 'cos': cos, 'tan': tan,
            'asin': asin, 'acos': acos, 'atan': atan,
            'sinh': sinh, 'cosh': cosh, 'tanh': tanh,
            'sqrt': sqrt, 'pow': pow, 'exp': exp,
            'log': log, 'log10': log10,
            'pi': pi, 'e': e,
            'radians': radians, 'degrees': degrees
        }
        
        exec(var + ' = 1', namespace)
        try:
            eval(expr, namespace)
            return True
        except Exception as e:
            return False
    
    @staticmethod
    def getMathNamespace():
        """Get namespace with math functions for eval()"""
        from math import sin, cos, tan, asin, acos, atan, sinh, cosh, tanh
        from math import sqrt, pow, exp, log, log10, pi, e, radians, degrees
        
        return {
            'sin': sin, 'cos': cos, 'tan': tan,
            'asin': asin, 'acos': acos, 'atan': atan,
            'sinh': sinh, 'cosh': cosh, 'tanh': tanh,
            'sqrt': sqrt, 'pow': pow, 'exp': exp,
            'log': log, 'log10': log10,
            'pi': pi, 'e': e,
            'radians': radians, 'degrees': degrees
        }
    
    def isInverted(expr):
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
        self.shapeSegCnt = params.mathFnResolution
        
        if(params.mathFnType == 'PARAMETRIC'):
            fn1 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric1)
            fn2 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric2)
            if(not MathFnDraw.testFn(fn1, 't')): return curvePts
            if(not MathFnDraw.testFn(fn2, 't')): return curvePts
            
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
                scaleFact = scaleFact / 25
            
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
            namespace = MathFnDraw.getMathNamespace()
            for step in range(intervals):
                try:
                    namespace['t'] = t
                    x = bbStart[idx0] + eval(fn1, namespace)
                    y = bbStart[idx1] + eval(fn2, namespace)
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
            
            if(inverted is None):
                return curvePts
            
            span = (bbEnd[idx1] - bbStart[idx1]) if(inverted) \
                else (bbEnd[idx0] - bbStart[idx0])
            
            intervals = int(self.shapeSegCnt * abs(span))
            if(intervals == 0): intervals = self.shapeSegCnt
            incr = span / intervals
            indep = bbStart[idx0] if(not inverted) else bbStart[idx1]
            
            namespace = MathFnDraw.getMathNamespace()
            for step in range(intervals):
                try:
                    if(inverted):
                        y = indep
                        namespace['y'] = indep
                    else:
                        x = indep
                        namespace['x'] = indep
                    dep = eval(expr, namespace)
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
        tools = bpy.context.workspace.tools
        tool = tools.from_space_view3d_mode('OBJECT', create = False)
        if not tool:
            tools.from_space_view3d_mode(GP_CONTEXT_MODE, create = False)
        
        # Local import to avoid circular dependency
        from ..tools.workspace_tools import FlexiDrawBezierTool
        
        if(tool.idname == FlexiDrawBezierTool.bl_idname and params.drawObjType == 'MATH'):
            
            if(params.mathFnType == 'PARAMETRIC'):
                fn1 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric1)
                fn2 = MathFnDraw.getEvaluatedExpr(params.drawMathFnParametric2)
                fnTxts = ['y = ' + (fn2 if MathFnDraw.testFn(fn2, 't') else INVALID), \
                    'x = '+ (fn1 if MathFnDraw.testFn(fn1, 't') else INVALID)]
            else:
                expr = MathFnDraw.getEvaluatedExpr(params.drawMathFn)
                inverted = MathFnDraw.isInverted(expr)
                if(inverted is None):
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
        if(MathFnDraw.mathFnItems is None or MathFnDraw.mathFnDirty):
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
            MathFnDraw.refreshDefaultParams()
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
        params.mathFnResolution = int(fnCurveRes)
        
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

# Math function operators
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
        if(fnName.strip() == ''):
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
            
            if(paramFn1.strip() == '' or paramFn2.strip() == ''):
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
            if(xyFn.strip() == ''):
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
        except:
            self.report({'ERROR'}, 'Error saving math function file')
        
        return {'FINISHED'}

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
