# bezier_utils/core/menus.py

import bpy
from bpy.types import Operator
from bpy.props import IntProperty
from .hotkeys import FTHotKeys
from ..constants import TOOL_TYPE_FLEXI_EDIT

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
        if(keyIds is not None and hotkeyId in keyIds):
            return FTMenu.idDataMap.get(hotkeyId)
        return None

    def procMenu(parent, context, event, outside):

        metakeys = parent.snapper.getMetakeys()
        evtType = event.type

        if(FTMenu.abandoned and evtType == 'ESC'):
            FTMenu.abandoned = False
            return True

        if(FTMenu.currMenuId is not None):
            params = bpy.context.window_manager.bezierToolkitParams
            opt = FTMenu.getCurrMenuSel()
            if(opt is not None or not evtType.startswith('TIMER')): # What's TIMER_REPORT?
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

        if(hkData is None): return False
        menuData = FTMenu.getMenuData(parent, hkData.id)
        if(menuData is None): return False
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
        if(menuData is not None):
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
    pass
