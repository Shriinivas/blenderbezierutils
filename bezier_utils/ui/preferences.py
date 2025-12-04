# bezier_utils/ui/preferences.py

import bpy
from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty, IntProperty, FloatProperty, StringProperty, FloatVectorProperty, EnumProperty
)
from .panel import updatePanel
from ..core.props import FTProps
from ..core.hotkeys import FTHotKeys
# BezierUtilsPreferences will be added here
class BezierUtilsPreferences(AddonPreferences):
    bl_idname = 'bezier_utils'


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

    showGuides: BoolProperty(
            name="Show Visual Guides", \
            description='Display visual guides: orientation axes, custom axis frame, pivot marker, constraint plane (Ctrl+Alt+H to toggle)', \
            default = True,
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
            col.label(text='Show Visual Guides (Ctrl+Alt+H):')
            col.prop(self, "showGuides", text = '')
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

