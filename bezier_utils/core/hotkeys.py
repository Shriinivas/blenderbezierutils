# bezier_utils/core/hotkeys.py

from ..constants import INVAL, TOOL_TYPES_FLEXI_ALL, TOOL_TYPES_FLEXI_DRAW_COMMON, \
    TOOL_TYPE_FLEXI_EDIT, TOOL_TYPE_FLEXI_GREASE

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
        self.exclTools = exclTools if(exclTools is not None) else set()
        self.inclTools = inclTools if(inclTools is not None) else TOOL_TYPES_FLEXI_ALL
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
    hkToggleGuides = 'hkToggleGuides'
    hkSwitchOut = 'hkSwitchOut'
    hkTweakPos = 'hkTweakPos'
    hkToggleDrwEd = 'hkToggleDrwEd'
    hkReorient = 'hkReorient'

    commonHotkeys = []
    commonHotkeys.append(FTHotKeyData(hkToggleKeyMap, 'Ctrl+Shift+H', \
        'Hide / Unhide Keymap', \
            'Hide / Unhide Keymap Displayed When Flexi Tool Is Active'))
    commonHotkeys.append(FTHotKeyData(hkToggleGuides, 'Alt+Shift+G', \
        'Hide / Unhide Visual Guides', \
            'Toggle visual guides: orientation axes, custom axis frame, pivot marker, constraint plane'))
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
        prefs = context.preferences.addons[__package__.split('.')[0]].preferences
        if prefs is None:
            FTHotKeys.updating = False
            return
        hmap = FTHotKeys.idDataMap
        # Checking entire map even if one key changed (TODO)
        for kId in hmap:
            if(reset): hmap[kId].key = hmap[kId].default
            # User pressed key
            prefKey = getattr(prefs, kId)
            combKey = ''
            for meta in FTHotKeys.metas:
                if(hasattr(prefs, kId + meta) and \
                    getattr(prefs, kId + meta)):
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

        # Local import
        from ..operators.modal_ops import ModalBaseFlexiOp
        ModalBaseFlexiOp.propsChanged()

    def getAvailKey(prefs):
        availKey = None
        oldKeys = set([m.upper() for m in FTHotKeys.metas])
        newKeys = set([getattr(prefs, kd.id) \
            for kd in FTHotKeys.snapHotkeysMeta])
        for availKey in (oldKeys - newKeys): break # Only one entry
        return availKey

    def initSnapMetaFromPref(context):
        prefs = context.preferences.addons[__package__.split('.')[0]].preferences
        if prefs is None:
            return
        for i, metakeyData in enumerate(FTHotKeys.snapHotkeysMeta):
            metaKeyId = metakeyData.id
            prefKeyMeta = getattr(prefs, metaKeyId)
            metakeyData.key = prefKeyMeta
            if(prefKeyMeta == 'KEY'):
                regKeyId = metaKeyId[0:metaKeyId.index('Meta')]
                regKeyData = FTHotKeys.idDataMap[regKeyId]
                prefKeyReg = getattr(prefs, regKeyId)
                regKeyData.key = prefKeyReg
        
        # Local import
        from ..operators.modal_ops import ModalBaseFlexiOp
        ModalBaseFlexiOp.propsChanged()

    # Validation for snap keys (format: EITHER meta key OR regular key)
    # (regular part validated by updateHotkeys)
    # UI Format: Drop-down with entries Ctrl, Alt, Shift, Keyboard and ...
    # a separate single key field activated with selection of 'Keyboard' in the dropdown
    def updateSnapMetaKeys(dummy, context):
        if(FTHotKeys.updating): return

        FTHotKeys.updating = True
        prefs = context.preferences.addons[__package__.split('.')[0]].preferences
        if prefs is None:
            FTHotKeys.updating = False
            return
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
        
        # Local import
        from ..operators.modal_ops import ModalBaseFlexiOp
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
            if(kcMap.get(i) is not None and kcMap[i] not in exclKeys):
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
