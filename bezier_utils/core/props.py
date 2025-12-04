# bezier_utils/core/props.py


class FTProps:
    propUpdating = False

    def updateProps(dummy, context):
        FTProps.updatePropsPrefs(context)

    def updatePropsPrefs(context, resetPrefs = False):
        if(FTProps.propUpdating): return
        FTProps.propUpdating = True
        try:
            # Local import to avoid circular dependency
            from ..operators.modal_ops import ModalBaseFlexiOp

            prefs = context.preferences.addons[__package__.split('.')[0]].preferences
            if prefs is None:
                FTProps.initDefault()
                FTProps.propUpdating = False
                return

            props = ['drawPtSize', 'lineWidth', 'axisLineWidth', 'editSubdivPtSize', \
            'greaseSubdivPtSize', 'colDrawSelSeg', 'colDrawNonHltSeg', 'colDrawHltSeg', \
            'colDrawMarker', 'colGreaseSelSeg', 'colGreaseNonHltSeg', 'colGreaseMarker', \
            'colHdlFree', 'colHdlVector', 'colHdlAligned', 'colHdlAuto', 'colSelTip', \
            'colHltTip', 'colBezPt', 'colHdlPtTip', 'colAdjBezTip', 'colEditSubdiv', \
            'colGreaseSubdiv', 'colGreaseBezPt', 'colKeymapText', 'colKeymapKey', 'snapDist', \
            'dispSnapInd', 'showGuides', 'snapPtSize', 'liveUpdate', 'dispCurveRes', \
            'showKeyMap', 'keyMapFontSize', 'keyMapLocX', 'keyMapLocY', 'keyMapNextToTool', \
            'defBevelFact', 'maxBevelFact', 'minBevelFact', 'bevelIncr', 'numpadEntry', \
            'mathFnTxtFontSize', 'colMathFnTxt']

            if(resetPrefs):
                FTProps.initDefault()
                for prop in props:
                    exec('prefs.' + prop +' = FTProps.' + prop)
                    # ~ setattr(prefs, prop, getattr(FTProps, prop))
            else:
                for prop in props:
                    exec('FTProps.' + prop +' = prefs.' + prop)
                    # ~ setattr(FTProps, prop, getattr(prefs, prop))
            FTProps.hdlColMap = {
                "FREE": FTProps.colHdlFree,
                "VECTOR": FTProps.colHdlVector,
                "ALIGNED": FTProps.colHdlAligned,
                "AUTO": FTProps.colHdlAuto,
            }

        except Exception as e:
            print("BezierUtils: Error fetching default sizes in Draw Bezier", e)
            FTProps.initDefault()

        # Local import
        from ..operators.modal_ops import ModalBaseFlexiOp
        ModalBaseFlexiOp.propsChanged()

        try: ModalBaseFlexiOp.opObj.refreshDisplaySelCurves()
        except: pass

        FTProps.propUpdating = False

    def initDefault():
        FTProps.drawPtSize = 5
        FTProps.lineWidth = 1.5
        FTProps.axisLineWidth = .25
        FTProps.editSubdivPtSize = 6
        FTProps.greaseSubdivPtSize = 4

        FTProps.colDrawSelSeg = (.6, .8, 1, 1)
        FTProps.colDrawNonHltSeg = (.1, .4, .6, 1)
        FTProps.colDrawHltSeg = (.2, .6, .9, 1)

        FTProps.colGreaseSelSeg = (0.2, .8, 0.2, 1)
        FTProps.colGreaseNonHltSeg = (0.2, .6, 0.2, 1)

        FTProps.colHdlFree = (.6, .05, .05, 1)
        FTProps.colHdlVector = (.4, .5, .2, 1)
        FTProps.colHdlAligned = (1, .3, .3, 1)
        FTProps.colHdlAuto = (.8, .5, .2, 1)

        FTProps.colDrawMarker = (.6, .8, 1, 1)
        FTProps.colGreaseMarker = (0.2, .8, 0.2, 1)

        FTProps.colSelTip = (.2, .7, .3, 1)
        FTProps.colHltTip = (.8, 1, .8, 1)
        FTProps.colBezPt = (1, 1, 0, 1)
        FTProps.colHdlPtTip = (.7, .7, 0, 1)
        FTProps.colAdjBezTip = (.1, .1, .1, 1)

        FTProps.colEditSubdiv = (.3, 0, 0, 1)

        FTProps.colGreaseSubdiv = (1, .3, 1, 1)
        FTProps.colGreaseBezPt = (1, .3, 1, 1)
        FTProps.colKeymapText = (1.0, 1.0, 1.0, 1.0)
        FTProps.colKeymapKey = (0.0, 1.0, 1.0, 1.0)

        FTProps.snapDist = 20
        FTProps.dispSnapInd = False
        FTProps.showGuides = True
        FTProps.showKeyMap = False
        FTProps.keyMapFontSize = 10
        FTProps.keyMapLocX = 10
        FTProps.keyMapLocY = 10
        FTProps.keyMapNextToTool = True
        FTProps.snapPtSize = 3
        FTProps.liveUpdate = False
        FTProps.dispCurveRes = .4
        FTProps.defBevelFact = 4
        FTProps.maxBevelFact = 15
        FTProps.minBevelFact = -15
        FTProps.bevelIncr = .5
        FTProps.numpadEntry = False

        FTProps.mathFnTxtFontSize = 20
        FTProps.colMathFnTxt = (0.6, 1.0, 0.03, 1.0)
        FTProps.hdlColMap = {
            "FREE": FTProps.colHdlFree,
            "VECTOR": FTProps.colHdlVector,
            "ALIGNED": FTProps.colHdlAligned,
            "AUTO": FTProps.colHdlAuto,
        }
# Initialize during registration, not at import time to avoid circular dependencies
# FTProps.initDefault()
