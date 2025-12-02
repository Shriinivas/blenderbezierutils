# bezier_utils/ui/preferences.py

from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty, IntProperty, FloatProperty, StringProperty
)
from .panel import updatePanel
from ..core.props import FTProps
# BezierUtilsPreferences will be added here
class BezierUtilsPreferences(AddonPreferences):
    bl_idname = __name__

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

    dispAxes: BoolProperty(
            name="Orientation / Origin Axis", \
            description='Display axes for selected orientation / origin', \
            default = False,
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
        name="Keymap Font Size",
        description="Font size for keymap display",
        default=12, min=8, max=32
    )
