# bezier_utils/__init__.py

bl_info = {
    "name": "Bezier Utilities",
    "author": "Shrinivas Kulkarni",
    "version": (1, 0, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Sidebar > Bezier Utilities",
    "description": "Utilities for Bezier curves",
    "warning": "",
    "doc_url": "",
    "category": "Curve",
}

# ruff: noqa: E402
# Import all modules
from . import constants  # noqa: F401
from .utils import math_utils, bezier_math, curve_utils, object_utils, view_utils, event_utils  # noqa: F401
from .core import props, hotkeys, menus, snap  # noqa: F401
from .drawing import primitives, math_fn  # noqa: F401
from .operators import simple_ops, modal_ops  # noqa: F401
from .tools import workspace_tools  # noqa: F401
from .ui import params, panel, preferences  # noqa: F401
from . import registration

def register():
    registration.register()
    print("Bezier Utilities: Registered")

def unregister():
    registration.unregister()
    print("Bezier Utilities: Unregistered")

if __name__ == "__main__":
    register()
