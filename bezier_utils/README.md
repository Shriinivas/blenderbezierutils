# Bezier Utils - Modular Refactored Version

**Status: Work in Progress - Testing Phase**

This is a modular refactored version of the original `blenderbezierutils.py` addon. The monolithic file has been restructured into a proper Python package with organized modules for better maintainability and development.

## ⚠️ Important Notice

This refactored version is currently **undergoing testing**. While it loads successfully and all operators are registered, comprehensive functional testing is still in progress.

**For production use, please continue using the stable version:** [`blenderbezierutils.py`](../blenderbezierutils.py)

## Installation

1. Copy the entire `bezier_utils` directory to your Blender addons folder:
   - **Linux**: `~/.config/blender/{version}/scripts/addons/`
   - **macOS**: `~/Library/Application Support/Blender/{version}/scripts/addons/`
   - **Windows**: `%APPDATA%\Blender Foundation\Blender\{version}\scripts\addons\`

2. Open Blender and go to Edit → Preferences → Add-ons
3. Search for "Bezier Utilities"
4. Enable the addon by checking the checkbox

## Package Structure

```
bezier_utils/
├── __init__.py           # Main addon entry point
├── constants.py          # Global constants and enums
├── registration.py       # Addon registration logic
│
├── core/                 # Core systems
│   ├── props.py         # Property definitions
│   ├── hotkeys.py       # Hotkey management
│   ├── menus.py         # Menu definitions
│   └── snap.py          # Snapping system
│
├── utils/               # Utility functions
│   ├── math_utils.py    # Mathematical utilities
│   ├── bezier_math.py   # Bezier curve mathematics
│   ├── curve_utils.py   # Curve manipulation utilities
│   ├── object_utils.py  # Object utilities
│   ├── view_utils.py    # Viewport utilities
│   └── event_utils.py   # Event handling utilities
│
├── drawing/             # Drawing systems
│   ├── primitives.py    # Primitive shape drawing
│   └── math_fn.py       # Mathematical function drawing
│
├── operators/           # Blender operators
│   ├── simple_ops.py    # Simple utility operators
│   └── modal_ops.py     # Modal interactive operators
│
├── tools/               # Workspace tools
│   └── workspace_tools.py
│
└── ui/                  # User interface
    ├── panel.py         # UI panels
    ├── params.py        # Tool parameters
    └── preferences.py   # Addon preferences
```

## Features

All features from the original `blenderbezierutils.py` are preserved:

- **Flexi Draw Bezier Tool** - Interactive Bezier curve drawing
- **Flexi Edit Bezier Tool** - Interactive curve editing
- **Flexi Grease Bezier Tool** - Grease pencil Bezier curves
- **Comprehensive snapping and locking framework**
- **Utility operators** for curve manipulation
- **Primitive shape drawing** (Rectangle, Ellipse, Polygon, Star)

## Development

This modular structure makes it easier to:
- Navigate and understand the codebase
- Add new features
- Fix bugs
- Maintain code quality
- Collaborate with other developers

## Testing Status

- ✅ Addon loads without errors
- ✅ All operators register successfully
- ✅ Code passes linting (ruff)
- ⏳ Functional testing in progress

## Reporting Issues

If you encounter any issues with this refactored version, please report them on the [GitHub Issues page](https://github.com/Shriinivas/blenderbezierutils/issues) with the tag `[refactored-version]`.

## License

Same as the original addon. See main repository for details.
