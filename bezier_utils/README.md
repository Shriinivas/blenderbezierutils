# Bezier Utilities - Blender Addon

**Version:** 1.0.0-beta

This is the modular version of the Bezier Utilities addon for Blender.

## Installation

### For Users

Download the latest `bezier_utils.zip` from [Releases](https://github.com/Shriinivas/blenderbezierutils/releases) and install via Blender Preferences → Add-ons → Install from File.

### For Developers

1. Copy or symlink the entire `bezier_utils` directory to your Blender addons folder:
   - **Linux**: `~/.config/blender/{version}/scripts/addons/`
   - **macOS**: `~/Library/Application Support/Blender/{version}/scripts/addons/`
   - **Windows**: `%APPDATA%\Blender Foundation\Blender\{version}\scripts\addons\`

2. Restart Blender or reload scripts (F3 → "Reload Scripts")

3. Enable in: Edit → Preferences → Add-ons → Search "Bezier Utilities"

## Documentation

See the main [README](../README.md) for full documentation, features, and usage instructions.

## Development

This modular structure makes it easier to navigate, maintain, and extend the codebase.

### Linting
```bash
ruff check .
```

See [CLAUDE.md](CLAUDE.md) for detailed architecture and development guidelines.

## Reporting Issues

Report issues on the [GitHub Issues page](https://github.com/Shriinivas/blenderbezierutils/issues).
