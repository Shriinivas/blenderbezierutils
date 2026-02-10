# blenderbezierutils (Noisysundae Fork)

**Disclaimer:** This is a quick and dirty edit of blenderbezierutils. I edited parts related to the curve-to-SVG exporter and added options that optimize the output SVG content for web graphics, especially icons. I have made a pull request with the hope that these changes get merged into the base repository, so I might not be actively maintaining this fork.

## Changes

* **NEW:** Copy SVG to Clipboard (curves only)
	* Same options as "Export to SVG", except file path
	* For seamless transfer with other softwares that support SVG data copy/paste
* SVG export file view (curves only)
	* Added the following options:
		* **View Attribute:** Select attribute(s) to specify bounding box, "[width](https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/width) / [height](https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/height)" (legacy) or "[viewBox](https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/viewBox)"
		* **ID Attribute:** Toggle "id" attribute insertion (legacy: always)
			* **Use Object Names as Ids:** If "ID Attribute" is enabled, switch between using sequence numbers and object names for id attributes
				* Names are used AS IS, so make sure to use naming patterns suitable for the attribute values (e.g. not containing double quotes).
		* **Style Attribute:** Toggle "style" attribute insertion (legacy: always)
			* Enabling this reveals the legacy color pick options.
		* **Max Path Value Precision:** Maximum precision of the coordinate values (legacy: infinite)
			* Defaults to 2. Adjust it as needed for a trade-off between precision and output SVG file size.
		* **Path Value Repeating Number Threshold:** Rounds a coordinate value to a precision where the first repeating group of 0s or 9s is found
			* Defaults to 3.
			* Without this, result of an integer coordinate value will sometime contain excess fractional values (e.g. 748.8 → 748.8000000000001, 323.41 → 323.40999999999997) due to floating point errors.
			* Also optimizes output file size.
		* **Use Relative Positions:** Makes use of relative path commands
			* Enabled by default.
			* In most cases, enabling this optimizes the output file size since relative position values occupy less digits than absolute ones, especially on large reference bounding boxes, or shapes with high amount of control points.
	* If there is one, the **Export View** option now defaults to the main camera in the active scene.
	* Changed default **Clip View** option value and style colors
	* Now excludes curve objects disabled in renders
* Removed XML tag in the output SVG file.
* Added curve shape key mix support.
	* This means shape keys can be used to manipulate output, useful for SVG animations.
* Output SVG paths for splines with **Cyclic U** enabled (closed shapes) are now ended with "z" path commands.
	* To eliminate duplicate first control points, as well as properly connect closed shapes.

Below is the content from the original README.

---

# Blender Add-on with Bézier Utility Operations

**Add-on Version: 1.1**

This add-on contains several tools and utility ops for working with Bézier curves.
**Supported Blender Versions:** 4.2 LTS, 4.3+, 5.0

> [!NOTE]
> This version is optimized for **Blender 5.0**. If you experience issues with earlier Blender versions (4.2, 4.3), please download [v1.0.0-beta](https://github.com/Shriinivas/blenderbezierutils/releases/tag/v1.0.0-beta) from the releases page.

**Video Tutorials:** https://www.youtube.com/playlist?list=PLxsh4i5F_h9G6QFoPzKvBRMayz8533fSW

---

# Installation

## Download and Install

1. Download the latest `bezier_utils.zip` from [Releases](https://github.com/Shriinivas/blenderbezierutils/releases/latest)
2. Open Blender and select Edit → Preferences
3. Click "Add-ons" tab and then "Install Add-on from File"
4. Select the downloaded `bezier_utils.zip` file
5. Check the 'Bézier Utilities' option in the add-ons dialog

After installation, a new tab 'Bézier Utilities' appears in Object and Edit modes on the properties panel. Three new tools appear on the toolshelf:
- **Flexi Draw Bézier** (Object Mode)
- **Flexi Edit Bézier** (Object Mode)
- **Flexi Grease Bézier** (GP Draw Mode)

---

# Overview

The Flexi Draw, Flexi Edit, and Flexi Grease tools are interactive tools for drawing and editing Bézier curves with comprehensive snapping and transform options.

## Flexi Draw Bézier Tool

Available in object mode, this tool allows drawing Bézier curves by manipulating control points.

### Drawing Bézier Curve

Activate the tool by clicking Flexi Draw Bézier on the toolshelf and select Bézier Curve from Shape Type dropdown. Click LMB on the starting point. Click and drag LMB on the end point to adjust curvature. Continue drawing subsequent segments. Double-click or hit Enter/Space to convert to curve object. Auto-close with Shift+Space or Shift+Enter.

![Draw Demo](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawdemo.gif)

**Repositioning the Bézier Point:**
While dragging to set handle location, press G (configurable) to grab the Bézier point and reposition it. All snapping options available. Press G again to release grab.

**Resetting Handle:**
Reset the handle for a straight line after curved segment by pressing Shift+R (configurable).

**Dissociating Handle:**
While adjusting handle, press V to change handle type to free, allowing cusp nodes.

**Undo:**
Press Backspace to undo one segment. Escape removes entire curve. After finishing, undo with Ctrl+Z.

### Drawing Primitive Shapes

![Shape Options](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawshapeoptions.png)

Select shape from Shape Type dropdown: Rectangle, Ellipse/Circle, Polygon, Star, Math Function. Click starting and end points. Hold Shift for equal height/width (perfect circle/square). Optionally snap to grid.

![Draw Primitives](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawprim.gif)

**Adjusting Segments/Sides:**
Use mouse wheel or +/- keys on numpad. Directly type value in toolbar edit box.

**Drawing from Center or Corner:**
Choose option in Drawing Mode dropdown.

**Sweep (Ellipse/Circle, Polygon, Star):**
Enter angle in toolbar or change incrementally with left/right arrow keys. Max: 360°.

**Starting Angle (Ellipse/Circle):**
Enter angle in toolbar or change with up/down arrow keys. Changes ellipse tilt.

**Copy Object Properties:**
Select existing curve object in toolbar to apply its material, dimension, bevel properties to new curves.

### Drawing Math Functions

Select **Math Function** from Shape Type dropdown to plot mathematical equations. Enter equation (e.g., `sin(x)`, `x**2`), set resolution and range. Supports XY and parametric functions with standard math operations. Functions can be saved and loaded for reuse.

## Flexi Edit Bézier Tool

![Edit Demo](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/editdemo.gif)

Available in object mode for editing Bézier curves.

**Edit Curve and Move Handles:**
Moving mouse highlights segments. Click to activate segment and show handles. Drag any curve point to edit. Drag handles and endpoints. Release to apply changes.

**Grab Edit Point:**
Double-click to grab point. Point follows mouse without dragging. Release with next single click.

**Adding Vertex:**
Hold Ctrl and click on curve. Shift+Ctrl for aligned handles. Alt+Ctrl for vector handles.

**Deleting Vertex:**
Select endpoint (dark green) and press Del. Del on handle point aligns it.

Toggle between Draw/Edit with **E**.
Toggle handle visibility with **H**.

**Subdivide Segments Uniformly:**
Select segments (Shift for multiple). Press W. Use mouse wheel or +/- to change subdivisions. Spacebar/Enter to confirm.

**Align Handle:**
Select handle point and press K to align with opposite handle.

## Flexi Grease Bézier Tool

![Grease Demo](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/greasedemo.gif)

Available in Grease Pencil Draw mode. Draw Bézier curves that convert to grease pencil strokes. All snapping/locking options available. Adjust stroke resolution with mouse wheel or +/-. Toggle subdivision visibility with H.

---

# Transform System

## Transform Presets

Quick-access buttons at the top of the toolbar for common configurations:

- **Free** - Global orientation, cursor pivot, standard workflow
- **Continue** - Reference orientation, offsets from previous point for smooth continuation
- **Axis** - Custom axis orientation and snapping for isometric/angled work
- **Surface** - Face-aligned orientation for drawing on mesh surfaces

## Pivot Point

Sets where transformation axes are centered (where RGB orientation lines originate).

**Options:**
- **3D Cursor** - Pivot at cursor location
- **Global Origin** - Pivot at world origin (0,0,0)
- **Active Object** - Pivot at active object location
- **Custom Axis Start** - Pivot at custom axis start point
- **Active Object Face** - Pivot at face center under cursor

## Offset Reference

Sets reference point for numeric input and angle snapping (separate from pivot point).

**Options:**
- **From Pivot** - Calculate offsets from selected pivot point
  *Example: X:5 means 5 units from pivot along X*
- **From Previous Point** - Calculate offsets from last drawn point
  *Example: X:5 means 5 units from previous point (useful for continuous drawing)*

---

# Snapping & Locking Framework

Comprehensive snapping and locking options common to all Flexi tools.

## Header Options

### Constraining Axes Dropdown
Select axis (X, Y, Z) or axis-pair (XY, XZ, YZ) to constrain drawing/editing to parallel axis or plane. Interpretation depends on Snapping Orientation.

### Snap to Plane Checkbox
Available when Constraining Axes has plane selection. Snaps point to plane of selected Pivot Point.

### Snapping Orientation Dropdown

**Options:**
- **Global Axes** - World space orientation
- **Reference Line** - Along previous segment/handle
- **Custom Axes** - Along custom axis with full XYZ coordinate frame
- **Active Object** - Object local space
- **Active Object Face** - Face normal under cursor
- **View** - Current viewport view axes

Affects constraining plane/axis and reference axis for angle snapping.

### Axis Scale Dropdown

- **Default** - Standard unit scale
- **Reference** - Scale based on reference line length (10 units = full length)
- **Axis** - Scale based on custom axis length (10 units = full length)

## Keyboard Input

Directly enter position values by typing after starting segment or grabbing edit point. First number sets movement along first free axis. Press Tab for next axis. Values relative to Offset Reference point along Snapping Orientation axes.

**Tweaking Location:**
After moving with mouse, press P to tweak via keyboard (e.g., round off coordinates). Point no longer moves with mouse.

**Polar Coordinates:**
When constrained to plane, press P twice to enter polar form (radius, angle).

## Reference Line

Context-dependent meaning:
- **Drawing:** Previous segment; point is previous endpoint
- **Editing Bézier Point:** Line to other endpoint of segment
- **Editing Handle:** Handle line itself; point is Bézier point
- **Moving Segment Point:** Point is original location before move

---

# Custom Axis

User-defined line serving multiple purposes.

## Creating Custom Axis

1. Set Pivot Point dropdown to 'Custom Axis Start'
2. Right-click starting point
3. Move pointer to end point
4. Right-click to confirm (or Enter/Space)
5. Type numeric coordinates during definition (e.g., X:10, Y:0)
6. Add snap points along axis with mouse wheel or +/-

**Snapping During Axis Definition:**
- Grid snapping (Ctrl)
- Vertex snapping (Alt)
- Angle snapping (Shift)
- Combined snapping (e.g., Ctrl+Shift)

## Custom Axis Uses

- **Snapping Orientation** - Align drawing to custom axis with full coordinate frame visualization (XYZ axes)
- **Pivot Point** - Use axis start as transformation center
- **Custom Scale** - Scale numeric input to axis length
- **Custom Snapping Points** - Snap to points defined along axis

---

# Visual Guides

Toggle visual feedback elements with **Ctrl+Alt+H** (configurable in preferences).

**Displayed Elements:**
1. **Orientation Axes** - Red (X), Green (Y), Blue (Z) lines from pivot point
2. **Custom Axis Frame** - Full 3D coordinate frame when custom axis active
3. **Pivot Point Marker** - Orange dot at transformation center
4. **Constraint Plane** - Outline showing active plane constraint
5. **Reference Line** - Highlighted line from offset reference to current point

---

# Hotkey Snapping Options

Active for point being drawn/edited:

- **Ctrl** - Snap to grid (1.0 unit intervals, aligned to integer world coordinates)
- **Shift** - Restrict angle to fixed increments (0°, 45°, 90°, etc.)
- **Alt** - Snap to Bézier points, vertices, faces, pivot point

After starting drawing, if Snapping Orientation or Pivot Point is 'Selected Object Face', orientation/origin locks to face under cursor. Press **U** to reposition to new face.

By default, snapping to endpoints joins new curve to existing curves. Hold Ctrl while ending curve to keep separate.

## Grid Overlay

When constraining to planes (XY, XZ, YZ), a visual grid overlay appears with 1.0 unit spacing. The grid aligns to integer world coordinates in Global mode. Ctrl-snapping snaps to grid intersections. Grid visibility is controlled by the guide toggle (**Ctrl+Alt+H**).

## Locking Options

- **X, Y, Z** - Lock segment to corresponding axis
- **Shift+X/Y/Z** - Lock to plane (e.g., Shift+Z = XY plane)
- **Esc** - Exit lock mode

Combine snapping and locking (e.g., Ctrl+Shift for grid + angle restriction).

---

# Bézier Utilities Panel

Utility operations arranged in collapsible sections on properties panel.

## Curve Operations

**Intersect / Boolean mode toggle**

### Intersect Curves
Finds intersection points between curves. Creates point objects or splits curves at intersections. Supports self-intersection.

**Options:**
- **Action** - Create points, insert vertices, or cut curves
- **Only Non-active** - Process only non-active curves
- **Self Intersection** - Find intersections within same curve

### Boolean Operations
Combines closed curves using boolean operations.

**Operations:**
- **Union** - Combine curves into single shape
- **Difference** - Subtract second curve from first
- **Intersection** - Keep only overlapping regions

**Requirements:** Works with closed (cyclic), coplanar curves.

## Split Curves

- **Separate Splines** - Create individual objects from each spline. Works with shape keys. Objects placed in separate collection.

- **Separate Segments** - Create individual objects from every segment. Works with shape keys.

- **Separate Points** - Create point objects from Bézier endpoints. Useful for snapping experiments.

## Join Curves

Joins selected objects at endpoints.
- **Join With Straight Segments** - Join with straight lines regardless of handle types
- **Join Optimized** - Next curve chosen by shortest distance, direction reversed if needed
- **Merge Distance** - Maximum distance for merging endpoints

## Align to Face

Aligns curve object to mesh face. Select mesh object, select curve, click button. Curve aligned to face under cursor.

**Options:**
- **Align Origin** - Position curve origin (Face Center, BBox Center, Original)
- **Align Location** - Move curve location to face center

## Select Objects in Collection

Select objects in active object's collection.

**Options:**
- **Interval** - Select every Nth object
- **Invert Selection** - Invert current selection within collection

Combine with split/separate ops for powerful workflows.

## Convert Curve to Mesh

> [!WARNING]
> The advanced meshing algorithms (Smart, Medial Axis, Grid TFI, Polar, Rectangle, Polygon, QMorph) are **experimental** and may not work correctly for all curve shapes. Use with caution.

Converts curve to quad mesh with intelligent meshing algorithms. Curve made 2D, all splines made cyclic.

**Fill Type Options:**
- **Smart** - Automatic shape detection with optimized mesher selection
- **Medial Axis** - Medial axis-based meshing with Grid Fill TFI
- **Grid TFI** - Direct Transfinite Interpolation grid
- **Polar** - O-Grid topology for circles/ellipses (pole-free)
- **Rectangle** - Pure quad grid for rectangular shapes (no poles)
- **Polygon** - Radial meshing for regular polygons (pole-free)
- **QMorph** - Q-Morph algorithm for convex freeform shapes
- **Quad** - Remesh modifier-based quad mesh
- **Quadriflow** - Quadriflow remeshing
- **Grid** - Grid-based remeshing
- **Offset** - Offset ring meshing
- **Triangulated** - Standard triangulated faces

**Additional Options:**
- **Remesh Depth** - Subdivision depth for quad mesh
- **Unsubdivide** - Reduce poly count after conversion
- **Resolution** - Mesh density for triangulated faces
- **Fill Detail** - Detail level for advanced meshers
- **Offset Size** - Offset ring size for offset meshing

## Set Handle Type

Set handle type (Vector, Aligned, Free, Auto) for all points in selected curves.

## Remove Duplicate Vertices

Removes coinciding vertices. If duplicate end vertices coincide with start, spline marked cyclic.

**Options:**
- **Proximity** - Maximum distance for considering vertices duplicate

## Set Curve Colors

Set viewport display colors for curves. Colors drawn on top of Blender curve objects. Toggle with **Apply Curve Color** button.

## Other Tools

- **Smart 2D Project** - Flatten 3D curves onto a best-fit plane while preserving shape. Uses SVD to compute optimal projection plane.
- **Export to SVG** - Export curves to SVG format. Preserves bezier data with adaptive recursive subdivision for accuracy.
- **Paste Length** - Match curve lengths to active curve. Scale unchanged.
- **Close Splines** - Mark non-cyclic splines as cyclic
- **Close with Straight Segment** - Close with straight line segment
- **Open Splines** - Open cyclic splines

## Edit Mode Tools

- **Mark Start Vertex** - Marks start vertex of closed splines. Important for curves with shape keys as changing start vertex may distort them.

---

# Configurable Entities

![Configurable Options](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/configitems.png)

Configure via Edit → Preferences → Add-ons → Bezier Utilities:

**Dimensions:**
- Draw Line Thickness, Handle Point Size, Marker Sizes, Axis Line Thickness, Snap Point Size

**Colors:**
- Selected/Adjacent Segments, Handle Tips, Bézier Points, Highlighted Points, Subdivision Markers

**Visual Guides:**
- Toggle visibility of orientation axes, custom axis frame, pivot marker, constraint plane

**Keymap:**
- Customize all hotkeys including toggle guides hotkey (Ctrl+Alt+H)

---

# Video Tutorials & Demos

- **Overview of Flexi Draw Bézier Tool:** https://youtu.be/C9PXp0XHgYQ
- **Overview of Flexi Edit Bézier Tool:** https://youtu.be/enQHGmluQIw
- **Overview of Snapping & Locking Framework:** https://youtu.be/VQCXZbOq47s
- **Overview of Flexi Grease Bézier & Uniform Subdiv:** https://youtu.be/4XrjpWwLU4M
- **Demo of Flexi Ellipse:** https://youtu.be/t7eVWP8gxeE
- **All Bézier Toolkit Videos:** https://www.youtube.com/playlist?list=PLxsh4i5F_h9G6QFoPzKvBRMayz8533fSW

---

# Credits

The curve editing functionality is adapted from Inkscape's edit curve tool:
https://gitlab.com/inkscape/inkscape/blob/master/src/ui/tool/curve-drag-point.cpp

The ellipse/circle tool uses python-converted a2c js function (Copyright © 2013-2015 by Vitaly Puzrin):
https://github.com/fontello/svgpath

Algorithms inspired by StackOverflow answers are referenced in the code.

---

# Known Issues

- In Flexi Draw Bézier, area under toolshelf and properties panel excluded from drawing. Hide these to maximize drawing area.

# Limitations

- Snapping doesn't work for curves with modifiers (intended behavior).

Exercise caution in production as all conditions not extensively tested. Report bugs on YouTube video comments or GitHub issues.
