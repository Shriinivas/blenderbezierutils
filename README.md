# Blender Add-on with Bézier Utility Operations
<b>Add-on Version: 0.9.96 </b>

This add-on contains several tools and utility ops for working with Bézier curves. <br>
Supported Blender Versions: <b>2.8+, 3.0+</b> <br>
Video Tutorials: </b> https://www.youtube.com/playlist?list=PLxsh4i5F_h9G6QFoPzKvBRMayz8533fSW<br>

# Installation
### *** Please do not download the github zip, you need only the blenderbezierutils.py file
- Download blenderbezierutils.py (save the file from <a href = 'https://raw.githubusercontent.com/Shriinivas/blenderbezierutils/master/blenderbezierutils.py'>this location</a> in a folder on disk.). For Blender versions below 2.93, download <a href='https://raw.githubusercontent.com/Shriinivas/blenderbezierutils/master/blenderbezierutils_2.8.py'> blenderbezierutils_2.8.py </a>.
- Open Blender and select Edit->Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bézier Utilities' option in the add-ons dialog

After installation, a new tab: 'Bézier Utilities' is displayed in Object and Edit modes on 'Active Tool and Workspace settings' area on the properties panel. <br/>
There will also appear two new buttons - <b>Flexi Draw Bézier</b> and <b>Flexi Edit Bézier</b> - in Object Mode and a <b>Flexi Grease Bézier</b> button in GP Draw Mode on the toolshelf. 

# Overview
The tools Flexi Draw , Flexi Edit and Flexi Grease are interactive tools that allow drawing and editing Bézier curves.

## Flexi Draw Bézier Tool
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). It allows drawing Bézier curves by manipulating the control points.<br>
### Drawing Bézier Curve:
To draw the curve, activate the tool by clicking the Flexi Draw Bézier tool on the toolshelf and select Bézier Curve from Shape Type drop-down. Click the LMB on the starting point of the curve. Then click and drag LMB on the end point to adjust the curvature. You can continue drawing subsequent segments in this fashion. Double clicking or hitting enter or space will convert the drawing to a curve object. You can auto-close the curve by pressing Shift+Space or Shift+Enter<br>
<p align="center"><img src="https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawdemo.gif" alt="Demo"/></p><br/>
<b>Repositioning the Bézier point:</b> <br>
At the time of dragging the LMB to set the handle location, you can grab the Bézier point to reposition it by pressing G (configirable). All the snapping options are available for setting the handle location and repositioning the Bézier point. Press G again to release the grab. <br>
<b>Resetting Handle:</b> <br>
You can reset the handle if you need to draw a straight line after a curved segment by pressing the hot key Shift+G (configurable).<br/>
<b>Dissociating Handle:</b> <br>
While adjusting the handle by dragging the mouse pointer, you can change the handle type to free by pressing the hot-key V. This allows creating cusp nodes, even while creating the curve.<br/>
<b>Undo:</b><br/> 
While drawing, you can undo one segment at a time by pressing backspace. Pressing escape removes the entire curve. After the drawing is finished, the curve creation can be undone by pressing ctrl-Z  <br><br>

### Drawing Primitive Shapes:
<p align="center"><img src="https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawshapeoptions.png" alt="Draw Options"/></p><br/>
Other than Bézier curve you can also create primitive shapes with Bézier segments. To do this, select the appropriate primitive shape from the Shape Type drop-down in the tollbar at the top. The shapes currently available are 1) Rectangle 2) Ellipse / Circle 3) Polygon 4) Star. Click the starting point and end point of the shape. To draw shapes with equal height and width (e.g. perfect circle or square), hold down shift key. You can optionally snap to the grid.
The shapes are all 2d, so for drawing in the 3d perspective view, you can either choose a constraining plane or set the orientation to view. In perspective view, without any constraining plain, the depth value of the current mouse location is applied to all the points.<br>
<p align="center"><img src="https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/drawprim.gif" alt="Demo"/></p><br/>
<b>Adjusting number of segments (Ellipse) or Sides (Polygon / Star):</b> <br>
To adjust the segment count or sides use the mouse wheel or + or - keys on numpad. You may also directly type the value in the edit box on toolbar.<br>
<b>Drawing from Center or from Corner:</b> <br>
You can start the drawing either from the corner (Bounding Box) or from the center by choosing appropriate option in the Drawing Mode drop-down.<br>
<b>Sweep (Available in Ellipse / Circle, Polygon and Star):</b> <br>
The sweep angle can be entered directly in the edit box of the toolbar or changed incrementally using left or right arrow keys while drawing the shape. The maximum angle value is 360, which makes the shape complete.<br>
<b>Starting Angle (Available in Ellipse / Circle):</b> <br>
The starting angle can be entered directly in the edit box of the toolbar or changed incrementally using up or down arrow keys while drawing the shape. This angle changes the tilt of the ellipse being drawn.<br>
<b>Copy Object Properties:</b> <br>
Select an existing curve object in the Copy Object Properties option in the toolbar to apply its properties like material, dimension, bevel depth, bevel object on to the curves / shapes drawn with the Flexi Bézier Tool. If a mesh object is selected, only its material is applied to the curve / shape drawn.<br>

## Flexi Edit Bézier Tool
![Demo](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/editdemo.gif)<br>
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). With it You can 1) edit a Bézier curve by dragging a point on the curve 2) Move Segment endpoints and manipulate handles 3) Add or delete a vertex at any arbitrary location on the curve <br><br>
<b>Edit Curve and move handles and end points: </b> When the tool is activated, moving a mouse cursor in the 3d viewport would highlight the individual curve segments under the mouse cursor. Clicking on a segment will make the segment active and it's handles will be visible. You can drag any point on the curve to edit it. Releasing the button will apply the changes to the curve. Also when the handles are visible, bringing the mouse cursor in the vicinity of any of the handle points will highlight that point (bright green), indicating the mouse click will operate on it. You can move the handles and segment points by dragging the mouse pointer. Releasing the mouse button makes the changes permanent.<br><br>
<b>Grab the Edit Point:</b>
The point being edited can be grabbed by double clicking it. Once grabbed the point will move along the mouse, without having to drag the pointer. The grab is released on the next single click.
<b>Adding a vertex: </b> Hold down control and click the mouse on any location on curve to add a vertex at that position. If you also hold down shift along with control the added point will have aligned handles. The handles will be of type vector if alt and control are held down while pressing the mouse button.<br><br>
<b>Deleting a vertex or handle point: </b> 
Select any end point (the selected point is marked in dark green) and press del to delete it. Pressing del when a handle point is selected will align it with the other point of the segment. <br><br>
You can toggle between Flexi Draw and Flexi Edit by pressing <b>e</b>.<br>
Press <b>h</b> to toggle the visibility of the selected segment handles.<br><br>
<b>Subdivide Segments Uniformly: </b> 
You can also subdivide the selected segments uniformly. To initiate the subdivision op, first select the segments (hold down shift to select multiple segments). Then press w. Now there will appear a subdivision marker at the middle of each selected segment. You can increase or decrease the number of subdivisions by scrolling the mouse wheel or pressing + or - keys. Press Spacebar or Enter to confirm the subdiv operation.<br><br>
<b>Align Handle:</b>
To align the handle with the opposite handle of the same end point, select the handle point and press K. This way you can quickly smooth out the sharp corners.

## Flexi Grease Bézier Tool
![Demo](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/greasedemo.gif)<br>
This tool will appear on the toolshelf in Grease Pencil Draw mode. You can draw Bézier curves just as you would draw with the Flexi Draw tool. After confirming the drawing is converted to grease pencil strokes. All the snapping and locking options of the Flexi Draw are available here also.<br>
Additionaly, you can increase of decrease the resolution of the stroke using the mouse wheel or pressing + or - keys.<br>
The subdivision point visibility can be toggled by pressing h key.<br>

## Configurable Entities
![Configurable Options](https://github.com/Shriinivas/etc/blob/master/blenderbezierutils/illustrations/configitems.png)<br>
Values of a number of entities are user configurable via Add-ons dialog (from Preferences->Add-on Menu). Some of these are:
- Bézier Toolkit Panel Tab<br>
<b>Dimensions</b>
- Draw Line Thickness
- Handle Point Size
- Uniform Subdiv Point Size
- Flexi Grease Resolution Point Size
- Draw Marker Size
- Axis Line Thickness
- Snap Point size<br>
<b>Colors</b>
- Selected and Adjacent Segment
- Handle Tips & Bézier Points
- Highlighted Points
- Subdivision & Resolution Markers

## Snapping & Locking Framework
The Framework provides comprihensive snapping and locking options common to all three Flexi tools.
After activating one of the Flexi tools the header portion displayes the following new options:

<b>1) Constraining Axes Dropdown:</b>
Selecting an axis or axis-pair from this dropdown will constrain the point being drawn or edited to the an axis or plane that is parallel to the selection. The actual axis or plane will be determined by click location in 3d space. The interpretation of Axes X, Y, Z will differ based on the selection in the Snapping Orientation (see below).

<b>2) Snap to Plane Checkbox:</b>
This option is available when the Costrain Axes drop down has an axis-pair (plane) selection. When checked the point being drawn or edited will be snapped to the plane of the point selected in the Snap Origin dialog (see below). 
  
<b>3) Snapping Orientation Dropdown:</b>
The options in the dropdown are:
- Global Axes: Orientation along the global axes
- Reference Line: Orientation along the Reference Line (see explanation of Reference Line below)
- Custom Axes: Orientation along Custom Axis (see explanation of Custom Axis below)
- Active Object: Orientation along local space of active object
- Active Object Face: Orientation along the normal of the active object face under mouse pointer
- View: Orientation along current viewport view axes

The Orientation affects the constraining plane and axis, as well as the reference axis for snapping to angle increment. For example, if the constraining axis are XY and the selected option in the Orienation dropdown is Active object, the point will be constrained to the XY plane of the active object local space, if it is available. 

<b>4) Snapping Origin Dropdown:</b>
The options in the dropdown are:
- Global Origin
- 3d Cursor Location
- Custom Axis Start
- Reference Line Point
- Active Object Location
- Active Object Face

The distance values are calculated / interpretated based on the selection in this dialog. Additionally, the snapping plane (when Snap to Plane option is selected) is the plane containing the Snapping Origin point.

<b>5) Axis Scale Checkbox: </b>
When this option is selected, the entries made via keyboard (see below) are interpreted in terms of the scale of the Custom Axis or the Reference Line. 10 Units on the scale represent the total length of the axis. 

<b> Keyboard Input: </b>
It's now possible to directly enter the position values of the point being drawn or edited via keyboard. To set the next point location, start typing a number after starting a new segment (draw) or grabbing the edit point (edit). This number will be the movement along the first free of the free axes. User can enter values for the next axis by pressing tab. The values entered are with respect to the current Snapping Origin and along the Snapping Orientation Axes.
  
<b>Tweaking the Location via Keyboard: </b>
After moving the point by mouse, user can further tweak the value via keyboard (e.g. round off the coordinates). To initiate tweaking press P, whereupon the point no more movable by mouse. Now you can type the coordinates just as in case of entering through keyboard.

<b>Entering Location value in The Form of Polar Coordinates: </b>
When the drawing / editing is constrained to a single plane (the Constrain Axes dropdown has an axis-pair selection), you can tweak the polar coordinate values. To do this press P twice, after the first time you can enter cartesian coordinates and after the second, it's possible to enter the coordinates in polar form. Please note all the values are interpreted based on the selections in the Snapping Orientation and Snapping Origin dropdown.

<b>Reference Line and Reference Line Point: </b>
Reference Line takes on different connotations based on the context. While drawing (Flexi Draw and Flexi Grease), the Reference Line is the segment previous to the one being drawn and the Reference Line Point is the immediate previous point. While editing, for a Bézier Point, the Reference Line is the line joining the current Bézier point with the other one of the segment (which is also the Reference Line Point). In case of handles, the Reference Line is the handle itself and the Reference Line point is the Bézier point in the handle. While moving a point on the segment, the Reference Line point is the location of the point before it is moved.

<b>Custom Axis:</b>
Custom Axis is a user defined line, that serves multiple purpose. To create a Custom Axis, make sure the selection in Snapping Origin dropdown is 'Custom Axis Origin' and rightclick the starting point, move the pointer to the end point and rightclick once again. You can snap to grid or bezier point location while creating the Custom Axis. Additionally, it's possible to define custom snapping points along the Custom Axis, using mouse wheel (or plus or minus keys).<br>
The Custom Axis can be used to define the Snappig Orientation, Snapping Origin, Custom scale and Custom snapping points.

<b>Hotkey Snapping Options (Active for the Point being Drawn / Edited):</b> 
- Holding down <b>ctrl</b> while moving mouse will snap the point or handle to the<b> grid.</b> 
- By holding down <b>shift</b> key the <b>angle</b> of the segment / handle being drawn / edited will be restricted to fixed values (0, 45, 90 etc). The reference axis for determining angle increment is the first free axis based on the selected Snapping Orientation.<br>
- Holding down <b>alt</b> key will snap the point being drawn / edited to the a) <b>Bézier points</b> of the splines within all the curve objects in the view b) <b>vertices</b> of the selected objects (if there are fewer than 1000 vertices) c) <b>the face</b> under the mouse pointer of the selected objects (if there are fewer than 1000 faces). d) <b>snapping Origin</b> <br>
By default, after the drawing is started or the segment is selected for editing, and if the Snapping Orientation or Snapping Origin is 'Selected Object Face', the orientation / orign will be locked to the normal / center of the face under the mouse pointer. During the drawing / editing operation, user can make the tool reposition the orientation / origin to the new face by pressing <b>U</b>.<br>

By default, snapping to the end points joins the new curve to the curve(s) it is being snapped to. You can hold down <b>ctrl</b> while ending the curve (by double click or space or return key) to keep the curve separate. The curve can also be separated from the snapped curves be pressing ctrl-Z after confirming. <br>

The snapping gets adapted to the viewport zoom level. <br>

<b>Locking Options:</b>
- Pressing <b>X</b>, <b>Y</b>, <b>Z</b> while the curve is being drawn will lock the segment to the corresponding axis. 
- Pressing shift together with one of these buttons will lock the segment to the axes other than the one denoted by the button (e.g. <b>shift+Z</b> - lock to <b>XY plane</b>).
- Press escape to get out of the lock mode.

Snapping and locking can be combined together. So user can hold down both control and shift to snap to grid as well as restrict the angle. Likewise, user can press shift-Z to lock to XY plane and then hold shift while moving the mouse to restrict the angle between the segment end points.

When Constrain Axes dropdown has an axis-pair selection (lock to plane), pressing a single axis key will allow users to draw lines parallel to the axis (or slide the point along this line).

# Other Tools
The utility ops are arranged in a collapsible Panel, grouped according to the functionality type. 
Here is a brief overview of a few of the ops in this add-on:
- Separate Bézier Splines: Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bézier Segments: Create individual objects out of every segment within the Curve. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bézier Points: Create individual objects out of the Bézier end point of each segment within the selected Curves. The newly created point objects do not inherit the shape keys of the original curve objects. New objects are put in a separate collection. The points can be used for snapping with the Flexi Bézier Tool for interesting experimentation :)
- Select Objects in collection: Allows selection of objects belonging to the collection of the active object. It's possible to select alternate objects or objects at fixed interval (based on the order in collection), as well as invert the selection. This can be combined with the split / separate ops to work on the newly created segment / spline objects.
- Close with Straight Segment: Closes all (non-cyclic) splines within the selected curves with straight segments.
- Remove Duplicate Curve Vertices: Removes vertices at the same location. If there are duplicate end vertices coinciding with the start node the spline is marked cyclic.
- Convert to Mesh: Converts the curve to a mesh with quad faces. The curve is first made 2d and all its splines cyclic. This op basically applies remesh modifier to make a quad mesh. Users can optionally check unsubdivide option to reduce the polygon count further.
- Join Bézier Curves: Joins the selected objects at their end points. If the 'Join optimized' option is unchecked the curves are joined in their order in the collection (alphabetical order of their names). If it is checked, the next curve to join is chosen based on its distance (shortest) from the current curve; the curve direction is reversed if needed. With the 'Join With Straight Segments' option, the curve objects are joined with straight line segments, regardless of the end point handle types.
- Set Curve Colors: This tool allows users to set the display color of the selected curves in viewport. The colored curves are drawn on top of the Blender curve objects. There is a button - Apply Curve Color - to toggle the curve coloring.
- Paste Length: Makes the length of all the selected curve objects the same as that of active curve object. The scale remains unchanged.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. (This tool is the same as the one included in the <a href ='https://github.com/Shriinivas/assignshapekey'>Assign Shape Keys</a> add-on.)<br> 


<b>Video Tutorials & Demos:</b><br>
<b>Overview of Flexi Draw Bézier Tool:</b> https://youtu.be/C9PXp0XHgYQ<br>
<b>Overview of Flexi Edit Bézier Tool:</b> https://youtu.be/enQHGmluQIw<br>
<b>Overview of Snapping & Locking Framework:</b> https://youtu.be/VQCXZbOq47s<br>
<b>Overview of Flexi Grease Bézier & Uniform Subdiv Op:</b> https://youtu.be/4XrjpWwLU4M<br>
<b>Demo of Flexi Ellipse: </b>https://youtu.be/t7eVWP8gxeE<br>
<b>All Bézier Toolkit videos: </b> https://www.youtube.com/playlist?list=PLxsh4i5F_h9G6QFoPzKvBRMayz8533fSW<br>

# Credits
The functionality of editing curve by grabbing a point on it, is adapted from Inkscape edit curve tool: https://gitlab.com/inkscape/inkscape/blob/master/src/ui/tool/curve-drag-point.cpp <br>
I am grateful to the authors of the module for making this great piece of code accessible to everyone.<br><br>

The add-on script includes python converted a2c js function (Copyright (C) 2013-2015 by Vitaly Puzrin) located at the github repository: https://github.com/fontello/svgpath. It is used in Draw Ellipse / Circle tool. Heartfelt thanks to the authors for this amazing function!<br><br>

Algorithms that were inspired by the answers on stackoverflow are mentioned with the corresponding links in the code. I would also like to thank the users who provided this very useful information, sometimes even with working sample code.<br><br>

# Known Issues
- In Flexi Draw Bézier, the part of the area under toolshelf and properties panel is excluded from drawing. Hide these elements to maximize the drawing area.<br>

# Limitations
- In Flexi Draw Bézier, snapping does not work for curves with modifiers. This is the intended functionality.<br>

In general, exercise caution when using this add-on in production, since all possible conditions have not been extensively tested.<br>
You may report bug as comment on the youtube videos or on the issues page here on Github. I will try and fix them as soon as I can.
