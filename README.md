
# Blender Add-on with Bezier Utility Operations
<b>Add-on Version: 0.9.54 </b>

This add-on contains several tools for working with Bezier curves. <br>
Supported Blender Version: <b>2.8</b> <br>

# Installation
- Download blenderbezierutils.py (save the file from <a href = 'https://raw.githubusercontent.com/Shriinivas/blenderbezierutils/master/blenderbezierutils.py'>this location</a> in a folder on disk.)
- Open Blender and select File->User Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bezier Utilities' option in the add-ons dialog

After installation, a new 'Bezier Utilities' tab is displayed in object mode on 'Active Tool and Workspace settings' tab on the properties panel. There will also appear two new buttons - <b>Flexi Draw Bezier</b> and <b>Flexi Edit Bezier</b> - in Object Mode and a <b>Flexi Grease Bezier</b> button in Draw Mode on the toolshelf. 

# Overview
The tools are arranged in a collapsible Panel, grouped according to the functionality type. 
Here is a brief overview of a few of the tools in this add-on:
- Separate Bezier Splines: Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bezier Segments: Create individual objects out of every segment within the Curve. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bezier Points: Create individual objects out of the Bezier end point of each segment within the selected Curves. The newly created point objects do not inherit the shape keys of the original curve objects. New objects are put in a separate collection. The points can be used for snapping with the Flexi Bezier Tool for interesting experimentation :)
- Select Objects in collection: Allows selection of objects belonging to the collection of the active object. It's possible to select alternate objects or objects at fixed interval (based on the order in collection), as well as invert the selection. This can be combined with the split / separate ops to work on the newly created segment / spline objects.
- Close with Straight Segment: Closes all (non-cyclic) splines within the selected curves with straight segments.
- Remove Duplicate Curve Vertices: Removes vertices at the same location. If there are duplicate end vertices coinciding with the start node the spline is marked cyclic.
- Convert to Mesh: Converts the curve to a mesh with quad faces. The curve is first made 2d and all its splines cyclic. This op basically applies remesh modifier to make a quad mesh. Users can optionally check unsubdivide option to reduce the polygon count further.
- Join Bezier Curves: Joins the selected objects at their end points. If the 'Join optimized' option is unchecked the curves are joined in their order in the collection (alphabetical order of their names). If it is checked, the next curve to join is chosen based on its distance (shortest) from the current curve; the curve direction is reversed if needed. With the 'Join With Straight Segments' option, the curve objects are joined with straight line segments, regardless of the end point handle types.
- Paste Length: Makes the length of all the selected curve objects the same as that of active curve object. The scale remains unchanged.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. (This tool is the same as the one included in the <a href ='https://github.com/Shriinivas/assignshapekey'>Assign Shape Keys</a> add-on.)<br> 

# Flexi Draw Bezier Tool
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/drawdemo.gif)<br>
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). It allows drawing Bezier curves by manipulating the control points.<br>
<b>Drawing the curve: </b> To draw the curve activate the tool by clicking the Flexi Draw Bezier tool on the toolbar. Click the LMB on the starting point of the curve. Then click and drag LMB on the end point to adjust the curvature. You can continue drawing subsequent segments in this fashion. Double clicking or hitting enter or space will convert the drawing to a curve object. <br>
<b>Repositioning the Bezier point:</b> <br>
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/grabrepos2o.gif)<br>
At the time of dragging the LMB to set the handle location, you can grab the Bezier point to reposition it by pressing g. All the snapping options are available for setting the handle location and repositioning the Bezier point. Press g again to release the grab. <br>
<b>Undo:</b> While drawing, you can undo one segment at a time by pressing backspace. Pressing escape removes the entire curve. After the drawing is finished, the curve op can be undone by pressing ctrl-Z  <br><br>

# Flexi Edit Bezier Tool
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/editdemo.gif)<br>
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). With it You can 1) edit a Bezier curve by dragging a point on the curve 2) Move Segment endpoints and manipulate handles 3) Add or delete a vertex at any arbitrary location on the curve <br><br>
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

# Flexi Grease Bezier Tool
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/greasedemo.gif)<br>
This tool will appear on the toolshelf in Grease Pencil Draw mode. You can draw Bezier curves just as you would draw with the Flexi Draw tool. After confirming the drawing is converted to grease pencil strokes. All the snapping and locking options of the Flexi Draw are available here also.<br>
Additionaly, you can increase of decrease the resolution of the stroke using the mouse wheel or pressing + or - keys.<br>
The subdivision point visibility can be toggled by pressing h key.<br>

# Configurable Entities
![Configurable Options](https://github.com/Shriinivas/blenderbezierutils/blob/master/configitems.png)<br>
Values of the following entities are user configurable via Add-ons dialog (from Preferences->Add-on Menu).
- Bezier Toolkit Panel Tab 
- Draw Line Thickness
- Handle Point Size
- Draw Marker Size
- Axis Line Thickness
- Snap Point size<br>

# Snapping & Locking Framework
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
- View: Orientation along current viewport view axes

The Orientation affects the constraining plane and axis, as well as the reference axis for snapping to angle increment. For example, if the constraining axis are XY and the selected option in the Orienation dropdown is Active object, the point will be constrained to the XY plane of the active object local space, if it is available. 

<b>4) Snapping Origin Dropdown:</b>
The options in the dropdown are:
- Global Origin
- 3d Cursor Location
- Custom Axis Start
- Reference Line Point
- Active Object Location

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
Reference Line takes on different connotations based on the context. While drawing (Flexi Draw and Flexi Grease), the Reference Line is the segment previous to the one being drawn and the Reference Line Point is the immediate previous point. While editing, for a Bezier Point, the Reference Line is the line joining the current Bezier point with the other one of the segment (which is also the Reference Line Point). In case of handles, the Reference Line is the handle itself and the Reference Line point is the Bezier point in the handle. While moving a point on the segment, the Reference Line point is the location of the point before it is moved.

<b>Custom Axis:</b>
Custom Axis is a user defined line, that serves multiple purpose. To create a Custom Axis, make sure the selection in Snapping Origin dropdown is 'Custom Axis Origin' and rightclick the starting point, move the pointer to the end point and rightclick once again. You can snap to grid or bezier point location while creating the Custom Axis. Additionally, it's possible to define custom snapping points along the Custom Axis, using mouse wheel (or plus or minus keys).<br>
The Custom Axis can be used to define the Snappig Orientation, Snapping Origin, Custom scale and Custom snapping points.

<b>Hotkey Snapping Options (Active for the Point being Drawn / Edited):</b> 
- Holding down ctrl while moving mouse will snap the point or handle to the grid. 
- By holding down shift key the angle of the segment / handle being drawn / edited will be restricted to fixed values (0, 45, 90 etc). The reference axis for determining angle increment is the first free axis based on the selected Snapping Orientation.<br>
- Holding down alt key will snap the point being drawn / edited to the a) Bezier points of the splines within all the curve objects in the view b) vertices of the active object (if there are fewer than 500 vertices) c) face of the active object under the mouse pointer (if there are fewer than 50 faces). By default snapping to the end points joins the new curve to the curve(s) it is being snapped to. You can hold down ctrl while ending the curve (by double click or space or return key) to keep the curve separate. The curve can also be separated from the snapped curves be pressing ctrl-Z after confirming. <br>

The snapping gets adapted to the viewport zoom level. <br>

<b>Locking Options:</b>
- Pressing X, Y, Z while the curve is being drawn will lock the segment to the corresponding axis. 
- Pressing shift together with one of these buttons will lock the segment to the axes other than the one denoted by the button (e.g. shift+Z - lock to XY plane).
- Press escape to get out of the lock mode.

Snapping and locking can be combined together. So user can hold down both control and shift to snap to grid as well as restrict the angle. Likewise, user can press shift-Z to lock to XY plane and then hold shift while moving the mouse to restrict the angle between the segment end points.


When Constrain Axes dropdown has an axis-pair selection (lock to plane), pressing a single axis key will allow users to draw lines parallel to the axis (or slide the point along this line).

<b>Video Tutorials & Demos:</b><br>
<b>Overview of Flexi Draw Bezier Tool:</b> https://youtu.be/C9PXp0XHgYQ<br>
<b>Overview of Flexi Edit Bezier Tool:</b> https://youtu.be/enQHGmluQIw<br>
<b>Overview of Snapping & Locking Framework:</b> https://youtu.be/VQCXZbOq47s<br>
<b>Overview of Flexi Grease Bezier & Uniform Subdiv Op:</b> https://youtu.be/4XrjpWwLU4M<br>
<b>All Bezier Toolkit videos: </b> https://www.youtube.com/playlist?list=PLxsh4i5F_h9G6QFoPzKvBRMayz8533fSW<br>

# Credits
The functionality of editing curve by grabbing a point on it is adapted from Inkscape edit curve tool: https://gitlab.com/inkscape/inkscape/blob/master/src/ui/tool/curve-drag-point.cpp <br>
I am grateful to the authors of the module for making this great piece of code accessible to everyone.

# Known Issues
- In Flexi Draw Bezier, the part of the area under toolshelf and properties panel is excluded from drawing. Hide these elements to maximize the drawing area.<br>

# Limitations
- In Flexi Draw Bezier, snapping does not work for curves with modifiers. This is the intended functionality.<br>

In general, exercise caution when using this add-on in production, since all possible conditions have not been extensively tested.<br>
You may report bug as comment on the youtube videos or on the issues page here on Github. I will try and fix them as soon as I can.
