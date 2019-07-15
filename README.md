
# Blender Add-on with Bezier Utility Operations
This add-on contains several tools for working with Bezier curves. <br>
Supported Blender Version: 2.8 Beta <br>

Here is a brief overview of a few of the tools in this add-on:
- Separate Bezier Splines (object mode): Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys. New objects are put in a separate collection.
- Split Bezier Segments (object mode): Create individual objects out of every segment within the Curve. This also works with Curves with shape keys. New objects are put in a separate collection.
- Select Objects in collection: Allows selection of objects belonging to the collection of the active object. It's possible to select alternate objects or objects at fixed interval (based on the order in collection) as well as invert the selection. This can be combined with the split / separate ops to work on the newly created segment / spline objects.
- Close with Straight Segment: Closes all (non-cyclic) splines within the selected curves with straight segments.
- Remove Duplicate Curve Vertices: Removes vertices at the same location. If there are duplicate end vertices coinciding with the start node the spline is marked cyclic.
- Convert to Mesh: Converts the curve to a mesh with quad faces. The curve is first made 2d and all its splines cyclic. This op basically applies remesh modifier to make a quad mesh. Users can optionally check unsubdivide option to reduce the polygon count further.
- Join Bezier Curves. Joins the selected objects at their end points. If the 'Join optimized' option is unchecked the curves are joined in their order in the collection (alphabetical order of their names). If it is checked, the next curve to join is chosen based on its distance (shortest) from the current curve; the curve direction is reversed if needed. With the 'Join With Straight Segments' option, the curve objects are joined with straight line segments, regardless of the end point handles.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. <br>

# Flexi Bezier Tool
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/drawdemo.gif)<br>
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). It allows drawing Bezier curves by manipulating the control points.<br>
<b>Drawing the curve: </b> To draw the curve activate the tool by clicking the Flexi Bezier tool on the toolbar. Click the LMB on the starting point of the curve. Then click and drag LMB on the end point to adjust the curvature. You can continue drawing subsequent segments in this fashion. Double clicking or hitting enter will convert the drawing to a curve object. To cancel press escape. <br>
<b>Undo:</b> You can undo one segment with backspace and entire curve with control+z.<br>
<b>Snapping Options :</b> Holding down ctrl while moving mouse will snap the point or handle to the grid. By holding down shift key the angle between the new point (or handle) and the previous one will be restricted to fixed values (0, 45, 90 etc). To snap to the end points of the open splines within all the curve objects in the view hold down alt key. <br>
The snapping gets adapted to the viewport zoom level. <br>

<b>Demo Video:</b> https://youtu.be/Wo-RzVI05po<br>
<b>Overview of Flexi Bezier Tool:</b> https://youtu.be/C9PXp0XHgYQ
  
# Installation
- Download blenderbezierutils.py
- Open Blender and select File->User Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bezier Utilities' option in the add-ons dialog

After installation, a new 'Bezier Utilities' tab is displayed in object mode on 'Active Tool and Workspace settings' tab on the properties panel.

# Known Issues
A couple of known (hopefully minor) issues I am still working on:
- Starting a new blend while the Flexi Bezier is active causes the tool to stop working. The workaround is to select any other tool from the toolbar before starting a new file (or restart Blender to make the tool operational again :)
- Clicking on the Flexi Bezier button while it is already active starts the curve from under the button on the viewport. In this case just press escape to undo the unwanted curve. The same thing may happen while clicking on the Blender menu.<br>
- The button probably won't work after Blender restart. This is likely an issue with Blender and I have opened a bug report. Till it gets fixed, to activate the tool again in new session, just disable and enable the add-on from preferences menu (thanks Nic for this work-around).
- Snapping does not work for curves with modifiers.

In general, exercise caution when using this add-on in production, since all possible conditions have not been extensively tested.<br>
You may report bug as comment on the youtube videos or on the issues page here on Github. I will try and fix them as soon as I can.
