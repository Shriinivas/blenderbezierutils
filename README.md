
# Blender Add-on with Bezier Utility Operations
<b>Add-on Version: 0.585 </b>

This add-on contains several tools for working with Bezier curves. <br>
Supported Blender Version: 2.8 Beta <br>

The tools are arranged in a collapsible Panel, grouped according to the functionality type. 
Here is a brief overview of a few of the tools in this add-on:
- Separate Bezier Splines: Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bezier Segments: Create individual objects out of every segment within the Curve. This also works with Curves with shape keys. New objects are put in a separate collection.
- Separate Bezier Points: Create individual objects out of the Bezier end point of each segment within the selected Curves. The newly created point objects do not inherit the shape keys of the original curve objects. New objects are put in a separate collection. The points can be used for snapping with the Flexi Bezier Tool for interesting experimentation :)
- Select Objects in collection: Allows selection of objects belonging to the collection of the active object. It's possible to select alternate objects or objects at fixed interval (based on the order in collection), as well as invert the selection. This can be combined with the split / separate ops to work on the newly created segment / spline objects.
- Close with Straight Segment: Closes all (non-cyclic) splines within the selected curves with straight segments.
- Remove Duplicate Curve Vertices: Removes vertices at the same location. If there are duplicate end vertices coinciding with the start node the spline is marked cyclic.
- Convert to Mesh: Converts the curve to a mesh with quad faces. The curve is first made 2d and all its splines cyclic. This op basically applies remesh modifier to make a quad mesh. Users can optionally check unsubdivide option to reduce the polygon count further.
- Join Bezier Curves. Joins the selected objects at their end points. If the 'Join optimized' option is unchecked the curves are joined in their order in the collection (alphabetical order of their names). If it is checked, the next curve to join is chosen based on its distance (shortest) from the current curve; the curve direction is reversed if needed. With the 'Join With Straight Segments' option, the curve objects are joined with straight line segments, regardless of the end point handle types.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. (This tool is the same as the one included in the <a href ='https://github.com/Shriinivas/assignshapekey'>Assign Shape Keys</a> add-on.)<br> 

# Flexi Bezier Tool
![Demo](https://github.com/Shriinivas/blenderbezierutils/blob/master/drawdemo.gif)<br>
This tool is available in object mode via a new button on the toolshelf (short cut to toggle the toolshelf - t). It allows drawing Bezier curves by manipulating the control points.<br>
<b>Drawing the curve: </b> To draw the curve activate the tool by clicking the Flexi Bezier tool on the toolbar. Click the LMB on the starting point of the curve. Then click and drag LMB on the end point to adjust the curvature. You can continue drawing subsequent segments in this fashion. Double clicking or hitting enter or space will convert the drawing to a curve object. <br>
<b>Undo:</b> While drawing, you can undo one segment at a time by pressing backspace. Pressing escape removes the entire curve. After the drawing is finished, the curve op can be undone by pressing ctrl-Z  <br><br>
<b>Snapping Options</b> 
- Holding down ctrl while moving mouse will snap the point or handle to the grid. 
- By holding down shift key the angle between the new point (or handle) and the previous one will be restricted to fixed values (0, 45, 90 etc). 
- To snap to the end points of the open splines within all the curve objects in the view hold down alt key. By default snapping to the end points joins the new curve to the curve(s) it is being snapped to. You can hold down ctrl while ending the curve (by double click or space or return key) to keep the curve separate. The curve can also be separated from the snapped curves be pressing ctrl-Z after confirming. <br>

The snapping gets adapted to the viewport zoom level. <br>

<b>Locking Options</b>
- Pressing X, Y, Z while the curve is being drawn will lock the segment to the corresponding axis. 
- Pressing shift together with one of these buttons will lock the segment to the axes other than the one denoted by the button (e.g. shift+Z - lock to XY plane). This option can be very useful to draw planar segments in perspective view.
- Press U to get out of the lock mode. 

Snapping and locking can be combined together. So user can hold down both control and shift to snap to grid as well as restrict the angle. Likewise, user can press shift-Z to lock to XY plane and then hold shift while moving the mouse to restrict the angle between the segment end points.

<b>Demo Video:</b> https://youtu.be/Wo-RzVI05po<br>
<b>Overview of Flexi Bezier Tool:</b> https://youtu.be/C9PXp0XHgYQ
  
# Installation
- Download blenderbezierutils.py
- Open Blender and select File->User Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bezier Utilities' option in the add-ons dialog

After installation, a new 'Bezier Utilities' tab is displayed in object mode on 'Active Tool and Workspace settings' tab on the properties panel. There will also appear a new button - Flexi Bezier Tool - on the toolshelf.

# Known Issues
A few known issues I am still working on:
- <strike>Starting a new blend while the Flexi Bezier is active causes the tool to stop working. The workaround is to select any other tool from the toolbar before starting a new file.</strike> This is fixed in version 0.51.
- <strike>Clicking on the Flexi Bezier button while it is already active starts the curve from under the button on the viewport. In this case just press escape to undo the unwanted curve. The same thing may happen while clicking on the Blender menu.</strike> This is fixed in version 0.55.<br>
- The Flexi Bezier button probably won't work after Blender restart. This is likely an issue with Blender and I have opened a bug report. Till it gets fixed, to activate the tool again in new session, just disable and enable the add-on from preferences menu (thanks Nic for this work-around).<br>
- To enable snapping, sometimes it may be required to reactivate the Flexi Bezier tool (by clicking on some other tool and clicking back on Flexi Bezier button).<br>

# Limitations
- Snapping does not work for curves with modifiers. This is intended functionality.<br>
- The part of the area under toolshelf and properties panel is excluded from drawing. Hide these elements to maximize the drawing area.<br>

In general, exercise caution when using this add-on in production, since all possible conditions have not been extensively tested.<br>
You may report bug as comment on the youtube videos or on the issues page here on Github. I will try and fix them as soon as I can.
