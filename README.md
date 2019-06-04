# Blender Add-on with Bezier Utility Operations
This add-on contains the following Bezier Utilities:
- Separate Bezier Splines (object mode): Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys. New objects are put in a separate collection.
- Split Bezier Segments (object mode): Create individual objects out of every segment within the Curve. This also works with Curves with shape keys. New objects are put in a separate collection.
- Select Objects in collection: Allows selection of objects belonging to the collection of the active object. It's possible to select alternate objects or objects at fixed interval (based on the order in collection) as well as invert the selection. This can be combined with the split / separate ops to work on the newly created segment / spline objects.
- Close with Straight Segment: Closes all (non-cyclic) splines within the selected curves with straight segments.
- Remove Duplicate Curve Vertices: Removes vertices at the same location. If there are duplicate end vertices coinciding with the start node the spline is marked cyclic.
- Convert to Mesh: Converts the 2d closed curve to a mesh with quad faces. This op basically applies solidify and remesh modifiers (with the specified depth) and removes the faces not belonging to the plane of the 2d curve.
- Join Bezier Curves. Joins the selected objects at their end points. If the 'Join optimized' option is unchecked the curves are joined in their order in the collection (alphabetical order of their names). If it is checked, the next curve to join is chosen based on its distance (shortest) from the current curve; the curve direction is reversed if needed. With the 'Join With Straight Segments' option, the curve objects are joined with straight line segments, regardless of the end point handles.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. <br>

Supported Blender Version: 2.8 Beta (Build dated after May 19, 2019)

<b>Demo Video: https://youtu.be/Wo-RzVI05po</b>

# Installation
- download blenderbezierutils.py
- Open Blender and select File->User Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bezier Utilities' option in the add-ons dialog

After installation, a new 'Bezier Utilities' tab is displayed in object mode on 'Active Tool and Workspace settings' tab on the properties panel.

# Limitations 
Exercise caution when using this add-on in production as it's in alpha stage
You may report bug as a comment on the youtube video or on the issues page here
