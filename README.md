# Blender Add-on with Bezier Utility Operations
This add-on contains the following Bezier Utilities:
- Separate Bezier Splines (object mode): Create individual objects out of the splines of the selected Curve objects. Only affects curves with multiple splines. This also works with Curves with shape keys.
- Split Bezier Segments (object mode): Create individual objects out of every segment within the Curve. This also works with Curves with shape keys.
- Mark Start Vertex (edit mode): Marks the start vertex of a closed (cyclic) splines of the selected curve objects. If the curves have shape keys, they may get distorted with change in the start vertex. <br>

Supported Blender Version: 2.8 Beta (Build dated after May 19, 2019)

# Installation
- download blenderbezierutils.py
- Open Blender and select File->User Preferences
- Click install Add-ons tab and then Install Add-on from File
- Select the downloaded file
- Check the 'Bezier Utilities' option in the add-ons dialog

After installation, a new 'Bezier Utilities' tab is displayed in object mode on 'Active Tool and Workspace settings' tab on the properties panel.

# Limitations 
Exercise caution when using this add-on in production as it's in alpha stage
