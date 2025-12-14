
import bpy
import math
import sys
import os

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

def test_hex():
    bpy.ops.wm.read_homefile(use_empty=True)
    try: 
        import bezier_utils
        import importlib
        importlib.reload(bezier_utils) # Reload to be sure
        bezier_utils.register() 
    except Exception as e: 
        print(f"Register failed: {e}")
        pass

    # Create Hexagon
    bpy.ops.curve.primitive_bezier_circle_add(radius=2, enter_editmode=True)
    bpy.ops.curve.subdivide(number_cuts=5) # Circle has 4 pts -> 24? No.
    # robust hex
    bpy.ops.curve.select_all(action='SELECT')
    bpy.ops.curve.delete(type='VERT')
    
    # Draw Hexagon manually
    bpy.ops.curve.primitive_bezier_circle_add(radius=2)
    c = bpy.context.active_object
    c.name = "Hexagon"
    # Make it a poly-like Bezier
    s = c.data.splines[0]
    # Default is BEZIER
    s.bezier_points.add(2) # Default 4 + 2 = 6
    s.use_cyclic_u = True
    
    import math
    radius = 2.0
    for i in range(6):
        angle = i * (2 * math.pi / 6)
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        pt = s.bezier_points[i]
        pt.co = (x, y, 0)
        pt.handle_left = (x, y, 0)
        pt.handle_right = (x, y, 0)
        pt.handle_left_type = 'VECTOR'
        pt.handle_right_type = 'VECTOR'
        
    c.select_set(True)
    bpy.context.view_layer.objects.active = c
    
    params = bpy.context.window_manager.bezierToolkitParams
    params.fillType = 'PAVING'
    params.fillDetail = 10 
    
    print(f"\n--- Running Paving check on {c.name} ---")
    try:
        bpy.ops.object.convert_2d_mesh()
        mesh = bpy.context.active_object
        print(f"Result Mesh: {mesh}")
        if mesh:
            print(f"Faces: {len(mesh.data.polygons)}")
            for p in mesh.data.polygons:
                print(f"Face {p.index}: {len(p.vertices)} verts")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_hex()
