import bpy
import math
import sys
import os

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

def test_star():
    bpy.ops.wm.read_homefile(use_empty=True)
    try:
        import bezier_utils
        import importlib
        importlib.reload(bezier_utils)
        bezier_utils.register()
    except Exception as e:
        print(f"Register failed: {e}")
        pass

    # Create 6-pointed star
    bpy.ops.curve.primitive_bezier_circle_add(radius=2)
    c = bpy.context.active_object
    c.name = "Star6"
    s = c.data.splines[0]

    # Create 12 points (6 tips + 6 valleys)
    s.bezier_points.add(8)  # 4 default + 8 = 12
    s.use_cyclic_u = True

    # Star parameters
    outer_radius = 10.0
    inner_radius = 3.0
    n_points = 6

    for i in range(12):
        angle = i * (2 * math.pi / 12)
        # Alternate between outer (tips) and inner (valleys) radius
        radius = outer_radius if i % 2 == 0 else inner_radius
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
    params.fillDetail = 3  # Low detail for easier debugging

    print(f"\n=== Testing Sector-based TFI on 6-pointed star ===")
    print(f"Outer radius: {outer_radius}, Inner radius: {inner_radius}")
    print(f"Expected: 6 tips, 6 valleys, inner polygon with 12 vertices")
    print(f"Fill detail: {params.fillDetail}\n")

    try:
        bpy.ops.object.convert_2d_mesh()
        mesh = bpy.context.active_object
        print(f"\n=== Result ===")
        print(f"Mesh: {mesh.name if mesh else 'None'}")
        if mesh and mesh.type == 'MESH':
            tris = sum(1 for p in mesh.data.polygons if len(p.vertices) == 3)
            quads = sum(1 for p in mesh.data.polygons if len(p.vertices) == 4)
            ngons = sum(1 for p in mesh.data.polygons if len(p.vertices) > 4)
            print(f"Total faces: {len(mesh.data.polygons)}")
            print(f"Triangles: {tris}")
            print(f"Quads: {quads}")
            print(f"N-gons (>4 verts): {ngons}")

            # Show first few ngons if any
            if ngons > 0:
                print("\nN-gon details:")
                count = 0
                for p in mesh.data.polygons:
                    if len(p.vertices) > 4:
                        print(f"  Face {p.index}: {len(p.vertices)} vertices")
                        count += 1
                        if count >= 5:
                            print(f"  ... and {ngons - 5} more")
                            break
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_star()
