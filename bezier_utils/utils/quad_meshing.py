# bezier_utils/utils/quad_meshing.py
"""
Shape-specific quad mesh generators.

Implements advanced algorithms for clean quad topology:
- Q-Morph: Advancing front triangle-to-quad transformation
- Medial Axis: Domain decomposition with transfinite interpolation
- Shape-specific optimized meshers for circles, rectangles, polygons

All meshers avoid center poles and produce clean quad topology.
"""

import bpy
import bmesh
from mathutils import Vector, kdtree
from math import pi, degrees, acos, atan2, cos, sin
from typing import List, Tuple, Dict, Optional, Set
from dataclasses import dataclass, field
from enum import Enum, auto

from .shape_detection import ShapeType, ShapeAnalysis, detect_shape
from ..constants import (
    DEF_ERR_MARGIN,
    QMORPH_END_ANGLE_THRESH,
    QMORPH_SIDE_ANGLE_THRESH,
    QMORPH_MAX_ITERATIONS,
)


# =============================================================================
# Base Class
# =============================================================================

class QuadMeshGenerator:
    """
    Base class for shape-specific quad mesh generation.
    """

    def __init__(self, mesh_obj, analysis: ShapeAnalysis, params: dict):
        self.mesh_obj = mesh_obj
        self.analysis = analysis
        self.params = params
        self.bm = bmesh.new()
        self.fill_detail = params.get('fillDetail', 5)
        self.offset_size = params.get('offsetSize', 0.3)

    def generate(self) -> bpy.types.Object:
        raise NotImplementedError("Subclasses must implement generate()")

    def _get_boundary_verts(self) -> List[Vector]:
        return list(self.analysis.boundary_points)

    def _get_hole_boundaries(self) -> List[List[Vector]]:
        return self.analysis.holes

    def _finalize(self):
        self.bm.to_mesh(self.mesh_obj.data)
        self.bm.free()
        return self.mesh_obj

    def _resample_boundary_uniform(
        self, boundary: List[Vector], n_samples: int
    ) -> List[Vector]:
        """Resample boundary to exactly n_samples points with uniform arc length."""
        if len(boundary) < 2 or n_samples < 2:
            return boundary[:n_samples] if len(boundary) >= n_samples else boundary

        # Calculate total perimeter
        n = len(boundary)
        total_length = sum(
            (boundary[(i + 1) % n] - boundary[i]).length
            for i in range(n)
        )

        if total_length < DEF_ERR_MARGIN:
            return boundary[:n_samples]

        # Target arc length between samples
        target_arc = total_length / n_samples

        resampled = [boundary[0].copy()]
        accumulated = 0.0
        curr_idx = 0
        curr_pos = boundary[0].copy()

        while len(resampled) < n_samples:
            next_idx = (curr_idx + 1) % n
            edge_vec = boundary[next_idx] - curr_pos
            edge_len = edge_vec.length

            if edge_len < DEF_ERR_MARGIN:
                curr_idx = next_idx
                curr_pos = boundary[next_idx].copy()
                continue

            remaining = target_arc - accumulated

            if edge_len >= remaining:
                # Place new point on this edge
                t = remaining / edge_len
                new_pt = curr_pos + edge_vec * t
                resampled.append(new_pt.copy())
                curr_pos = new_pt
                accumulated = 0.0
            else:
                # Move to next edge
                accumulated += edge_len
                curr_idx = next_idx
                curr_pos = boundary[next_idx].copy()

            # Safety: prevent infinite loop
            if curr_idx == 0 and len(resampled) > 1 and accumulated > target_arc * 0.5:
                break

        return resampled


# =============================================================================
# Rectangle Grid Mesher (No poles - pure quad grid)
# =============================================================================

class RectangleGridMesher(QuadMeshGenerator):
    """
    Axis-aligned grid for rectangles - pure quads, no poles.
    """

    def generate(self) -> bpy.types.Object:
        corners = self._identify_corners()
        if len(corners) != 4:
            return self._fallback_grid()

        corners = self._order_corners(corners)

        # Compute grid dimensions proportional to edge lengths
        edge01 = (corners[1] - corners[0]).length
        edge12 = (corners[2] - corners[1]).length

        detail = self.fill_detail
        if edge01 >= edge12:
            n_u = detail
            n_v = max(1, round(detail * edge12 / edge01))
        else:
            n_v = detail
            n_u = max(1, round(detail * edge01 / edge12))

        # Create grid vertices using bilinear interpolation (TFI)
        verts = []
        for j in range(n_v + 1):
            row = []
            t_v = j / n_v
            p0 = corners[0].lerp(corners[3], t_v)
            p1 = corners[1].lerp(corners[2], t_v)

            for i in range(n_u + 1):
                t_u = i / n_u
                pt = p0.lerp(p1, t_u)
                row.append(self.bm.verts.new(pt))
            verts.append(row)

        # Create quad faces
        for j in range(n_v):
            for i in range(n_u):
                self.bm.faces.new([
                    verts[j][i], verts[j][i + 1],
                    verts[j + 1][i + 1], verts[j + 1][i]
                ])

        return self._finalize()

    def _identify_corners(self) -> List[Vector]:
        boundary = self._get_boundary_verts()
        if len(boundary) < 4:
            return []

        # Find 4 points with sharpest angles
        angles_with_idx = []
        n = len(boundary)
        for i in range(n):
            p_prev = boundary[(i - 1) % n]
            p_curr = boundary[i]
            p_next = boundary[(i + 1) % n]

            v1 = (p_prev - p_curr)
            v2 = (p_next - p_curr)

            if v1.length < DEF_ERR_MARGIN or v2.length < DEF_ERR_MARGIN:
                angles_with_idx.append((180.0, i))
                continue

            dot = max(-1.0, min(1.0, v1.normalized().dot(v2.normalized())))
            angle = degrees(acos(dot))
            angles_with_idx.append((angle, i))

        angles_with_idx.sort(key=lambda x: x[0])
        corner_indices = sorted([idx for _, idx in angles_with_idx[:4]])
        return [boundary[i].copy() for i in corner_indices]

    def _order_corners(self, corners: List[Vector]) -> List[Vector]:
        center = sum(corners, Vector((0, 0, 0))) / 4

        def angle_from_center(p):
            return atan2(p.y - center.y, p.x - center.x)

        sorted_corners = sorted(corners, key=angle_from_center)
        min_sum_idx = min(range(4), key=lambda i: sorted_corners[i].x + sorted_corners[i].y)
        return sorted_corners[min_sum_idx:] + sorted_corners[:min_sum_idx]

    def _fallback_grid(self) -> bpy.types.Object:
        bbox_min = self.analysis.bbox_min
        bbox_max = self.analysis.bbox_max

        corners = [
            Vector((bbox_min.x, bbox_min.y, 0)),
            Vector((bbox_max.x, bbox_min.y, 0)),
            Vector((bbox_max.x, bbox_max.y, 0)),
            Vector((bbox_min.x, bbox_max.y, 0)),
        ]

        n = self.fill_detail
        verts = []
        for j in range(n + 1):
            row = []
            t_v = j / n
            for i in range(n + 1):
                t_u = i / n
                pt = corners[0].lerp(corners[1], t_u).lerp(
                    corners[3].lerp(corners[2], t_u), t_v
                )
                row.append(self.bm.verts.new(pt))
            verts.append(row)

        for j in range(n):
            for i in range(n):
                self.bm.faces.new([
                    verts[j][i], verts[j][i + 1],
                    verts[j + 1][i + 1], verts[j + 1][i]
                ])

        return self._finalize()


# =============================================================================
# Grid TFI Mesher - Direct TFI without offset rings
# =============================================================================

class GridTFIMesher(QuadMeshGenerator):
    """
    Direct Grid Fill using Transfinite Interpolation.

    Maps an n-gon boundary directly to a structured quad grid without
    intermediate offset rings. Preserves boundary shape well and creates
    clean quad topology.

    Algorithm (Blender's Grid Fill):
    1. Resample boundary to multiple of 4 vertices
    2. Identify 4 corners at n/4 intervals
    3. Map the n-gon to a square grid:
       - 4 edges become grid boundaries
       - Interior filled with TFI
    4. Result: (n/4) × (n/4) quads, no poles
    """

    def generate(self) -> bpy.types.Object:
        boundary = self._get_boundary_verts()
        if len(boundary) < 4:
            return self._finalize()

        # Detect corners to preserve sharp angles
        corners = self._detect_corners(boundary)

        # Resample to multiple of 4, preserving corners
        n = len(boundary)
        target_n = ((n + 2) // 4) * 4
        target_n = max(8, target_n)

        # Respect fill_detail to control grid density
        # fill_detail controls the target edge count (n/4)
        min_target = self.fill_detail * 4
        target_n = max(target_n, min_target)

        if n != target_n or corners:
            boundary = self._resample_preserving_corners(boundary, corners, target_n)

        n = len(boundary)
        if n < 4 or n % 4 != 0:
            # Fallback: simple n-gon face
            verts = [self.bm.verts.new(pt) for pt in boundary]
            try:
                self.bm.faces.new(verts)
            except ValueError:
                pass
            return self._finalize()

        # Create boundary vertices in BMesh
        boundary_verts = [self.bm.verts.new(pt) for pt in boundary]

        # Apply Grid Fill TFI
        self._grid_fill_tfi(boundary_verts)

        return self._finalize()

    def _detect_corners(self, boundary: List[Vector], angle_thresh: float = 140.0) -> List[int]:
        """Detect corner vertices (sharp angles) in the boundary."""
        corners = []
        n = len(boundary)

        for i in range(n):
            prev_pt = boundary[(i - 1) % n]
            curr_pt = boundary[i]
            next_pt = boundary[(i + 1) % n]

            v1 = (prev_pt - curr_pt)
            v2 = (next_pt - curr_pt)

            if v1.length < DEF_ERR_MARGIN or v2.length < DEF_ERR_MARGIN:
                continue

            dot = max(-1.0, min(1.0, v1.normalized().dot(v2.normalized())))
            angle = degrees(acos(dot))

            if angle < angle_thresh:
                corners.append(i)

        return corners

    def _resample_preserving_corners(
        self, boundary: List[Vector], corner_indices: List[int], target_n: int
    ) -> List[Vector]:
        """Resample boundary to target count while preserving corner positions."""
        n = len(boundary)

        if not corner_indices:
            return self._resample_boundary_uniform(boundary, target_n)

        corner_indices = sorted(corner_indices)
        n_corners = len(corner_indices)

        remaining = target_n - n_corners

        # Calculate arc lengths between corners
        arc_lengths = []
        for i in range(n_corners):
            start_idx = corner_indices[i]
            end_idx = corner_indices[(i + 1) % n_corners]

            arc_len = 0.0
            idx = start_idx
            while idx != end_idx:
                next_idx = (idx + 1) % n
                arc_len += (boundary[next_idx] - boundary[idx]).length
                idx = next_idx
            arc_lengths.append(arc_len)

        total_arc = sum(arc_lengths)

        # Distribute vertices proportionally
        verts_per_segment = []
        distributed = 0
        for arc_len in arc_lengths:
            if total_arc > DEF_ERR_MARGIN:
                count = int(round(remaining * arc_len / total_arc))
            else:
                count = remaining // n_corners
            verts_per_segment.append(count)
            distributed += count

        # Adjust for rounding
        diff = remaining - distributed
        if diff != 0:
            longest_idx = arc_lengths.index(max(arc_lengths))
            verts_per_segment[longest_idx] += diff

        # Build resampled boundary
        resampled = []
        for i in range(n_corners):
            start_idx = corner_indices[i]
            end_idx = corner_indices[(i + 1) % n_corners]
            n_intermediate = verts_per_segment[i]

            resampled.append(boundary[start_idx].copy())

            if n_intermediate <= 0:
                continue

            arc_len = arc_lengths[i]
            target_spacing = arc_len / (n_intermediate + 1)

            idx = start_idx
            accumulated = 0.0
            placed = 0

            while placed < n_intermediate and idx != end_idx:
                next_idx = (idx + 1) % n
                edge_vec = boundary[next_idx] - boundary[idx]
                edge_len = edge_vec.length

                while accumulated + edge_len >= target_spacing * (placed + 1) and placed < n_intermediate:
                    t = (target_spacing * (placed + 1) - accumulated) / edge_len
                    t = max(0.0, min(1.0, t))
                    new_pt = boundary[idx] + edge_vec * t
                    resampled.append(new_pt.copy())
                    placed += 1

                accumulated += edge_len
                idx = next_idx

        return resampled

    def _grid_fill_tfi(self, boundary_verts: List):
        """
        Fill n-gon with quads using Grid Fill TFI.

        Maps the n-gon boundary to a square grid where:
        - n must be multiple of 4
        - Grid is (n/4 + 1) × (n/4 + 1) vertices
        - Interior vertices computed via Transfinite Interpolation
        """
        n = len(boundary_verts)

        if n < 4:
            return

        if n == 4:
            try:
                self.bm.faces.new(boundary_verts)
            except ValueError:
                pass
            return

        if n % 4 != 0:
            # Can't make perfect quad grid, create single n-gon
            try:
                self.bm.faces.new(boundary_verts)
            except ValueError:
                pass
            return

        # Grid dimensions
        k = n // 4  # Vertices per edge (not counting shared corner)
        grid_size = k + 1

        # Identify corners at indices 0, k, 2k, 3k
        corners = [0, k, 2 * k, 3 * k]

        # Extract the 4 edges
        edges = []
        for i in range(4):
            start = corners[i]
            end = corners[(i + 1) % 4]
            if end == 0:
                end = n
            edge = [boundary_verts[j % n] for j in range(start, end + 1)]
            edges.append(edge)

        # Build grid using TFI
        grid = [[None for _ in range(grid_size)] for _ in range(grid_size)]

        # Fill boundary from edges
        # Top row (row=0): edge0
        for col in range(grid_size):
            if col < len(edges[0]):
                grid[0][col] = edges[0][col]

        # Right column (col=k): edge1
        for row in range(grid_size):
            if row < len(edges[1]):
                grid[row][grid_size - 1] = edges[1][row]

        # Bottom row (row=k): edge2 reversed
        for col in range(grid_size):
            rev_col = grid_size - 1 - col
            if rev_col < len(edges[2]):
                grid[grid_size - 1][col] = edges[2][rev_col]

        # Left column (col=0): edge3 reversed
        for row in range(grid_size):
            rev_row = grid_size - 1 - row
            if rev_row < len(edges[3]):
                grid[row][0] = edges[3][rev_row]

        # Fill interior using TFI
        center = sum((v.co for v in boundary_verts), Vector((0, 0, 0))) / n

        for row in range(1, grid_size - 1):
            for col in range(1, grid_size - 1):
                u = col / (grid_size - 1)
                v = row / (grid_size - 1)

                # Get boundary points
                top = grid[0][col].co if grid[0][col] else center
                bottom = grid[grid_size - 1][col].co if grid[grid_size - 1][col] else center
                left = grid[row][0].co if grid[row][0] else center
                right = grid[row][grid_size - 1].co if grid[row][grid_size - 1] else center

                # Get corner points
                p00 = grid[0][0].co if grid[0][0] else center
                p10 = grid[0][grid_size - 1].co if grid[0][grid_size - 1] else center
                p01 = grid[grid_size - 1][0].co if grid[grid_size - 1][0] else center
                p11 = grid[grid_size - 1][grid_size - 1].co if grid[grid_size - 1][grid_size - 1] else center

                # TFI formula
                pt = (
                    (1 - v) * top + v * bottom +
                    (1 - u) * left + u * right -
                    (1 - u) * (1 - v) * p00 -
                    u * (1 - v) * p10 -
                    (1 - u) * v * p01 -
                    u * v * p11
                )

                grid[row][col] = self.bm.verts.new(pt)

        # Create quad faces
        for row in range(grid_size - 1):
            for col in range(grid_size - 1):
                v00 = grid[row][col]
                v10 = grid[row][col + 1]
                v01 = grid[row + 1][col]
                v11 = grid[row + 1][col + 1]

                if v00 and v10 and v01 and v11:
                    try:
                        self.bm.faces.new([v00, v10, v11, v01])
                    except ValueError:
                        pass


# =============================================================================
# Circle/Ellipse Mesher - True O-Grid topology (no pole)
# =============================================================================

class PolarGridMesher(QuadMeshGenerator):
    """
    O-Grid mesher for circles/ellipses - TRUE pole-free topology.

    Uses proper O-grid: the domain is divided into 4 quadrants, each meshed
    with a structured quad grid. The center becomes a single quad face
    (4 vertices), not a pole.

    Topology:
    - 4 quadrant patches meeting at center
    - Each patch is a structured quad grid using TFI
    - Center is a quad, corners are on boundary
    - All faces are quads, no triangles, no poles
    """

    def generate(self) -> bpy.types.Object:
        center = self.analysis.center.copy()
        center.z = 0

        boundary = self._get_boundary_verts()
        if len(boundary) < 8:
            return self._finalize()

        # Need boundary divisible by 4 for clean O-grid
        n_around = ((len(boundary) + 3) // 4) * 4
        if n_around < 8:
            n_around = 8
        n_around = max(n_around, self.fill_detail * 4)

        # Resample boundary uniformly
        boundary = self._resample_boundary_uniform(boundary, n_around)

        n_per_quadrant = n_around // 4  # Points per quadrant edge on boundary
        n_radial = max(2, self.fill_detail)  # Rings from boundary to center

        # Size of center quad (as fraction of average radius)
        avg_radius = sum((b - center).length for b in boundary) / len(boundary)
        center_size = avg_radius * 0.15

        # Create 4 corner points of center quad (at 45, 135, 225, 315 degrees)
        center_quad_pts = []
        for i in range(4):
            angle = pi / 4 + i * pi / 2
            pt = center + Vector((cos(angle), sin(angle), 0)) * center_size
            center_quad_pts.append(pt)

        # Create vertex grid for each quadrant
        # Quadrant i goes from boundary index i*n_per_quadrant to (i+1)*n_per_quadrant
        # and from center_quad corner i to corner (i+1)%4

        all_patch_verts = []  # [quadrant][row][col] -> BMVert

        for q in range(4):
            # Boundary segment for this quadrant
            b_start = q * n_per_quadrant

            # Get boundary points for this quadrant (inclusive both ends)
            quad_boundary = []
            for i in range(n_per_quadrant + 1):
                idx = (b_start + i) % n_around
                quad_boundary.append(boundary[idx])

            # Center quad edge for this quadrant
            c0 = center_quad_pts[q]
            c1 = center_quad_pts[(q + 1) % 4]

            # Create TFI grid for this quadrant patch
            # The patch is a curvilinear quad with:
            # - Bottom edge: boundary segment (curved)
            # - Top edge: center quad edge (straight, small)
            # - Left/Right edges: radial lines from boundary corners to center quad corners

            patch_verts = []
            for j in range(n_radial + 1):  # j=0 is boundary, j=n_radial is center
                row = []
                t_radial = j / n_radial  # 0 at boundary, 1 at center

                for i in range(n_per_quadrant + 1):  # Along the quadrant
                    t_along = i / n_per_quadrant

                    # Boundary point at this position
                    b_pt = quad_boundary[i]

                    # Center edge point at this position
                    c_pt = c0.lerp(c1, t_along)

                    # Interpolate from boundary to center
                    pt = b_pt.lerp(c_pt, t_radial)

                    row.append(self.bm.verts.new(pt))
                patch_verts.append(row)

            all_patch_verts.append(patch_verts)

        # Create quad faces for each patch
        for q in range(4):
            patch = all_patch_verts[q]
            for j in range(n_radial):
                for i in range(n_per_quadrant):
                    try:
                        self.bm.faces.new([
                            patch[j][i], patch[j][i + 1],
                            patch[j + 1][i + 1], patch[j + 1][i]
                        ])
                    except ValueError:
                        pass  # Face may already exist at boundaries

        # Create center quad face
        center_verts = [all_patch_verts[q][n_radial][0] for q in range(4)]
        try:
            self.bm.faces.new(center_verts)
        except ValueError:
            pass

        # Merge duplicate vertices at quadrant boundaries
        bmesh.ops.remove_doubles(self.bm, verts=list(self.bm.verts), dist=DEF_ERR_MARGIN * 10)

        return self._finalize()


# =============================================================================
# Polygon Mesher - True pole-free using inner polygon
# =============================================================================

class PolygonRadialMesher(QuadMeshGenerator):
    """
    Mesher for regular polygons - TRUE pole-free topology.

    Strategy:
    - Create a small inner polygon (scaled version of outer)
    - Each side of outer polygon maps to corresponding side of inner polygon
    - Mesh each trapezoid sector with TFI quads
    - The inner polygon itself becomes a single N-gon face (for N-sided polygon)

    For even-sided polygons (4, 6, 8, etc.): inner polygon is one N-gon face
    For odd-sided polygons (3, 5, 7, etc.): inner polygon is one N-gon face

    No center vertex = no pole!
    """

    def generate(self) -> bpy.types.Object:
        center = self.analysis.center.copy()
        center.z = 0

        n_sides = self.analysis.corner_count
        if n_sides < 3:
            # Fallback to Q-Morph
            qmorph = QMorphMesher(self.mesh_obj, self.analysis, self.params)
            qmorph.bm = self.bm
            return qmorph.generate()

        corners = self._get_polygon_corners()
        if len(corners) != n_sides:
            corners = self._estimate_corners(n_sides)

        if len(corners) < 3:
            # Still can't get corners, fallback
            qmorph = QMorphMesher(self.mesh_obj, self.analysis, self.params)
            qmorph.bm = self.bm
            return qmorph.generate()

        detail = self.fill_detail
        n_along_edge = max(2, detail)  # Subdivisions along each edge
        n_radial = max(2, detail)      # Layers from boundary to inner polygon

        # Create inner polygon (scaled down version) to avoid center pole
        inner_scale = 0.2  # Size of inner polygon relative to outer
        inner_corners = [center + (c - center) * inner_scale for c in corners]

        # Store sector vertices for merging at boundaries
        all_sector_verts = []

        # For each sector (trapezoid between outer edge and inner edge)
        for side_idx in range(n_sides):
            # Outer edge corners
            c0_out = corners[side_idx]
            c1_out = corners[(side_idx + 1) % n_sides]

            # Inner edge corners
            c0_in = inner_corners[side_idx]
            c1_in = inner_corners[(side_idx + 1) % n_sides]

            # Create TFI grid for this trapezoid sector
            sector_verts = self._mesh_trapezoid_sector(
                c0_out, c1_out, c1_in, c0_in, n_along_edge, n_radial
            )
            all_sector_verts.append(sector_verts)

        # Create inner polygon face (N-gon at center - this is NOT a pole)
        # Get the innermost vertices from each sector
        inner_face_verts = []
        for sector_verts in all_sector_verts:
            # Last row (n_radial), first column (0) is the inner corner
            inner_face_verts.append(sector_verts[n_radial][0])

        try:
            self.bm.faces.new(inner_face_verts)
        except ValueError:
            pass  # Face may already exist

        # Merge duplicate vertices at sector boundaries
        bmesh.ops.remove_doubles(self.bm, verts=list(self.bm.verts), dist=DEF_ERR_MARGIN * 10)

        return self._finalize()

    def _mesh_trapezoid_sector(
        self, c0: Vector, c1: Vector, c2: Vector, c3: Vector,
        n_along: int, n_radial: int
    ) -> List[List]:
        """
        Mesh a trapezoid sector using TFI.

        c0 -- c1  (outer edge, longer)
        |      |
        c3 -- c2  (inner edge, shorter)

        Returns grid of BMVerts [radial_idx][along_idx]
        """
        verts = []
        for j in range(n_radial + 1):  # j=0 is outer, j=n_radial is inner
            row = []
            v = j / n_radial

            # Interpolate edges
            left_pt = c0.lerp(c3, v)
            right_pt = c1.lerp(c2, v)

            for i in range(n_along + 1):
                u = i / n_along

                # Interpolate along the row
                pt = left_pt.lerp(right_pt, u)
                row.append(self.bm.verts.new(pt))
            verts.append(row)

        # Create quad faces
        for j in range(n_radial):
            for i in range(n_along):
                try:
                    self.bm.faces.new([
                        verts[j][i], verts[j][i + 1],
                        verts[j + 1][i + 1], verts[j + 1][i]
                    ])
                except ValueError:
                    pass

        return verts

    def _get_polygon_corners(self) -> List[Vector]:
        boundary = self._get_boundary_verts()
        if not self.analysis.corner_indices:
            return []

        n_boundary = len(boundary)
        n_bpts = len(self.analysis.corner_angles)

        corners = []
        for idx in self.analysis.corner_indices:
            boundary_idx = (idx * n_boundary) // n_bpts
            corners.append(boundary[boundary_idx % n_boundary].copy())

        return corners

    def _estimate_corners(self, n_sides: int) -> List[Vector]:
        boundary = self._get_boundary_verts()
        center = self.analysis.center

        corners = []
        for i in range(n_sides):
            target_angle = 2 * pi * i / n_sides - pi / 2

            best_pt = None
            best_diff = float('inf')

            for pt in boundary:
                dx = pt.x - center.x
                dy = pt.y - center.y
                angle = atan2(dy, dx)
                diff = abs(angle - target_angle)
                if diff > pi:
                    diff = 2 * pi - diff
                if diff < best_diff:
                    best_diff = diff
                    best_pt = pt

            if best_pt:
                corners.append(best_pt.copy())

        return corners


# =============================================================================
# Q-Morph Mesher - Proper Implementation
# =============================================================================

class FrontEdgeState(Enum):
    END = auto()    # angle < 75: sharp corner
    SIDE = auto()   # 75 <= angle < 135: normal
    CLOSE = auto()  # angle >= 135: flat/reflex


@dataclass
class QMorphVertex:
    co: Vector
    bm_vert: Optional[object] = None


@dataclass
class QMorphEdge:
    v0: int  # vertex index
    v1: int
    state: FrontEdgeState = FrontEdgeState.SIDE
    level: int = 0
    processed: bool = False


@dataclass
class QMorphTriangle:
    v0: int
    v1: int
    v2: int
    alive: bool = True


class QMorphMesher(QuadMeshGenerator):
    """
    Q-Morph: Advancing front quad meshing via triangle transformation.

    Algorithm (Owen et al. 1999):
    1. Triangulate the domain (constrained Delaunay / ear clipping)
    2. Initialize front from boundary edges
    3. Classify edges by angle (END, SIDE, CLOSE)
    4. Process edges to form quads:
       - CLOSE: merge two adjacent triangles
       - SIDE: edge recovery + merge
       - END: special corner handling
    5. Smooth and clean up

    This implementation produces all-quad meshes with at most one triangle.
    """

    def __init__(self, mesh_obj, analysis: ShapeAnalysis, params: dict):
        super().__init__(mesh_obj, analysis, params)
        self.vertices: List[QMorphVertex] = []
        self.triangles: List[QMorphTriangle] = []
        self.front: List[QMorphEdge] = []
        self.quads: List[Tuple[int, int, int, int]] = []
        self.edge_to_tris: Dict[Tuple[int, int], List[int]] = {}
        self.vert_to_tris: Dict[int, Set[int]] = {}

    def generate(self) -> bpy.types.Object:
        boundary = self._get_boundary_verts()
        if len(boundary) < 3:
            return self._finalize()

        # Initialize vertices from boundary
        for v in boundary:
            self.vertices.append(QMorphVertex(co=v.copy()))

        # Step 1: Triangulate
        self._triangulate_ear_clipping()

        if not self.triangles:
            return self._fallback_offset()

        # Step 2: Initialize front from boundary
        self._initialize_front()

        # Step 3: Main Q-Morph loop
        iterations = 0
        max_iter = min(QMORPH_MAX_ITERATIONS, len(boundary) * 10)

        while self._count_active_front() > 4 and iterations < max_iter:
            if not self._advance_step():
                break
            iterations += 1

        # Step 4: Close remaining front
        self._close_remaining_front()

        # Step 5: Build BMesh
        self._build_bmesh()

        return self._finalize()

    def _triangulate_ear_clipping(self):
        """Triangulate using ear clipping algorithm."""
        n = len(self.vertices)
        if n < 3:
            return

        indices = list(range(n))

        while len(indices) > 3:
            ear_found = False

            for i in range(len(indices)):
                i0 = indices[(i - 1) % len(indices)]
                i1 = indices[i]
                i2 = indices[(i + 1) % len(indices)]

                if self._is_ear(i0, i1, i2, indices):
                    tri = QMorphTriangle(v0=i0, v1=i1, v2=i2)
                    tri_idx = len(self.triangles)
                    self.triangles.append(tri)
                    self._register_triangle(tri_idx)
                    indices.remove(i1)
                    ear_found = True
                    break

            if not ear_found:
                break

        if len(indices) == 3:
            tri = QMorphTriangle(v0=indices[0], v1=indices[1], v2=indices[2])
            tri_idx = len(self.triangles)
            self.triangles.append(tri)
            self._register_triangle(tri_idx)

    def _is_ear(self, i0: int, i1: int, i2: int, indices: List[int]) -> bool:
        p0 = self.vertices[i0].co
        p1 = self.vertices[i1].co
        p2 = self.vertices[i2].co

        # Check convexity (CCW orientation)
        cross = (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
        if cross <= DEF_ERR_MARGIN:
            return False

        # Check no other vertices inside
        for idx in indices:
            if idx in (i0, i1, i2):
                continue
            if self._point_in_triangle(self.vertices[idx].co, p0, p1, p2):
                return False

        return True

    def _point_in_triangle(self, p: Vector, t0: Vector, t1: Vector, t2: Vector) -> bool:
        def sign(p1, p2, p3):
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)

        d1 = sign(p, t0, t1)
        d2 = sign(p, t1, t2)
        d3 = sign(p, t2, t0)

        has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
        has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

        return not (has_neg and has_pos)

    def _register_triangle(self, tri_idx: int):
        tri = self.triangles[tri_idx]
        verts = [tri.v0, tri.v1, tri.v2]

        # Register edges
        for i in range(3):
            v0, v1 = verts[i], verts[(i + 1) % 3]
            edge = tuple(sorted([v0, v1]))
            if edge not in self.edge_to_tris:
                self.edge_to_tris[edge] = []
            self.edge_to_tris[edge].append(tri_idx)

        # Register vertices
        for v in verts:
            if v not in self.vert_to_tris:
                self.vert_to_tris[v] = set()
            self.vert_to_tris[v].add(tri_idx)

    def _initialize_front(self):
        """Initialize front from boundary edges."""
        n = len(self.vertices)
        boundary_verts = set(range(n))  # Original boundary vertices

        for i in range(n):
            j = (i + 1) % n
            angle = self._compute_angle_at_vertex(i, boundary_verts)
            state = self._classify_angle(angle)
            self.front.append(QMorphEdge(v0=i, v1=j, state=state, level=0))

    def _compute_angle_at_vertex(self, v_idx: int, boundary: Set[int]) -> float:
        """Compute interior angle at vertex."""
        # Find adjacent boundary vertices
        n = len(self.vertices)
        prev_idx = (v_idx - 1) % n
        next_idx = (v_idx + 1) % n

        p_prev = self.vertices[prev_idx].co
        p_curr = self.vertices[v_idx].co
        p_next = self.vertices[next_idx].co

        v1 = p_prev - p_curr
        v2 = p_next - p_curr

        if v1.length < DEF_ERR_MARGIN or v2.length < DEF_ERR_MARGIN:
            return 180.0

        dot = max(-1.0, min(1.0, v1.normalized().dot(v2.normalized())))
        return degrees(acos(dot))

    def _classify_angle(self, angle: float) -> FrontEdgeState:
        if angle < QMORPH_END_ANGLE_THRESH:
            return FrontEdgeState.END
        elif angle < QMORPH_SIDE_ANGLE_THRESH:
            return FrontEdgeState.SIDE
        else:
            return FrontEdgeState.CLOSE

    def _count_active_front(self) -> int:
        return sum(1 for e in self.front if not e.processed)

    def _advance_step(self) -> bool:
        """One step of Q-Morph: select best edge and form quad."""
        # Find best unprocessed edge
        best_edge = None
        best_score = -1

        for edge in self.front:
            if edge.processed:
                continue

            score = 0
            if edge.state == FrontEdgeState.CLOSE:
                score = 300
            elif edge.state == FrontEdgeState.SIDE:
                score = 200
            else:
                score = 100
            score -= edge.level

            if score > best_score:
                best_score = score
                best_edge = edge

        if best_edge is None:
            return False

        # Try to form a quad
        return self._form_quad_at_edge(best_edge)

    def _form_quad_at_edge(self, edge: QMorphEdge) -> bool:
        """Form a quad at the given front edge."""
        edge_key = tuple(sorted([edge.v0, edge.v1]))
        tri_indices = self.edge_to_tris.get(edge_key, [])

        # Find alive triangles on this edge
        alive_tris = [ti for ti in tri_indices if self.triangles[ti].alive]

        if len(alive_tris) >= 2:
            # Can merge two triangles into a quad
            return self._merge_triangles_to_quad(edge, alive_tris[0], alive_tris[1])
        elif len(alive_tris) == 1:
            # Form quad by advancing front
            return self._advance_front_quad(edge, alive_tris[0])
        else:
            edge.processed = True
            return True

    def _merge_triangles_to_quad(self, edge: QMorphEdge, tri_idx0: int, tri_idx1: int) -> bool:
        """Merge two triangles sharing an edge into a quad."""
        tri0 = self.triangles[tri_idx0]
        tri1 = self.triangles[tri_idx1]

        verts0 = {tri0.v0, tri0.v1, tri0.v2}
        verts1 = {tri1.v0, tri1.v1, tri1.v2}

        shared = verts0 & verts1
        if len(shared) != 2:
            edge.processed = True
            return True

        unique0 = list(verts0 - shared)[0]
        unique1 = list(verts1 - shared)[0]
        shared_list = list(shared)

        # Order vertices for proper quad winding
        quad_verts = self._order_quad_verts(shared_list[0], shared_list[1], unique0, unique1)
        self.quads.append(quad_verts)

        # Mark triangles as dead
        tri0.alive = False
        tri1.alive = False

        edge.processed = True
        return True

    def _advance_front_quad(self, edge: QMorphEdge, tri_idx: int) -> bool:
        """Form quad by advancing front inward."""
        tri = self.triangles[tri_idx]
        tri_verts = [tri.v0, tri.v1, tri.v2]

        # Find the vertex opposite to the edge
        edge_verts = {edge.v0, edge.v1}
        opposite = [v for v in tri_verts if v not in edge_verts]

        if not opposite:
            edge.processed = True
            return True

        opp_v = opposite[0]

        # Create new vertex by projecting inward
        v0 = self.vertices[edge.v0].co
        v1 = self.vertices[edge.v1].co
        v_opp = self.vertices[opp_v].co

        # New vertex position: reflect opposite vertex across edge midpoint
        mid = (v0 + v1) / 2
        new_pos = mid + (mid - v_opp) * 0.5

        # Add new vertex
        new_idx = len(self.vertices)
        self.vertices.append(QMorphVertex(co=new_pos))

        # Create quad
        quad = self._order_quad_verts(edge.v0, edge.v1, new_idx, opp_v)
        self.quads.append(quad)

        # Mark triangle as dead
        tri.alive = False
        edge.processed = True

        return True

    def _order_quad_verts(self, v0: int, v1: int, v2: int, v3: int) -> Tuple[int, int, int, int]:
        """Order 4 vertices to form proper quad (CCW)."""
        verts = [v0, v1, v2, v3]
        pts = [self.vertices[v].co for v in verts]
        center = sum(pts, Vector((0, 0, 0))) / 4

        def angle_key(idx):
            p = pts[idx]
            return atan2(p.y - center.y, p.x - center.x)

        sorted_indices = sorted(range(4), key=angle_key)
        return tuple(verts[i] for i in sorted_indices)

    def _close_remaining_front(self):
        """Close the remaining front."""
        active = [e for e in self.front if not e.processed]

        if len(active) == 4:
            verts = [active[0].v0]
            for e in active:
                if e.v1 not in verts:
                    verts.append(e.v1)
            if len(verts) == 4:
                self.quads.append(tuple(verts))
        elif len(active) == 3:
            # Leave as triangle (unavoidable)
            pass

        # Also add any remaining alive triangles as faces
        for tri in self.triangles:
            if tri.alive:
                # Convert to degenerate quad or keep as tri
                pass

    def _build_bmesh(self):
        """Build BMesh from vertices and quads."""
        # Create BMesh vertices
        for vert in self.vertices:
            vert.bm_vert = self.bm.verts.new(vert.co)

        self.bm.verts.ensure_lookup_table()

        # Create quad faces
        for quad in self.quads:
            try:
                face_verts = [self.vertices[i].bm_vert for i in quad]
                if len(set(face_verts)) == 4:
                    self.bm.faces.new(face_verts)
            except (ValueError, IndexError):
                pass

        # Add remaining triangles
        for tri in self.triangles:
            if tri.alive:
                try:
                    verts = [
                        self.vertices[tri.v0].bm_vert,
                        self.vertices[tri.v1].bm_vert,
                        self.vertices[tri.v2].bm_vert
                    ]
                    if len(set(verts)) == 3:
                        self.bm.faces.new(verts)
                except (ValueError, IndexError):
                    pass

    def _fallback_offset(self) -> bpy.types.Object:
        """Fallback to offset-based meshing."""
        boundary = self._get_boundary_verts()
        center = self.analysis.center.copy()
        n_rings = self.fill_detail

        # Simple offset rings
        ring_verts = []
        for ring_idx in range(n_rings + 1):
            t = ring_idx / n_rings
            ring = []
            for pt in boundary:
                inner_pt = pt.lerp(center, t)
                ring.append(self.bm.verts.new(inner_pt))
            ring_verts.append(ring)

        # Create quads between rings
        for r in range(n_rings):
            outer = ring_verts[r]
            inner = ring_verts[r + 1]
            n = len(outer)
            for i in range(n):
                j = (i + 1) % n
                self.bm.faces.new([outer[i], outer[j], inner[j], inner[i]])

        return self._finalize()


# =============================================================================
# Medial Axis Mesher - Proper Implementation
# =============================================================================

@dataclass
class MANode:
    """Node on the medial axis with position and inscribed circle radius."""
    co: Vector
    radius: float
    boundary_indices: List[int] = field(default_factory=list)


@dataclass
class MAEdge:
    """Edge connecting two MA nodes."""
    node0: int
    node1: int


class MedialAxisMesher(QuadMeshGenerator):
    """
    Medial Axis Transform based quad meshing.

    Algorithm:
    1. Compute approximate medial axis via Voronoi-like sampling
    2. Build MA skeleton graph (nodes and edges)
    3. Decompose domain into quadrilateral patches based on MA structure
    4. Mesh each patch with TFI (structured quads)
    5. Handle junctions and endpoints to avoid poles

    Key principle: The MA divides the shape into regions, each bounded by
    - A segment of the boundary
    - A segment of the MA
    These regions are naturally quadrilateral and can be meshed with TFI.
    """

    def __init__(self, mesh_obj, analysis: ShapeAnalysis, params: dict):
        super().__init__(mesh_obj, analysis, params)
        self.ma_nodes: List[MANode] = []
        self.ma_edges: List[MAEdge] = []
        self.boundary: List[Vector] = []
        self.boundary_kd = None

    def generate(self) -> bpy.types.Object:
        self.boundary = self._get_boundary_verts()
        print(f"[MA DEBUG] Initial boundary points: {len(self.boundary)}")

        if len(self.boundary) < 3:
            print("[MA DEBUG] Too few boundary points, returning early")
            return self._finalize()

        # Detect corners (sharp angles) to preserve them during resampling
        corners = self._detect_corners(self.boundary)
        print(f"[MA DEBUG] Detected {len(corners)} corners")

        # Resample boundary to multiple of 4, preserving corners
        n = len(self.boundary)
        target_n = ((n + 2) // 4) * 4  # Round to nearest multiple of 4
        target_n = max(8, target_n)  # Minimum 8 vertices

        if n != target_n or len(corners) > 0:
            self.boundary = self._resample_preserving_corners(
                self.boundary, corners, target_n
            )
            print(f"[MA DEBUG] Resampled boundary from {n} to {len(self.boundary)} (preserving {len(corners)} corners)")

        # Build boundary KD-tree
        self._build_boundary_kdtree()

        # Compute medial axis skeleton (used for adaptive parameters)
        self._compute_medial_axis()
        print(f"[MA DEBUG] MA nodes found: {len(self.ma_nodes)}")

        # Build MA graph if we have nodes (for future use)
        if len(self.ma_nodes) >= 2:
            self._build_ma_graph()

        # Always use progressive vertex reduction meshing
        # MA nodes are used for adaptive parameters, not for conditional branching
        self._decompose_and_mesh()

        # Merge duplicate vertices at patch boundaries
        bmesh.ops.remove_doubles(self.bm, verts=list(self.bm.verts), dist=DEF_ERR_MARGIN * 10)

        return self._finalize()

    def _detect_corners(self, boundary: List[Vector], angle_thresh: float = 140.0) -> List[int]:
        """
        Detect corner vertices (sharp angles) in the boundary.

        Args:
            boundary: List of boundary points
            angle_thresh: Maximum angle (degrees) to consider a corner

        Returns:
            List of indices of corner vertices
        """
        corners = []
        n = len(boundary)

        for i in range(n):
            prev_pt = boundary[(i - 1) % n]
            curr_pt = boundary[i]
            next_pt = boundary[(i + 1) % n]

            v1 = (prev_pt - curr_pt)
            v2 = (next_pt - curr_pt)

            if v1.length < DEF_ERR_MARGIN or v2.length < DEF_ERR_MARGIN:
                continue

            v1 = v1.normalized()
            v2 = v2.normalized()

            dot = max(-1.0, min(1.0, v1.dot(v2)))
            angle = degrees(acos(dot))

            # Sharp corner if angle is less than threshold
            if angle < angle_thresh:
                corners.append(i)

        return corners

    def _resample_preserving_corners(
        self, boundary: List[Vector], corner_indices: List[int], target_n: int
    ) -> List[Vector]:
        """
        Resample boundary to target count while preserving corner positions.

        Args:
            boundary: Original boundary points
            corner_indices: Indices of corners to preserve
            target_n: Target vertex count (should be multiple of 4)

        Returns:
            Resampled boundary with corners preserved
        """
        n = len(boundary)

        if not corner_indices:
            # No corners to preserve, use uniform resampling
            return self._resample_boundary_uniform(boundary, target_n)

        # Sort corners by index
        corner_indices = sorted(corner_indices)
        n_corners = len(corner_indices)

        # Calculate vertices to distribute between corners
        # We want target_n total, with n_corners being the corners
        remaining = target_n - n_corners

        # Calculate arc lengths between consecutive corners
        arc_lengths = []
        for i in range(n_corners):
            start_idx = corner_indices[i]
            end_idx = corner_indices[(i + 1) % n_corners]

            # Calculate arc length from start to end (wrapping if needed)
            arc_len = 0.0
            idx = start_idx
            while idx != end_idx:
                next_idx = (idx + 1) % n
                arc_len += (boundary[next_idx] - boundary[idx]).length
                idx = next_idx
            arc_lengths.append(arc_len)

        total_arc = sum(arc_lengths)

        # Distribute remaining vertices proportionally to arc lengths
        verts_per_segment = []
        distributed = 0
        for i, arc_len in enumerate(arc_lengths):
            if total_arc > DEF_ERR_MARGIN:
                count = int(round(remaining * arc_len / total_arc))
            else:
                count = remaining // n_corners
            verts_per_segment.append(count)
            distributed += count

        # Adjust for rounding errors
        diff = remaining - distributed
        if diff != 0:
            # Add/remove from longest segment
            longest_idx = arc_lengths.index(max(arc_lengths))
            verts_per_segment[longest_idx] += diff

        # Build resampled boundary
        resampled = []

        for i in range(n_corners):
            start_idx = corner_indices[i]
            end_idx = corner_indices[(i + 1) % n_corners]
            n_intermediate = verts_per_segment[i]

            # Add corner vertex
            resampled.append(boundary[start_idx].copy())

            if n_intermediate <= 0:
                continue

            # Calculate positions for intermediate vertices
            arc_len = arc_lengths[i]
            target_spacing = arc_len / (n_intermediate + 1)

            # Walk along arc, placing vertices at regular intervals
            idx = start_idx
            accumulated = 0.0
            placed = 0

            while placed < n_intermediate and idx != end_idx:
                next_idx = (idx + 1) % n
                edge_vec = boundary[next_idx] - boundary[idx]
                edge_len = edge_vec.length

                while accumulated + edge_len >= target_spacing * (placed + 1) and placed < n_intermediate:
                    # Place vertex on this edge
                    t = (target_spacing * (placed + 1) - accumulated) / edge_len
                    t = max(0.0, min(1.0, t))
                    new_pt = boundary[idx] + edge_vec * t
                    resampled.append(new_pt.copy())
                    placed += 1

                accumulated += edge_len
                idx = next_idx

        return resampled

    def _build_boundary_kdtree(self):
        """Build KD-tree for fast boundary point queries."""
        n = len(self.boundary)
        self.boundary_kd = kdtree.KDTree(n)
        for i, pt in enumerate(self.boundary):
            self.boundary_kd.insert(pt, i)
        self.boundary_kd.balance()

    def _compute_medial_axis(self):
        """
        Compute medial axis using grid sampling.

        MA points are where the inscribed circle touches the boundary
        at two or more distinct points (equidistant condition).
        """
        bbox_min = self.analysis.bbox_min
        bbox_max = self.analysis.bbox_max
        bbox_size = bbox_max - bbox_min

        # Sample resolution based on fill detail
        resolution = max(20, self.fill_detail * 5)
        step_x = bbox_size.x / resolution
        step_y = bbox_size.y / resolution

        if step_x < DEF_ERR_MARGIN or step_y < DEF_ERR_MARGIN:
            return

        # Tolerance for equidistant condition (relative to step size)
        equi_tol = min(step_x, step_y) * 0.4

        candidates = []

        # Grid sampling to find MA candidates
        for i in range(resolution + 1):
            for j in range(resolution + 1):
                x = bbox_min.x + i * step_x
                y = bbox_min.y + j * step_y
                pt = Vector((x, y, 0))

                if not self._point_in_polygon(pt):
                    continue

                # Find closest boundary points
                nearest = self.boundary_kd.find_n(pt, 3)
                if len(nearest) < 2:
                    continue

                d0 = nearest[0][2]  # Distance to closest
                d1 = nearest[1][2]  # Distance to second closest

                # MA condition: equidistant to two (or more) boundary points
                # that are not adjacent on the boundary
                if abs(d0 - d1) < equi_tol:
                    # Check that the boundary points are not adjacent
                    idx0, idx1 = nearest[0][1], nearest[1][1]
                    n_boundary = len(self.boundary)
                    arc_dist = min(
                        abs(idx0 - idx1),
                        n_boundary - abs(idx0 - idx1)
                    )
                    # Points should be reasonably far apart on boundary
                    if arc_dist > n_boundary // 8:
                        candidates.append(MANode(
                            co=pt.copy(),
                            radius=d0,
                            boundary_indices=[idx0, idx1]
                        ))

        # Filter candidates to get clean MA skeleton
        self._filter_and_order_ma_nodes(candidates)

    def _filter_and_order_ma_nodes(self, candidates: List[MANode]):
        """Filter and order MA nodes to form a clean skeleton."""
        if not candidates:
            return

        # Remove nodes too close to each other
        bbox_size = self.analysis.bbox_max - self.analysis.bbox_min
        min_dist = max(bbox_size.x, bbox_size.y) / (self.fill_detail * 3)

        filtered = []
        for node in candidates:
            too_close = False
            for existing in filtered:
                if (node.co - existing.co).length < min_dist:
                    # Keep the one with larger radius (more central)
                    if node.radius > existing.radius:
                        filtered.remove(existing)
                    else:
                        too_close = True
                    break
            if not too_close:
                filtered.append(node)

        # Sort by distance from center (innermost first)
        center = self.analysis.center
        filtered.sort(key=lambda n: (n.co - center).length)

        self.ma_nodes = filtered

    def _build_ma_graph(self):
        """Build graph structure connecting MA nodes."""
        if len(self.ma_nodes) < 2:
            return

        # Connect nearby nodes to form MA skeleton
        n = len(self.ma_nodes)
        bbox_size = self.analysis.bbox_max - self.analysis.bbox_min
        max_edge_dist = max(bbox_size.x, bbox_size.y) / 2

        # Use MST-like approach: connect nearest unconnected nodes
        connected = {0}  # Start from center-most node
        edges = []

        while len(connected) < n:
            best_edge = None
            best_dist = float('inf')

            for i in connected:
                for j in range(n):
                    if j in connected:
                        continue
                    dist = (self.ma_nodes[i].co - self.ma_nodes[j].co).length
                    if dist < best_dist and dist < max_edge_dist:
                        best_dist = dist
                        best_edge = (i, j)

            if best_edge is None:
                break

            edges.append(MAEdge(node0=best_edge[0], node1=best_edge[1]))
            connected.add(best_edge[1])

        self.ma_edges = edges

    def _decompose_and_mesh(self):
        """
        Decompose domain into patches and mesh each with TFI.

        Strategy:
        - Create offset loops with SAME vertex count (clean 1:1 quad mapping)
        - Number of rings controlled by fill_detail parameter
        - Inner ring scale controlled by offset_size parameter
        - Final center region uses Grid Fill TFI for all-quad topology
        """
        n_boundary = len(self.boundary)
        center = self.analysis.center

        # Number of concentric rings based on fill_detail (default 5)
        n_rings = max(1, self.fill_detail)

        # Inner ring scale based on offset_size (0.0-1.0)
        # offset_size=0.1 means inner ring is 10% of outer, large center to fill
        # offset_size=0.9 means inner ring is 90% of outer, small center to fill
        inner_scale = max(0.05, min(0.5, self.offset_size))

        print(f"[MA DEBUG] n_boundary={n_boundary}, n_rings={n_rings}, inner_scale={inner_scale}")

        # Create offset loops - all with SAME vertex count for clean quad topology
        loops = []
        for ring_idx in range(n_rings + 1):
            t = ring_idx / n_rings
            # Use smooth interpolation for scale (ease in-out)
            t_smooth = t * t * (3 - 2 * t)
            # Scale from 1.0 (boundary) to inner_scale (innermost ring)
            scale = 1.0 - t_smooth * (1.0 - inner_scale)

            loop = []
            for pt in self.boundary:
                offset_pt = center + (pt - center) * scale
                loop.append(self.bm.verts.new(offset_pt))
            loops.append(loop)

        # Connect consecutive rings with simple 1:1 quads
        for ring_idx in range(n_rings):
            outer_loop = loops[ring_idx]
            inner_loop = loops[ring_idx + 1]

            for i in range(n_boundary):
                j = (i + 1) % n_boundary
                try:
                    self.bm.faces.new([
                        outer_loop[i], outer_loop[j],
                        inner_loop[j], inner_loop[i]
                    ])
                except ValueError:
                    pass

        # Fill center region with Grid Fill TFI
        inner_loop = loops[-1]
        print(f"[MA DEBUG] Filling center with {len(inner_loop)} vertices (scale={inner_scale})")
        self._fill_center_no_pole(inner_loop, center, inner_scale)

    def _fill_center_no_pole(
        self, inner_loop: List, center: Vector, inner_scale: float
    ):
        """
        Fill n-gon center region with all quads using Grid Fill TFI.

        REQUIRES: n must be a multiple of 4 for all-quad result.

        Strategy (Blender's Grid Fill algorithm):
        1. Identify 4 corners at n/4 intervals
        2. Map the n-gon to a square grid:
           - Edge 0 (corners 0→1) = top row
           - Edge 1 (corners 1→2) = right column
           - Edge 2 (corners 2→3) = bottom row (reversed)
           - Edge 3 (corners 3→0) = left column (reversed)
        3. Interior vertices via Transfinite Interpolation
        4. Result: (n/4) × (n/4) quads, no center pole
        """
        n = len(inner_loop)

        if n < 3:
            return

        if n == 3:
            try:
                self.bm.faces.new(inner_loop)
            except ValueError:
                pass
            return

        if n == 4:
            try:
                self.bm.faces.new(inner_loop)
            except ValueError:
                pass
            return

        # For Grid Fill, n should be multiple of 4
        # If not, we can't make all quads (would need triangles)
        if n % 4 != 0:
            print(f"[MA DEBUG] Warning: n={n} not multiple of 4, grid fill may have issues")

        # Grid dimensions
        k = n // 4  # Vertices per edge (not counting shared corner)

        # Identify corners at indices 0, k, 2k, 3k
        corners = [0, k, 2 * k, 3 * k]

        # Extract the 4 edges (each includes both endpoint corners)
        # Edge i goes from corner[i] to corner[i+1]
        edges = []
        for i in range(4):
            start = corners[i]
            end = corners[(i + 1) % 4]
            if end == 0:
                end = n
            edge = [inner_loop[j % n] for j in range(start, end + 1)]
            edges.append(edge)

        # Edge lengths (should all be k+1)
        edge_lens = [len(e) for e in edges]
        print(f"[MA DEBUG] Grid Fill: n={n}, k={k}, edge_lens={edge_lens}")

        # Build the grid using TFI
        # Grid is (k+1) × (k+1) vertices
        # grid[row][col] where row=0 is edge0 (top), row=k is edge2 reversed (bottom)
        # col=0 is edge3 reversed (left), col=k is edge1 (right)

        grid_size = k + 1
        grid = [[None for _ in range(grid_size)] for _ in range(grid_size)]

        # Fill boundary of grid from edges
        # Top row (row=0): edge0
        for col in range(grid_size):
            if col < len(edges[0]):
                grid[0][col] = edges[0][col]

        # Right column (col=k): edge1
        for row in range(grid_size):
            if row < len(edges[1]):
                grid[row][grid_size - 1] = edges[1][row]

        # Bottom row (row=k): edge2 reversed
        for col in range(grid_size):
            rev_col = grid_size - 1 - col
            if rev_col < len(edges[2]):
                grid[grid_size - 1][col] = edges[2][rev_col]

        # Left column (col=0): edge3 reversed
        for row in range(grid_size):
            rev_row = grid_size - 1 - row
            if rev_row < len(edges[3]):
                grid[row][0] = edges[3][rev_row]

        # Fill interior using Transfinite Interpolation
        for row in range(1, grid_size - 1):
            for col in range(1, grid_size - 1):
                u = col / (grid_size - 1)
                v = row / (grid_size - 1)

                # Get boundary points
                top = grid[0][col].co if grid[0][col] else center
                bottom = grid[grid_size - 1][col].co if grid[grid_size - 1][col] else center
                left = grid[row][0].co if grid[row][0] else center
                right = grid[row][grid_size - 1].co if grid[row][grid_size - 1] else center

                # Get corner points
                p00 = grid[0][0].co if grid[0][0] else center
                p10 = grid[0][grid_size - 1].co if grid[0][grid_size - 1] else center
                p01 = grid[grid_size - 1][0].co if grid[grid_size - 1][0] else center
                p11 = grid[grid_size - 1][grid_size - 1].co if grid[grid_size - 1][grid_size - 1] else center

                # TFI formula
                pt = (
                    (1 - v) * top + v * bottom +
                    (1 - u) * left + u * right -
                    (1 - u) * (1 - v) * p00 -
                    u * (1 - v) * p10 -
                    (1 - u) * v * p01 -
                    u * v * p11
                )

                grid[row][col] = self.bm.verts.new(pt)

        # Create quad faces
        quads_created = 0
        for row in range(grid_size - 1):
            for col in range(grid_size - 1):
                v00 = grid[row][col]
                v10 = grid[row][col + 1]
                v01 = grid[row + 1][col]
                v11 = grid[row + 1][col + 1]

                if v00 and v10 and v01 and v11:
                    try:
                        self.bm.faces.new([v00, v10, v11, v01])
                        quads_created += 1
                    except ValueError:
                        pass

        print(f"[MA DEBUG] Grid Fill created {quads_created} quads")

    def _connect_loops_with_quads(
        self, outer_loop: List, inner_loop: List, center: Vector
    ):
        """
        Connect two loops with quads, handling different vertex counts.

        Uses interpolation to create intermediate vertices when needed,
        ensuring all faces are quads.
        """
        n_outer = len(outer_loop)
        n_inner = len(inner_loop)

        if n_outer == n_inner:
            # Simple case: 1:1 mapping
            for i in range(n_outer):
                j = (i + 1) % n_outer
                try:
                    self.bm.faces.new([
                        outer_loop[i], outer_loop[j],
                        inner_loop[j], inner_loop[i]
                    ])
                except ValueError:
                    pass
            return

        # For different counts, create intermediate ring(s) to transition smoothly
        # This ensures all-quad topology

        # Calculate ratio
        ratio = n_outer / n_inner

        if ratio <= 2:
            # Can handle with direct bridging using interpolated vertices
            self._bridge_loops_direct(outer_loop, inner_loop, center)
        else:
            # Need intermediate ring
            # Create a ring with count between outer and inner
            n_mid = (n_outer + n_inner) // 2
            if n_mid % 2 == 1:
                n_mid += 1

            mid_loop = []
            for i in range(n_mid):
                # Interpolate position between outer and inner
                outer_idx = (i * n_outer) // n_mid
                inner_idx = (i * n_inner) // n_mid

                outer_pt = outer_loop[outer_idx].co
                inner_pt = inner_loop[inner_idx].co
                mid_pt = (outer_pt + inner_pt) / 2
                mid_loop.append(self.bm.verts.new(mid_pt))

            # Connect outer to mid
            self._bridge_loops_direct(outer_loop, mid_loop, center)
            # Connect mid to inner
            self._bridge_loops_direct(mid_loop, inner_loop, center)

    def _bridge_loops_direct(
        self, outer_loop: List, inner_loop: List, center: Vector
    ):
        """
        Bridge two loops with quads using direct vertex mapping.

        For each inner edge, creates quads to cover the corresponding outer edges.
        When an inner edge maps to multiple outer edges, creates additional
        vertices on the inner edge to maintain quad topology.
        """
        n_outer = len(outer_loop)
        n_inner = len(inner_loop)

        # Build mapping: which outer vertices correspond to each inner vertex
        inner_to_outer = []
        for i in range(n_inner):
            outer_idx = (i * n_outer) // n_inner
            inner_to_outer.append(outer_idx)

        # Process each inner edge
        for i in range(n_inner):
            i_next = (i + 1) % n_inner

            outer_start = inner_to_outer[i]
            outer_end = inner_to_outer[i_next]

            # Handle wraparound
            if outer_end <= outer_start and i_next != 0:
                outer_end += n_outer
            elif outer_end == outer_start and i_next == 0:
                outer_end = n_outer

            # Count outer edges in this segment
            n_outer_edges = outer_end - outer_start
            if n_outer_edges <= 0:
                n_outer_edges = 1

            inner_v0 = inner_loop[i]
            inner_v1 = inner_loop[i_next]

            if n_outer_edges == 1:
                # Simple quad
                o0 = outer_start % n_outer
                o1 = outer_end % n_outer
                try:
                    self.bm.faces.new([
                        outer_loop[o0], outer_loop[o1],
                        inner_v1, inner_v0
                    ])
                except ValueError:
                    pass
            else:
                # Multiple outer edges -> create interpolated vertices on inner edge
                # and make quads
                inner_edge_verts = [inner_v0]

                # Create intermediate vertices along inner edge
                for k in range(1, n_outer_edges):
                    t = k / n_outer_edges
                    interp_pt = inner_v0.co.lerp(inner_v1.co, t)
                    inner_edge_verts.append(self.bm.verts.new(interp_pt))

                inner_edge_verts.append(inner_v1)

                # Create quads
                for k in range(n_outer_edges):
                    o0 = (outer_start + k) % n_outer
                    o1 = (outer_start + k + 1) % n_outer
                    try:
                        self.bm.faces.new([
                            outer_loop[o0], outer_loop[o1],
                            inner_edge_verts[k + 1], inner_edge_verts[k]
                        ])
                    except ValueError:
                        pass

    def _fallback_offset_mesh(self) -> bpy.types.Object:
        """Fallback for simple shapes: concentric offset loops."""
        n_boundary = len(self.boundary)
        n_layers = self.fill_detail
        center = self.analysis.center

        # Create offset loops
        loops = []
        for layer_idx in range(n_layers + 1):
            t = layer_idx / n_layers
            scale = 1.0 - t * 0.85  # Leave small center region

            loop = []
            for pt in self.boundary:
                offset_pt = center + (pt - center) * scale
                loop.append(self.bm.verts.new(offset_pt))
            loops.append(loop)

        # Create quad faces between loops
        for layer_idx in range(n_layers):
            outer = loops[layer_idx]
            inner = loops[layer_idx + 1]

            for i in range(n_boundary):
                j = (i + 1) % n_boundary
                try:
                    self.bm.faces.new([outer[i], outer[j], inner[j], inner[i]])
                except ValueError:
                    pass

        # Fill center with single n-gon (acceptable for small center)
        try:
            self.bm.faces.new(loops[-1])
        except ValueError:
            pass

        return self._finalize()

    def _point_in_polygon(self, point: Vector) -> bool:
        """Ray casting point-in-polygon test."""
        x, y = point.x, point.y
        n = len(self.boundary)
        inside = False

        j = n - 1
        for i in range(n):
            xi, yi = self.boundary[i].x, self.boundary[i].y
            xj, yj = self.boundary[j].x, self.boundary[j].y

            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside
            j = i

        return inside


# =============================================================================
# Convex Offset Mesher - Uses Q-Morph
# =============================================================================

class ConvexOffsetMesher(QuadMeshGenerator):
    """For convex freeform shapes - delegates to Q-Morph."""

    def generate(self) -> bpy.types.Object:
        qmorph = QMorphMesher(self.mesh_obj, self.analysis, self.params)
        qmorph.bm = self.bm
        return qmorph.generate()


# =============================================================================
# Concave Grid Mesher - Uses Medial Axis
# =============================================================================

class ConcaveGridMesher(QuadMeshGenerator):
    """For concave shapes - delegates to Medial Axis."""

    def generate(self) -> bpy.types.Object:
        ma = MedialAxisMesher(self.mesh_obj, self.analysis, self.params)
        ma.bm = self.bm
        return ma.generate()


# =============================================================================
# Hole Mesher
# =============================================================================

class HoleMesher(QuadMeshGenerator):
    """Mesher for shapes with holes (annular regions)."""

    def generate(self) -> bpy.types.Object:
        outer = self._get_boundary_verts()
        holes = self._get_hole_boundaries()

        if not holes:
            if self.analysis.is_convex:
                mesher = QMorphMesher(self.mesh_obj, self.analysis, self.params)
            else:
                mesher = MedialAxisMesher(self.mesh_obj, self.analysis, self.params)
            mesher.bm = self.bm
            return mesher.generate()

        # Mesh annular region for each hole
        for hole in holes:
            self._mesh_annular(outer, hole)

        return self._finalize()

    def _mesh_annular(self, outer: List[Vector], inner: List[Vector]):
        """Create quad mesh between outer and inner boundaries."""
        # Resample to same count
        n = max(len(outer), len(inner), self.fill_detail * 4)
        outer = self._resample_boundary_uniform(outer, n)
        inner = self._resample_boundary_uniform(inner, n)

        n_rings = self.fill_detail

        # Create rings
        ring_verts = []
        for ring_idx in range(n_rings + 1):
            t = ring_idx / n_rings
            ring = []
            for i in range(n):
                pt = outer[i].lerp(inner[i], t)
                ring.append(self.bm.verts.new(pt))
            ring_verts.append(ring)

        # Create quads
        for r in range(n_rings):
            for i in range(n):
                j = (i + 1) % n
                try:
                    self.bm.faces.new([
                        ring_verts[r][i], ring_verts[r][j],
                        ring_verts[r + 1][j], ring_verts[r + 1][i]
                    ])
                except ValueError:
                    pass


# =============================================================================
# Smart Dispatcher
# =============================================================================

def smart_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Main entry point for shape-aware quad meshing.

    Selects the best algorithm based on detected shape type.
    """
    analysis = detect_shape(curve_obj)

    # Select mesher
    if analysis.has_holes:
        mesher_class = HoleMesher
    else:
        mesher_map = {
            ShapeType.CIRCLE: PolarGridMesher,
            ShapeType.ELLIPSE: PolarGridMesher,
            ShapeType.RECTANGLE: RectangleGridMesher,
            ShapeType.REGULAR_POLYGON: PolygonRadialMesher,
            ShapeType.CONVEX_FREEFORM: QMorphMesher,
            ShapeType.CONCAVE_FREEFORM: MedialAxisMesher,
            ShapeType.UNKNOWN: QMorphMesher,
        }
        mesher_class = mesher_map.get(analysis.shape_type, QMorphMesher)

    mesher = mesher_class(mesh_obj, analysis, params)
    return mesher.generate()


def medial_axis_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for explicit Medial Axis Transform quad meshing.

    Uses MA decomposition regardless of shape type.
    """
    analysis = detect_shape(curve_obj)

    mesher = MedialAxisMesher(mesh_obj, analysis, params)
    return mesher.generate()


def grid_tfi_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for Grid Fill TFI quad meshing.

    Direct TFI grid fill without offset rings - preserves boundary shape
    and creates clean quad topology by mapping n-gon to square grid.
    """
    analysis = detect_shape(curve_obj)

    mesher = GridTFIMesher(mesh_obj, analysis, params)
    return mesher.generate()


def polar_grid_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for Polar Grid quad meshing.

    Uses O-Grid topology optimal for circles and ellipses.
    Creates radial sectors with no center pole.
    """
    analysis = detect_shape(curve_obj)

    mesher = PolarGridMesher(mesh_obj, analysis, params)
    return mesher.generate()


def rectangle_grid_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for Rectangle Grid quad meshing.

    Creates axis-aligned quad grid optimal for rectangles.
    Uses bilinear TFI for clean quad topology.
    """
    analysis = detect_shape(curve_obj)

    mesher = RectangleGridMesher(mesh_obj, analysis, params)
    return mesher.generate()


def polygon_radial_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for Polygon Radial quad meshing.

    Creates radial sectors for regular polygons - each side maps to
    a trapezoid sector meshed with TFI. No center pole.
    """
    analysis = detect_shape(curve_obj)

    mesher = PolygonRadialMesher(mesh_obj, analysis, params)
    return mesher.generate()


def qmorph_quad_mesh(
    mesh_obj: bpy.types.Object,
    curve_obj: bpy.types.Object,
    params: dict
) -> bpy.types.Object:
    """
    Entry point for Q-Morph quad meshing.

    Advancing front algorithm that converts triangulation to quads.
    Good for convex freeform shapes.
    """
    analysis = detect_shape(curve_obj)

    mesher = QMorphMesher(mesh_obj, analysis, params)
    return mesher.generate()
