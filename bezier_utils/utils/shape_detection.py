# bezier_utils/utils/shape_detection.py
"""
Shape detection and classification for Bezier curves.

Analyzes curve geometry to determine the best meshing algorithm:
- Circle/Ellipse: Polar grid
- Rectangle: Axis-aligned grid
- Regular Polygon: Radial sectors
- Convex Freeform: Q-Morph or offset
- Concave Freeform: Medial axis decomposition
"""

from enum import Enum, auto
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from math import pi, degrees, acos, sqrt
from mathutils import Vector

from ..constants import (
    SHAPE_CIRCULARITY_CIRCLE_THRESH,
    SHAPE_CIRCULARITY_ELLIPSE_THRESH,
    SHAPE_ANGLE_TOLERANCE,
    SHAPE_SIDE_VARIANCE_THRESH,
    SHAPE_ASPECT_RATIO_CIRCLE_MIN,
    SHAPE_ASPECT_RATIO_CIRCLE_MAX,
    SHAPE_DETECTION_SAMPLE_RES,
    SHAPE_MIN_CORNER_ANGLE,
    SHAPE_MAX_CORNER_ANGLE,
    DEF_ERR_MARGIN,
)


class ShapeType(Enum):
    """Classification of detected shape types."""
    CIRCLE = auto()
    ELLIPSE = auto()
    RECTANGLE = auto()
    REGULAR_POLYGON = auto()
    CONVEX_FREEFORM = auto()
    CONCAVE_FREEFORM = auto()
    UNKNOWN = auto()


@dataclass
class ShapeAnalysis:
    """Result of shape detection analysis."""
    shape_type: ShapeType
    confidence: float  # 0.0 to 1.0
    center: Vector
    bbox_min: Vector
    bbox_max: Vector

    # Closure
    is_closed: bool
    segment_count: int

    # Circularity metrics
    circularity: float  # 1.0 = perfect circle
    area: float
    perimeter: float

    # Corner/polygon metrics
    corner_count: int
    corner_angles: List[float] = field(default_factory=list)
    corner_indices: List[int] = field(default_factory=list)
    straight_segment_indices: List[int] = field(default_factory=list)

    # Side metrics (for polygons)
    side_lengths: List[float] = field(default_factory=list)
    side_length_variance: float = 0.0

    # Ellipse metrics
    major_axis: float = 0.0
    minor_axis: float = 0.0
    aspect_ratio: float = 1.0
    ellipse_center: Optional[Vector] = None
    ellipse_fit_error: float = 0.0

    # Convexity
    is_convex: bool = True

    # Holes (inner boundaries)
    holes: List[List[Vector]] = field(default_factory=list)
    has_holes: bool = False

    # Sampled boundary points (for meshing)
    boundary_points: List[Vector] = field(default_factory=list)


class ShapeDetector:
    """
    Analyzes Bezier curves to detect their geometric shape type.

    Uses a multi-metric approach:
    1. Segment analysis (straight vs curved)
    2. Corner detection via angle analysis
    3. Circularity measurement
    4. Ellipse fitting
    5. Convexity check
    """

    def __init__(self, curve_obj, spline_idx: int = 0):
        """
        Initialize detector with a Blender curve object.

        Args:
            curve_obj: Blender curve object (must be Bezier type)
            spline_idx: Index of spline to analyze (default: first)
        """
        self.curve = curve_obj
        self.spline_idx = spline_idx
        self.spline = curve_obj.data.splines[spline_idx]
        self.bpts = self.spline.bezier_points
        self.mw = curve_obj.matrix_world
        self._sampled_pts = None

    def analyze(self) -> ShapeAnalysis:
        """
        Perform full shape analysis.

        Returns:
            ShapeAnalysis with detected shape type and metrics
        """
        # Get sampled points for analysis
        sampled_pts = self._get_sampled_points(SHAPE_DETECTION_SAMPLE_RES)

        # Basic metrics
        is_closed = self.spline.use_cyclic_u
        segment_count = len(self.bpts) if is_closed else len(self.bpts) - 1

        # Bounding box
        bbox_min, bbox_max = self._compute_bounding_box(sampled_pts)

        # Center (centroid)
        center = self._compute_center(sampled_pts)

        # Area and perimeter
        area = self._compute_area(sampled_pts)
        perimeter = self._compute_perimeter(sampled_pts)

        # Circularity
        circularity = self._compute_circularity(area, perimeter)

        # Straight segments
        straight_indices = self._detect_straight_segments()

        # Corner analysis
        corner_angles, corner_indices = self._compute_corner_angles()

        # Side lengths (for polygon detection)
        side_lengths = self._compute_side_lengths(corner_indices, sampled_pts)
        side_variance = self._compute_side_variance(side_lengths)

        # Ellipse fitting
        major, minor, fit_error, ellipse_center = self._fit_ellipse(sampled_pts)
        aspect_ratio = major / minor if minor > DEF_ERR_MARGIN else 1.0

        # Convexity
        is_convex = self._check_convexity(sampled_pts)

        # Detect holes (other splines)
        holes, has_holes = self._detect_holes()

        # Build analysis result
        analysis = ShapeAnalysis(
            shape_type=ShapeType.UNKNOWN,  # Will be set by classification
            confidence=0.0,
            center=center,
            bbox_min=bbox_min,
            bbox_max=bbox_max,
            is_closed=is_closed,
            segment_count=segment_count,
            circularity=circularity,
            area=area,
            perimeter=perimeter,
            corner_count=len(corner_indices),
            corner_angles=corner_angles,
            corner_indices=corner_indices,
            straight_segment_indices=straight_indices,
            side_lengths=side_lengths,
            side_length_variance=side_variance,
            major_axis=major,
            minor_axis=minor,
            aspect_ratio=aspect_ratio,
            ellipse_center=ellipse_center,
            ellipse_fit_error=fit_error,
            is_convex=is_convex,
            holes=holes,
            has_holes=has_holes,
            boundary_points=sampled_pts,
        )

        # Classify shape
        shape_type, confidence = self._classify_shape(analysis)
        analysis.shape_type = shape_type
        analysis.confidence = confidence

        return analysis

    def _get_sampled_points(self, resolution: int) -> List[Vector]:
        """Sample points uniformly along the curve."""
        if self._sampled_pts is not None:
            return self._sampled_pts

        from mathutils.geometry import interpolate_bezier

        pts = []
        n_pts = len(self.bpts)
        is_closed = self.spline.use_cyclic_u
        n_segs = n_pts if is_closed else n_pts - 1

        # Samples per segment
        samples_per_seg = max(2, resolution // n_segs)

        for seg_idx in range(n_segs):
            i0 = seg_idx
            i1 = (seg_idx + 1) % n_pts

            p0 = self.mw @ self.bpts[i0].co
            p1 = self.mw @ self.bpts[i0].handle_right
            p2 = self.mw @ self.bpts[i1].handle_left
            p3 = self.mw @ self.bpts[i1].co

            seg_pts = interpolate_bezier(p0, p1, p2, p3, samples_per_seg)
            # Don't include last point (will be first of next segment)
            pts.extend(seg_pts[:-1])

        self._sampled_pts = pts
        return pts

    def _compute_bounding_box(
        self, pts: List[Vector]
    ) -> Tuple[Vector, Vector]:
        """Calculate axis-aligned bounding box."""
        if not pts:
            return Vector((0, 0, 0)), Vector((0, 0, 0))

        min_x = min(p.x for p in pts)
        max_x = max(p.x for p in pts)
        min_y = min(p.y for p in pts)
        max_y = max(p.y for p in pts)
        min_z = min(p.z for p in pts)
        max_z = max(p.z for p in pts)

        return Vector((min_x, min_y, min_z)), Vector((max_x, max_y, max_z))

    def _compute_center(self, pts: List[Vector]) -> Vector:
        """Compute geometric centroid."""
        if not pts:
            return Vector((0, 0, 0))
        return sum(pts, Vector((0, 0, 0))) / len(pts)

    def _compute_area(self, pts: List[Vector]) -> float:
        """Compute area using shoelace formula (2D, XY plane)."""
        if len(pts) < 3:
            return 0.0

        area = 0.0
        n = len(pts)
        for i in range(n):
            j = (i + 1) % n
            area += pts[i].x * pts[j].y
            area -= pts[j].x * pts[i].y

        return abs(area) / 2.0

    def _compute_perimeter(self, pts: List[Vector]) -> float:
        """Compute perimeter (total edge length)."""
        if len(pts) < 2:
            return 0.0

        perimeter = sum(
            (pts[i] - pts[(i + 1) % len(pts)]).length
            for i in range(len(pts))
        )
        return perimeter

    def _compute_circularity(self, area: float, perimeter: float) -> float:
        """
        Compute circularity metric.

        circularity = 4 * pi * area / perimeter^2
        Perfect circle = 1.0
        """
        if perimeter < DEF_ERR_MARGIN:
            return 0.0
        return (4.0 * pi * area) / (perimeter * perimeter)

    def _detect_straight_segments(self) -> List[int]:
        """Identify segments that are effectively straight lines."""
        straight_indices = []
        n_pts = len(self.bpts)
        is_closed = self.spline.use_cyclic_u
        n_segs = n_pts if is_closed else n_pts - 1

        for seg_idx in range(n_segs):
            i0 = seg_idx
            i1 = (seg_idx + 1) % n_pts

            bp0 = self.bpts[i0]
            bp1 = self.bpts[i1]

            # Check handle types
            if bp0.handle_right_type == 'VECTOR' and bp1.handle_left_type == 'VECTOR':
                straight_indices.append(seg_idx)
                continue

            # Check if handles are collinear with segment
            p0 = self.mw @ bp0.co
            h0 = self.mw @ bp0.handle_right
            h1 = self.mw @ bp1.handle_left
            p1 = self.mw @ bp1.co

            seg_dir = (p1 - p0).normalized()
            if seg_dir.length < DEF_ERR_MARGIN:
                continue

            # Check handle alignment
            h0_dir = (h0 - p0).normalized() if (h0 - p0).length > DEF_ERR_MARGIN else seg_dir
            h1_dir = (p1 - h1).normalized() if (p1 - h1).length > DEF_ERR_MARGIN else seg_dir

            if abs(seg_dir.dot(h0_dir)) > 0.999 and abs(seg_dir.dot(h1_dir)) > 0.999:
                straight_indices.append(seg_idx)

        return straight_indices

    def _compute_corner_angles(self) -> Tuple[List[float], List[int]]:
        """
        Calculate internal angles at each control point.

        Returns:
            Tuple of (all_angles, corner_indices)
            corner_indices contains indices where angle indicates a corner
        """
        angles = []
        corner_indices = []
        n_pts = len(self.bpts)
        is_closed = self.spline.use_cyclic_u

        for i in range(n_pts):
            if not is_closed and (i == 0 or i == n_pts - 1):
                angles.append(180.0)  # Endpoints get 180 (no corner)
                continue

            prev_idx = (i - 1) % n_pts
            next_idx = (i + 1) % n_pts

            p_prev = self.mw @ self.bpts[prev_idx].co
            p_curr = self.mw @ self.bpts[i].co
            p_next = self.mw @ self.bpts[next_idx].co

            v1 = (p_prev - p_curr)
            v2 = (p_next - p_curr)

            if v1.length < DEF_ERR_MARGIN or v2.length < DEF_ERR_MARGIN:
                angles.append(180.0)
                continue

            v1 = v1.normalized()
            v2 = v2.normalized()

            # Angle between vectors
            dot = max(-1.0, min(1.0, v1.dot(v2)))
            angle = degrees(acos(dot))
            angles.append(angle)

            # Check if this is a corner
            if SHAPE_MIN_CORNER_ANGLE < angle < SHAPE_MAX_CORNER_ANGLE:
                corner_indices.append(i)

        return angles, corner_indices

    def _compute_side_lengths(
        self, corner_indices: List[int], sampled_pts: List[Vector]
    ) -> List[float]:
        """Compute lengths between consecutive corners."""
        if len(corner_indices) < 2:
            return []

        is_closed = self.spline.use_cyclic_u

        side_lengths = []
        n_corners = len(corner_indices)

        for i in range(n_corners):
            if not is_closed and i == n_corners - 1:
                break

            i0 = corner_indices[i]
            i1 = corner_indices[(i + 1) % n_corners]

            p0 = self.mw @ self.bpts[i0].co
            p1 = self.mw @ self.bpts[i1].co

            # Use direct distance for now (could use arc length for curved sides)
            side_lengths.append((p1 - p0).length)

        return side_lengths

    def _compute_side_variance(self, side_lengths: List[float]) -> float:
        """Compute coefficient of variation for side lengths."""
        if len(side_lengths) < 2:
            return 0.0

        mean = sum(side_lengths) / len(side_lengths)
        if mean < DEF_ERR_MARGIN:
            return 0.0

        variance = sum((s - mean) ** 2 for s in side_lengths) / len(side_lengths)
        std_dev = sqrt(variance)

        return std_dev / mean  # Coefficient of variation

    def _fit_ellipse(
        self, pts: List[Vector]
    ) -> Tuple[float, float, float, Vector]:
        """
        Fit an ellipse to the sampled points.

        Uses algebraic fitting method.

        Returns:
            (major_axis, minor_axis, fit_error, center)
        """
        if len(pts) < 5:
            bbox_min, bbox_max = self._compute_bounding_box(pts)
            size = bbox_max - bbox_min
            return max(size.x, size.y) / 2, min(size.x, size.y) / 2, 1.0, self._compute_center(pts)

        # Simple ellipse fitting: use bounding box as approximation
        # For more accurate fitting, implement least-squares ellipse fit
        center = self._compute_center(pts)
        bbox_min, bbox_max = self._compute_bounding_box(pts)

        # Half-widths
        half_w = (bbox_max.x - bbox_min.x) / 2
        half_h = (bbox_max.y - bbox_min.y) / 2

        major = max(half_w, half_h)
        minor = min(half_w, half_h)

        # Compute fit error as average distance from ellipse
        fit_error = 0.0
        if major > DEF_ERR_MARGIN and minor > DEF_ERR_MARGIN:
            for pt in pts:
                # Normalized coordinates
                dx = (pt.x - center.x) / major
                dy = (pt.y - center.y) / minor
                # Distance from unit circle
                r = sqrt(dx * dx + dy * dy)
                fit_error += abs(r - 1.0)
            fit_error /= len(pts)

        return major, minor, fit_error, center

    def _check_convexity(self, pts: List[Vector]) -> bool:
        """Check if the shape is convex using cross product method."""
        if len(pts) < 3:
            return True

        n = len(pts)
        sign = None

        for i in range(n):
            p0 = pts[i]
            p1 = pts[(i + 1) % n]
            p2 = pts[(i + 2) % n]

            # 2D cross product (z component of 3D cross)
            v1 = p1 - p0
            v2 = p2 - p1
            cross_z = v1.x * v2.y - v1.y * v2.x

            if abs(cross_z) < DEF_ERR_MARGIN:
                continue  # Collinear points

            current_sign = cross_z > 0

            if sign is None:
                sign = current_sign
            elif sign != current_sign:
                return False  # Sign changed, not convex

        return True

    def _detect_holes(self) -> Tuple[List[List[Vector]], bool]:
        """
        Detect holes (inner boundaries) from other splines.

        A hole is a spline that is fully contained within the main spline.
        """
        holes = []
        curve_data = self.curve.data

        if len(curve_data.splines) <= 1:
            return holes, False

        # Get main boundary points
        main_pts = self._get_sampled_points(SHAPE_DETECTION_SAMPLE_RES)
        if not main_pts:
            return holes, False

        # Check each other spline
        for spline_idx, spline in enumerate(curve_data.splines):
            if spline_idx == self.spline_idx:
                continue

            if spline.type != 'BEZIER':
                continue

            if not spline.use_cyclic_u:
                continue  # Holes must be closed

            # Sample this spline
            hole_pts = self._sample_spline(spline)
            if not hole_pts:
                continue

            # Check if this spline is inside the main boundary
            # Use point-in-polygon test for first point
            if self._point_in_polygon(hole_pts[0], main_pts):
                holes.append(hole_pts)

        return holes, len(holes) > 0

    def _sample_spline(self, spline) -> List[Vector]:
        """Sample points from a spline."""
        from mathutils.geometry import interpolate_bezier

        pts = []
        bpts = spline.bezier_points
        n_pts = len(bpts)
        is_closed = spline.use_cyclic_u
        n_segs = n_pts if is_closed else n_pts - 1
        samples_per_seg = max(2, SHAPE_DETECTION_SAMPLE_RES // max(1, n_segs))

        for seg_idx in range(n_segs):
            i0 = seg_idx
            i1 = (seg_idx + 1) % n_pts

            p0 = self.mw @ bpts[i0].co
            p1 = self.mw @ bpts[i0].handle_right
            p2 = self.mw @ bpts[i1].handle_left
            p3 = self.mw @ bpts[i1].co

            seg_pts = interpolate_bezier(p0, p1, p2, p3, samples_per_seg)
            pts.extend(seg_pts[:-1])

        return pts

    def _point_in_polygon(self, point: Vector, polygon: List[Vector]) -> bool:
        """Ray casting point-in-polygon test."""
        x, y = point.x, point.y
        n = len(polygon)
        inside = False

        j = n - 1
        for i in range(n):
            xi, yi = polygon[i].x, polygon[i].y
            xj, yj = polygon[j].x, polygon[j].y

            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside
            j = i

        return inside

    def _classify_shape(self, analysis: ShapeAnalysis) -> Tuple[ShapeType, float]:
        """
        Apply classification rules to determine shape type.

        Decision tree:
        1. If has holes -> use shape type but mark has_holes
        2. If all segments straight and 4 corners at ~90 deg -> RECTANGLE
        3. If all segments straight and N corners at equal angles -> REGULAR_POLYGON
        4. If circularity > 0.95 and aspect_ratio ~1 -> CIRCLE
        5. If circularity > 0.85 -> ELLIPSE
        6. If convex -> CONVEX_FREEFORM
        7. Otherwise -> CONCAVE_FREEFORM
        """
        if not analysis.is_closed:
            return ShapeType.UNKNOWN, 0.0

        all_straight = (
            len(analysis.straight_segment_indices) == analysis.segment_count
        )
        n_corners = analysis.corner_count
        angles = analysis.corner_angles

        # Check for rectangle (4 corners at ~90 degrees, all straight)
        if all_straight and n_corners == 4:
            # Get angles at corners
            corner_angles = [angles[i] for i in analysis.corner_indices]
            angle_deviation = sum(abs(a - 90.0) for a in corner_angles) / 4.0
            if angle_deviation < SHAPE_ANGLE_TOLERANCE:
                confidence = 1.0 - (angle_deviation / 90.0)
                return ShapeType.RECTANGLE, confidence

        # Check for regular polygon (N corners at equal angles, all straight)
        if all_straight and n_corners >= 3:
            expected_angle = (n_corners - 2) * 180.0 / n_corners
            corner_angles = [angles[i] for i in analysis.corner_indices]
            angle_deviation = sum(abs(a - expected_angle) for a in corner_angles) / n_corners

            if (angle_deviation < SHAPE_ANGLE_TOLERANCE and
                    analysis.side_length_variance < SHAPE_SIDE_VARIANCE_THRESH):
                confidence = 1.0 - max(angle_deviation / 90.0, analysis.side_length_variance)
                return ShapeType.REGULAR_POLYGON, confidence

        # Check for circle (high circularity, aspect ratio ~1)
        if analysis.circularity > SHAPE_CIRCULARITY_CIRCLE_THRESH:
            if (SHAPE_ASPECT_RATIO_CIRCLE_MIN < analysis.aspect_ratio <
                    SHAPE_ASPECT_RATIO_CIRCLE_MAX):
                return ShapeType.CIRCLE, analysis.circularity

        # Check for ellipse (high circularity, non-unit aspect ratio)
        if analysis.circularity > SHAPE_CIRCULARITY_ELLIPSE_THRESH:
            return ShapeType.ELLIPSE, analysis.circularity

        # Check convexity for freeform classification
        if analysis.is_convex:
            return ShapeType.CONVEX_FREEFORM, 0.7

        return ShapeType.CONCAVE_FREEFORM, 0.5


def detect_shape(curve_obj, spline_idx: int = 0) -> ShapeAnalysis:
    """
    Convenience function to detect shape of a curve.

    Args:
        curve_obj: Blender curve object
        spline_idx: Index of spline to analyze

    Returns:
        ShapeAnalysis with detected shape and metrics
    """
    detector = ShapeDetector(curve_obj, spline_idx)
    return detector.analyze()
