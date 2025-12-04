# bezier_utils/utils/view_utils.py

import bpy
import gpu
from gpu_extras.batch import batch_for_shader
from mathutils import Vector, geometry
from math import log, pi
from bpy_extras.view3d_utils import (
    region_2d_to_vector_3d,
    region_2d_to_location_3d,
    location_3d_to_region_2d,
    region_2d_to_origin_3d,
)
from ..core.props import FTProps
from ..constants import LARGE_NO, MAX_NONSEL_CURVE_RES
from .bezier_math import isStraightSeg, getPtsAlongBezier2D, getLinesFromPts
from .object_utils import isBezier
from .curve_utils import moveSplineStart, createClipElem, getSVGPathElem
from .math_utils import toHexStr
from xml.dom import minidom
import random


def getGridSubdiv(space3d):
    return space3d.overlay.grid_subdivisions


def getUnit():
    return bpy.context.scene.unit_settings.length_unit


def getUnitSystem():
    return bpy.context.scene.unit_settings.system


def getUnitScale():
    fact = 3.28084 if (getUnitSystem() == "IMPERIAL") else 1
    return fact * bpy.context.scene.unit_settings.scale_length


def get3dLoc(region, rv3d, xy, vec=None):
    if vec is None:
        vec = region_2d_to_vector_3d(region, rv3d, xy)
    return region_2d_to_location_3d(region, rv3d, xy, vec)


# TODO: Rework after grid subdiv is enabled (in a version later than 2.8)
def getViewDistRounding(space3d, rv3d):
    viewDist = rv3d.view_distance * getUnitScale()
    gridDiv = getGridSubdiv(space3d)
    subFact = 1
    # TODO: Separate logic for 1
    if gridDiv == 1:
        gridDiv = 10
    elif gridDiv == 2:
        subFact = 5
    elif viewDist < 0.5:
        subFact = 2
    return int(log(viewDist, gridDiv)) - subFact


# Return axis-indices (x:0, y:1, z:2) of plane with closest orientation
# to view
def getClosestPlaneToView(rv3d):
    viewmat = rv3d.view_matrix
    trans, quat, scale = viewmat.decompose()
    tm = quat.to_matrix().to_4x4()
    viewnormal = tm.inverted() @ Vector((0, 0, 1))

    xynormal = [Vector((0, 0, 1)), [0, 1]]
    yznormal = [Vector((1, 0, 0)), [1, 2]]
    xznormal = [Vector((0, 1, 0)), [0, 2]]
    normals = [xynormal, yznormal, xznormal]
    minAngle = pi
    minIdx = None
    for idx, normal in enumerate(normals):
        rotDiff = viewnormal.rotation_difference(normal[0]).angle
        if rotDiff > pi / 2:
            rotDiff = pi - rotDiff
        if rotDiff < minAngle:
            minIdx = idx
            minAngle = rotDiff
    return normals[minIdx][1]


def getCoordFromLoc(region, rv3d, loc):
    coord = location_3d_to_region_2d(region, rv3d, loc)
    # return a unlocatable pt if None to avoid errors
    return coord if (coord is not None) else Vector((9000, 9000))


# To be called only from 3d view
def getCurrAreaRegion(context):
    a, r = [
        (a, r)
        for a in bpy.context.screen.areas
        if a.type == "VIEW_3D"
        for r in a.regions
        if (r == context.region)
    ][0]
    return a, r


def isOutside(context, event, exclInRgns=True):
    x = event.mouse_region_x
    y = event.mouse_region_y
    region = context.region

    if x < 0 or x > region.width or y < 0 or y > region.height:
        return True

    elif not exclInRgns:
        return False

    area, r = getCurrAreaRegion(context)

    for r in area.regions:
        if r == region:
            continue
        xR = r.x - region.x
        yR = r.y - region.y
        if x >= xR and y >= yR and x <= (xR + r.width) and y <= (yR + r.height):
            return True

    return False


def getPtProjOnPlane(region, rv3d, xy, p1, p2, p3, p4=None):
    vec = region_2d_to_vector_3d(region, rv3d, xy)
    orig = region_2d_to_origin_3d(region, rv3d, xy)

    pt = geometry.intersect_ray_tri(p1, p2, p3, vec, orig, False)  # p4 != None)
    # ~ if(not pt and p4):
    # ~ pt = geometry.intersect_ray_tri(p2, p4, p3, vec, orig, True)
    return pt


# find the location on 3d line p1-p2 if xy is already on 2d projection (in rv3d) of p1-p2
def getPtProjOnLine(region, rv3d, xy, p1, p2):
    # Just find a non-linear point (TODO: simpler way)
    pd1 = p2 - p1
    pd2 = Vector(sorted(pd1, key=lambda x: abs(x), reverse=True))
    maxIdx0 = [i for i in range(3) if abs(pd1[i]) == abs(pd2[0])][0]
    maxIdx1 = [i for i in range(3) if abs(pd1[i]) == abs(pd2[1])][0]
    pd = Vector()
    pd[maxIdx0] = -pd2[1]
    pd[maxIdx1] = pd2[0]
    p3 = p2 + pd

    # Raycast from 2d point onto the plane
    return getPtProjOnPlane(region, rv3d, xy[:2], p1, p2, p3)


def getLineTransMatrices(pt0, pt1):
    diffV = (pt1 - pt0).normalized()

    # Calculate rotation from global X axis to custom axis
    # This gives a proper orthogonal coordinate frame by rotating the entire global frame
    global_x = Vector((1, 0, 0))
    quat = global_x.rotation_difference(diffV)

    # Convert quaternion to 4x4 matrix
    # invTm transforms from custom space to global space
    invTm = quat.to_matrix().to_4x4()
    tm = invTm.inverted_safe()

    return tm, invTm


def getWindowRegionIdx(area, regionIdx):  # For finding quad view index
    idx = 0
    for j, r in enumerate(area.regions):
        if j == regionIdx:
            return idx
        if r.type == "WINDOW":
            idx += 1
    return None


def getAreaRegionIdxs(xy, exclInRgns=True):
    x, y = xy
    areas = [a for a in bpy.context.screen.areas]
    idxs = None
    for i, a in enumerate(areas):
        if a.type != "VIEW_3D":
            continue
        regions = [r for r in a.regions]
        for j, r in enumerate(regions):
            if x > r.x and x < r.x + r.width and y > r.y and y < r.y + r.height:
                if r.type == "WINDOW":
                    if not exclInRgns:
                        return [i, j]
                    idxs = [i, j]
                elif exclInRgns:
                    return None
    return idxs


def getAllAreaRegions():
    info = []
    areas = []
    i = 0

    areas = [a for a in bpy.context.screen.areas if (a.type == "VIEW_3D")]

    # bpy.context.screen doesn't work in case of Add-on Config window
    while len(areas) == 0 and i < len(bpy.data.screens):
        areas = [a for a in bpy.data.screens[i].areas if (a.type == "VIEW_3D")]
        i += 1

    for a in areas:
        regions = [r for r in a.regions if r.type == "WINDOW"]
        if len(a.spaces[0].region_quadviews) > 0:
            for i, r in enumerate(a.spaces[0].region_quadviews):
                info.append([a, regions[i], r])
        else:
            r = a.spaces[0].region_3d
            info.append([a, regions[0], r])
    return info


def getResetBatch(shader, btype):  # "LINES" or "POINTS"
    from gpu.types import GPUBatch

    return GPUBatch(type=btype, init=dict(pos=[], color=[]))
    # return batch_for_shader(shader, btype, {"pos": [], "color": []})


# From python template
def getFaceUnderMouse(obj, region, rv3d, xy, maxFaceCnt):
    if obj is None or obj.type != "MESH" or len(obj.data.polygons) > maxFaceCnt:
        return None, None, None
    viewVect = region_2d_to_vector_3d(region, rv3d, xy)
    rayOrig = region_2d_to_origin_3d(region, rv3d, xy)
    mw = obj.matrix_world
    invMw = mw.inverted_safe()

    rayTarget = rayOrig + viewVect
    rayOrigObj = invMw @ rayOrig
    rayTargetObj = invMw @ rayTarget
    rayDirObj = rayTargetObj - rayOrigObj

    success, location, normal, faceIdx = obj.ray_cast(rayOrigObj, rayDirObj)
    if success:
        return mw @ location, normal, faceIdx
    else:
        return None, None, None


def get2dBBox(obj, region, rv3d, precise=False):
    mw = obj.matrix_world
    if precise:
        co2ds = [
            getCoordFromLoc(region, rv3d, mw @ Vector(v.co)) for v in obj.data.vertices
        ]
    else:
        co2ds = [getCoordFromLoc(region, rv3d, mw @ Vector(b)) for b in obj.bound_box]
    minX = min(c[0] for c in co2ds)
    maxX = max(c[0] for c in co2ds)
    minY = min(c[1] for c in co2ds)
    maxY = max(c[1] for c in co2ds)

    return minX, minY, maxX, maxY


def isPtIn2dBBox(obj, region, rv3d, xy, extendBy=0, precise=False):
    minX, minY, maxX, maxY = get2dBBox(obj, region, rv3d, precise)
    if (
        xy[0] > (minX - extendBy)
        and xy[0] < (maxX + extendBy)
        and xy[1] > (minY - extendBy)
        and xy[1] < (maxY + extendBy)
    ):
        return True
    else:
        return False


def getSnappableObjs(region, rv3d, xy):
    objs = bpy.context.selected_objects
    if bpy.context.object is not None:
        objs.append(bpy.context.object)
    return [
        o
        for o in objs
        if (
            o.type == "MESH"
            and len(o.modifiers) == 0
            and isPtIn2dBBox(o, region, rv3d, xy)
        )
    ]


def getClosestEdgeLoc2d(obj, region, rv3d, xy, faceIdx=None):
    mw = obj.matrix_world
    minDist = LARGE_NO
    closestLoc = None
    pt = Vector(xy).to_3d()
    edgeWSCos = None
    edgeIdx = None
    closestIntersect = None
    vertPairs = (
        obj.data.polygons[faceIdx].edge_keys
        if faceIdx is not None
        else [e.vertices for e in obj.data.edges]
    )
    for i, vertPair in enumerate(vertPairs):
        co0 = mw @ obj.data.vertices[vertPair[0]].co
        co1 = mw @ obj.data.vertices[vertPair[1]].co

        pt0 = getCoordFromLoc(region, rv3d, co0).to_3d()
        pt1 = getCoordFromLoc(region, rv3d, co1).to_3d()
        intersect, percDist = geometry.intersect_point_line(pt, pt0, pt1)
        if percDist < 0:
            intersect = pt0
            percDist = 0
        elif percDist > 1:
            intersect = pt1
            percDist = 1
        dist = (intersect - pt).length
        if dist < minDist:
            minDist = dist
            closestIntersect = intersect
            edgeWSCos = [co0, co1]
            edgeIdx = i
    if edgeWSCos is not None:
        closestLoc = getPtProjOnLine(
            region, rv3d, closestIntersect, edgeWSCos[0], edgeWSCos[1]
        )

    return edgeIdx, edgeWSCos, closestLoc, minDist


def getSelFaceLoc(region, rv3d, xy, maxFaceCnt, objs=None, checkEdge=False):
    # Import locally to avoid circular dependency if FTProps is needed
    # But FTProps is in core/props.py. view_utils is in utils/.
    # If we need FTProps.snapDist, we might need to pass it as arg or import.
    # For now, let's assume we can import it or use a default.
    # Actually, the original code uses FTProps.snapDist.
    # I will import FTProps inside the function if needed, but FTProps is in core.
    # view_utils is in utils. utils should not depend on core.
    # I will use a constant or pass it as argument.
    # The original code: if(closestLoc != None and minDist < FTProps.snapDist):
    # I will use a hardcoded default or import constant.
    # FTProps.snapDist default is 20.
    SNAP_DIST = 20  # Fallback

    if objs is None:
        objs = getSnappableObjs(region, rv3d, xy)
    if len(objs) > maxFaceCnt:
        return None, None, None, None
    for obj in objs:
        loc, normal, faceIdx = getFaceUnderMouse(obj, region, rv3d, xy, maxFaceCnt)
        if loc is not None:
            if checkEdge:
                edgeIdx, edgeWSCos, closestLoc, minDist = getClosestEdgeLoc2d(
                    obj, region, rv3d, xy, faceIdx
                )
                # Try to get snapDist from context if possible, else default
                try:
                    snapDist = bpy.context.window_manager.bezierToolkitParams.snapDist
                except:
                    snapDist = SNAP_DIST

                if closestLoc is not None and minDist < snapDist:
                    return obj, closestLoc, normal, faceIdx
            return obj, loc, normal, faceIdx
    return None, None, None, None


class RegionMouseXYInfo:
    @staticmethod
    def getRegionMouseXYInfo(event, exclInRgns):
        xyScreen = [event.mouse_x, event.mouse_y]
        idxs = getAreaRegionIdxs(xyScreen, exclInRgns)
        if idxs is None:
            return None
        else:
            i, j = idxs
            area = bpy.context.screen.areas[i]
            region = area.regions[j]
            space3d = area.spaces[0]
            if len(space3d.region_quadviews) > 0:
                qIdx = getWindowRegionIdx(area, j)
                rv3d = space3d.region_quadviews[qIdx]
            else:
                rv3d = space3d.region_3d
            xy = [xyScreen[0] - region.x, xyScreen[1] - region.y]
            return RegionMouseXYInfo(area, space3d, region, rv3d, xy, xyScreen)

    def __init__(self, area, space3d, region, rv3d, xy, xyScreen):
        self.area = area
        self.space3d = space3d
        self.region = region
        self.rv3d = rv3d
        self.xy = xy
        self.xyScreen = xyScreen

    def __eq__(self, other):
        if other is None:
            return False
        return (
            self.area == other.area
            and self.region == other.region
            and self.rv3d == other.rv3d
        )


def getLineShades(lineCos, baseColor, start, end, mid=True):
    if len(lineCos) == 0:
        return [], []
    if len(lineCos) == 1:
        return lineCos[0], [baseColor]
    if mid:
        midPt = lineCos[0] + (lineCos[1] - lineCos[0]) / 2
    col1 = [start * c for c in baseColor]
    col2 = [end * c for c in baseColor]
    if mid:
        return [lineCos[0], midPt, midPt, lineCos[1]], [col1, col2, col2, col1]
    else:
        return [lineCos[0], lineCos[1]], [col1, col2]


# Display Info Classes
class SegDisplayInfo:
    def __init__(self, segPts, segColor):
        self.segPts = segPts
        self.segColor = segColor


def getLineShades(lineCos, baseColor, start, end, mid=True):
    if len(lineCos) == 0:
        return [], []
    if len(lineCos) == 1:
        return lineCos[0], [baseColor]
    if mid:
        midPt = lineCos[0] + (lineCos[1] - lineCos[0]) / 2
    col1 = [start * c for c in baseColor]
    col2 = [end * c for c in baseColor]
    if mid:
        return [lineCos[0], midPt, midPt, lineCos[1]], [col1, col2, col2, col1]
    else:
        return [lineCos[0], lineCos[1]], [col1, col2]


class BGLDrawInfo:
    def __init__(self, size, color, pts):
        self.size = size
        self.color = color
        self.pts = pts


class BGLDrawInfoLine(BGLDrawInfo):
    def __init__(
        self, size, color, pts, gradientStart=None, gradientEnd=None, mid=True
    ):
        super(BGLDrawInfoLine, self).__init__(size, color, pts)
        self.gradientStart = gradientStart
        self.gradientEnd = gradientEnd
        self.mid = mid


class BGLDrawMgr:
    def __init__(self, shader):
        self.lineInfoMap = {}
        self.ptInfoMap = {}
        self.shader = shader

    def addLineInfo(
        self, infoId, size, color, pts, gradientStart=None, gradientEnd=None, mid=True
    ):
        self.lineInfoMap[infoId] = BGLDrawInfoLine(
            size, color, pts, gradientStart, gradientEnd, mid
        )

    def addPtInfo(self, infoId, size, color, pts):
        self.ptInfoMap[infoId] = BGLDrawInfo(size, color, pts)

    def redraw(self):
        lineInfos = sorted(self.lineInfoMap.values(), key=lambda x: (x.size))
        pos = []
        col = []
        batches = []
        for i, info in enumerate(lineInfos):
            if i == 0 or info.size != lineInfos[i - 1].size:
                if i > 0:
                    # bgl.glLineWidth(lineInfos[i-1].size)
                    gpu.state.line_width_set(lineInfos[i - 1].size)
                    batch = batch_for_shader(
                        self.shader, "LINES", {"pos": pos, "color": col}
                    )
                    batch.draw(self.shader)
                pos = []
                col = []

            if info.gradientEnd is not None and info.gradientStart is not None:
                if len(info.pts) != 2 and len(info.color) != 1:
                    raise ValueError(
                        "Exactly two " + "coordinates and one color for gradient line"
                    )
                linePos, lineCols = getLineShades(
                    info.pts,
                    info.color[0],
                    info.gradientStart,
                    info.gradientEnd,
                    info.mid,
                )
            else:
                if len(info.pts) == 0:
                    continue
                linePos = info.pts[:]
                lineCols = info.color[:]
                diff = len(linePos) - len(lineCols)
                if diff >= 0:
                    for j in range(diff):
                        lineCols.append(info.color[-1])
                else:
                    for j in range(-diff):
                        lineCols.pop()
            pos += linePos
            col += lineCols

        if len(pos) > 0:
            # bgl.glLineWidth(lineInfos[-1].size)
            gpu.state.line_width_set(lineInfos[-1].size)
            batch = batch_for_shader(self.shader, "LINES", {"pos": pos, "color": col})
            batch.draw(self.shader)

        ptInfos = sorted(self.ptInfoMap.values(), key=lambda x: (x.size))
        for i, info in enumerate(ptInfos):
            if i == 0 or info.size != ptInfos[i - 1].size:
                if i > 0:
                    # bgl.glPointSize(ptInfos[i-1].size)
                    gpu.state.point_size_set(ptInfos[i - 1].size)
                    batch = batch_for_shader(
                        self.shader, "POINTS", {"pos": pos, "color": col}
                    )
                    batch.draw(self.shader)
                pos = []
                col = []

            if len(info.pts) == 0:
                continue
            ptCols = info.color[:]
            diff = len(info.pts) - len(info.color)
            if diff >= 0:
                for j in range(diff):
                    ptCols.append(info.color[-1])
            else:
                for j in range(-diff):
                    ptCols.pop()

            pos += info.pts[:]
            col += ptCols

        if len(pos) > 0:
            # bgl.glPointSize(ptInfos[-1].size)
            gpu.state.point_size_set(ptInfos[-1].size)
            batch = batch_for_shader(self.shader, "POINTS", {"pos": pos, "color": col})
            batch.draw(self.shader)

    def resetLineInfo(self, infoId):
        drawInfo = self.lineInfoMap.get(infoId)
        if drawInfo is not None:
            drawInfo.pts = []

    def resetPtInfo(self, infoId):
        drawInfo = self.ptInfoMap.get(infoId)
        if drawInfo is not None:
            drawInfo.pts = []

    def reset(self):
        for key in list(self.ptInfoMap.keys()):
            self.ptInfoMap[key].pts = []
        for key in list(self.lineInfoMap.keys()):
            self.lineInfoMap[key].pts = []


# Return line batch for bezier line segments and handles and point batch for handle tips
def updateBezierBatches(
    bglDrawMgr, segDispInfos, bptDispInfos, areaRegionInfo, defHdlType="ALIGNED"
):
    lineCos = []  # segment is also made up of lines
    lineColors = []
    for i, info in enumerate(segDispInfos):
        segPts = info.segPts
        if isStraightSeg(segPts):
            lineCos += [segPts[0][1], segPts[1][1]]
            lineColors += [info.segColor, info.segColor]
        else:
            pts = getPtsAlongBezier2D(
                segPts,
                areaRegionInfo,
                FTProps.dispCurveRes,
                getCoordFromLoc,
                maxRes=MAX_NONSEL_CURVE_RES,
            )
            segLineCos = getLinesFromPts(pts)
            lineCos += segLineCos
            lineColors += [info.segColor for j in range(0, len(segLineCos))]

    tipCos = []
    tipColors = []
    for i, info in enumerate(bptDispInfos):
        pt = info.pt
        for hn in info.handleNos:
            lineCos += [pt[hn], pt[hn + 1]]

            if len(pt) < 5:
                htype = defHdlType  # For Draw
            else:
                htype = pt[3 + hn]

            lineColors += [
                FTProps.hdlColMap[htype],
                FTProps.hdlColMap[htype],
            ]

        # Re-arrange tips so handles are on top of Bezier point
        tc = info.tipColors
        tc = [tc[1], tc[0], tc[2]]
        pt = [pt[1], pt[0], pt[2]]
        for j, tipColor in enumerate(tc):
            if tipColor is not None:
                tipCos.append(pt[j])
                tipColors.append(tipColor)

    bglDrawMgr.addLineInfo("bezLineBatch", FTProps.lineWidth, lineColors, lineCos)
    bglDrawMgr.addPtInfo("bezTipBatch", FTProps.drawPtSize, tipColors, tipCos)


def resetToolbarTool():
    win = bpy.context.window
    scr = win.screen
    areas3d = [area for area in scr.areas if area.type == "VIEW_3D"]
    override = {"window": win, "screen": scr, "scene": bpy.context.scene}
    for a in areas3d:
        override["area"] = a
        regions = [region for region in a.regions if region.type == "WINDOW"]
        for r in regions:
            override["region"] = r
            with bpy.context.temp_override(**override):
                bpy.ops.wm.tool_set_by_index()


def updateMetaBtns(caller, event, keymap=None):
    if keymap is None:
        keymap = {
            "LEFT_SHIFT": "shift",
            "RIGHT_SHIFT": "shift",
            "LEFT_CTRL": "ctrl",
            "RIGHT_CTRL": "ctrl",
            "LEFT_ALT": "alt",
            "RIGHT_ALT": "alt",
        }

    var = keymap.get(event.type)

    if var is not None:
        expr = "caller." + var + " = "
        if event.value == "PRESS":
            exec(expr + "True")
        if event.value == "RELEASE":
            exec(expr + "False")
        return True

    return False


unitMap = {"FEET": "'", "METERS": "m"}


def showSnapToPlane(params):
    return (
        params.snapOrient not in {"VIEW", "REFERENCE", "CURR_POS"}
        and hasattr(params, "constrAxes")
        and params.constrAxes.startswith("shift")
    )


class MarkerController:
    drawHandlerRef = None
    defPointSize = 6
    ptColor = (0, 0.8, 0.8, 1)

    def createSMMap(self, context):
        objs = context.selected_objects
        smMap = {}
        for curve in objs:
            if not isBezier(curve):
                continue

            smMap[curve.name] = {}
            mw = curve.matrix_world
            for splineIdx, spline in enumerate(curve.data.splines):
                if not spline.use_cyclic_u:
                    continue

                # initialize to the curr start vert co and idx
                smMap[curve.name][splineIdx] = [
                    mw @ curve.data.splines[splineIdx].bezier_points[0].co,
                    0,
                ]

                for pt in spline.bezier_points:
                    pt.select_control_point = False

            if len(smMap[curve.name]) == 0:
                del smMap[curve.name]

        return smMap

    def createBatch(self, context):
        positions = [s[0] for cn in self.smMap.values() for s in cn.values()]
        colors = [MarkerController.ptColor for i in range(0, len(positions))]

        self.batch = batch_for_shader(
            self.shader, "POINTS", {"pos": positions, "color": colors}
        )

        if context.area:
            context.area.tag_redraw()

    def drawHandler(self):
        # bgl.glPointSize(MarkerController.defPointSize)
        gpu.state.point_size_set(MarkerController.defPointSize)
        self.batch.draw(self.shader)

    def removeMarkers(self, context):
        if MarkerController.drawHandlerRef is not None:
            bpy.types.SpaceView3D.draw_handler_remove(
                MarkerController.drawHandlerRef, "WINDOW"
            )

            if context.area and hasattr(context.space_data, "region_3d"):
                context.area.tag_redraw()

            MarkerController.drawHandlerRef = None
        self.deselectAll()

    def __init__(self, context):
        self.smMap = self.createSMMap(context)
        self.shader = gpu.shader.from_builtin("FLAT_COLOR")
        # ~ self.shader.bind()

        try:
            markerSize = (
                context.preferences.addons.get("bezier_utils").preferences.markerSize
                if context.preferences.addons.get("bezier_utils")
                else 7
            )
            MarkerController.defPointSize = markerSize
        except Exception as e:
            # ~ print("BezierUtils: Fetching marker size", e)
            MarkerController.defPointSize = 6

        MarkerController.drawHandlerRef = bpy.types.SpaceView3D.draw_handler_add(
            self.drawHandler, (), "WINDOW", "POST_VIEW"
        )

        self.createBatch(context)

    def saveStartVerts(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            spMap = self.smMap[curveName]

            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                if markerInfo[1] != 0:
                    loc, idx = markerInfo[0], markerInfo[1]
                    moveSplineStart(curve, splineIdx, idx)

    def updateSMMap(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            spMap = self.smMap[curveName]
            mw = curve.matrix_world

            for splineIdx in spMap.keys():
                markerInfo = spMap[splineIdx]
                loc, idx = markerInfo[0], markerInfo[1]
                pts = curve.data.splines[splineIdx].bezier_points

                selIdxs = [x for x in range(0, len(pts)) if pts[x].select_control_point]

                selIdx = selIdxs[0] if (len(selIdxs) > 0) else idx
                co = mw @ pts[selIdx].co
                self.smMap[curveName][splineIdx] = [co, selIdx]

    def deselectAll(self):
        for curveName in self.smMap.keys():
            curve = bpy.data.objects[curveName]
            for spline in curve.data.splines:
                for pt in spline.bezier_points:
                    pt.select_control_point = False

    def getSpaces3D(context):
        areas3d = [
            area for area in context.window.screen.areas if area.type == "VIEW_3D"
        ]

        return [s for a in areas3d for s in a.spaces if s.type == "VIEW_3D"]

    def hideHandles(context):
        states = []
        spaces = MarkerController.getSpaces3D(context)
        for s in spaces:
            if hasattr(s.overlay, "show_curve_handles"):
                states.append(s.overlay.show_curve_handles)
                s.overlay.show_curve_handles = False
            elif hasattr(s.overlay, "display_handle"):  # 2.90
                states.append(s.overlay.display_handle)
                s.overlay.display_handle = "NONE"
        return states

    def resetShowHandleState(context, handleStates):
        spaces = MarkerController.getSpaces3D(context)
        for i, s in enumerate(spaces):
            if hasattr(s.overlay, "show_curve_handles"):
                s.overlay.show_curve_handles = handleStates[i]
            elif hasattr(s.overlay, "display_handle"):  # 2.90
                s.overlay.display_handle = handleStates[i]


class FTHotKeyData:
    pass


def world_to_camera_view(scene, obj, coord):
    co_local = obj.matrix_world.normalized().inverted() @ coord
    z = -co_local.z

    camera = scene.camera.data
    frame = [-v for v in camera.view_frame(scene=scene)[:3]]
    if camera.type != "ORTHO":
        if z == 0.0:
            return Vector([0.5, 0.5, 0.0])
        frame = [v / (v.z / z) for v in frame]

    min_x, max_x = frame[1].x, frame[2].x
    min_y, max_y = frame[0].y, frame[1].y

    x = (co_local.x - min_x) / (max_x - min_x)
    y = (co_local.y - min_y) / (max_y - min_y)

    return Vector([x, y, z])


# Round to logarithmic scale .1, 0, 10, 100 etc.
# (47.538, -1) -> 47.5; (47.538, 0) -> 48.0; (47.538, 1) -> 50.0; (47.538, 2) -> 0,
# TODO: Rework after grid subdiv is enabled (in a version later than 2.8)
def roundedVect(space3d, vect, rounding, axes):
    rounding += 1
    subdiv = getGridSubdiv(space3d)
    # TODO: Separate logic for 1
    if subdiv == 1:
        subdiv = 10
    fact = ((subdiv**rounding) / subdiv) / getUnitScale()
    retVect = vect.copy()
    # ~ Vector([round(vect[i] / fact) * fact for i in axes])
    for i in axes:
        retVect[i] = round(vect[i] / fact) * fact
    return retVect


def getSVGPt(co, docW, docH, camera=None, region=None, rv3d=None):
    if camera is not None:
        scene = bpy.context.scene
        xy = world_to_camera_view(scene, camera, co)
        return complex(xy[0] * docW, docH - (xy[1] * docH))
    elif region is not None and rv3d is not None:
        xy = getCoordFromLoc(region, rv3d, co)
        return complex(xy[0], docH - xy[1])


def exportSVG(
    context,
    filepath,
    exportView,
    clipView,
    lineWidth,
    lineColorOpts,
    lineColor,
    fillColorOpts,
    fillColor,
):
    svgXML = '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
    clipElemId = "BBoxClipElem"

    if lineColorOpts == "PICK":
        lineCol, lineAlpha = toHexStr(lineColor)

    if fillColorOpts == "PICK":
        fillCol, fillAlpha = toHexStr(fillColor)

    if exportView == "ACTIVE_VIEW":
        area = context.area
        if area.type != "VIEW_3D":
            area = [a for a in bpy.context.screen.areas if a.type == "VIEW_3D"][0]
        region = [r for r in area.regions if r.type == "WINDOW"][0]
        space3d = area.spaces[0]
        if len(space3d.region_quadviews) > 0:
            rv3d = space3d.region_quadviews[3]
        else:
            rv3d = space3d.region_3d
        camera = None
        docW = region.width
        docH = region.height
    else:
        region = None
        rv3d = None
        camera = bpy.data.objects[exportView]
        docW = bpy.context.scene.render.resolution_x
        docH = bpy.context.scene.render.resolution_y

    doc = minidom.parseString(svgXML)
    svgElem = doc.documentElement

    svgElem.setAttribute("width", str(docW))
    svgElem.setAttribute("height", str(docH))

    if clipView:
        createClipElem(doc, svgElem, docW, docH, clipElemId)

    idx = 0
    for o in bpy.context.scene.objects:
        mw = o.matrix_world
        if isBezier(o) and o.visible_get():
            path = []
            filledPath = []
            for spline in o.data.splines:
                part = []
                bpts = spline.bezier_points
                for i in range(1, len(bpts)):
                    prevBezierPt = bpts[i - 1]
                    pt = bpts[i]
                    seg = [
                        prevBezierPt.co,
                        prevBezierPt.handle_right,
                        pt.handle_left,
                        pt.co,
                    ]
                    part.append(
                        [
                            getSVGPt(mw @ co, docW, docH, camera, region, rv3d)
                            for co in seg
                        ]
                    )

                if spline.use_cyclic_u:
                    seg = [
                        bpts[-1].co,
                        bpts[-1].handle_right,
                        bpts[0].handle_left,
                        bpts[0].co,
                    ]
                    part.append(
                        [
                            getSVGPt(mw @ co, docW, docH, camera, region, rv3d)
                            for co in seg
                        ]
                    )

                if len(part) > 0:
                    if (
                        spline.use_cyclic_u
                        and o.data.dimensions == "2D"
                        and o.data.fill_mode != "NONE"
                    ):
                        filledPath.append(part)
                    else:
                        path.append(part)

            for p in [path, filledPath]:
                if len(p) == 0:
                    continue

                if lineColorOpts == "RANDOM":
                    lineColor = [random.random() for i in range(3)] + [1]
                    lineCol, lineAlpha = toHexStr(lineColor)

                if p == path:
                    fc, fa = None, None
                elif fillColorOpts == "RANDOM":
                    fillColor = [random.random() for i in range(3)] + [1]
                    fc, fa = toHexStr(fillColor)
                else:
                    fc, fa = fillCol, fillAlpha

                svgPathElem = getSVGPathElem(
                    doc,
                    docW,
                    docH,
                    p,
                    idx,
                    lineWidth,
                    lineCol,
                    lineAlpha,
                    fc,
                    fa,
                    clipView,
                    clipElemId,
                )
                if svgPathElem is not None:
                    svgElem.appendChild(svgPathElem)
                    idx += 1

    doc.writexml(open(filepath, "w"))
