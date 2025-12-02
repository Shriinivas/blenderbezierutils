# bezier_utils/utils/curve_utils.py

import bpy
import bmesh
from math import ceil
from mathutils import Vector, Matrix
from xml.dom import minidom
from ..constants import DEF_ERR_MARGIN
from .math_utils import vectCmpWithMargin, floatCmpWithMargin, getBBox, toHexStr
from .bezier_math import getInterpolatedVertsCo
from .bezier_math import (
    getTForPt,
    getPartialSeg,
    getInterpBezierPts,
    getSegLen,
    getIntersectPts,
)
from .object_utils import isBezier, safeRemoveObj, copyObjAttr
from .view_utils import getCoordFromLoc, world_to_camera_view


def copyBezierPt(src, target, freeHandles=None, srcMw=Matrix(), invDestMW=Matrix()):
    target.handle_left_type = "FREE"
    target.handle_right_type = "FREE"

    target.co = invDestMW @ (srcMw @ src.co)
    target.handle_left = invDestMW @ (srcMw @ src.handle_left)
    target.handle_right = invDestMW @ (srcMw @ src.handle_right)

    if freeHandles is None or not freeHandles[0]:
        target.handle_left_type = src.handle_left_type
    if freeHandles is None or not freeHandles[1]:
        target.handle_right_type = src.handle_right_type


def createSplineForSeg(curveData, bezierPts):
    spline = curveData.splines.new("BEZIER")
    spline.bezier_points.add(len(bezierPts) - 1)
    spline.use_cyclic_u = False

    for i, pt in enumerate(bezierPts):
        if i == 0:
            freeHandles = [False, True]
        elif i == len(bezierPts) - 1:
            freeHandles = [True, False]
        else:
            freeHandles = None
        copyBezierPt(pt, spline.bezier_points[i], freeHandles=freeHandles)

    return spline


def createSpline(curveData, srcSpline, excludePtIdxs={}):
    spline = curveData.splines.new("BEZIER")
    spline.bezier_points.add(len(srcSpline.bezier_points) - len(excludePtIdxs) - 1)
    spline.use_cyclic_u = srcSpline.use_cyclic_u

    ptIdx = 0
    for i in range(0, len(srcSpline.bezier_points)):
        if i not in excludePtIdxs:
            copyBezierPt(
                srcSpline.bezier_points[i],
                spline.bezier_points[ptIdx],
                freeHandles=None,
            )
            ptIdx += 1

    return spline


def createSkeletalCurve(obj, collections):
    objCopy = obj.copy()
    objCopy.name = obj.name
    dataCopy = obj.data.copy()
    dataCopy.splines.clear()
    objCopy.data = dataCopy
    # Duplicate the action (animation data)
    if obj.animation_data and obj.animation_data.action:
        objCopy.animation_data_clear()  # Clear any existing animation data links
        objCopy.animation_data_create()  # Create a new animation data block
        action_copy = obj.animation_data.action.action.copy()  # Copy the action
        objCopy.animation_data.action = action_copy  # Assign the copied action

    for coll in collections:
        coll.objects.link(objCopy)

    return objCopy


def createObjFromPts(
    curvePts, dimensions="3D", collection=None, closed=False, calcHdlTypes=True
):
    data = bpy.data.curves.new("BezierCurve", "CURVE")
    data.dimensions = dimensions
    obj = bpy.data.objects.new("BezierCurve", data)
    # ~ collection = context.collection
    if collection is None:
        collection = bpy.context.scene.collection
    collection.objects.link(obj)
    # ~ obj.location = context.scene.cursor.location

    # ~ depsgraph = context.evaluated_depsgraph_get()
    # ~ depsgraph.update()

    # ~ invM = obj.matrix_world.inverted_safe()

    spline = data.splines.new("BEZIER")
    spline.use_cyclic_u = False

    if vectCmpWithMargin(curvePts[0][1], curvePts[-1][1]):
        curvePts[0][0] = curvePts[-1][0]
        spline.use_cyclic_u = True
        curvePts.pop()

    if closed:
        spline.use_cyclic_u = True

    spline.bezier_points.add(len(curvePts) - 1)
    prevPt = None
    for i, pt in enumerate(curvePts):
        currPt = spline.bezier_points[i]
        currPt.co = pt[1]
        currPt.handle_right = pt[2]
        if not calcHdlTypes and len(pt) > 3:
            currPt.handle_right_type = pt[3]
            currPt.handle_left_type = pt[4]
        elif (
            prevPt is not None
            and prevPt.handle_right == prevPt.co
            and pt[0] == pt[1]
            and currPt.co != prevPt.co
        ):  # straight line
            if prevPt.handle_left_type != "VECTOR":
                prevPt.handle_left_type = "FREE"
            prevPt.handle_right_type = "VECTOR"
            currPt.handle_right_type = "FREE"
            currPt.handle_left_type = "VECTOR"
        else:
            currPt.handle_left_type = "FREE"
            currPt.handle_right_type = "FREE"
            currPt.handle_left = pt[0]
            ldiffV = pt[1] - pt[0]
            rdiffV = pt[2] - pt[1]
            if vectCmpWithMargin(ldiffV, rdiffV) and not floatCmpWithMargin(
                ldiffV.length, 0
            ):
                currPt.handle_left_type = "ALIGNED"
                currPt.handle_right_type = "ALIGNED"
        prevPt = currPt

    bpts = spline.bezier_points
    if (
        spline.use_cyclic_u
        and vectCmpWithMargin(bpts[-1].handle_right, bpts[-1].co)
        and vectCmpWithMargin(bpts[0].handle_left, bpts[0].co)
    ):
        if bpts[-1].handle_left_type != "VECTOR":
            bpts[-1].handle_left_type = "FREE"
        bpts[-1].handle_right_type = "VECTOR"
        if bpts[0].handle_right_type != "VECTOR":
            bpts[0].handle_right_type = "FREE"
        bpts[0].handle_left_type = "VECTOR"

    return obj


def removeShapeKeys(obj):
    if obj.data.shape_keys is None:
        return

    keyblocks = reversed(obj.data.shape_keys.key_blocks)
    for sk in keyblocks:
        obj.shape_key_remove(sk)


def getShapeKeyInfo(obj):
    keyData = []
    keyNames = []

    if obj.data.shape_keys is not None:
        keyblocks = obj.data.shape_keys.key_blocks
        for key in keyblocks:
            keyData.append(
                [
                    [d.handle_left.copy(), d.co.copy(), d.handle_right.copy()]
                    for d in key.data
                ]
            )
            keyNames.append(key.name)

    return keyNames, keyData


def updateShapeKeyData(obj, keyData, keyNames, startIdx, cnt=None, add=False):
    if obj.data.shape_keys is None and not add:
        return

    currIdx = obj.active_shape_key_index
    if not add:
        removeShapeKeys(obj)
    if cnt is None:
        cnt = len(keyData[0])

    for i, name in enumerate(keyNames):
        key = obj.shape_key_add(name=name)
        for j in range(0, cnt):
            keyIdx = j + startIdx
            key.data[j].handle_left = keyData[i][keyIdx][0].copy()
            key.data[j].co = keyData[i][keyIdx][1].copy()
            key.data[j].handle_right = keyData[i][keyIdx][2].copy()

    obj.active_shape_key_index = currIdx


def getLastSegIdx(obj, splineIdx):
    spline = obj.data.splines[splineIdx]
    ptCnt = len(spline.bezier_points)
    return ptCnt - 1 if (spline.use_cyclic_u) else ptCnt - 2


def addLastSeg(spline):
    if spline.use_cyclic_u:
        lt0 = spline.bezier_points[0].handle_left_type
        rt0 = spline.bezier_points[0].handle_right_type
        pt = spline.bezier_points[0]
        pt.handle_left_type = "FREE"
        pt.handle_right_type = "FREE"
        spline.use_cyclic_u = False
        spline.bezier_points.add(1)
        copyObjAttr(spline.bezier_points[0], spline.bezier_points[-1])
        spline.bezier_points[0].handle_left_type = lt0
        spline.bezier_points[-1].handle_right_type = rt0


def moveSplineStart(obj, splineIdx, idx):
    pts = obj.data.splines[splineIdx].bezier_points
    cnt = len(pts)

    ptCopy = [
        [
            p.co.copy(),
            p.handle_right.copy(),
            p.handle_left.copy(),
            p.handle_right_type,
            p.handle_left_type,
        ]
        for p in pts
    ]

    for i, pt in enumerate(pts):
        srcIdx = (idx + i) % cnt
        p = ptCopy[srcIdx]

        pt.handle_left_type = "FREE"
        pt.handle_right_type = "FREE"
        pt.co = p[0]
        pt.handle_right = p[1]
        pt.handle_left = p[2]
        pt.handle_right_type = p[3]
        pt.handle_left_type = p[4]


def joinCurves(curves):
    obj = curves[0]
    invMW = obj.matrix_world.inverted_safe()
    for curve in curves[1:]:
        mw = curve.matrix_world
        for spline in curve.data.splines:
            newSpline = obj.data.splines.new("BEZIER")
            copyObjAttr(spline, newSpline)
            newSpline.bezier_points.add(len(spline.bezier_points) - 1)
            for i, pt in enumerate(spline.bezier_points):
                copyObjAttr(pt, newSpline.bezier_points[i], invDestMW=invMW, mw=mw)
        safeRemoveObj(curve)
    return obj


# Insert spline at location insertIdx, duplicated from existing spline at
# location srcSplineIdx and remove points with indices in removePtIdxs from new spline
def insertSpline(obj, srcSplineIdx, insertIdx, removePtIdxs):
    srcSpline = obj.data.splines[srcSplineIdx]
    # Appended at end
    createSpline(obj.data, srcSpline, removePtIdxs)
    splineCnt = len(obj.data.splines)
    nextIdx = insertIdx
    for idx in range(nextIdx, splineCnt - 1):
        srcSpline = obj.data.splines[nextIdx]
        createSpline(obj.data, srcSpline)
        obj.data.splines.remove(srcSpline)


def removeBezierPts(obj, splineIdx, removePtIdxs):
    oldSpline = obj.data.splines[splineIdx]
    bpts = oldSpline.bezier_points
    if min(removePtIdxs) >= len(bpts):
        return

    if len(set(range(len(bpts))) - set(removePtIdxs)) == 0:
        obj.data.splines.remove(oldSpline)
        if len(obj.data.splines) == 0:
            safeRemoveObj(obj)
        return

    insertSpline(obj, splineIdx, splineIdx, removePtIdxs)
    obj.data.splines.remove(obj.data.splines[splineIdx + 1])


def getAdjIdx(obj, splineIdx, startIdx, offset=1, ptCnt=None):
    spline = obj.data.splines[splineIdx]
    if ptCnt is None:
        ptCnt = len(spline.bezier_points)
    if not spline.use_cyclic_u and (
        (startIdx + offset) >= ptCnt or (startIdx + offset) < 0
    ):
        return None
    return (ptCnt + startIdx + offset) % ptCnt  # add ptCnt for negative offset


# Returns a tuple with first value indicating change in spline index (-1, 0, 1)
# and second indicating shift in seg index (negative) due to removal
def removeBezierSeg(obj, splineIdx, segIdx):
    nextIdx = getAdjIdx(obj, splineIdx, segIdx)
    if nextIdx is None:
        return
    spline = obj.data.splines[splineIdx]
    bpts = spline.bezier_points
    ptCnt = len(bpts)
    lastSegIdx = getLastSegIdx(obj, splineIdx)
    splineIdxIncr = 0
    segIdxIncr = 0
    if ptCnt <= 2:
        removeBezierPts(obj, splineIdx, {segIdx, nextIdx})
        # Spline removed by above call
        splineIdxIncr = -1
    else:
        bpt = obj.data.splines[splineIdx].bezier_points[segIdx]
        bpt.handle_right_type = "FREE"
        bpt.handle_left_type = "FREE"
        nextIdx = getAdjIdx(obj, splineIdx, segIdx)
        bpt = obj.data.splines[splineIdx].bezier_points[nextIdx]
        bpt.handle_right_type = "FREE"
        bpt.handle_left_type = "FREE"
        if spline.use_cyclic_u:
            spline.use_cyclic_u = False
            if segIdx != lastSegIdx:
                moveSplineStart(obj, splineIdx, getAdjIdx(obj, splineIdx, segIdx))
                segIdxIncr = -(segIdx + 1)
        else:
            if segIdx == lastSegIdx:
                removeBezierPts(obj, splineIdx, {lastSegIdx + 1})
            elif segIdx == 0:
                removeBezierPts(obj, splineIdx, {0})
                segIdxIncr = -1
            else:
                insertSpline(obj, splineIdx, splineIdx, set(range(segIdx + 1, ptCnt)))
                removeBezierPts(obj, splineIdx + 1, range(segIdx + 1))
                splineIdxIncr = 1
                segIdxIncr = -(segIdx + 1)
    return splineIdxIncr, segIdxIncr


def insertBezierPts(obj, splineIdx, startIdx, cos, handleType, margin=DEF_ERR_MARGIN):
    spline = obj.data.splines[splineIdx]
    bpts = spline.bezier_points

    nextIdx = getAdjIdx(obj, splineIdx, startIdx)

    firstPt = bpts[startIdx]
    nextPt = bpts[nextIdx]

    if firstPt.handle_right_type == "AUTO":
        firstPt.handle_left_type = "ALIGNED"
        firstPt.handle_right_type = "ALIGNED"
    if nextPt.handle_left_type == "AUTO":
        nextPt.handle_left_type = "ALIGNED"
        nextPt.handle_right_type = "ALIGNED"

    fhdl = firstPt.handle_right_type
    nhdl = nextPt.handle_left_type

    firstPt.handle_right_type = "FREE"
    nextPt.handle_left_type = "FREE"

    ptCnt = len(bpts)
    addCnt = len(cos)

    bpts.add(addCnt)
    nextIdx = startIdx + 1

    for i in range(0, (ptCnt - nextIdx)):
        idx = ptCnt - i - 1  # reversed
        offsetIdx = idx + addCnt
        copyObjAttr(bpts[idx], bpts[offsetIdx])

    endIdx = getAdjIdx(obj, splineIdx, nextIdx, addCnt)
    firstPt = bpts[startIdx]
    nextPt = bpts[endIdx]

    prevPt = firstPt
    for i, pt in enumerate(bpts[nextIdx : nextIdx + addCnt]):
        pt.handle_left_type = "FREE"
        pt.handle_right_type = "FREE"

        co = cos[i]
        seg = [prevPt.co, prevPt.handle_right, nextPt.handle_left, nextPt.co]
        t = getTForPt(seg, co, margin)
        ctrlPts0 = getPartialSeg(seg, 0, t)
        ctrlPts1 = getPartialSeg(seg, t, 1)

        segPt = [ctrlPts0[2], ctrlPts1[0], ctrlPts1[1]]

        prevRight = ctrlPts0[1]
        nextLeft = ctrlPts1[2]

        pt.handle_left = segPt[0]
        pt.co = segPt[1]
        pt.handle_right = segPt[2]
        pt.handle_left_type = handleType
        pt.handle_right_type = handleType
        prevPt.handle_right = prevRight

        prevPt = pt
        nextPt.handle_left = nextLeft

    firstPt.handle_right_type = fhdl
    nextPt.handle_left_type = nhdl


# Change position of bezier points according to new matrix_world
def changeMW(obj, newMW):
    invMW = newMW.inverted_safe()
    for spline in obj.data.splines:
        for pt in spline.bezier_points:
            pt.co = invMW @ (obj.mw @ pt.co)
            pt.handle_left = invMW @ (obj.mw @ pt.handle_left)
            pt.handle_right = invMW @ (obj.mw @ pt.handle_right)
    obj.matrix_world = newMW


# Return map in the form of objName->[splineIdx, [startPt, endPt]]
# Remove the invalid keys (if any) from it.
def updateCurveEndPtMap(endPtMap, addObjNames=None, removeObjNames=None):
    invalOs = set()
    if addObjNames is None:
        addObjNames = [o.name for o in bpy.context.scene.objects]
        invalOs = endPtMap.keys() - set(addObjNames)  # In case of redo

    if removeObjNames is not None:
        invalOs.union(set(removeObjNames))

    for o in invalOs:
        del endPtMap[o]

    for objName in addObjNames:
        obj = bpy.context.scene.objects.get(objName)
        if obj is not None and isBezier(obj) and obj.visible_get():
            endPtMap[objName] = []
            mw = obj.matrix_world
            for i, s in enumerate(obj.data.splines):
                pts = [mw @ pt.co for pt in s.bezier_points]
                endPtMap[objName].append([i, pts])


# reverseCurve
def reverseCurve(curve):
    cp = curve.data.copy()
    curve.data.splines.clear()
    for s in reversed(cp.splines):
        ns = curve.data.splines.new("BEZIER")
        copyObjAttr(s, ns)
        ns.bezier_points.add(len(s.bezier_points) - 1)
        for i, p in enumerate(reversed(s.bezier_points)):
            copyObjAttr(p, ns.bezier_points[i])
            ns.bezier_points[i].handle_left_type = "FREE"
            ns.bezier_points[i].handle_right_type = "FREE"
            ns.bezier_points[i].handle_left = p.handle_right
            ns.bezier_points[i].handle_right = p.handle_left

            ns.bezier_points[i].handle_left_type = p.handle_right_type
            ns.bezier_points[i].handle_right_type = p.handle_left_type
    bpy.data.curves.remove(cp)


# splitCurveSelPts
def splitCurveSelPts(selPtMap, newColl=True):
    changeCnt = 0
    newObjs = []

    if len(selPtMap) == 0:
        return newObjs, changeCnt

    for obj in selPtMap.keys():
        splinePtMap = selPtMap.get(obj)

        if (
            len(obj.data.splines) == 1
            and len(obj.data.splines[0].bezier_points) <= 2
            and not obj.data.splines[0].use_cyclic_u
        ) or len(splinePtMap) == 0:
            continue

        keyNames, keyData = getShapeKeyInfo(obj)
        collections = obj.users_collection

        if newColl:
            objGrp = bpy.data.collections.new(obj.name)
            parentColls = [objGrp]
        else:
            parentColls = collections

        splineCnt = len(obj.data.splines)

        endSplineIdx = splineCnt - 1
        if endSplineIdx not in splinePtMap.keys():
            splinePtMap[endSplineIdx] = [
                len(obj.data.splines[endSplineIdx].bezier_points) - 1
            ]

        splineIdxs = sorted(splinePtMap.keys())

        lastSplineIdx = -1
        objCopy = createSkeletalCurve(obj, parentColls)
        newObjs.append(objCopy)
        for i in splineIdxs:
            for j in range(lastSplineIdx + 1, i):
                srcSpline = obj.data.splines[j]
                createSpline(objCopy.data, srcSpline)
                # ~ updateShapeKeyData(objCopy, keyData, keyNames, skStart, ptCnt)
            srcSpline = obj.data.splines[i]
            selPtIdxs = sorted(splinePtMap[i])

            if len(selPtIdxs) == 0:
                createSpline(objCopy.data, srcSpline)
            else:
                bpts = srcSpline.bezier_points
                cyclic = srcSpline.use_cyclic_u
                if cyclic:
                    firstIdx = selPtIdxs[0]
                    moveSplineStart(obj, i, firstIdx)
                    selPtIdxs = [getAdjIdx(obj, i, s, -firstIdx) for s in selPtIdxs]
                    addLastSeg(srcSpline)
                if len(selPtIdxs) > 0 and selPtIdxs[0] == 0:
                    selPtIdxs.pop(0)
                if len(selPtIdxs) > 0 and selPtIdxs[-1] == len(bpts) - 1:
                    selPtIdxs.pop(-1)
                bpts = srcSpline.bezier_points

                if len(selPtIdxs) == 0:
                    segBpts = bpts[: len(bpts)]
                    createSplineForSeg(objCopy.data, segBpts)
                else:
                    lastSegIdx = 0
                    bpts = srcSpline.bezier_points
                    for j in selPtIdxs:
                        segBpts = bpts[lastSegIdx : j + 1]
                        createSplineForSeg(objCopy.data, segBpts)
                        # ~ updateShapeKeyData(objCopy, keyData, keyNames, \
                        # ~ len(newObjs), 2)
                        objCopy = createSkeletalCurve(obj, parentColls)
                        newObjs.append(objCopy)
                        lastSegIdx = j
                    if j != len(bpts) - 1:
                        createSplineForSeg(objCopy.data, bpts[j:])

            lastSplineIdx = i

        if len(objCopy.data.splines) == 0:
            newObjs.remove(objCopy)
            safeRemoveObj(objCopy)

        if newColl:
            for collection in collections:
                collection.children.link(objGrp)

        safeRemoveObj(obj)
        changeCnt += 1

    for obj in newObjs:
        obj.data.splines.active = obj.data.splines[0]

    return newObjs, changeCnt


# splitCurve
def splitCurve(selObjs, split, newColl=True):
    changeCnt = 0
    newObjs = []

    if len(selObjs) == 0:
        return newObjs, changeCnt

    for obj in selObjs:
        if not isBezier(obj) or len(obj.data.splines) == 0:
            continue

        if len(obj.data.splines) == 1:
            if split == "spline":
                newObjs.append(obj)
                continue
            if split == "seg" and len(obj.data.splines[0].bezier_points) <= 2:
                newObjs.append(obj)
                continue
            if split == "point" and len(obj.data.splines[0].bezier_points) == 1:
                newObjs.append(obj)
                continue

        keyNames, keyData = getShapeKeyInfo(obj)
        collections = obj.users_collection

        if newColl:
            objGrp = bpy.data.collections.new(obj.name)
            parentColls = [objGrp]
        else:
            parentColls = collections

        segCnt = 0

        for i, spline in enumerate(obj.data.splines):
            if split == "seg" or split == "point":
                ptLen = len(spline.bezier_points)
                if split == "seg":
                    ptLen -= 1

                for j in range(0, ptLen):
                    objCopy = createSkeletalCurve(obj, parentColls)
                    if split == "seg":
                        createSplineForSeg(
                            objCopy.data, spline.bezier_points[j : j + 2]
                        )
                        updateShapeKeyData(objCopy, keyData, keyNames, len(newObjs), 2)
                    else:  # (split == 'point')
                        mw = obj.matrix_world.copy()
                        newSpline = objCopy.data.splines.new("BEZIER")
                        newPtCo = mw @ spline.bezier_points[j].co.copy()
                        newWM = Matrix()
                        newWM.translation = newPtCo
                        objCopy.matrix_world = newWM
                        copyObjAttr(
                            spline.bezier_points[j],
                            newSpline.bezier_points[0],
                            newWM.inverted_safe(),
                            mw,
                        )

                        # No point having shapekeys (pun intended :)
                        removeShapeKeys(objCopy)

                    newObjs.append(objCopy)

                if split == "seg" and spline.use_cyclic_u:
                    objCopy = createSkeletalCurve(obj, parentColls)
                    createSplineForSeg(
                        objCopy.data,
                        [spline.bezier_points[-1], spline.bezier_points[0]],
                    )
                    updateShapeKeyData(objCopy, keyData, keyNames, -1, 2)
                    newObjs.append(objCopy)

            else:  # split == 'spline'
                objCopy = createSkeletalCurve(obj, parentColls)
                createSpline(objCopy.data, spline)
                currSegCnt = len(objCopy.data.splines[0].bezier_points)
                updateShapeKeyData(objCopy, keyData, keyNames, segCnt, currSegCnt)
                newObjs.append(objCopy)
                segCnt += currSegCnt

        if newColl:
            for collection in collections:
                collection.children.link(objGrp)

        safeRemoveObj(obj)
        changeCnt += 1

    for obj in newObjs:
        obj.data.splines.active = obj.data.splines[0]

    return newObjs, changeCnt

    # getClosestCurve


def getClosestCurve(srcMW, pt, curves, minDist=9e99):
    closestCurve = None
    for i, curve in enumerate(curves):
        mw = curve.matrix_world
        addLastSeg(curve.data.splines[0])

        start = curve.data.splines[0].bezier_points[0]
        end = curve.data.splines[-1].bezier_points[-1]
        dist = ((mw @ start.co) - (srcMW @ pt)).length
        if dist < minDist:
            minDist = dist
            closestCurve = curve
        dist = ((mw @ end.co) - (srcMW @ pt)).length
        if dist < minDist:
            minDist = dist
            reverseCurve(curve)
            closestCurve = curve
    return closestCurve, minDist


# getCurvesArrangedByDist
def getCurvesArrangedByDist(curves):
    idMap = {c.name: c for c in curves}
    orderedCurves = [curves[0].name]
    nextCurve = curves[0]
    remainingCurves = curves[1:]

    # Arrange in order
    while len(remainingCurves) > 0:
        addLastSeg(nextCurve.data.splines[-1])

        srcMW = nextCurve.matrix_world

        ncEnd = nextCurve.data.splines[-1].bezier_points[-1]
        closestCurve, dist = getClosestCurve(srcMW, ncEnd.co, remainingCurves)

        # Check the start also for the first curve
        if len(orderedCurves) == 1:
            ncStart = nextCurve.data.splines[0].bezier_points[0]
            closestCurve2, dist2 = getClosestCurve(
                srcMW, ncStart.co, remainingCurves, dist
            )
            if closestCurve2 is not None:
                reverseCurve(nextCurve)
                closestCurve = closestCurve2

        orderedCurves.append(closestCurve.name)
        nextCurve = closestCurve
        remainingCurves.remove(closestCurve)
    return [idMap[cn] for cn in orderedCurves]


def joinSegs(curves, optimized, straight, srcCurve=None, margin=DEF_ERR_MARGIN):
    if len(curves) == 0:
        return None
    if len(curves) == 1:
        return curves[0]

    if optimized:
        curves = getCurvesArrangedByDist(curves)

    firstCurve = curves[0]

    if srcCurve is None:
        srcCurve = firstCurve

    elif srcCurve != firstCurve:
        srcCurveData = srcCurve.data.copy()
        changeMW(firstCurve, srcCurve.matrix_world)
        srcCurve.data = firstCurve.data

    srcMW = srcCurve.matrix_world
    invSrcMW = srcMW.inverted_safe()
    newCurveData = srcCurve.data

    for curve in curves[1:]:
        if curve == srcCurve:
            curveData = srcCurveData
        else:
            curveData = curve.data

        mw = curve.matrix_world

        currSpline = newCurveData.splines[-1]
        nextSpline = curveData.splines[0]

        addLastSeg(currSpline)
        addLastSeg(nextSpline)

        currBezierPt = currSpline.bezier_points[-1]
        nextBezierPt = nextSpline.bezier_points[0]

        # Don't add new point if the last one and the current one are the 'same'

        # === curve_helpers ===
        # Don't add new point if the last one and the current one are the 'same'
        if vectCmpWithMargin(srcMW @ currBezierPt.co, mw @ nextBezierPt.co, margin):
            currBezierPt.handle_right_type = nextBezierPt.handle_right_type
            if currBezierPt.handle_right_type != "VECTOR":
                currBezierPt.handle_right_type = "FREE"
            currBezierPt.handle_right = invSrcMW @ (mw @ nextBezierPt.handle_right)
            ptIdx = 1
        else:
            ptIdx = 0

        if straight and ptIdx == 0:
            currBezierPt.handle_left_type = "FREE"
            currBezierPt.handle_right_type = "VECTOR"
            # ~ currBezierPt.handle_right = currBezierPt.co

        for i in range(ptIdx, len(nextSpline.bezier_points)):
            if (i == len(nextSpline.bezier_points) - 1) and vectCmpWithMargin(
                mw @ nextSpline.bezier_points[i].co,
                srcMW @ currSpline.bezier_points[0].co,
                margin,
            ):
                currSpline.bezier_points[0].handle_left_type = "FREE"
                currSpline.bezier_points[0].handle_left = invSrcMW @ (
                    mw @ nextSpline.bezier_points[i].handle_left
                )

                currSpline.use_cyclic_u = True
                break

            currSpline.bezier_points.add(1)
            currBezierPt = currSpline.bezier_points[-1]
            copyObjAttr(nextSpline.bezier_points[i], currBezierPt, invSrcMW, mw)

            if straight and i == 0:
                currBezierPt.handle_right_type = "FREE"
                currBezierPt.handle_left_type = "VECTOR"
                # ~ currBezierPt.handle_left = currBezierPt.co

        # Simply add the remaining splines
        for spline in curveData.splines[1:]:
            newSpline = newCurveData.splines.new("BEZIER")
            copyObjAttr(spline, newSpline)
            for i, pt in enumerate(spline.bezier_points):
                if i > 0:
                    newSpline.bezier_points.add(1)
                copyObjAttr(pt, newSpline.bezier_points[-1], invSrcMW, mw)

        if curve != srcCurve:
            safeRemoveObj(curve)

    if firstCurve != srcCurve:
        safeRemoveObj(firstCurve)

    return srcCurve


def removeDupliVert(curve, margin):
    newCurveData = curve.data.copy()
    newCurveData.splines.clear()
    dupliFound = False
    for spline in curve.data.splines:
        newCurveData.splines.new("BEZIER")
        currSpline = newCurveData.splines[-1]
        copyObjAttr(spline, currSpline)

        if len(spline.bezier_points) == 1:
            copyObjAttr(spline.bezier_points[0], currSpline.bezier_points[0])
            continue

        cmpPts = spline.bezier_points[:]
        pt0 = cmpPts[0]
        while vectCmpWithMargin(cmpPts[-1].co, pt0.co, margin) and len(cmpPts) > 1:
            endPt = cmpPts.pop()
            pt0.handle_left_type = "FREE"
            pt0.handle_right_type = "FREE"
            pt0.handle_left = endPt.handle_left
            currSpline.use_cyclic_u = True
            dupliFound = True

        prevPt = None
        for pt in cmpPts:
            if prevPt is not None and vectCmpWithMargin(prevPt.co, pt.co, margin):
                copyObjAttr(pt, currSpline.bezier_points[-1])
                dupliFound = True
            else:
                if prevPt is not None:
                    currSpline.bezier_points.add(1)
                copyObjAttr(pt, currSpline.bezier_points[-1])
            prevPt = pt

    if dupliFound:
        curve.data = newCurveData
    else:
        bpy.data.curves.remove(newCurveData)


def convertToFace(curve, remeshRes, perSeg, fillType, optimized):
    bptData = getBptData(curve, local=True)
    bm = bmesh.new()
    splineLens = [spline.calc_length() for spline in curve.data.splines]
    maxSplineLen = max(splineLens)
    centers = []
    normals = []

    for splineIdx, spline in enumerate(curve.data.splines):
        bpts = spline.bezier_points
        verts = []
        addLastVert = spline.use_cyclic_u or fillType != "NONE"
        if not perSeg and remeshRes > 0:
            segPts = [bptData[splineIdx][x] for x in range(len(bpts))]
            if addLastVert:
                segPts.append(segPts[0])
            numSegs = int(remeshRes * splineLens[splineIdx] / maxSplineLen)
            if numSegs <= 2:
                vertCos = [bpts[0].co, bpts[-1].co]
            else:
                pts = getInterpBezierPts(segPts, subdivPerUnit=100, segLens=None)
                vertCos = getInterpolatedVertsCo(pts, numSegs)
            for co in vertCos:
                verts.append(bm.verts.new(co))
        else:
            if remeshRes > 0:
                segPtPairs = [
                    getSegPtsInSpline(bptData, splineIdx, ptIdx, addLastVert)
                    for ptIdx in range(len(bpts))
                ]
                if not addLastVert:
                    segPtPairs.pop()
                segs = [
                    [segPts[0][1], segPts[0][2], segPts[1][0], segPts[1][1]]
                    for segPts in segPtPairs
                ]
                segLens = [getSegLen(seg) for seg in segs]
                maxLen = max(segLens)
                for i, segPts in enumerate(segPtPairs):
                    if optimized and isStraightSeg(segPts):
                        verts.append(bm.verts.new(segPts[0][1]))
                        verts.append(bm.verts.new(segPts[1][1]))
                    else:
                        pts = getInterpBezierPts(
                            segPts, subdivPerUnit=100, segLens=[segLens[i]]
                        )
                        numSegs = ceil(remeshRes * segLens[i] / maxLen)
                        vertCos = getInterpolatedVertsCo(pts, numSegs)
                        for j, co in enumerate(vertCos):
                            if j == 0 and i > 0:
                                continue
                            verts.append(bm.verts.new(co))
            else:
                for ptIdx, bpt in enumerate(bpts):
                    verts.append(bm.verts.new(bpts[ptIdx].co))
                if addLastVert:
                    verts.append(bm.verts.new(bpts[0].co))
        if len(verts) < 2:
            pass
        elif len(verts) == 2:
            bm.edges.new(verts)
        else:
            if spline.use_cyclic_u:
                bm.verts.remove(verts[-1])
                verts.pop()

            vertCos = [v.co for v in verts]
            center = Vector(
                [sum(vertCos[i][j] for i in range(len(vertCos))) for j in range(3)]
            ) / len(vertCos)
            normal = geometry.normal(vertCos)

            centers.append(center)
            normals.append(normal)

            if fillType == "NGON":
                bm.faces.new(verts)

            elif fillType == "NONE":
                for i in range(1, len(verts)):
                    bm.edges.new([verts[i - 1], verts[i]])
                if spline.use_cyclic_u:
                    bm.edges.new([verts[-1], verts[0]])
            elif fillType == "FAN":
                centerVert = bm.verts.new(center)
                for i in range(1, len(verts)):
                    bm.faces.new([centerVert, verts[i - 1], verts[i]])
                if spline.use_cyclic_u:
                    bm.faces.new([centerVert, verts[-1], verts[0]])
    cnt = len(centers)
    if cnt > 0:
        center = (
            Vector([sum(centers[i][j] for i in range(cnt)) for j in range(3)]) / cnt
        )
        normal = (
            Vector([sum(normals[i][j] for i in range(cnt)) for j in range(3)]) / cnt
        )
    else:
        center = None
        normal = None
    m = bpy.data.meshes.new(curve.data.name)
    bm.to_mesh(m)
    meshObj = bpy.data.objects.new(curve.name, m)
    collections = curve.users_collection
    for c in collections:
        c.objects.link(meshObj)

    return meshObj, center, normal


def convertToMesh(curve):
    mt = curve.to_mesh()  # Can't be used directly
    bm = bmesh.new()
    bm.from_mesh(mt)
    m = bpy.data.meshes.new(curve.data.name)
    bm.to_mesh(m)
    meshObj = bpy.data.objects.new(curve.name, m)
    collections = curve.users_collection
    for c in collections:
        c.objects.link(meshObj)
    return meshObj


def applyMeshModifiers(meshObj, remeshDepth):
    bpy.context.view_layer.objects.active = meshObj
    normal = geometry.normal([v.co for v in meshObj.data.vertices])
    normal = Vector([round(c, 5) for c in normal])
    if vectCmpWithMargin(normal, Vector()):
        return meshObj

    planeVert = Vector([round(c, 5) for c in meshObj.data.vertices[0].co])
    mod = meshObj.modifiers.new("mod", type="SOLIDIFY")
    mod.thickness = 20
    bpy.ops.object.modifier_apply(modifier=mod.name)

    mod = meshObj.modifiers.new("mod", type="REMESH")
    mod.octree_depth = remeshDepth
    mod.use_remove_disconnected = False

    bpy.ops.object.modifier_apply(modifier=mod.name)

    bm = bmesh.new()
    bm.from_mesh(meshObj.data)
    bm.verts.ensure_lookup_table()
    toRemove = []
    for i, v in enumerate(bm.verts):
        co = Vector([round(c, 5) for c in v.co])
        if (
            abs(geometry.distance_point_to_plane(co, planeVert, normal))
            > DEF_ERR_MARGIN
        ):
            toRemove.append(v)
    for v in toRemove:
        bm.verts.remove(v)
    bm.to_mesh(meshObj.data)


def unsubdivideObj(meshObj):
    bm = bmesh.new()
    bm.from_object(meshObj, bpy.context.evaluated_depsgraph_get())
    bmesh.ops.unsubdivide(bm, verts=bm.verts)
    bm.to_mesh(meshObj.data)


def pasteLength(src, dests):
    tmp = bpy.data.curves.new("t", "CURVE")
    ts = tmp.splines.new("BEZIER")
    ts.bezier_points.add(1)
    mw = src.matrix_world
    srcLen = sum(getSplineLenTmpObj(ts, s, mw) for s in src.data.splines)

    for c in dests:
        mw = c.matrix_world
        destLen = sum(getSplineLenTmpObj(ts, s, mw) for s in c.data.splines)
        fact = srcLen / destLen
        for s in c.data.splines:
            lts = []
            rts = []
            for pt in s.bezier_points:
                lts.append(pt.handle_left_type)
                rts.append(pt.handle_right_type)
                pt.handle_left_type = "FREE"
                pt.handle_right_type = "FREE"
            for pt in s.bezier_points:
                pt.co = fact * pt.co
                pt.handle_left = fact * pt.handle_left
                pt.handle_right = fact * pt.handle_right
            for i, pt in enumerate(s.bezier_points):
                pt.handle_left_type = lts[i]
                pt.handle_right_type = rts[i]
    bpy.data.curves.remove(tmp)


def intersectCurves(curves, action, firstActive, margin, rounding):
    allIntersectCos, intersectMap = getCurveIntersectPts(
        curves, firstActive, margin, rounding
    )

    if action in {"MARK_EMPTY", "MARK_POINT"}:
        newObjs = []
        collection = bpy.data.collections.new("Intersect Markers")
        bpy.context.scene.collection.children.link(collection)
        objName = "Marker"
        for co in allIntersectCos:
            if action == "MARK_EMPTY":
                obj = bpy.data.objects.new(objName, None)
            elif action == "MARK_POINT":
                curveData = bpy.data.curves.new(objName, "CURVE")
                curveData.splines.new("BEZIER")
                obj = bpy.data.objects.new(objName, curveData)
            obj.location = co
            newObjs.append(obj)
            collection.objects.link(obj)
        for obj in newObjs:
            obj.select_set(True)

    if action in {"INSERT_PT", "CUT"}:
        mapKeys = sorted(intersectMap.keys(), key=lambda x: (x[0], x[1], x[2]))
        selPtMap = {}
        prevCnt = 0
        prevCurveIdx = None
        prevSplineIdx = None
        for i, key in enumerate(mapKeys):
            intersectPts = intersectMap[key]

            curveIdx, splineIdx, segIdx = key
            if prevCurveIdx == curveIdx and prevSplineIdx == splineIdx:
                segIdx += prevCnt
            else:
                prevCnt = 0

            curve = curves[curveIdx]

            if firstActive and curve == curves[0]:
                continue

            mw = curve.matrix_world
            pts = curve.data.splines[splineIdx].bezier_points
            nextIdx = getAdjIdx(curve, splineIdx, segIdx)

            seg = [
                mw @ pts[segIdx].co,
                mw @ pts[segIdx].handle_right,
                mw @ pts[nextIdx].handle_left,
                mw @ pts[nextIdx].co,
            ]

            sortedCos = getCosSortedByT(seg, intersectPts, margin)
            sortedCos = removeDupliCos(sortedCos, margin)
            insertCos = sortedCos.copy()
            start = seg[0].freeze()
            end = seg[1].freeze()
            startIdxIncr = 1
            endIdxIncr = 1
            if start in insertCos:
                startIdxIncr = 0  # Keep in split list
                insertCos.remove(start)  # Remove from insert list
            if end in insertCos:
                insertCos.remove(end)  # Remove from insert list
                endIdxIncr = 2  # Keep in split list
            for j, co in enumerate(insertCos):
                insertCos[j] = mw.inverted_safe() @ co
            insertBezierPts(curve, splineIdx, segIdx, insertCos, "FREE", margin)
            prevCurveIdx = curveIdx
            prevSplineIdx = splineIdx
            prevCnt += len(insertCos)

            if action == "CUT":
                if selPtMap.get(curve) is None:
                    selPtMap[curve] = {}
                if selPtMap[curve].get(splineIdx) is None:
                    selPtMap[curve][splineIdx] = []
                selPtMap[curve][splineIdx] += list(
                    range(segIdx + startIdxIncr, segIdx + len(sortedCos) + endIdxIncr)
                )

        if action == "CUT":
            for curve in list(selPtMap.keys()):
                splineIdxs = sorted(selPtMap[curve].keys())
                if len(curve.data.splines) > 1:
                    newObjs, changeCnt = splitCurve(
                        [curve], "spline", curve.users_collection
                    )
                    splineIdxs = list(selPtMap[curve].keys())
                    for idx in splineIdxs:
                        newCurve = newObjs[idx]
                        selPtMap[newCurve] = {}
                        selPtMap[newCurve][0] = selPtMap[curve][idx]
                        selPtMap[curve].pop(idx)
                        if len(selPtMap[curve]) == 0:
                            selPtMap.pop(curve)
            splitCurveSelPts(selPtMap)


def getSVGPt(co, docW, docH, camera=None, region=None, rv3d=None):
    if camera is not None:
        scene = bpy.context.scene
        xy = world_to_camera_view(scene, camera, co)
        return complex(xy[0] * docW, docH - (xy[1] * docH))
    elif region is not None and rv3d is not None:
        xy = getCoordFromLoc(region, rv3d, co)
        return complex(xy[0], docH - xy[1])


def getPathD(path):
    curve = ""

    for i, part in enumerate(path):
        comps = []
        for j, segment in enumerate(part):
            if j == 0:
                comps.append("M {},{} C".format(segment[0].real, segment[0].imag))
            args = (
                segment[1].real,
                segment[1].imag,
                segment[2].real,
                segment[2].imag,
                segment[3].real,
                segment[3].imag,
            )
            comps.append("{},{} {},{} {},{}".format(*args))
        curve += " ".join(comps)

    return curve


def getPathBBox(path):
    minX, minY, maxX, maxY = [None, None, None, None]
    for part in path:
        for seg in part:
            seg3d = [(seg[i].real, seg[i].imag, 0) for i in range(len(seg))]
            leftBotFront, rgtTopBack = getBBox(seg3d)
            if minX is None or leftBotFront[0] < minX:
                minX = leftBotFront[0]
            if minY is None or leftBotFront[1] < minY:
                minY = leftBotFront[1]
            if maxX is None or rgtTopBack[0] > maxX:
                maxX = rgtTopBack[0]
            if maxY is None or rgtTopBack[1] > maxY:
                maxY = rgtTopBack[1]
    return minX, minY, maxX, maxY


def createClipElem(doc, svgElem, docW, docH, clipElemId):
    elem = doc.createElement("defs")
    svgElem.appendChild(elem)
    clipElem = doc.createElement("clipPath")
    clipElem.setAttribute("clipPathUnits", "userSpaceOnUse")
    clipElem.setAttribute("id", clipElemId)
    elem.appendChild(clipElem)
    rectElem = doc.createElement("rect")
    rectElem.setAttribute("x", "0")
    rectElem.setAttribute("y", "0")
    rectElem.setAttribute("width", str(docW))
    rectElem.setAttribute("height", str(docH))
    clipElem.appendChild(rectElem)


def getSVGPathElem(
    doc,
    docW,
    docH,
    path,
    idx,
    lineWidth,
    lineCol,
    lineAlpha,
    fillCol,
    fillAlpha,
    clipView,
    clipElemId,
):
    idPrefix = "id"
    style = {
        "opacity": "1",
        "stroke": "#000000",
        "stroke-width": "1",
        "fill": "none",
        "stroke-linecap": "round",
        "stroke-linejoin": "miter",
        "stroke-miterlimit": "4",
    }

    clipped = False
    if clipView:
        minX, minY, maxX, maxY = getPathBBox(path)
        if maxX < 0 or maxY < 0 or minX > docW or minY > docH:
            return None

        if minX < 0 or minY < 0 or maxX > docW or maxY > docH:
            clipped = True

    elem = doc.createElement("path")
    elem.setAttribute("id", idPrefix + str(idx).zfill(3))
    elem.setAttribute("d", getPathD(path))
    style["stroke-width"] = str(lineWidth)
    style["stroke"] = "#" + lineCol
    style["opacity"] = lineAlpha
    if fillCol is not None:
        style["fill"] = "#" + fillCol
        style["opacity"] = fillAlpha  # Overwrite
    styleStr = ";".join([k + ":" + style[k] for k in style])
    elem.setAttribute("style", styleStr)

    if clipped:
        elem.setAttribute("clip-path", "url(#" + clipElemId + ")")

    return elem


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


###################### Operators ######################


def getBptData(
    obj,
    withShapeKey=True,
    shapeKeyIdx=None,
    fromMix=True,
    updateDeps=False,
    local=False,
):
    # Less readable but more convenient than class
    # Format: [handle_left, co, handle_right, handle_left_type, handle_right_type]
    worldSpaceData = []
    mw = Matrix() if local else obj.matrix_world

    keydata = None
    dataIdx = 0
    shapeKey = obj.active_shape_key
    tmpsk = None
    if withShapeKey and shapeKey is not None:
        if shapeKeyIdx is None:
            shapeKeyIdx = obj.active_shape_key_index
        if fromMix:
            if not obj.data.shape_keys.use_relative:
                val = obj.data.shape_keys.eval_time
            else:
                val = obj.data.shape_keys.key_blocks[obj.active_shape_key_index].value

            if floatCmpWithMargin(val, 0):
                keyBlock = obj.data.shape_keys.key_blocks[0]
            else:
                tmpsk = obj.shape_key_add(name="tmp", from_mix=True)
                keyBlock = obj.data.shape_keys.key_blocks[tmpsk.name]
        else:
            keyBlock = obj.data.shape_keys.key_blocks[shapeKeyIdx]
        keydata = keyBlock.data

    for spline in obj.data.splines:
        pts = []
        for pt in spline.bezier_points:
            lt, rt = pt.handle_left_type, pt.handle_right_type
            if keydata is not None:
                pt = keydata[dataIdx]
                dataIdx += 1

            pts.append([mw @ pt.handle_left, mw @ pt.co, mw @ pt.handle_right, lt, rt])
        worldSpaceData.append(pts)
    if tmpsk is not None:
        obj.shape_key_remove(tmpsk)
        obj.active_shape_key_index = shapeKeyIdx
        if updateDeps:
            bpy.context.evaluated_depsgraph_get().update()
    return worldSpaceData


def getBezierDataForSeg(
    obj,
    splineIdx,
    segIdx,
    withShapeKey=True,
    shapeKeyIdx=None,
    fromMix=True,
    updateDeps=False,
):
    wsData = getBptData(obj, withShapeKey, shapeKeyIdx, fromMix, updateDeps)
    pt0 = wsData[splineIdx][segIdx]
    segEndIdx = getAdjIdx(obj, splineIdx, segIdx)
    if segEndIdx is None:
        return []
    pt1 = wsData[splineIdx][segEndIdx]
    return [pt0, pt1]

def getRoundedSplineSegs(mw, spline, reverse = False, rounding = 5):

    def getRoundedVect(mw, co, rounding):
        return Vector((round(x, rounding) for x in (mw @ co)))

    bpts = spline.bezier_points
    if(reverse): bpts = reversed(bpts)
    segs = []
    for i in range(1, len(bpts)):
        segPts = [bpts[i-1].co, bpts[i-1].handle_right, bpts[i].handle_left, bpts[i].co]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    if(spline.use_cyclic_u):
        segPts = [bpts[-1].co, bpts[-1].handle_right, bpts[0].handle_left, bpts[0].co]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    return segs

# Because there is some discrepancy between this and getSegLen
# This seems to be more accurate
def getSegLenTmpObj(tmpSpline, bpts, mw = Matrix()):
    tmpSpline.bezier_points[0].co = mw @ bpts[0].co
    tmpSpline.bezier_points[0].handle_right = mw @ bpts[0].handle_right
    tmpSpline.bezier_points[1].handle_left = mw @ bpts[1].handle_left
    tmpSpline.bezier_points[1].co = mw @ bpts[1].co
    return tmpSpline.calc_length()

def getSplineLenTmpObj(tmpSpline, spline, mw):
    l = 0
    bpts = spline.bezier_points
    l += sum(getSegLenTmpObj(tmpSpline, bpts[i:i+2], mw) for i in range(len(bpts) -1))
    if(spline.use_cyclic_u): l += getSegLenTmpObj(tmpSpline, [bpts[-1], bpts[0]], mw)
    return l

def getSplineIntersectPts(curves, splineInfos, firstActive, margin, rounding):
    soln = []
    solnRounded = set()
    allIntersectCos = []
    intersectMap = {} # (curveIdx, splineIdx) -> [ptIdx, ptIdx...]

    for i in range(len(splineInfos)):
        curveIdx0, splineIdx0 = splineInfos[i]
        spline0 = curves[curveIdx0].data.splines[splineIdx0]
        mw0 = curves[curveIdx0].matrix_world
        segs0 = getRoundedSplineSegs(mw0, spline0, rounding=rounding)

        for j in range(i + 1, len(splineInfos)):
            curveIdx1, splineIdx1 = splineInfos[j]
            if(firstActive and curveIdx0 == 0 and curveIdx1 == 0): continue

            spline1 = curves[curveIdx1].data.splines[splineIdx1]
            mw1 = curves[curveIdx1].matrix_world
            segs1 = getRoundedSplineSegs(mw1, spline1, rounding=rounding)

            for k, seg0 in enumerate(segs0):
                for l, seg1 in enumerate(segs1):
                    if(getIntersectPts(seg0, seg1, soln, solnRounded, 0, margin, rounding)):
                        pass

    for co in soln:
        allIntersectCos.append(co)

    return allIntersectCos, intersectMap

def getCurveIntersectPts(curves, firstActive, margin, rounding):
    splineInfos = [(x, y) for x in range(len(curves)) \
        for y in range(len(curves[x].data.splines))]

    return getSplineIntersectPts(curves, splineInfos, firstActive, margin, rounding)
