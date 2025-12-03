# bezier_utils/utils/curve_utils.py

import bpy
import bmesh
from math import ceil
from mathutils import Vector, Matrix
from ..constants import DEF_ERR_MARGIN
from .math_utils import vectCmpWithMargin, floatCmpWithMargin, toHexStr
from .bezier_math import getInterpolatedVertsCo, getBBox
from .bezier_math import (
    getTForPt,
    getPartialSeg,
    getInterpBezierPts,
    getSegLen,
    getIntersectPts,
    getPtFromT,
)
from .object_utils import isBezier, safeRemoveObj, copyObjAttr


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


def intersectCurves(curves, action, firstActive, margin, rounding, selfIntersect=False):
    allIntersectCos, intersectMap = getCurveIntersectPts(
        curves, firstActive, margin, rounding, selfIntersect
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


def getRoundedSplineSegs(mw, spline, reverse=False, rounding=5):
    def getRoundedVect(mw, co, rounding):
        return Vector((round(x, rounding) for x in (mw @ co)))

    bpts = spline.bezier_points
    if reverse:
        bpts = reversed(bpts)
    segs = []
    for i in range(1, len(bpts)):
        segPts = [
            bpts[i - 1].co,
            bpts[i - 1].handle_right,
            bpts[i].handle_left,
            bpts[i].co,
        ]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    if spline.use_cyclic_u:
        segPts = [bpts[-1].co, bpts[-1].handle_right, bpts[0].handle_left, bpts[0].co]
        segs.append([getRoundedVect(mw, pt, rounding) for pt in segPts])
    return segs


# Because there is some discrepancy between this and getSegLen
# This seems to be more accurate
def getSegLenTmpObj(tmpSpline, bpts, mw=Matrix()):
    tmpSpline.bezier_points[0].co = mw @ bpts[0].co
    tmpSpline.bezier_points[0].handle_right = mw @ bpts[0].handle_right
    tmpSpline.bezier_points[1].handle_left = mw @ bpts[1].handle_left
    tmpSpline.bezier_points[1].co = mw @ bpts[1].co
    return tmpSpline.calc_length()


def getSplineLenTmpObj(tmpSpline, spline, mw):
    l = 0
    bpts = spline.bezier_points
    l += sum(
        getSegLenTmpObj(tmpSpline, bpts[i : i + 2], mw) for i in range(len(bpts) - 1)
    )
    if spline.use_cyclic_u:
        l += getSegLenTmpObj(tmpSpline, [bpts[-1], bpts[0]], mw)
    return l


def getSplineIntersectPts(
    curves, splineInfos, firstActive, margin, rounding, selfIntersect=False
):
    segPairMap = {}
    for i, c0Info in enumerate(splineInfos):
        idxCurve0, idxSpline0 = c0Info
        if firstActive and idxCurve0 != splineInfos[0][0]:
            break
        c0 = curves[idxCurve0]
        spline0 = c0.data.splines[idxSpline0]

        # Self-intersection: check segments within same spline
        if selfIntersect:
            segs0 = getRoundedSplineSegs(c0.matrix_world, spline0)
            numSegs = len(segs0)
            for idxSeg0 in range(numSegs):
                # Skip adjacent segments (they share endpoints)
                for idxSeg1 in range(idxSeg0 + 2, numSegs):
                    # For cyclic splines, skip last-to-first pair
                    if spline0.use_cyclic_u and idxSeg0 == 0 and idxSeg1 == numSegs - 1:
                        continue
                    key = ((idxCurve0, idxSpline0), (idxCurve0, idxSpline0))
                    if key not in segPairMap:
                        segPairMap[key] = []
                    segPairMap[key].append(
                        ((idxSeg0, segs0[idxSeg0]), (idxSeg1, segs0[idxSeg1]))
                    )

        # Cross-spline intersections
        for j in range(i + 1, len(splineInfos)):
            idxCurve1, idxSpline1 = splineInfos[j]
            c1 = curves[idxCurve1]
            spline1 = c1.data.splines[idxSpline1]

            segPairMap[((idxCurve0, idxSpline0), (idxCurve1, idxSpline1))] = [
                ((idxSeg0, s0), (idxSeg1, s1))
                for idxSeg0, s0 in enumerate(
                    getRoundedSplineSegs(c0.matrix_world, spline0)
                )
                for idxSeg1, s1 in enumerate(
                    getRoundedSplineSegs(c1.matrix_world, spline1)
                )
            ]

    intersectMap = {}
    allIntersectCos = []
    for key in segPairMap:
        (idxCurve0, idxSpline0), (idxCurve1, idxSpline1) = key
        segPairInfo = segPairMap[key]
        for info in segPairInfo:
            solnRounded = set()
            soln = []
            (idxSeg0, seg0), (idxSeg1, seg1) = info
            ret = getIntersectPts(
                seg0,
                seg1,
                soln,
                solnRounded,
                recurs=0,
                margin=margin,
                rounding=rounding,
            )
            if ret:
                extKey = (idxCurve0, idxSpline0, idxSeg0)
                if intersectMap.get(extKey) is None:
                    intersectMap[extKey] = []
                intersectMap[extKey] += soln

                # For self-intersection, both segments are on same spline
                if not (idxCurve0 == idxCurve1 and idxSpline0 == idxSpline1):
                    extKey = (idxCurve1, idxSpline1, idxSeg1)
                    if intersectMap.get(extKey) is None:
                        intersectMap[extKey] = []
                    intersectMap[extKey] += soln
                else:
                    # Self-intersection: add to second segment too
                    extKey = (idxCurve1, idxSpline1, idxSeg1)
                    if intersectMap.get(extKey) is None:
                        intersectMap[extKey] = []
                    intersectMap[extKey] += soln

                allIntersectCos += soln

    return allIntersectCos, intersectMap


def getCurveIntersectPts(curves, firstActive, margin, rounding, selfIntersect=False):
    splineInfos = [
        (x, y) for x in range(len(curves)) for y in range(len(curves[x].data.splines))
    ]

    return getSplineIntersectPts(
        curves, splineInfos, firstActive, margin, rounding, selfIntersect
    )


def getCosSortedByT(seg, cos, margin):
    from .bezier_math import getTForPt

    coInfo = set()

    for co in cos:
        t = getTForPt(seg, co, margin)
        if t >= 1:
            coInfo.add(((seg[3]).freeze(), 1))
        elif t <= 0:
            coInfo.add(((seg[0]).freeze(), 0))
        else:
            coInfo.add((co.freeze(), t))

    return [inf[0] for inf in sorted(coInfo, key=lambda x: x[1])]


def removeDupliCos(sortedCos, margin):
    from .math_utils import vectCmpWithMargin

    prevCo = sortedCos[0]
    newCos = [prevCo]
    for i in range(1, len(sortedCos)):
        co = sortedCos[i]
        if not vectCmpWithMargin(co, prevCo, margin):
            newCos.append(co)
        prevCo = co
    return newCos


# ============== Boolean Operations ==============


def isPointInsideCurve(pt, curve, margin=DEF_ERR_MARGIN):
    """Ray casting algorithm to check if point is inside a closed curve (2D, XY plane)."""
    mw = curve.matrix_world
    crossings = 0

    for spline in curve.data.splines:
        if not spline.use_cyclic_u:
            continue
        segs = getRoundedSplineSegs(mw, spline)

        for seg in segs:
            # Sample segment and count ray crossings (ray goes in +X direction)
            samples = 50  # Increased for better accuracy
            for i in range(samples):
                t0, t1 = i / samples, (i + 1) / samples
                p0 = getPtFromT(seg[0], seg[1], seg[2], seg[3], t0)
                p1 = getPtFromT(seg[0], seg[1], seg[2], seg[3], t1)

                # Check if edge crosses the horizontal ray from pt going right
                if (p0.y > pt.y) != (p1.y > pt.y):
                    x_intersect = p0.x + (pt.y - p0.y) * (p1.x - p0.x) / (p1.y - p0.y)
                    if pt.x < x_intersect:
                        crossings += 1

    return crossings % 2 == 1


def getSegmentMidpoint(curve, splineIdx, segIdx):
    """Get midpoint of a bezier segment."""
    mw = curve.matrix_world
    spline = curve.data.splines[splineIdx]
    bpts = spline.bezier_points
    nextIdx = getAdjIdx(curve, splineIdx, segIdx)

    seg = [
        mw @ bpts[segIdx].co,
        mw @ bpts[segIdx].handle_right,
        mw @ bpts[nextIdx].handle_left,
        mw @ bpts[nextIdx].co,
    ]
    return getPtFromT(seg[0], seg[1], seg[2], seg[3], 0.5)


def booleanCurves(curves, operation, margin, rounding):
    """
    Perform boolean operation on multiple closed bezier curves.
    operation: 'UNION', 'DIFFERENCE', or 'INTERSECTION'
    Returns new curve object(s).

    For multiple curves:
    - UNION: Combines all curves
    - INTERSECTION: Keeps only area common to all curves
    - DIFFERENCE: Subtracts all other curves from the first (active) curve
    """
    if len(curves) < 2:
        return []

    # Check all curves are closed
    if not all(s.use_cyclic_u for c in curves for s in c.data.splines):
        return []

    # For 2 curves, use direct operation
    if len(curves) == 2:
        return booleanCurvesPair(curves[0], curves[1], operation, margin, rounding)

    # For UNION with multiple curves, we need to union each new curve with ALL existing results
    if operation == "UNION":
        return booleanUnionMultiple(curves, margin, rounding)

    # For INTERSECTION and DIFFERENCE, apply iteratively
    result = booleanCurvesPair(curves[0], curves[1], operation, margin, rounding)

    for i in range(2, len(curves)):
        if not result:
            break

        # Merge multiple results into one curve for next operation
        if len(result) > 1:
            workingCurve = mergeSplinesToCurve(result)
            for obj in result:
                safeRemoveObj(obj)
        else:
            workingCurve = result[0]

        newResult = booleanCurvesPair(
            workingCurve, curves[i], operation, margin, rounding
        )
        safeRemoveObj(workingCurve)
        result = newResult

    return result


def booleanUnionMultiple(curves, margin, rounding):
    """
    Union multiple curves properly by trying to merge each new curve with existing results.
    """
    # Start with first curve as a copy
    results = [duplicateCurve(curves[0])]

    for i in range(1, len(curves)):
        newCurve = curves[i]
        mergedWith = -1

        # Try to union newCurve with each existing result
        for j, existing in enumerate(results):
            unionResult = booleanCurvesPair(
                existing, newCurve, "UNION", margin, rounding
            )

            if len(unionResult) == 1:
                # Successfully merged
                results[j] = unionResult[0]
                safeRemoveObj(existing)
                mergedWith = j
                break
            else:
                # Disjoint - clean up the copies returned
                for obj in unionResult:
                    safeRemoveObj(obj)

        if mergedWith < 0:
            # newCurve didn't merge with any existing result, add it as new
            results.append(duplicateCurve(newCurve))
        else:
            # After merging, try to merge the new combined result with other results
            # This handles cases where A and C don't intersect, but (A∪B) intersects C
            changed = True
            while changed and len(results) > 1:
                changed = False
                for j in range(len(results)):
                    if j == mergedWith:
                        continue
                    unionResult = booleanCurvesPair(
                        results[mergedWith], results[j], "UNION", margin, rounding
                    )
                    if len(unionResult) == 1:
                        # Merged two results
                        safeRemoveObj(results[mergedWith])
                        safeRemoveObj(results[j])
                        results[mergedWith] = unionResult[0]
                        results.pop(j)
                        changed = True
                        break
                    else:
                        for obj in unionResult:
                            safeRemoveObj(obj)

    # Merge all remaining disjoint results into single curve object
    if len(results) > 1:
        merged = mergeSplinesToCurve(results)
        for obj in results:
            safeRemoveObj(obj)
        return [merged]

    return results


def duplicateCurve(curve):
    """Create a duplicate of a curve object."""
    newData = curve.data.copy()
    newObj = bpy.data.objects.new(curve.name + "_copy", newData)
    newObj.matrix_world = curve.matrix_world.copy()
    for coll in curve.users_collection:
        coll.objects.link(newObj)
    return newObj


def mergeSplinesToCurve(curveObjs):
    """Merge splines from multiple curve objects into one."""
    if not curveObjs:
        return None

    base = curveObjs[0]
    resultData = bpy.data.curves.new("BoolResult", "CURVE")
    resultData.dimensions = base.data.dimensions
    resultObj = bpy.data.objects.new("BoolResult", resultData)
    resultObj.matrix_world = base.matrix_world.copy()

    for coll in base.users_collection:
        coll.objects.link(resultObj)

    invMW = resultObj.matrix_world.inverted_safe()

    for curveObj in curveObjs:
        mw = curveObj.matrix_world
        for spline in curveObj.data.splines:
            newSpline = resultData.splines.new("BEZIER")
            bpts = spline.bezier_points
            if len(bpts) > 1:
                newSpline.bezier_points.add(len(bpts) - 1)
            for i, bp in enumerate(bpts):
                newBp = newSpline.bezier_points[i]
                newBp.co = invMW @ (mw @ bp.co)
                newBp.handle_left = invMW @ (mw @ bp.handle_left)
                newBp.handle_right = invMW @ (mw @ bp.handle_right)
                newBp.handle_left_type = "FREE"
                newBp.handle_right_type = "FREE"
            newSpline.use_cyclic_u = spline.use_cyclic_u

    return resultObj


def booleanCurvesPair(curveA, curveB, operation, margin, rounding):
    """
    Perform boolean operation on two closed bezier curves.
    curveA is the base curve, curveB is applied to it.
    For DIFFERENCE: result = curveA - curveB (curveB is subtracted from curveA)
    Returns new curve object(s).
    """
    # Work on copies to avoid modifying originals
    workA = duplicateCurve(curveA)
    workB = duplicateCurve(curveB)
    curves = [workA, workB]

    # Get intersections and cut both curves
    allIntersectCos, intersectMap = getCurveIntersectPts(
        curves, False, margin, rounding
    )

    if len(allIntersectCos) < 2:
        # No valid intersections - handle special cases
        midA = getSegmentMidpoint(workA, 0, 0)
        midB = getSegmentMidpoint(workB, 0, 0)
        aInsideB = isPointInsideCurve(midA, workB, margin)
        bInsideA = isPointInsideCurve(midB, workA, margin)

        # Clean up work copies we won't use
        if operation == "UNION":
            if aInsideB:
                safeRemoveObj(workA)
                return [workB]  # A inside B, keep B
            elif bInsideA:
                safeRemoveObj(workB)
                return [workA]  # B inside A, keep A
            else:
                return [workA, workB]  # Disjoint, keep both
        elif operation == "INTERSECTION":
            if aInsideB:
                safeRemoveObj(workB)
                return [workA]  # A inside B, keep A
            elif bInsideA:
                safeRemoveObj(workA)
                return [workB]  # B inside A, keep B
            else:
                safeRemoveObj(workA)
                safeRemoveObj(workB)
                return []  # Disjoint, no intersection
        else:  # DIFFERENCE: A - B
            if aInsideB:
                safeRemoveObj(workA)
                safeRemoveObj(workB)
                return []  # A completely inside B, nothing left
            elif bInsideA:
                return [workA, workB]  # B inside A, A with hole
            else:
                safeRemoveObj(workB)
                return [workA]  # Disjoint, keep A unchanged

    # Insert points at intersections on both curves
    mapKeys = sorted(intersectMap.keys(), key=lambda x: (x[0], x[1], x[2]))
    prevCnt = {0: 0, 1: 0}
    prevSplineIdx = {0: None, 1: None}

    for key in mapKeys:
        intersectPts = intersectMap[key]
        curveIdx, splineIdx, segIdx = key
        curve = curves[curveIdx]

        if prevSplineIdx[curveIdx] == splineIdx:
            segIdx += prevCnt[curveIdx]
        else:
            prevCnt[curveIdx] = 0

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

        insertCos = [
            mw.inverted_safe() @ co
            for co in sortedCos
            if not vectCmpWithMargin(co, seg[0], margin)
            and not vectCmpWithMargin(co, seg[3], margin)
        ]

        if insertCos:
            insertBezierPts(curve, splineIdx, segIdx, insertCos, "FREE", margin)
            prevCnt[curveIdx] += len(insertCos)

        prevSplineIdx[curveIdx] = splineIdx

    # Build result based on operation
    result = buildBooleanResult(workA, workB, operation, allIntersectCos, margin)

    # Clean up work copies
    safeRemoveObj(workA)
    safeRemoveObj(workB)

    return result


def buildBooleanResult(curveA, curveB, operation, intersectCos, margin):
    """Build the result curve from boolean operation."""
    collections = curveA.users_collection

    # Create result curve
    resultData = bpy.data.curves.new("BoolResult", "CURVE")
    resultData.dimensions = curveA.data.dimensions
    resultObj = bpy.data.objects.new("BoolResult", resultData)

    for coll in collections:
        coll.objects.link(resultObj)

    resultObj.matrix_world = curveA.matrix_world.copy()
    invMW = resultObj.matrix_world.inverted_safe()

    # Collect segments from both curves, marking inside/outside
    segmentsA = collectCurveSegments(curveA, curveB, margin)
    segmentsB = collectCurveSegments(curveB, curveA, margin)

    # Select segments based on operation
    if operation == "UNION":
        # Keep segments of A outside B, and segments of B outside A
        keepA = [s for s in segmentsA if not s["inside"]]
        keepB = [s for s in segmentsB if not s["inside"]]
    elif operation == "INTERSECTION":
        # Keep segments of A inside B, and segments of B inside A
        keepA = [s for s in segmentsA if s["inside"]]
        keepB = [s for s in segmentsB if s["inside"]]
    else:  # DIFFERENCE
        # Keep segments of A outside B, and segments of B inside A (reversed)
        keepA = [s for s in segmentsA if not s["inside"]]
        keepB = [s for s in segmentsB if s["inside"]]
        for s in keepB:
            s["reversed"] = True

    # Build splines from kept segments
    allSegs = keepA + keepB
    if allSegs:
        buildSplinesFromSegments(resultData, allSegs, invMW, margin)

    if len(resultData.splines) == 0:
        safeRemoveObj(resultObj)
        return []

    return [resultObj]


def collectCurveSegments(curve, otherCurve, margin):
    """Collect all segments from curve with inside/outside info relative to otherCurve."""
    segments = []
    mw = curve.matrix_world

    for splineIdx, spline in enumerate(curve.data.splines):
        bpts = spline.bezier_points
        numSegs = len(bpts) if spline.use_cyclic_u else len(bpts) - 1

        for segIdx in range(numSegs):
            nextIdx = (segIdx + 1) % len(bpts)

            # Check multiple points along segment for robust inside detection
            seg = [
                mw @ bpts[segIdx].co,
                mw @ bpts[segIdx].handle_right,
                mw @ bpts[nextIdx].handle_left,
                mw @ bpts[nextIdx].co,
            ]

            # Sample at 0.25, 0.5, 0.75 and use majority vote
            insideCount = 0
            for t in [0.25, 0.5, 0.75]:
                pt = getPtFromT(seg[0], seg[1], seg[2], seg[3], t)
                if isPointInsideCurve(pt, otherCurve, margin):
                    insideCount += 1
            inside = insideCount >= 2  # Majority vote

            segments.append(
                {
                    "curve": curve,
                    "splineIdx": splineIdx,
                    "segIdx": segIdx,
                    "start": seg[0],
                    "hr": seg[1],
                    "hl": seg[2],
                    "end": seg[3],
                    "inside": inside,
                    "reversed": False,
                }
            )

    return segments


def getSegEndpoints(seg):
    """Get effective start/end points and handles for a segment, accounting for reversal."""
    if seg["reversed"]:
        return seg["end"], seg["hl"], seg["hr"], seg["start"]
    return seg["start"], seg["hr"], seg["hl"], seg["end"]


def buildSplinesFromSegments(curveData, segments, invMW, margin):
    """Build splines from collected segments, attempting to join connected ones."""
    if not segments:
        return

    used = [False] * len(segments)

    for i, seg in enumerate(segments):
        if used[i]:
            continue

        # Start new spline
        spline = curveData.splines.new("BEZIER")
        chain = [seg]
        used[i] = True

        # Try to extend chain in both directions
        changed = True
        while changed:
            changed = False
            _, _, _, chainEnd = getSegEndpoints(chain[-1])
            chainStart, _, _, _ = getSegEndpoints(chain[0])

            for j, s in enumerate(segments):
                if used[j]:
                    continue
                sStart, _, _, sEnd = getSegEndpoints(s)

                # Try to append to end
                if vectCmpWithMargin(chainEnd, sStart, margin):
                    chain.append(s)
                    used[j] = True
                    changed = True
                    break
                # Try to prepend to beginning
                if vectCmpWithMargin(sEnd, chainStart, margin):
                    chain.insert(0, s)
                    used[j] = True
                    changed = True
                    break

        # Check if closed
        _, _, _, chainEnd = getSegEndpoints(chain[-1])
        chainStart, _, _, _ = getSegEndpoints(chain[0])

        # Use a relaxed margin for closure check to handle precision issues
        # The boolean operation involves multiple matrix transformations which can accumulate error
        closureMargin = margin * 10
        isClosed = vectCmpWithMargin(chainEnd, chainStart, closureMargin)

        # Build spline points - each segment contributes its start point
        for idx, s in enumerate(chain):
            if idx > 0:
                spline.bezier_points.add(1)
            bp = spline.bezier_points[-1]
            bp.handle_left_type = "FREE"
            bp.handle_right_type = "FREE"

            start, hrOut, hlIn, end = getSegEndpoints(s)
            bp.co = invMW @ start
            bp.handle_right = invMW @ hrOut
            if idx == 0:
                bp.handle_left = invMW @ start  # placeholder for first point

        # Set handle_left for each point from previous segment's incoming handle
        for idx in range(1, len(spline.bezier_points)):
            _, _, hlIn, _ = getSegEndpoints(chain[idx - 1])
            spline.bezier_points[idx].handle_left = invMW @ hlIn

        # Handle closing or final point
        if isClosed:
            # Set first point's handle_left from last segment's incoming handle
            _, _, hlIn, _ = getSegEndpoints(chain[-1])
            spline.bezier_points[0].handle_left = invMW @ hlIn
        else:
            # Add final endpoint
            _, _, hlIn, end = getSegEndpoints(chain[-1])
            spline.bezier_points.add(1)
            bp = spline.bezier_points[-1]
            bp.handle_left_type = "FREE"
            bp.handle_right_type = "FREE"
            bp.co = invMW @ end
            bp.handle_left = invMW @ hlIn
            bp.handle_right = invMW @ end

        spline.use_cyclic_u = isClosed
