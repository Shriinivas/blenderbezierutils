# bezier_utils/utils/bezier_math.py

from mathutils import Vector
from math import sqrt, sin, cos, acos, tan, pi
from ..constants import DEF_ERR_MARGIN, LARGE_NO
from .math_utils import floatCmpWithMargin, vectCmpWithMargin, getBBoxCenter

def getPtFromT(p0, p1, p2, p3, t):
    c = (1 - t)
    pt = (c ** 3) * p0 + 3 * (c ** 2) * t * p1 + \
        3 * c * (t ** 2) * p2 + (t ** 3) * p3
    return pt

def getTangentAtT(p0, p1, p2, p3, t):
    c = (1 - t)
    tangent = -3 * (c * c) * p0 + 3 * c * c * p1 - 6 * t * c * p1 - \
        3 * t * t * p2 + 6 * t * c * p2 + 3 * t * t * p3
    return tangent

# iterative brute force, not optimized, some iterations maybe redundant
def getTsForPt(p0, p1, p2, p3, co, coIdx, tolerance = 0.000001, maxItr = 1000):
    ts = set()
    # check t from start to end and end to start
    for T in [1., 0.]:
        # check clockwise as well as anticlockwise
        for dirn in [1, -1]:
            t = T
            t2 = 1
            rhs = getPtFromT(p0, p1, p2, p3, t)[coIdx]
            error = rhs - co
            i = 0

            while(abs(error) > tolerance and i < maxItr):
                t2 /= 2
                if(dirn * error < 0):
                    t += t2
                else:
                    t -= t2
                rhs = getPtFromT(p0, p1, p2, p3, t)[coIdx]
                error = rhs - co

                i += 1

            if(i < maxItr and t >= 0 and t <= 1):
                ts.add(round(t, 3))
    return ts

#TODO: There may be a more efficient approach, but this seems foolproof
def getTForPt(curve, testPt, tolerance = .000001):
    minLen = LARGE_NO
    retT = None
    for coIdx in range(0, 3):
        ts = getTsForPt(curve[0], curve[1], curve[2], curve[3], \
            testPt[coIdx], coIdx, tolerance)
        for t in ts:
            pt = getPtFromT(curve[0], curve[1], curve[2], curve[3], t)
            pLen = (testPt - pt).length
            if(pLen < minLen):
                minLen = pLen
                retT = t
    return retT

def getLinesFromPts(pts):
    positions = []
    for i, pt in enumerate(pts):
        positions.append(pt)
        if(i > 0 and i < (len(pts)-1)):
            positions.append(pt)
    return positions

#see https://stackoverflow.com/questions/878862/drawing-part-of-a-b%c3%a9zier-curve-by-reusing-a-basic-b%c3%a9zier-curve-function/879213#879213
def getPartialSeg(seg, t0, t1):
    pts = [seg[0], seg[1], seg[2], seg[3]]

    if(t0 > t1):
        tt = t1
        t1 = t0
        t0 = tt

    u0 = 1.0 - t0
    u1 = 1.0 - t1

    qa = [pts[0][i]*u0*u0 + pts[1][i]*2*t0*u0 + pts[2][i]*t0*t0 for i in range(0, 3)]
    qb = [pts[0][i]*u1*u1 + pts[1][i]*2*t1*u1 + pts[2][i]*t1*t1 for i in range(0, 3)]
    qc = [pts[1][i]*u0*u0 + pts[2][i]*2*t0*u0 + pts[3][i]*t0*t0 for i in range(0, 3)]
    qd = [pts[1][i]*u1*u1 + pts[2][i]*2*t1*u1 + pts[3][i]*t1*t1 for i in range(0, 3)]

    pta = Vector([qa[i]*u0 + qc[i]*t0 for i in range(0, 3)])
    ptb = Vector([qa[i]*u1 + qc[i]*t1 for i in range(0, 3)])
    ptc = Vector([qb[i]*u0 + qd[i]*t0 for i in range(0, 3)])
    ptd = Vector([qb[i]*u1 + qd[i]*t1 for i in range(0, 3)])

    return [pta, ptb, ptc, ptd]

def getSegLen(pts, error = DEF_ERR_MARGIN, start = None, end = None, t1 = 0, t2 = 1):
    if(start is None): start = pts[0]
    if(end is None): end = pts[-1]

    t1_5 = (t1 + t2)/2
    mid = getPtFromT(*pts, t1_5)
    l = (end - start).length
    l2 = (mid - start).length + (end - mid).length
    if (l2 - l > error):
        return (getSegLen(pts, error, start, mid, t1, t1_5) +
                getSegLen(pts, error, mid, end, t1_5, t2))
    return l2

def getInterpolatedVertsCo(curvePts, numDivs):
    # Can be calculated only once
    curveLength = sum((curvePts[i] - curvePts[i-1]).length
        for i in range(1, len(curvePts)))

    if(floatCmpWithMargin(curveLength, 0)):
        return [curvePts[0]] * numDivs

    segLen = curveLength / numDivs
    vertCos = [curvePts[0]]

    actualLen = 0
    vertIdx = 0

    for i in range(1, numDivs):
        co = None
        targetLen = i * segLen

        while(not floatCmpWithMargin(actualLen, targetLen)
            and actualLen < targetLen):

            vertCo = curvePts[vertIdx]
            vertIdx += 1
            nextVertCo = curvePts[vertIdx]
            actualLen += (nextVertCo - vertCo).length

        if(floatCmpWithMargin(actualLen, targetLen)):
            co = curvePts[vertIdx]

        else:   #interpolate
            diff = actualLen - targetLen
            co = (nextVertCo - (nextVertCo - vertCo) * \
                (diff/(nextVertCo - vertCo).length))

            #Revert to last pt
            vertIdx -= 1
            actualLen -= (nextVertCo - vertCo).length
        vertCos.append(co)

    # ~ if(not vectCmpWithMargin(curvePts[0], curvePts[-1])):
    vertCos.append(curvePts[-1])

    return vertCos

# Arc conversion helpers (from SVG arc to bezier)
TAU = pi * 2

def unit_vector_angle(ux, uy, vx, vy):
    if(ux * vy - uy * vx < 0):
        sign = -1
    else:
        sign = 1
    
    dot  = ux * vx + uy * vy
    
    if (round(dot, 3) >=  1.0):
        dot =  1.0
    
    if (round(dot, 3) <= -1.0):
        dot = -1.0
    
    return sign * acos(dot)

def get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi):
    x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2
    y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2
    
    rx_sq  =  rx * rx
    ry_sq  =  ry * ry
    x1p_sq = x1p * x1p
    y1p_sq = y1p * y1p
    
    radicant = (rx_sq * ry_sq) - (rx_sq * y1p_sq) - (ry_sq * x1p_sq)
    
    if (radicant < 0):
        radicant = 0
    
    radicant /=   (rx_sq * y1p_sq) + (ry_sq * x1p_sq)
    factor = 1
    if(fa == fs):
        factor = -1
    radicant = sqrt(radicant) * factor
    
    cxp = radicant *  rx/ry * y1p
    cyp = radicant * -ry/rx * x1p
    
    cx = cos_phi*cxp - sin_phi*cyp + (x1+x2)/2
    cy = sin_phi*cxp + cos_phi*cyp + (y1+y2)/2
    
    v1x =  (x1p - cxp) / rx
    v1y =  (y1p - cyp) / ry
    v2x = (-x1p - cxp) / rx
    v2y = (-y1p - cyp) / ry
    
    theta1 = unit_vector_angle(1, 0, v1x, v1y)
    delta_theta = unit_vector_angle(v1x, v1y, v2x, v2y)
    
    if (fs == 0 and delta_theta > 0):
        delta_theta -= TAU
    
    if (fs == 1 and delta_theta < 0):
        delta_theta += TAU
    
    return [ cx, cy, theta1, delta_theta ]

def approximate_unit_arc(theta1, delta_theta):
    alpha = 4.0/3 * tan(delta_theta/4)
    
    x1 = cos(theta1)
    y1 = sin(theta1)
    x2 = cos(theta1 + delta_theta)
    y2 = sin(theta1 + delta_theta)
    
    return [ x1, y1, x1 - y1*alpha, y1 + x1*alpha, x2 + y2*alpha, y2 - x2*alpha, x2, y2 ]

def getInterpBezierPts(segPts, subdivPerUnit, segLens = None, maxRes = None):
    from mathutils import geometry
    
    if(len(segPts) < 2):
        return []
    
    curvePts = []
    for i in range(1, len(segPts)):
        seg = [segPts[i-1][1], segPts[i-1][2], segPts[i][0], segPts[i][1]]
        if(segLens is not None and len(segLens) > (i-1)):
            res = int(segLens[i-1] * subdivPerUnit)
        else:
            res = int(getSegLen(seg) * subdivPerUnit)
        if(res < 2): res = 2
        if(maxRes is not None and res > maxRes): res = maxRes
        curvePts += geometry.interpolate_bezier(*seg, res)
    
    return curvePts

def getMappedList(result, rx, ry, sin_phi, cos_phi, cc):
    mappedList = []
    for elem in result:
        curve = []
        for i in range(0, len(elem), 2):
            x = elem[i + 0]
            y = elem[i + 1]
            
            x *= rx
            y *= ry
            
            xp = cos_phi*x - sin_phi*y
            yp = sin_phi*x + cos_phi*y
            
            elem[i + 0] = xp + cc[0]
            elem[i + 1] = yp + cc[1]
            curve.append(complex(elem[i + 0], elem[i + 1]))
        mappedList.append(curve)
    return mappedList

def a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi, noSegs):
    sin_phi = sin(phi * TAU / 360)
    cos_phi = cos(phi * TAU / 360)
    
    x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2
    y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2
    
    if (x1p == 0 and y1p == 0):
        return []
    
    if (rx == 0 or ry == 0):
        return []
    
    rx = abs(rx)
    ry = abs(ry)
    
    lmbd = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry)
    if (lmbd > 1):
        rx *= sqrt(lmbd)
        ry *= sqrt(lmbd)
    
    cc = get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi)
    
    result = []
    theta1 = cc[2]
    delta_theta = cc[3]
    
    segments = noSegs
    delta_theta /= segments
    
    for i in range(0, segments):
        result.append(approximate_unit_arc(theta1, delta_theta))
        theta1 += delta_theta
    
    return getMappedList(result, rx, ry, sin_phi, cos_phi, cc)

# 3D conversion helpers
def get3DVector(cmplx, axisIdxs, z):
    vElems = [None] * 3
    vElems[axisIdxs[0]] = cmplx.real
    vElems[axisIdxs[1]] = cmplx.imag
    vElems[axisIdxs[2]] = z
    return Vector(vElems)

def getSegsForArc(start, radius, sweep, end, noSegs, axisIdxs, z):
    x1, y1 = start.real, start.imag
    x2, y2 = end.real, end.imag
    fa = 0
    fs = sweep
    rx, ry = radius.real, radius.imag
    phi = 0
    curvesPts = a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi, noSegs)
    newSegs = []
    for curvePts in curvesPts:
        newSegs.append([get3DVector(curvePts[0], axisIdxs, z), get3DVector(curvePts[1], axisIdxs, z), \
            get3DVector(curvePts[2], axisIdxs, z), get3DVector(curvePts[3], axisIdxs, z)])
    
    return newSegs

def getWSDataForSegs(segs):
    from .math_utils import vectCmpWithMargin
    
    prevSeg = None
    wsData = []
    
    for j, seg in enumerate(segs):
        pt = seg[0]
        handleRight = seg[1]
        
        if(j == 0): handleLeft = pt
        else: handleLeft = prevSeg[2]
        
        ht = 'ALIGNED' if(vectCmpWithMargin(pt - handleLeft, handleRight - pt)) else 'FREE'
        wsData.append([handleLeft, pt, handleRight, ht, ht])
        prevSeg = seg
    
    if(prevSeg is not None): wsData.append([seg[2], seg[3], seg[3], 'FREE', 'FREE'])
    else: return []
    
    return wsData

def isStraightSeg(segPts):
    if(len(segPts) != 2): return False
    if((len(segPts[0]) == 5 or len(segPts[1]) == 5) and \
        segPts[0][4] == 'VECTOR' and segPts[1][3] == 'VECTOR'):
            return True
    if(vectCmpWithMargin((segPts[0][2]-segPts[0][1]).normalized(), \
        (segPts[1][1] - segPts[1][0]).normalized())):
        return True
    return False

def hasAlignedHandles(pt):
    if(len(pt) == 5 and 'ALIGNED' in {pt[3], pt[4]} and 'FREE' not in {pt[3], pt[4]}):
        return True
    diffV1 = pt[1] - pt[0]
    diffV2 = pt[2] - pt[1]
    if(vectCmpWithMargin(diffV1.normalized(), diffV2.normalized())):
        return True
    return False

# Used in functions where actual locs of pts on curve matter (like subdiv Bezier)
# (... kind of expensive)
def getPtsAlongBezier3D(segPts, rv3d, curveRes, minRes = 200):

    viewDist = rv3d.view_distance

    # The smaller the view dist (higher zoom level),
    # the higher the num of subdivisions
    curveRes = curveRes / viewDist

    if(curveRes < minRes): curveRes = minRes

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes)

# Used in functions where only visual resolution of curve matters (like draw Bezier)
# (... not so expensive)
# TODO: Calculate maxRes dynamically
def getPtsAlongBezier2D(segPts, areaRegionInfo, curveRes, getCoordFromLoc_func, maxRes = None):
    segLens = []
    for i in range(1, len(segPts)):
        seg = [segPts[i-1][1], segPts[i-1][2], segPts[i][0], segPts[i][1]]

        #TODO: A more optimized solution... (Called very frequently)
        segLen = 0
        for info in areaRegionInfo:
            seg2D = [getCoordFromLoc_func(info[1], info[2], loc) for loc in seg]
            sl = getSegLen(seg2D)
            if(sl > segLen):
                segLen = sl
        segLens.append(segLen)

    return getInterpBezierPts(segPts, subdivPerUnit = curveRes, \
        segLens = segLens, maxRes = maxRes)
def getIntersectPts(seg0, seg1, soln, solnRounded, recurs, margin, rounding, \
    maxRecurs = 100):
    overlap, bbox0, bbox1 = getBBoxOverlapInfo(seg0, seg1)
    if(overlap):
        if(recurs == maxRecurs):
            print('Maximum recursions in getIntersectPts!')
            return False

        if(all(abs(bbox0[0][i] - bbox0[1][i]) < margin for i in range(3))):
            center = getBBoxCenter(bbox0)
            roundedVect = Vector([round(x, rounding) for x in center]).freeze()
            if(roundedVect not in solnRounded):
                soln.append(center)
                solnRounded.add(roundedVect)
            return True
        elif(all(abs(bbox1[0][i] - bbox1[1][i]) < margin for i in range(3))):
            center = getBBoxCenter(bbox1)
            roundedVect = Vector([round(x, rounding) for x in center]).freeze()
            if(roundedVect not in solnRounded):
                soln.append(center)
                solnRounded.add(roundedVect)
            return True
        else:
            seg01 = getPartialSeg(seg0, t0 = 0, t1 = 0.5)
            seg02 = getPartialSeg(seg0, t0 = 0.5, t1 = 1)

            seg11 = getPartialSeg(seg1, t0 = 0, t1 = 0.5)
            seg12 = getPartialSeg(seg1, t0 = 0.5, t1 = 1)

            r0 = getIntersectPts(seg01, seg11, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r1 = getIntersectPts(seg01, seg12, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r2 = getIntersectPts(seg02, seg11, soln, solnRounded, \
                recurs + 1, margin, rounding)
            r3 = getIntersectPts(seg02, seg12, soln, solnRounded, \
                recurs + 1, margin, rounding)

            return any((r0, r1, r2, r3))
    return False

# https://stackoverflow.com/questions/24809978/calculating-the-bounding-box-of-cubic-bezier-curve
#(3 D - 9 C + 9 B - 3 A) t^2 + (6 A - 12 B + 6 C) t + 3 (B - A)
def getBBox(seg):
    A = seg[0]
    B = seg[1]
    C = seg[2]
    D = seg[3]

    leftBotFront = Vector([min([A[i], D[i]]) for i in range(0, 3)])
    rgtTopBack = Vector([max([A[i], D[i]]) for i in range(0, 3)])

    a = [3 * D[i] - 9 * C[i] + 9 * B[i] - 3 * A[i] for i in range(0, 3)]
    b = [6 * A[i] - 12 * B[i] + 6 * C[i] for i in range(0, 3)]
    c = [3 * (B[i] - A[i]) for i in range(0, 3)]

    solnsxyz = []
    for i in range(0, 3):
        solns = []
        if(a[i] == 0):
            if(b[i] == 0):
                solns.append(0)#Independent of t so lets take the starting pt
            else:
                solns.append(c[i] / b[i])
        else:
            rootFact = b[i] * b[i] - 4 * a[i] * c[i]
            if(rootFact >=0 ):
                #Two solutions with + and - sqrt
                solns.append((-b[i] + sqrt(rootFact)) / (2 * a[i]))
                solns.append((-b[i] - sqrt(rootFact)) / (2 * a[i]))
        solnsxyz.append(solns)

    for i, soln in enumerate(solnsxyz):
        for j, t in enumerate(soln):
            if(t <= 1 and t >= 0):
                co = getPtFromT(A[i], B[i], C[i], D[i], t)
                if(co < leftBotFront[i]): leftBotFront[i] = co
                if(co > rgtTopBack[i]): rgtTopBack[i] = co

    return leftBotFront, rgtTopBack

def getBBoxOverlapInfo(seg0, seg1):
    bbox0 = getBBox(seg0)
    max0 = [max(bbox0[i][axis] for i in range(2)) for axis in range(3)]
    min0 = [min(bbox0[i][axis] for i in range(2)) for axis in range(3)]
    bbox1 = getBBox(seg1)

    overlap = True
    if any(all(bbox1[i][j] < min0[j] for i in range(2)) for j in range(3)) or any(
        all(bbox1[i][j] > max0[j] for i in range(2)) for j in range(3)
    ):
        overlap = False

    return overlap, bbox0, bbox1
