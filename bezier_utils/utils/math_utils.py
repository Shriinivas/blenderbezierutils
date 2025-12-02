# bezier_utils/utils/math_utils.py

from ..constants import DEF_ERR_MARGIN
from mathutils import Vector


def floatCmpWithMargin(float1, float2, margin=DEF_ERR_MARGIN):
    return abs(float1 - float2) < margin


def vectCmpWithMargin(v1, v2, margin=DEF_ERR_MARGIN):
    return all(floatCmpWithMargin(v1[i], v2[i], margin) for i in range(0, len(v1)))


def getBBox(seg):
    s = seg
    bbox = [s[1], s[1], s[1], s[1], s[1], s[1]]
    for p in seg:
        for i, co in enumerate(p):
            if co < bbox[i * 2]:
                bbox[i * 2] = co
            if co > bbox[i * 2 + 1]:
                bbox[i * 2 + 1] = co
    return bbox


def toHexStr(rgba):
    ch = []
    for c in rgba[:3]:
        if c < 0.0031308:
            cc = 0.0 if c < 0.0 else c * 12.92
        else:
            cc = 1.055 * pow(c, 1.0 / 2.4) - 0.055
        ch.append(hex(max(min(int(cc * 255 + 0.5), 255), 0))[2:])
    return "".join(ch), str(rgba[-1])


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


def getBBoxCenter(bbox):  # bbox -> [leftBotFront, rightTopBack]
    return Vector(((bbox[0][i] + bbox[1][i]) / 2 for i in range(3)))
