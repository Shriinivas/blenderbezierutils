# bezier_utils/utils/math_utils.py

from ..constants import DEF_ERR_MARGIN
from mathutils import Vector
from math import sqrt


def floatCmpWithMargin(float1, float2, margin=DEF_ERR_MARGIN):
    return abs(float1 - float2) < margin


def vectCmpWithMargin(v1, v2, margin=DEF_ERR_MARGIN):
    return all(floatCmpWithMargin(v1[i], v2[i], margin) for i in range(0, len(v1)))

def toHexStr(rgba):
    ch = []
    for c in rgba[:3]:
        if c < 0.0031308:
            cc = 0.0 if c < 0.0 else c * 12.92
        else:
            cc = 1.055 * pow(c, 1.0 / 2.4) - 0.055
        ch.append(hex(max(min(int(cc * 255 + 0.5), 255), 0))[2:])
    return "".join(ch), str(rgba[-1])



def getBBoxCenter(bbox):  # bbox -> [leftBotFront, rightTopBack]
    return Vector(((bbox[0][i] + bbox[1][i]) / 2 for i in range(3)))
