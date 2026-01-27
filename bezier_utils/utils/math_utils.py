# bezier_utils/utils/math_utils.py

from ..constants import DEF_ERR_MARGIN
from mathutils import Vector, Matrix
from math import sqrt
import numpy as np


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


def get_best_fit_matrix(coords):
    if not coords:
        return Matrix.Identity(4)
        
    pts = np.array(coords)
    centroid = Vector(np.mean(pts, axis=0))
    centered = pts - np.mean(pts, axis=0)
    
    # Check if we have enough variation to compute SVD
    if np.allclose(centered, 0):
        return Matrix.Translation(centroid)

    try:
        _, _, vh = np.linalg.svd(centered)
    except np.linalg.LinAlgError:
        return Matrix.Translation(centroid)

    # 1. Get the Normal (last row of Vh corresponds to smallest singular value)
    normal = Vector(vh[2]).normalized()
    if normal.dot(Vector((0, 0, 1))) < 0:
        normal *= -1

    # 2. Stable Axes (Align Local X with World X)
    world_x = Vector((1, 0, 0))
    if abs(normal.dot(world_x)) > 0.9:
        world_x = Vector((0, 1, 0))

    local_x = (world_x - normal * world_x.dot(normal)).normalized()
    local_y = normal.cross(local_x).normalized()

    rot_matrix = Matrix((local_x, local_y, normal)).transposed()
    return Matrix.LocRotScale(centroid, rot_matrix, None)
