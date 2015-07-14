import numpy as np

sqrt = np.sqrt
sin = np.sin
cos = np.cos
sinh = np.sinh
cosh = np.cosh

def stumpff_c(z):
    if hasattr(z,'__iter__'):
        gt0 = z > 0
        lt0 = z < 0
        eq0 = z == 0

        zgt0 = z[gt0]
        zlt0 = z[lt0]

        res = np.empty(z.shape)
        res[gt0] = (1 - cos(sqrt(zgt0))) / zgt0
        res[lt0] = (1 - cosh(sqrt(-zlt0))) / zlt0
        res[eq0] = 1/2

        return res

    else:
        if z > 0:
            return (1 - cos(sqrt(z))) / z
        elif z < 0:
            return (1 - cosh(sqrt(-z))) / z
        else:
            return 1/2


def stumpff_s(z):
    if hasattr(z,'__iter__'):
        gt0 = z > 0
        lt0 = z < 0
        eq0 = z == 0

        zgt0 = z[gt0]
        zlt0 = z[lt0]

        sqrt_zgt0 = sqrt(zgt0)
        sqrt_neg_zlt0 = sqrt(-zlt0)

        res = np.empty(z.shape)
        res[gt0] = (sqrt_zgt0 - sin(sqrt_zgt0)) / sqrt_zgt0**3
        res[lt0] = (sinh(sqrt_neg_zlt0) - sqrt_neg_zlt0) / sqrt_neg_zlt0**3
        res[eq0] = 1/6

        return res

    else:
        if z > 0:
            sqrt_z = sqrt(z)
            return (sqrt_z - sin(sqrt_z)) / sqrt_z**3
        elif z < 0:
            sqrt_neg_z = sqrt(-z)
            return (sinh(sqrt_neg_z) - sqrt_neg_z) / sqrt_neg_z**3
        else:
            return 1/6

