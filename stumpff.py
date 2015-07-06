import numpy as np

π = np.pi
inf = np.inf
nan = np.nan

sqrt = np.sqrt
sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin
arccos = np.arccos
arctan = np.arctan
arctan2 = np.arctan2
sinh = np.sinh
cosh = np.cosh
arctanh = np.arctanh

def stumpff_s(z):
    if hasattr(z,'__iter__'):
        lt0 = z < 0
        gt0 = z > 0
        eq0 = z == 0

        zlt0 = z[lt0]
        zgt0 = z[gt0]

        sqrt_neg_zlt0 = sqrt(-zlt0)
        sqrt_zgt0 = sqrt(zgt0)

        res = np.empty(z.shape)
        res[lt0] = (sinh(sqrt_neg_zlt0) - sqrt_neg_zlt0) / sqrt_neg_zlt0**3
        res[gt0] = (sqrt_zgt0 - sin(sqrt_zgt0)) / sqrt_zgt0**3
        res[eq0] = 1/6

        return res

    else:
        if z < 0:
            sqrt_neg_z = sqrt(-z)
            return (sinh(sqrt_neg_z) - sqrt_neg_z) / sqrt_neg_z**3
        elif z > 0:
            sqrt_z = sqrt(z)
            return (sqrt_z - sin(sqrt_z)) / sqrt_z**3
        else:
            return 1/6

def stumpff_c(z):
    if hasattr(z,'__iter__'):
        lt0 = z < 0
        gt0 = z > 0
        eq0 = z == 0

        zlt0 = z[lt0]
        zgt0 = z[gt0]

        res = np.empty(z.shape)
        res[lt0] = (1 - cosh(sqrt(-zlt0))) / zlt0
        res[gt0] = (1 - cos(sqrt(zgt0))) / zgt0
        res[eq0] = 1/2

        return res

    else:
        if z < 0:
            return (1 - cosh(sqrt(-z))) / z
        elif z > 0:
            return (1 - cos(sqrt(z))) / z
        else:
            return 1/2

