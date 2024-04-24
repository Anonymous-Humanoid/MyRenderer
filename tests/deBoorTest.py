from typing import Optional

from myRenderer import KnotVector, PointLike, PointsLike


# Credit to "The NURBS Book", page 75
# PDF page 41 of: https://argos.vu/wp-content/uploads/2019/03/Springer-TheNURBSBook.pdf
# N_(i,p)[U](u) or N_(i,p)(u) or N_i(u)
def bSplineBasis(i: int,
                 p: int,
                 U: KnotVector,
                 u: float) \
                        -> float:
    m = len(U) - 1
    
    if ((i == 0 and u == U[0]) or
            (i == m - p - 1 and u == U[m])):
        return 1
    elif (u < U[i] or u >= U[i + p + 1]):
        return 0
    
    N = {}
    
    for j in range(p + 1):
        if (u >= U[i + j] and u < U[i + j + 1]):
            N[j] = 1
        else:
            N[j] = 0
    
    for k in range(1, p + 1):
        if N[0] == 0:
            saved = 0
        else:
            saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i])
        
        for j in range(p - k + 1):
            leftKnot, rightKnot = U[i + j + 1], U[i + j + k + 1]
            
            if N[j + 1] == 0:
                N[j], saved = saved, 0
            else:
                temp = N[j + 1] / (rightKnot - leftKnot)
                N[j] = saved + (rightKnot - u) * temp
                saved = (u - leftKnot) * temp
    
    return N[0]

# Doesn't implement the Cox-De Boor algorithm (i.e, is inefficient)
# B(u)
def bSpline(controlPts: PointsLike,
            U: KnotVector,
            u: float) \
                    -> PointLike:
    # Extra optimization for special cases
    if u == U[0]:
        return controlPts[0]
    elif u == U[-1]:
        return controlPts[-1]
    
    dim = len(controlPts[0])
    
    if u < U[0] or u > U[-1]:
        return tuple(0 for _ in range(dim))
    
    m, n = len(U) - 1, len(controlPts) - 1
    p = m - n - 1
    
    Ns = [bSplineBasis(i, p, U, u) for i in range(n + 1)]
    
    return tuple(sum(N * c[i] for c, N in zip(controlPts, Ns)) for i in range(dim))

# Implements the Cox-De Boor algorithm
# N_(i,p)[U](u) or N_(i,p)(u) or N_i(u)
def _bSplineBasis(i: int,
                  p: int,
                  U: KnotVector,
                  u: float) \
                        -> float:
    if p == 0:
        return 1 if U[i] <= u < U[i + 1] else 0
    
    denom1 = U[i + p] - U[i]
    denom2 = U[i + p + 1] - U[i + 1]

    result = 0
    if denom1 != 0:
        result += (u - U[i]) / denom1 * _bSplineBasis(i, p - 1, U, u)
    if denom2 != 0:
        result += (U[i + p + 1] - u) / denom2 * _bSplineBasis(i + 1, p - 1, U, u)

    return result

# B(u)
def _bSpline(controlPts: PointsLike,
             U: KnotVector,
             u: float,
             maxDegree: Optional[int] = None) \
                    -> PointLike:
    # Extra optimization for special cases
    if u == U[0]:
        return tuple(c for c in controlPts[0])
    elif u == U[-1]:
        return tuple(c for c in controlPts[-1])
    
    dim = len(controlPts[0])
    
    if u < U[0] or u > U[-1]:
        return tuple(0 for _ in range(dim))
    
    m, n = len(U) - 1, len(controlPts) - 1
    tempP = m - n - 1
    p = tempP if maxDegree is None else min(maxDegree, tempP)

    Ns = [_bSplineBasis(i, p, U, u) for i in range(n + 1)]
    
    return tuple(sum(N * c[i] for c, N in zip(controlPts, Ns)) for i in range(dim))

# Example usage:
if __name__ == '__main__':
    U = KnotVector({0: 3, 1: 1, 2: 1, 3: 1, 4: 3})
    # U.normalizeAll()

    maxDegree = 3
    controlPts = (
        (0, 0),
        (1, 2),
        (3, 4),
        (5, 0)
    )

    u = 3

    print(_bSpline(controlPts, U, u, maxDegree))
    print(bSpline(controlPts, U, u))
