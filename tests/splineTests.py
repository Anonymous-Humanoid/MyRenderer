'''
Created on Apr. 22, 2023

@author: Matthew
'''

# TODO Annotate individual methods with what test output they generate, if any 
''' Test implementation notes:
- 2D spline tests are built for Desmos (https://www.desmos.com/calculator)
- 3D spline tests are built for GeoGebra (https://www.geogebra.org/calculator)
- Spline tests of generalized dimension are also built for GeoGebra
- Tables of 2D points in generalized dimension function tests are built for Desmos
'''

from math import comb
from random import random, sample  # @UnusedImport
from typing import List, Optional

from myRenderer import SplineException, KnotVector, Point, Point2D, Point3D, PointsLike, Points2DLike, Points3DLike


# TODO TEST Bezier curve
# Doesn't implement De Casteljau's Algorithm: beware of numeric instability!
def bezierCurve(pts: PointsLike,
                _debug: bool = False) \
                        -> List[List[float]]:
    n, tot = len(pts[0]), len(pts) - 1
    
    for pt in pts:
        if len(pt) != n:
            raise SplineException('Points must be of equal dimension')
    
    c = [comb(tot, i) for i in range(tot + 1)]
    out = [[c[i] * pt[j] for i, pt in enumerate(pts)] for j in range(n)]
    
    # Testing code (parametric form)
    if _debug and n < 4: # Can't visualize higher than 3 dimensions
        components = [','.join([str(pt[i]) for pt in pts]) for i in range(n)]
        
        while len(components) < 3:
            components.append(','.join(['0' for _ in range(tot + 1)]))
        
        points = f'({{{"},{".join(components)}}})'
        
        eqs = []
        for coeffs in out:
            eq = []
            for i, coeff in enumerate(coeffs):
                eq.append(f'{coeff}(1-t)^({tot-i})t^({i})')
            eqs.append('+'.join(eq))
        
        while len(eqs) < 3:
            eqs.append('0')
        
        curve = (f'Curve({",".join(eqs)},t,0,1)'
                 .replace('+-', '-')
                 .replace('-+', '-')
                 .replace('(1-t)^(0)', '')
                 .replace('t^(0)', '')
                 .replace('(1-t)^(1)', '(1-t)')
                 .replace('t^(1)', 't'))
        
        print(f'{curve}\n{points}')
    
    return out

# TODO Bezier spline
def bezierSpline(pts: PointsLike,
                 _debug: bool = False):
    pass

# TODO Bezier surface
def bezierSurface(pts: PointsLike,
                  _debug: bool = False):
    pass

# Credit to "The NURBS Book", page 75
# PDF page 41 of: https://argos.vu/wp-content/uploads/2019/03/Springer-TheNURBSBook.pdf
# N_(i,p)[U](u) or N_(i,p)(u) or N_i(u) or N_i^m(u)
def bSplineBasis(i: int,
                 p: int,
                 U: KnotVector,
                 u: float,
                 weights: Optional[List[float]] = None,
                 _debug: bool = False) \
                        -> float: # TODO Use _debug parameter
    m = len(U) - 1
    
    # TODO TESTING
    if i < 0 or i + p + 1 >= len(U):
        return 0
    
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
    
    return N[0] if weights is None else weights[i] * N[0]

# TODO Add weights to B-spline functions

# Doesn't implement the Cox-De Boor algorithm (i.e, is inefficient)
# B(u)
def bSpline(controlPts: PointsLike,
            U: KnotVector,
            u: float,
            _debug: bool = False) \
                    -> Point: # TODO Use _debug parameter
    # Extra optimization for special cases
    if u == U[0]:
        return controlPts[0]
    elif u == U[-1]:
        return controlPts[-1]
    
    m, n = len(U) - 1, len(controlPts) - 1
    p = m - n - 1
    dim = len(controlPts[0])
    
    Ns = [bSplineBasis(i, p, U, u) for i in range(n + 1)]
    
    return tuple(sum(N * c[i] for c, N in zip(controlPts, Ns)) for i in range(dim))

# TODO Implement open, closed and clamped: B-splines and NURBS curves and surfaces, as well as T-splines

# N_(i,j)(u,v)
def bSplineSurfaceBasis(i: int,
                        j: int,
                        pu: int,
                        pv: int,
                        U: KnotVector,
                        V: KnotVector,
                        u: float,
                        v: float,
                        _debug: bool = False) \
                                -> float: # TODO Use _debug parameter
    return bSplineBasis(i, pu, U, u) * bSplineBasis(j, pv, V, v)

# B(u,v)
def bSplineSurface(controlPtArray: List[PointsLike],
                   U: KnotVector,
                   V: KnotVector,
                   u: float,
                   v: float,
                   _debug: bool = False) \
                        -> Point: # TODO Use _debug parameter
    mu, mv = len(U) - 1, len(V) - 1
    rows, cols = len(controlPtArray) - 1, len(controlPtArray[0]) - 1
    pu, pv = mu - cols - 1, mv - rows - 1
    dim = len(controlPtArray[0][0])
    
    Ns = [
        [
            bSplineSurfaceBasis(i, j, pu, pv, U, V, u, v)
            for i in range(rows + 1)
        ]
        for j in range(cols + 1)
    ]
    
    return tuple(
        sum(Ns[i][j] * controlPtArray[i][j][k]
                for i in range(rows + 1)
                for j in range(cols + 1))
                    for k in range(dim))

# R_(i,p)(u)
def nurbsBasis(i: int,
               p: float,
               U: KnotVector,
               u: float,
               weights: Optional[List[float]] = None,
               n: Optional[int] = None,
               _debug: bool = False) \
                    -> float: # TODO Use _debug parameter
    if n is None:
        if weights is None:
            raise SplineException("At most one of n and weights can be None")
        n = len(weights) - 1
    
    Ns = [bSplineBasis(j, p, U, u, weights) for j in range(n + 1)]
    
    return Ns[i] / sum(Ns)

# C(u)
def nurbs(controlPts: PointsLike,
          U: KnotVector,
          u: float,
          weights: Optional[List[float]] = None,
          _debug: bool = False) \
                -> Point: # TODO Use _debug parameter
    n = len(controlPts) - 1
    
    if weights is not None and len(weights) != n + 1:
        raise SplineException(f'NURBS curves require as many weights as control points. Expected {n + 1}, got: {len(weights)}')
    
    # Extra optimization for special cases
    if u == U[0]:
        return controlPts[0]
    elif u == U[-1]:
        return controlPts[-1]
    
    m = len(U) - 1
    p = m - n - 1
    dim = len(controlPts[0])
    
    Ns = [nurbsBasis(i, p, U, u, weights) for i in range(n + 1)]
    
    return tuple(sum(N * c[i] for c, N in zip(controlPts, Ns)) for i in range(dim))

# R_(i,j)(u,v)
def nurbsSurfaceBasis(i: int,
                      j: int,
                      pu: int,
                      pv: int,
                      U: KnotVector,
                      V: KnotVector,
                      u: float,
                      v: float,
                      weights: Optional[List[List[float]]] = None,
                      nu: Optional[int] = None,
                      nv: Optional[int] = None,
                      _debug: bool = False) \
                            -> float: # TODO Use _debug parameter
    if nu is None:
        if weights is None:
            raise SplineException('Only one of nu and weights can be None')
        nu = len(weights) - 1
    if nv is None:
        if weights is None:
            raise SplineException('Only one of nv and weights can be None')
        nv = len(weights[0]) - 1
    
    if weights is None:
        Ns = [
            [bSplineSurfaceBasis(_i, _j, pu, pv, U, V, u, v) for _i in range(nv + 1)]
            for _j in range(nu + 1)
        ]
    else:
        Ns = [
            [weights[_i][_j] * bSplineSurfaceBasis(_i, _j, pu, pv, U, V, u, v) for _i in range(nv + 1)]
            for _j in range(nu + 1)
        ]
    
    return Ns[i][j] / sum(c for cs in Ns for c in cs)

# C(u,v)
def nurbsSurface(controlPtArray: List[PointsLike],
                 U: KnotVector,
                 V: KnotVector,
                 u: float,
                 v: float,
                 weights: Optional[List[List[float]]] = None,
                 _debug: bool = False) \
                        -> Point: # TODO Use _debug parameter
    mu, mv = len(U) - 1, len(V) - 1
    rows, cols = len(controlPtArray) - 1, len(controlPtArray[0]) - 1
    pu, pv = mu - cols - 1, mv - rows - 1
    dim = len(controlPtArray[0][0])
    
    Rs = [
        [
            nurbsSurfaceBasis(i, j, pu, pv, U, V, u, v, weights)
            for i in range(rows + 1)
        ]
        for j in range(cols + 1)
    ]
    
    return tuple(
        sum(Rs[i][j] * controlPtArray[i][j][k]
                for i in range(rows + 1)
                for j in range(cols + 1))
                    for k in range(dim))


if __name__ == '__main__':
    print('Running main...')
    
    # Fixed knot vector index testing
    '''
    U = KnotVector([0, 1, 2, 3, 4, 5, 6])
    U.normalizeAll()
    
    for i, j in zip(range(len(U)), [0.0000001, 1.5, 2.5, 3.5, 4.5, 5.5, 6.9999999]):
        assert(U[i] == U[j])
    
    print('Testing successful!')
    '''
    
    # Randomized 3D spline testing code
    '''
    size = 5
    samples = [float(i) for i in range(-size, size + 1)]
    pts = list(zip(
        sample(samples, size),
        sample(samples, size),
        sample(samples, size)))
    # normals = list(zip(
    #     sample(samples, size),
    #     sample(samples, size),
    #     sample(samples, size)))
    # t = random()
    # pts.sort()
    
    bezierCurve(pts, _debug=True)
    '''
    
    # Fixed 2D B-spline testing code
    '''
    # controlPts = [(0, 1), (0.25, 0.75), (0.5, 0.5), (0.75, 0.25), (1, 0)]
    # controlPts = [(0, 0), (1, 2), (3, 0), (4, 2)]
    # controlPts = [(5, 5), (1, 1), (2, 0), (3, 1), (4, 0)]
    # controlPts = [(5, 5), (1, 1), (2, 0), (3, 1), (4, 0), (6, 2), (1, -5)]
    # controlPts = [(0, 0), (0, 1), (1, 1), (1, 0)]
    # controlPts = [(-3.5, -3.5), (-4.5, -1), (-3.5, 2), (-0.5, 3.5), (3.5, 2), (4.5, -2), (1, -3.5)]
    # controlPts = [(-1, 0), (0, 1), (1, 0)]
    # controlPts = [(167, 339), (69, 57)]
    controlPts = [(341, 116), (167, 339), (69, 57)]
    
    # U = KnotVector([0, 0, 0, 0, 0.14, 0.28, 0.42, 0.57, 0.71, 0.85, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 0, 0.3, 0.5, 0.5, 0.6, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 0, 1/5, 2/5, 3/5, 4/5, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 0, 0.1, 0.25, 0.5, 0.6, 0.85, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 0, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    U = KnotVector([0, 0, 0, 0, 1, 1, 1, 1])
    # U = KnotVector([0, 0, 0, 1, 1, 1])
    
    step = (U[-1] - U[0]) / (3 * len(U))
    knots = U.toUnnormalizedVector()
    interpolation = []
    mathematicaCode = []
    u = U[0]
    
    while u < U[-1]:
        pt = bSpline(controlPts, U, u, _debug=False)
        interpolation.append(f'{pt[0] :f},{pt[1] :f}')
        mathematicaCode.append(f'F[{u}]')
        u += step
    
    u = U[-1]
    pt = bSpline(controlPts, U, u, _debug=False)
    interpolation.append(f'{pt[0] :f},{pt[1] :f}')
    mathematicaCode.append(f'F[{u}]')
    
    interpolationStr = '\n'.join(interpolation)
    controlPtStrs = [f'{{{pt[0]},{pt[1]}}}' for pt in controlPts]
    controlPtStr = f'{{{",".join(controlPtStrs)}}}'
    
    # S[x_, i_] := DecimalForm[BSplineFunction[controlPtStr][x][[i]]]
    # F[u_] := Print[ToString@StringForm["``,``", S[u, 1], S[u, 2]]]
    mathematicaCode.insert(0, f'F[u_]:=Print[ToString@StringForm["``,``",S[u, 1],S[u, 2]]]')
    mathematicaCode.insert(0, f'S[x_,i_]:=DecimalForm[BSplineFunction[{controlPtStr}][x][[i]]]')
    
    mathematicaStr = ';'.join(mathematicaCode) # @UnusedVariable
    controlPtsStr = '\n'.join([f'{pt[0] :f},{pt[1] :f}' for pt in controlPts])
    
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnots:\n{knots}\n\nControl points:\n{controlPtsStr}')
    # print('\n\nMathematica test code:\n{mathematicaStr}')
    '''
    
    # TODO Fixed 3D B-spline surface testing code
    '''
    xs = (-3.5, -4.5, -3.5, -0.5, 3.5, 4.5, 1)
    ys = (-3.5, -1, 2, 3.5, 2, -2, -3.5)
    controlPtArray = [[(x, y, z) for x, y in zip(xs, ys)] for z in xs]
    
    U = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    V = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    
    # controlPtArray = [[(-1, 1, 0), (0, 1, 0), (1, 1, 0)],
    #                   [(-1, 0, 0), (0, 0, 1), (1, 0, 0)],
    #                   [(-1, -1, 0), (0, -1, 0), (1, -1, 0)]]
    #
    # U = KnotVector([0, 0, 0, 1, 1, 1])
    # V = KnotVector([0, 0, 0, 1, 1, 1])
    
    stepU = (U[-1] - U[0]) / (3 * len(U))
    stepV = (V[-1] - V[0]) / (3 * len(V))
    
    knotsU, knotsV = U.toUnnormalizedVector(), V.toUnnormalizedVector()
    interpolation = []
    
    u = knotsU[0]
    while u < V[-1]:
        v = knotsV[0]
        while v < V[-1]:
            pt = bSplineSurface(controlPtArray, U, V, u, v)
            interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
            v += stepV
        
        v = V[-1]
        pt = bSplineSurface(controlPtArray, U, V, u, v)
        interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
        u += stepU
    
    u, v = knotsU[-1], knotsV[0]
    while v < knotsV[-1]:
        pt = bSplineSurface(controlPtArray, U, V, u, v)
        interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
        v += stepV
    
    v = knotsV[-1]
    pt = bSplineSurface(controlPtArray, U, V, u, v)
    interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
    
    interpolationStr = '{' + ('\n,'.join(interpolation)) + '}'
    controlPtsStr = '{' + (','.join([f'({pt[0] :f},{pt[1] :f},{pt[2] :f})' for controlPts in controlPtArray for pt in controlPts])) + '}'
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnots (u):\n{knotsU}\n\nKnots (v):\n{knotsV}\n\nControl points:\n{controlPtsStr}')
    '''
    
    # Fixed 2D NURBS curve test code
    '''
    # controlPts = [(-3.5, -3.5), (-4.5, -1), (-3.5, 2), (-0.5, 3.5), (3.5, 2), (4.5, -2), (1, -3.5)]
    controlPts = [(-1, 0), (0, 1), (1, 0)]
    
    # U = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    U = KnotVector([0, 0, 0, 1, 1, 1])
    
    # weights = [1, 1, 1, 1, 1, 1, 1]
    # weights = [3, 3, 3, 3, 3, 3, 3]
    # weights = [1, 2, 3, 4, 3, 2, 1]
    # weights = [5, 1, 1, 1, 1, 1, 1]
    # weights = [1, 1, 1, 5, 1, 1, 1]
    # weights = [1, 1, 1]
    weights = [1, 3, 1]
    
    step = (U[-1] - U[0]) / (3 * len(U))
    knots = U.toUnnormalizedVector()
    interpolation = []
    u = U[0]
    
    while u < U[-1]:
        pt = nurbs(controlPts, U, u, weights, _debug=False)
        interpolation.append(f'{pt[0] :f},{pt[1] :f}')
        u += step
    
    u = U[-1]
    pt = nurbs(controlPts, U, u, weights, _debug=False)
    interpolation.append(f'{pt[0] :f},{pt[1] :f}')
    
    interpolationStr = '\n'.join(interpolation)
    controlPtsStr = '\n'.join([f'{pt[0] :f},{pt[1] :f}' for pt in controlPts])
    
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnots:\n{knots}\n\nControl points:\n{controlPtsStr}')
    '''
    
    # Fixed 3D NURBS surface testing code
    '''
    # xs = (-3.5, -4.5, -3.5, -0.5, 3.5, 4.5, 1)
    # ys = (-3.5, -1, 2, 3.5, 2, -2, -3.5)
    # controlPtArray = [[(x, y, z) for x, y in zip(xs, ys)] for z in xs]
    # weights = [[1, 1, 1],
    #            [1, 2, 1],
    #            [1, 1, 1]]
    #
    # U = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    # V = KnotVector([0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
    
    controlPtArray = [[(-1, 1, 0), (0, 1, 0), (1, 1, 0)],
                      [(-1, 0, 0), (0, 0, 2), (1, 0, 0)],
                      [(-1, -1, 0), (0, -1, 0), (1, -1, 0)]]
    weights = [[1, 1, 1],
               [1, 3, 1],
               [1, 1, 1]]
    
    U = KnotVector([0, 0, 0, 1, 1, 1])
    V = KnotVector([0, 0, 0, 1, 1, 1])
    
    stepU = (U[-1] - U[0]) / (3 * len(U))
    stepV = (V[-1] - V[0]) / (3 * len(V))
    
    knotsU, knotsV = U.toUnnormalizedVector(), V.toUnnormalizedVector()
    interpolation = []
    
    u = knotsU[0]
    while u < V[-1]:
        v = knotsV[0]
        while v < V[-1]:
            pt = nurbsSurface(controlPtArray, U, V, u, v, weights)
            interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
            v += stepV
        
        v = V[-1]
        pt = nurbsSurface(controlPtArray, U, V, u, v, weights)
        interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
        u += stepU
    
    u, v = knotsU[-1], knotsV[0]
    while v < knotsV[-1]:
        pt = nurbsSurface(controlPtArray, U, V, u, v, weights)
        interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
        v += stepV
    
    v = knotsV[-1]
    pt = nurbsSurface(controlPtArray, U, V, u, v, weights)
    interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
    
    interpolationStr = '{' + (','.join(interpolation)) + '}'
    controlPtsStr = '{' + (','.join([f'({pt[0] :f},{pt[1] :f},{pt[2] :f})' for controlPts in controlPtArray for pt in controlPts])) + '}'
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnots (u):\n{knotsU}\n\nKnots (v):\n{knotsV}\n\nControl points:\n{controlPtsStr}\n\nWeights:\n{weights}')
    '''
    