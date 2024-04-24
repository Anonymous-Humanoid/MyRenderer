'''
Created on Nov. 10, 2023

@author: Matthew
'''

from math import sqrt
from typing import List, Optional

from myRenderer import bSplineSurfaceBasis, KnotVector, Points3DLike, SplineException


# B_i(u,v) or B_i^m(u,v)
def tSplineBasis(U: KnotVector,
                 V: KnotVector,
                 u: float,
                 v: float,
                 pu: Optional[int] = None,
                 pv: Optional[int] = None) \
                        -> float:
    if pu is None:
        pu = len(U) - 2
    if pv is None:
        pv = len(V) - 2
    
    return bSplineSurfaceBasis(0, 0, pu, pv, U, V, u, v)

# S(u,v)
def tSpline(controlPts: Points3DLike,
            Us: List[KnotVector],
            Vs: List[KnotVector],
            u: float,
            v: float,
            weights: Optional[List[float]] = None,
            _debug: bool = False): # TODO Use _debug parameter
    n = len(controlPts) - 1
    
    if weights is not None and len(weights) - 1 != n:
        raise SplineException(f'T-splines require as many weights as control points. Expected {n + 1}, got: {len(weights)}')
    elif len(Us) - 1 != n:
        raise SplineException(f'T-splines require as many U knots as control points. Expected {n + 1}, got: {len(Us)}')
    elif len(Vs) - 1 != n:
        raise SplineException(f'T-splines require as many V knots as control points. Expected {n + 1}, got: {len(Vs)}')
    
    blends = [weights[i] * tSplineBasis(U, V, u, v)
              for i, U, V in zip(range(len(controlPts)), Us, Vs)]
    
    x, y, z, denom = 0, 0, 0, 0
    for pt, blend in zip(controlPts, blends):
        x += pt[0] * blend
        y += pt[1] * blend
        z += pt[2] * blend
        denom += blend
    
    x /= denom
    y /= denom
    z /= denom
    
    return x, y, z


# Fixed T-spline testing
if __name__ == '__main__':
    # Simple wave
    '''
    controlPts = [(-1, 1, 0),  (0, 1, 0),  (1, 1, 0),
                  (-1, 0, 0),  (0, 0, 2),  (1, 0, 0),
                  (-1, -1, 0), (0, -1, 0), (1, -1, 0)]
    
    weights = [1, 1, 1,
               1, -3, 1,
               1, 1, 1]
    
    Us = [
        KnotVector([0, 0, 0.25, 0.5, 0.75]), KnotVector([0, 0.25, 0.5, 0.75, 1]), KnotVector([0.25, 0.5, 0.75, 1, 1]),
        KnotVector([0, 0, 0.25, 0.5, 0.75]), KnotVector([0, 0.25, 0.5, 0.75, 1]), KnotVector([0.25, 0.5, 0.75, 1, 1]),
        KnotVector([0, 0, 0.25, 0.5, 0.75]), KnotVector([0, 0.25, 0.5, 0.75, 1]), KnotVector([0.25, 0.5, 0.75, 1, 1])
    ]
    
    Vs = [
        KnotVector([0, 0, 0.25, 0.5, 0.75]), KnotVector([0, 0, 0.25, 0.5, 0.75]), KnotVector([0, 0, 0.25, 0.5, 0.75]),
        KnotVector([0, 0.25, 0.5, 0.75, 1]), KnotVector([0, 0.25, 0.5, 0.75, 1]), KnotVector([0, 0.25, 0.5, 0.75, 1]),
        KnotVector([0.25, 0.5, 0.75, 1, 1]), KnotVector([0.25, 0.5, 0.75, 1, 1]), KnotVector([0.25, 0.5, 0.75, 1, 1]),
    ]
    '''
    
    # TODO Simple box
    '''
    controlPts = [
                      (-1,  1, -1), (-1, -1, -1),
                      ( 1,  1, -1), ( 1, -1, -1),
        ( 1,  1, -1), ( 1,  1,  1), ( 1, -1,  1), ( 1, -1, -1), ( 1,  1, -1),
        (-1,  1, -1), (-1,  1,  1), (-1, -1,  1), (-1, -1, -1), (-1,  1, -1),
                      (-1,  1, -1), (-1, -1, -1),
    ]
    
    weights = [
           1, 1,
           1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
           1, 1,
    ]
    
    Us = [
                                           KnotVector([0, 0, 1/3, 2/3, 1]),     KnotVector([0, 1/3, 2/3, 1, 1]),
                                           KnotVector([0, 0, 1/3, 2/3, 1]),     KnotVector([0, 1/3, 2/3, 1, 1]),
        KnotVector([0, 0, 1/6, 2/6, 3/6]), KnotVector([0, 1/6, 2/6, 3/6, 4/6]), KnotVector([1/6, 2/6, 3/6, 4/6, 5/6]), KnotVector([2/6, 3/6, 4/6, 5/6, 1]), KnotVector([3/6, 4/6, 5/6, 1, 1]),
        KnotVector([0, 0, 1/6, 2/6, 3/6]), KnotVector([0, 1/6, 2/6, 3/6, 4/6]), KnotVector([1/6, 2/6, 3/6, 4/6, 5/6]), KnotVector([2/6, 3/6, 4/6, 5/6, 1]), KnotVector([3/6, 4/6, 5/6, 1, 1]),
                                           KnotVector([0, 0, 1/3, 2/3, 1]),     KnotVector([0, 1/3, 2/3, 1, 1]),
    ]
    
    Vs = [
                                         KnotVector([3/6, 4/6, 5/6, 1, 1]),     KnotVector([3/6, 4/6, 5/6, 6/6, 1]),
                                         KnotVector([2/6, 3/6, 4/6, 5/6, 1]),   KnotVector([2/6, 3/6, 4/6, 5/6, 1]),
        KnotVector([0, 0, 2/3, 1, 1]),   KnotVector([1/6, 2/6, 3/6, 4/6, 5/6]), KnotVector([1/6, 2/6, 3/6, 4/6, 5/6]), KnotVector([0, 0, 2/3, 1, 1]),   KnotVector([0, 0, 2/3, 1, 1]),
        KnotVector([0, 0, 1/3, 2/3, 1]), KnotVector([0, 1/6, 2/6, 3/6, 4/6]),   KnotVector([0, 1/6, 2/6, 3/6, 4/6]),   KnotVector([0, 0, 1/3, 2/3, 1]), KnotVector([0, 0, 1/3, 2/3, 1]),
                                         KnotVector([0, 0, 1/6, 2/6, 3/6]),     KnotVector([0, 0, 1/6, 2/6, 3/6]),
    ]
    '''
    
    # TODO Simple sphere
    '''
    controlPts = [
        ( 0,  1,  0),
        ( 1,  0,  0),
        ( 0, -1,  0),
        (-1,  0,  0),
        ( 0,  0, -1),
        ( 0,  0,  1)
    ]
    
    weights = [
        1, 1, 1, 1, 1, 1
    ]
    
    Vs = [
        KnotVector([0, 0, 1/6, 2/6, 3/6]),
        KnotVector([0, 1/6, 2/6, 3/6, 4/6]),
        KnotVector([1/6, 2/6, 3/6, 4/6, 5/6]),
        KnotVector([2/6, 3/6, 4/6, 5/6, 1]),
    ]
    
    Us = [
        KnotVector([0, 0, 1/2, 1, 1]),
        KnotVector([0, 0, 1/2, 1, 1]),
        KnotVector([0, 0, 1/2, 1, 1]),
        KnotVector([0, 0, 1/2, 1, 1]),
        KnotVector([0, 0, 1/2, 1, 1]),
    ]
    '''
    
    invsqrt2 = 1 / sqrt(2)
    
    controlPts = [
        (-invsqrt2, invsqrt2, -1),
        (invsqrt2, -invsqrt2, 1),
        (-invsqrt2, 1 + invsqrt2, 0),
        (invsqrt2, -1 - invsqrt2, 0),
        (-1 - invsqrt2, invsqrt2, 0),
        (1 + invsqrt2, -invsqrt2, 0)
    ]
    
    weights = [
        1, 1, 1, 1, 1, 1
    ]
    
    Us = [
        KnotVector([0, 1, 2, 3, 3]),
        KnotVector([0, 0, 1, 2, 3]),
        KnotVector([0, 0, 1, 2, 2]),
        KnotVector([0, 0, 1, 2, 2]),
        KnotVector([0, 0, 1, 2, 3]),
        KnotVector([0, 1, 2, 3, 3])
    ]
    
    Vs = [
        KnotVector([0, 0, 1, 2, 3]),
        KnotVector([0, 1, 2, 3, 3]),
        KnotVector([0, 1, 2, 3, 3]),
        KnotVector([0, 0, 1, 2, 3]),
        KnotVector([0, 0, 1, 2, 2]),
        KnotVector([0, 0, 1, 2, 2])
    ]
    
    for U in Us:
        U.normalizeAll()
    for V in Vs:
        V.normalizeAll()
    
    stepU = 1 / 30
    stepV = 1 / 30
    
    interpolation = []
    knot0, knotn = 0, 1
    
    u = knot0 + stepU
    while u < knotn:
        # print(u) # todo TESTING
        v = knot0 + stepV
        while v < knotn:
            pt = tSpline(controlPts, Us, Vs, u, v, weights)
            interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
            v += stepV
        
        # v = knotn
        # pt = tSpline(controlPts, Us, Vs, u, v, weights)
        # interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
        u += stepU
    
    # u, v = knotn, knot0
    # while v < 1:
    #     pt = tSpline(controlPts, Us, Vs, u, v, weights)
    #     interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
    #     v += stepV
    
    # v = knotn
    # pt = tSpline(controlPts, Us, Vs, u, v, weights)
    # interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
    
    interpolationStr = '{' + (','.join(interpolation)) + '}'
    controlPtsStr = '{' + (','.join([f'({pt[0] :f},{pt[1] :f},{pt[2] :f})' for pt in controlPts])) + '}'
    UsStr = '[' + (', '.join([str(U) for U in Us])) + ']'
    VsStr = '[' + (', '.join([str(V) for V in Vs])) + ']'
    
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnot Vectors (u):\n{UsStr}\n\nKnot Vectors (v):\n{VsStr}\n\nControl points:\n{controlPtsStr}\n\nWeights:\n{weights}')
    