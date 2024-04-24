'''
Created on Oct. 16, 2023

@author: Matthew
'''

''' Credit to:
https://people.computing.clemson.edu/~dhouse/courses/405/notes/splines.pdf
https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Explicit_definition
https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node4.html
https://argos.vu/wp-content/uploads/2019/03/Springer-TheNURBSBook.pdf
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS/NURBS-property.html
(section 3.1) https://web.archive.org/web/20110725162753/http://alice.loria.fr/publications/papers/2006/SGP_Splines/SGP06-fin-electronic_mesh_to_spline.pdf
(equation 15.9) https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=1000&context=facpub
https://indico.ictp.it/event/a12191/session/35/contribution/20/material/0/0.pdf
'''

''' Additional reading:
https://www.ibiblio.org/e-notes/Splines/basis.html
https://www.cs.utexas.edu/~theshark/courses/cs354/lectures/cs354-16.pdf
https://math.stackexchange.com/a/421572/714488
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/
https://web.archive.org/web/20120714062727/http://www.tsplines.com/educationportal.html
https://en.wikipedia.org/wiki/De_Boor%27s_algorithm#Example_implementation
https://github.com/caadxyz/DeBoorAlgorithmNurbs
https://people.engr.tamu.edu/schaefer/research/slides/LRNonUniformTalk.pdf
https://www.inf.ed.ac.uk/teaching/courses/cg/lectures/lect15.pdf
https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/advgraphnotes.html
https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/
https://mirror.its.dal.ca/ctan/graphics/pstricks/contrib/pst-bspline/pst-bspline-doc.pdf
https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline
https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2020/huang/
https://www.researchgate.net/publication/234827617_T-splines_and_T-NURCCs
https://web.archive.org/web/20120713095920/http://www.tsplines.com/technicalpapers.html
https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
https://zingl.github.io/Bresenham.pdf
https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2009/c9brown/
'''

''' Surface notes:
- A surface is obtained as the tensor product of two curves,
  and as such, requires twice as many parameters as a single curve.
- This concept can be generalized to higher dimensions,
  however, we do not implement this functionality beyond three dimensions.
'''

# TODO Migrate to Point data structure

from . import KnotVector, Matrix, PointLike, PointsLike, Point2DLike, Point3DLike, Points2DLike, Points3DLike, SplineException
import numpy as np
from typing import List, Optional, Self, Sequence, Tuple

def _toVector(values: Sequence) -> Matrix:
    return Matrix([[i] for i in values])

def polynomialSpline2D(pts: Points2DLike) -> Matrix:
    degree = len(pts) - 1
    
    # Solving system of equations
    A = Matrix([[x ** i for i in range(degree + 1)] for x, _y in pts])
    B = _toVector([y for _x, y in pts])
    X = np.linalg.solve(A, B)
    
    return X

def polynomialSpline3D(pts: Points3DLike) -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = polynomialSpline2D(xys)
    XZ = polynomialSpline2D(xzs)
    
    return XY, XZ

def _fill(inner: List[float],
          leading: int,
          trailing: int,
          blocks: int) \
            -> List[float]:
    if leading < 0:
        raise SplineException('Cannot have less than 0 leading zeroes')
    elif trailing < 0:
        raise SplineException('Cannot have less than 0 trailing zeroes')
    elif blocks < 0:
        raise SplineException('Cannot have blocks of less than 0 zeroes')
    
    start = [0] * blocks * leading
    end = [0] * blocks * trailing
    
    return start + inner + end

def piecewisePolynomialSpline2D(pts: Points2DLike,
                                startSlope: float,
                                endSlope: float) \
                                    -> Matrix:
    degree = len(pts) - 1
    a = []
    b = [startSlope]
    
    # Start slope s0
    x0 = pts[0][0]
    inner = [0, 1, x0, 2 * x0**2]
    a.append(_fill(inner, 0, degree - 1, 4))
    
    for i, (xi, yi) in enumerate(pts[:-1]):
        # Left c0 continuity
        inner = [1, xi, xi**2, xi**3]
        a.append(_fill(inner, i, degree - i - 1, 4))
        b.append(yi)
        
        # Right c0 continuity
        xi_1, yi_1 = pts[i + 1]
        inner = [1, xi_1, xi_1**2, xi_1**3]
        a.append(_fill(inner, i, degree - i - 1, 4))
        b.append(yi_1)
        
        if degree - i - 2 >= 0:
            # c1 continuity
            inner = [0, 1, 2 * xi_1, 3 * xi_1**2,
                     0, -1, -2 * xi_1, -3 * xi_1**2]
            a.append(_fill(inner, i, degree - i - 2, 4))
            b.append(0)
            
            # c2 continuity
            inner = [0, 0, 1, 3 * xi_1,
                     0, 0, -1, -3 * xi_1]
            a.append(_fill(inner, i, degree - i - 2, 4))
            b.append(0)
    
    # End slope sn
    xn = pts[-1][0]
    inner = [0, 1, xn, 2 * xn**2]
    a.append(_fill(inner, degree - 1, 0, 4))
    b.append(endSlope)
    
    # Solving system of equations
    A = Matrix(a)
    B = _toVector(b)
    X = np.linalg.solve(A, B)
    
    return X

def piecewisePolynomialSpline3D(pts: Points3DLike,
                                startSlope: Point2DLike,
                                endSlope: Point2DLike) \
                                    -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = piecewisePolynomialSpline2D(xys, startSlope[0], endSlope[0])
    XZ = piecewisePolynomialSpline2D(xzs, startSlope[1], endSlope[1])
    
    return XY, XZ

def parameterizedSpline2D(pts: Points2DLike,
                          startSlope: float,
                          endSlope: float) \
                                -> Matrix:
    degree = len(pts) - 1
    a = []
    b = []
    
    # Start slope s0
    inner = [1, 0, 0]
    a.append(_fill(inner, 0, degree - 1, 3))
    b.append(startSlope)
    
    for i, (xi, yi) in enumerate(pts[:-1]):
        # Right c0 continuity (left is ai = yi, which is assumed)
        xi_1, yi_1 = pts[i + 1]
        inner = [1, 1, 1]
        a.append(_fill(inner, i, degree - i - 1, 3))
        b.append(yi_1 - yi)
        
        if degree - i - 2 >= 0:
            # c1 continuity
            xi_2 = pts[i + 2][0]
            inner = [1, 2, 3,
                     -(xi_1 - xi) / (xi_2 - xi_1), 0, 0]
            a.append(_fill(inner, i, degree - i - 2, 3))
            b.append(0)
            
            # c2 continuity
            inner = [0, 1, 3,
                     0, -((xi_1 - xi) / (xi_2 - xi_1))**2, 0]
            a.append(_fill(inner, i, degree - i - 2, 3))
            b.append(0)
    
    # End slope sn
    inner = [1, 2, 3]
    a.append(_fill(inner, degree - 1, 0, 3))
    b.append(endSlope)
    
    # Solving system of equations
    A = Matrix(a)
    B = _toVector(b)
    X = np.linalg.solve(A, B)
    
    return X

def parameterizedSpline3D(pts: Points3DLike,
                          startSlope: Point2DLike,
                          endSlope: Point2DLike) \
                                -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = parameterizedSpline2D(xys, startSlope[0], endSlope[0])
    XZ = parameterizedSpline2D(xzs, startSlope[1], endSlope[1])
    
    return XY, XZ

def _spaceSplineAxis(pts: Sequence[float],
                     closed: bool = False,
                     startSlope: float = 1,
                     endSlope: float = 1) \
                        -> Matrix:
    degree = len(pts) - 1
    a, b = [], []
    
    if closed:
        # Endpoint slope
        leading = [1, 0, 0]
        inner = [0] * 3 * (degree - 2)
        trailing = [-1, -2, -3]
        slope = leading + inner + trailing
        a.append(slope)
        b.append(0)
    else:
        # Start slope s0
        inner = [1, 0, 0]
        a.append(_fill(inner, 0, degree - 1, 3))
        b.append(startSlope)
    
    for i, ui in enumerate(pts[:-1]):
        # Right c0 continuity
        ui_1 = pts[i + 1]
        inner = [1, 1, 1]
        a.append(_fill(inner, i, degree - i - 1, 3))
        b.append(ui_1 - ui)
        
        if degree - i - 2 >= 0:
            # c1 continuity
            inner = [1, 2, 3,
                     -1, 0, 0]
            a.append(_fill(inner, i, degree - i - 2, 3))
            b.append(0)
            
            # c2 continuity
            inner = [0, 1, 3,
                     0, -1, 0]
            a.append(_fill(inner, i, degree - i - 2, 3))
            b.append(0)
    
    if closed:
        # Endpoint curvature
        leading = [0, 1, 0]
        inner = [0] * 3 * (degree - 2)
        trailing = [0, -1, -3]
        curvature = leading + inner + trailing
        a.append(curvature)
        b.append(0)
    else:
        # End slope sn
        inner = [1, 2, 3]
        a.append(_fill(inner, degree - 1, 0, 3))
        b.append(endSlope)
    
    # Solving system of equations
    A = Matrix(a)
    B = _toVector(b)
    U = np.linalg.solve(A, B)
    
    return U

def spaceSpline2D(pts: Points2DLike,
                  startSlope: Optional[float] = None,
                  endSlope: Optional[float] = None,
                  closed: bool = False) \
                        -> Tuple[Matrix, Matrix]:
    if closed and not (startSlope is None or endSlope is None):
        raise SplineException('Closed space splines cannot have start and end slopes specified')
    elif (not closed) and (startSlope is None or endSlope is None):
        raise SplineException('Open space splines must have start and end slopes specified')
    
    if closed:
        # Endpoint c0 continuity
        pts.append(pts[0])
    
    xs = [xi for xi, _yi in pts]
    ys = [yi for _xi, yi in pts]
    X = _spaceSplineAxis(xs, closed)
    Y = _spaceSplineAxis(ys, closed, startSlope, endSlope)
    
    if closed:
        # Removing appended element
        del pts[-1]
    
    return X, Y

def spaceSpline3D(pts: Points3DLike,
                  startSlope: Optional[Point2DLike] = None,
                  endSlope: Optional[Point2DLike] = None,
                  closed: bool = False) \
                        -> Tuple[Matrix, Matrix, Matrix]:
    if closed and not (startSlope is None or endSlope is None):
        raise SplineException('Closed space splines cannot have start and end slopes specified')
    elif (not closed) and (startSlope is None or endSlope is None):
        raise SplineException('Open space splines must have start and end slopes specified')
    
    if closed:
        # Endpoint c0 continuity
        pts.append(pts[0])
    
    xs = [xi for xi, _yi, _zi in pts]
    ys = [yi for _xi, yi, _zi in pts]
    zs = [zi for _xi, _yi, zi in pts]
    
    X = _spaceSplineAxis(xs, closed)
    Y = _spaceSplineAxis(ys, closed,
                         None if startSlope is None else startSlope[0],
                         None if endSlope is None else endSlope[0])
    Z = _spaceSplineAxis(zs, closed,
                         None if startSlope is None else startSlope[1],
                         None if endSlope is None else endSlope[1])
    
    if closed:
        # Removing appended element
        del pts[-1]
    
    return X, Y, Z

def hermiteCubicSplineSegment2D(pt1: Point2DLike,
                                pt2: Point2DLike,
                                vec1: Point2DLike,
                                vec2: Point2DLike) \
                                        -> Tuple[Matrix, Matrix]:
    (x1, y1), (x2, y2), (vx1, vy1), (vx2, vy2) = pt1, pt2, vec1, vec2
    
    # Solving system of equations
    invA = Matrix([[1, 0, 0, 0],
                   [0, 0, 1, 0],
                   [-3, 3, -2, -1],
                   [2, -2, 1, 1]])
    # B can theoretically be optimized for updating the spline based on changes in only vectors or control points
    B = Matrix([[x1, y1],
                [x2, y2],
                [vx1, vy1],
                [vx2, vy2]])
    U = invA * B
    X, Y = U[:,0], U[:,1]
    
    return X, Y

def hermiteCubicSplineSegment3D(pt1: Point3DLike,
                                pt2: Point3DLike,
                                vec1: Point3DLike,
                                vec2: Point3DLike) \
                                        -> Tuple[Matrix, Matrix, Matrix]:
    (x1, y1, z1), (x2, y2, z2), (vx1, vy1, vz1), (vx2, vy2, vz2) = pt1, pt2, vec1, vec2
    
    # Solving system of equations
    invA = Matrix([[1, 0, 0, 0],
                   [0, 0, 1, 0],
                   [-3, 3, -2, -1],
                   [2, -2, 1, 1]])
    # B can theoretically be optimized for updating the spline based on changes in only vectors or control points
    B = Matrix([[x1, y1, z1],
                [x2, y2, z2],
                [vx1, vy1, vz1],
                [vx2, vy2, vz2]])
    U = invA * B
    X, Y, Z = U[:,0], U[:,1], U[:,2]
    
    return X, Y, Z

def hermiteCubicSpline2D(pts: Points2DLike,
                         vecs: Points2DLike) \
                                -> Tuple[List[Matrix], List[Matrix]]:
    if len(pts) != len(vecs):
        raise SplineException('Must have equal quantities of control points and vectors in a hermite cubic spline')
    
    Xs, Ys = [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], vecs[:-1])):
        pt2, vec2 = pts[i + 1], vecs[i + 1]
        X, Y = hermiteCubicSplineSegment2D(pt1, pt2, vec1, vec2)
        Xs.append(X)
        Ys.append(Y)
    
    return Xs, Ys

def hermiteCubicSpline3D(pts: Points3DLike,
                         vecs: Points3DLike) \
                                -> Tuple[List[Matrix], List[Matrix], List[Matrix]]:
    if len(pts) != len(vecs):
        raise SplineException('Must have equal quantities of control points and vectors in a hermite cubic spline')
    
    Xs, Ys, Zs = [], [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], vecs[:-1])):
        pt2, vec2 = pts[i + 1], vecs[i + 1]
        X, Y, Z = hermiteCubicSplineSegment3D(pt1, pt2, vec1, vec2)
        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
    
    return Xs, Ys, Zs

def catmullRomSplineSegment2D(pt1: Point2DLike,
                              pt2: Point2DLike,
                              normal1: Point2DLike,
                              normal2: Point2DLike) \
                                    -> Tuple[Matrix, Matrix]:
    return cardinalSplineSegment2D(pt1, pt2, normal1, normal2, 0)

def catmullRomSplineSegment3D(pt1: Point3DLike,
                              pt2: Point3DLike,
                              normal1: Point3DLike,
                              normal2: Point3DLike) \
                                    -> Tuple[Matrix, Matrix]:
    return cardinalSplineSegment3D(pt1, pt2, normal1, normal2, 0)

def catmullRomSpline2D(pts: Points2DLike,
                       normals: Points2DLike) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    return cardinalSpline2D(pts, normals, 0)

def catmullRomSpline3D(pts: Points3DLike,
                       normals: Points3DLike) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    return cardinalSpline3D(pts, normals, 0)

def cardinalSplineSegment2D(pt1: Point2DLike,
                            pt2: Point2DLike,
                            normal1: Point2DLike,
                            normal2: Point2DLike,
                            t: float) \
                                    -> Tuple[Matrix, Matrix]:
    p = (1 - t) / 2
    vec1 = (p * (pt2[0] - normal1[0]),
            p * (pt2[1] - normal1[1]))
    vec2 = (p * (normal2[0] - pt1[0]),
            p * (normal2[1] - pt1[1]))
    return hermiteCubicSplineSegment2D(pt1, pt2, vec1, vec2)

def cardinalSplineSegment3D(pt1: Point3DLike,
                            pt2: Point3DLike,
                            normal1: Point3DLike,
                            normal2: Point3DLike,
                            t: float) \
                                    -> Tuple[Matrix, Matrix, Matrix]:
    p = (1 - t) / 2
    vec1 = (p * (pt2[0] - normal1[0]),
            p * (pt2[1] - normal1[1]),
            p * (pt2[2] - normal1[2]))
    vec2 = (p * (normal2[0] - pt1[0]),
            p * (normal2[1] - pt1[1]),
            p * (normal2[2] - pt1[2]))
    return hermiteCubicSplineSegment3D(pt1, pt2, vec1, vec2)

def cardinalSpline2D(pts: Points2DLike,
                     normals: Points2DLike,
                     t: float) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    if len(normals) != len(pts):
        raise SplineException('Cardinal splines require equal quantities of points and normals')
    
    Xs, Ys = [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], normals[:-1])):
        pt2, vec2 = pts[i + 1], normals[i + 1]
        X, Y = cardinalSplineSegment2D(pt1, pt2, vec1, vec2, t)
        Xs.append(X)
        Ys.append(Y)
    
    return Xs, Ys

def cardinalSpline3D(pts: Points3DLike,
                     normals: Points3DLike,
                     t: float) \
                            -> Tuple[List[Matrix], List[Matrix], List[Matrix]]:
    if len(normals) != len(pts):
        raise SplineException('Cardinal splines require equal quantities of points and normals')
    
    Xs, Ys, Zs = [], [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], normals[:-1])):
        pt2, vec2 = pts[i + 1], normals[i + 1]
        X, Y, Z = cardinalSplineSegment3D(pt1, pt2, vec1, vec2, t)
        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
    
    return Xs, Ys, Zs

# TODO Use DefaultWeights class in place of future weights: Optional[List[float]] = None
class DefaultWeights:
    def __init__(self, length: int) -> Self:
        self.__length = length
    
    def __getitem__(self, _) -> int:
        return 1
    
    def __len__(self) -> int:
        return self.__length

# Implements the Cox-De Boor algorithm
# N_(i,p)[U](u) or N_(i,p)(u) or N_i(u)
def bSplineBasis(i: int,
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
        result += (u - U[i]) / denom1 * bSplineBasis(i, p - 1, U, u)
    if denom2 != 0:
        result += (U[i + p + 1] - u) / denom2 * bSplineBasis(i + 1, p - 1, U, u)

    return result

# B(u)
def bSpline(controlPts: PointsLike,
            U: KnotVector,
            u: float,
            maxDegree: Optional[int] = None) \
                    -> PointLike:
    # Extra optimization for special cases
    # NOTE: Techincally should be a zero vector for u < 0 or u > 1,
    # but due to floating point errors in checking equality, we don't care
    if u <= 0:
        return tuple(c for c in controlPts[0])
    elif u >= 1:
        return tuple(c for c in controlPts[-1])
    
    dim = len(controlPts[0])
    
    m, n = len(U) - 1, len(controlPts) - 1
    tempP = m - n - 1
    p = tempP if maxDegree is None else min(maxDegree, tempP)

    Ns = [bSplineBasis(i, p, U, u) for i in range(n + 1)]
    
    return tuple(sum(N * c[i] for c, N in zip(controlPts, Ns)) for i in range(dim))

# TODO Take advantage of the maxDegree parameter of the bspline function

# N_(i,j)(u,v)
def bSplineSurfaceBasis(i: int,
                        j: int,
                        pu: int,
                        pv: int,
                        U: KnotVector,
                        V: KnotVector,
                        u: float,
                        v: float) \
                                -> float:
    return bSplineBasis(i, pu, U, u) * bSplineBasis(j, pv, V, v)

# B(u,v) or B[U,V](s,t)
def bSplineSurface(controlPtArray: List[PointsLike],
                   U: KnotVector,
                   V: KnotVector,
                   u: float,
                   v: float) \
                        -> PointLike:
    dim = len(controlPtArray[0][0])
    
    mu, mv = len(U) - 1, len(V) - 1
    rows, cols = len(controlPtArray) - 1, len(controlPtArray[0]) - 1
    pu, pv = mu - cols - 1, mv - rows - 1
    
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
               weights: List[float]) \
                    -> float:
    n = len(weights) - 1
    
    Ns = [bSplineBasis(j, p, U, u) * weights[j] for j in range(n + 1)]
    
    return Ns[i] / sum(Ns)

# C(u)
def nurbs(controlPts: PointsLike,
          U: KnotVector,
          u: float,
          weights: List[float]) \
                -> PointLike:
    n = len(controlPts) - 1
    
    if len(weights) != n + 1:
        raise SplineException(f'NURBS curves require as many weights as control points. Expected {n + 1}, got: {len(weights)}')
    
    # Extra optimization for special cases
    # NOTE: Techincally should be a zero vector for u < 0 or u > 1,
    # but due to floating point errors in checking equality, we don't care
    if u <= 0:
        return controlPts[0]
    elif u >= 1:
        return controlPts[-1]
    
    dim = len(controlPts[0])
    
    m = len(U) - 1
    p = m - n - 1
    
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
                      weights: List[List[float]]) \
                            -> float:
    rows, cols = len(weights) - 1, len(weights[0]) - 1
    
    Ns = [
        [bSplineSurfaceBasis(_i, _j, pu, pv, U, V, u, v) * weights[_i][_j] for _i in range(cols + 1)]
        for _j in range(rows + 1)
    ]
    
    return Ns[i][j] / sum(c for cs in Ns for c in cs)

# C(u,v)
def nurbsSurface(controlPtArray: List[PointsLike],
                   U: KnotVector,
                   V: KnotVector,
                   u: float,
                   v: float,
                   weights: List[List[float]]) \
                        -> PointLike:
    # Extra optimization for special cases
    dim = len(controlPtArray[0][0])
    
    mu, mv = len(U) - 1, len(V) - 1
    rows, cols = len(controlPtArray) - 1, len(controlPtArray[0]) - 1
    pu, pv = mu - cols - 1, mv - rows - 1
    
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
            weights: Optional[List[float]] = None):
    n = len(controlPts) - 1
    
    if weights is not None and len(weights) - 1 != n:
        raise SplineException(f'T-splines require as many weights as control points. Expected {n + 1}, got: {len(weights)}')
    elif len(Us) - 1 != n:
        raise SplineException(f'T-splines require as many U knots as control points. Expected {n + 1}, got: {len(Us)}')
    elif len(Vs) - 1 != n:
        raise SplineException(f'T-splines require as many V knots as control points. Expected {n + 1}, got: {len(Vs)}')
    
    if weights is None:
        weights = DefaultWeights(len(controlPts))

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

# TODO Test and add Beziers
