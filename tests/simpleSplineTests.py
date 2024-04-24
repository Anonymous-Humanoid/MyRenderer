'''
Created on Nov. 1, 2023

@author: Matthew
'''

from random import random, sample  # @UnusedImport
from typing import List, Optional, Sequence, Tuple

from myRenderer import Matrix, PointLike, Point2DLike, Point3DLike, PointsLike, Points2DLike, Points3DLike, SplineException

import numpy as np


# TODO Unused
# Evaluates A + Bx + Cx^2 + ... = A + x(B + x(C + ...))
def evalPolynomial(coefficients: Sequence[float],
                   x: float) -> float:
    out = 0
    for coefficient in coefficients[::-1]:
        out = coefficient + x * out
    return out

def _toVector(values: Sequence) -> Matrix:
    return Matrix([[i] for i in values])

def polynomialSpline2D(pts: Points2DLike,
                       _debug: bool = False) -> Matrix:
    degree = len(pts) - 1
    
    # Solving system of equations
    A = Matrix([[x ** i for i in range(degree + 1)] for x, _y in pts])
    B = _toVector([y for _x, y in pts])
    X = np.linalg.solve(A, B)
    
    # Testing code
    if _debug:
        coefficients = [i[0] for i in np.array(X)]
        eq = ('+'.join([(f'{coefficients[i]:.15f}x^{{{i}}}') for i in range(len(pts))])
                 .replace('x^{0}', '')
                 .replace('^{1}', ''))
        xs = [x for x, _y in pts]
        ys = [y for _x, y in pts]
        xMin = min(xs)
        xMax = max(xs)
        
        print(f'f(x)={eq}\\left\\{{{xMin}\le x\le {xMax}\\right\\}}'
              .replace('+-', '-')
              .replace('-+', '-'))
        print(f'X={xs}\n(X,f(X))\n(X,{ys})'
              .replace(' ', ''))
    
    return X

def polynomialSpline3D(pts: Points3DLike,
                       _debug: bool = False) -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = polynomialSpline2D(xys)
    XZ = polynomialSpline2D(xzs)
    
    # Testing code
    if _debug:
        xycoefficients = [i[0] for i in np.array(XY)]
        xzcoefficients = [i[0] for i in np.array(XZ)]
        xyeq = ('+'.join([(f'{xycoefficients[i]:.15f}x^({i})') for i in range(len(pts))])
                  .replace('x^(0)', '')
                  .replace('^(1)', ''))
        xzeq = ('+'.join([(f'{xzcoefficients[i]:.15f}x^({i})') for i in range(len(pts))])
                  .replace('x^(0)', '')
                  .replace('^(1)', ''))
        xs = [x for  x, _y, _z in pts]
        ys = [y for _x,  y, _z in pts]
        zs = [z for _x, _y,  z in pts]
        xMin = min(xs)
        xMax = max(xs)
        
        print(f'f=Curve(x,{xyeq},{xzeq},x,{xMin},{xMax}))'
              .replace('+-', '-')
              .replace('-+', '-'))
        print(f'X={xs}\nf(X)\n(X,{ys},{zs})'
              .replace('[', '{')
              .replace(']', '}')
              .replace(' ', ''))
    
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
                                endSlope: float,
                                _debug: bool = False) \
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
    
    # Testing code
    if _debug:
        polynomials = []
        for i, row in enumerate(np.array(X)):
            xi = row[0]
            if i % 4 == 0:
                polynomials.append([])
            polynomials[-1].append(xi)
        polynomialStrs = [f'{coefs[0]}+{coefs[1]}x+{coefs[2]}x^{{2}}+{coefs[3]}x^{{3}}'
                          for coefs in polynomials]
        xs = []
        ys = []
        for xi, yi, in pts:
            xs.append(xi)
            ys.append(yi)
        
        conditions = []
        for i, (xi, _yi) in enumerate(pts[:-1]):
            xi_1, _yi_1 = pts[i + 1]
            conditions.append(f'{xi}\le x\le {xi_1}')
        
        graph = '\\left\\{'
        for polynomial, condition in zip(polynomialStrs, conditions):
            graph += f'{condition}:{polynomial},'
        graph = ((graph[:-1] + '\\right\\}')
                 .replace('+-', '-')
                 .replace('-+', '-')
                 if len(graph) > 1 else '\\left\\{\\right\\}')
        
        print(graph)
        print(f'({xs},{ys})'.replace(' ', ''))
    
    return X

def piecewisePolynomialSpline3D(pts: Points3DLike,
                                startSlope: Point2DLike,
                                endSlope: Point2DLike,
                                _debug: bool = False) \
                                    -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = piecewisePolynomialSpline2D(xys, startSlope[0], endSlope[0])
    XZ = piecewisePolynomialSpline2D(xzs, startSlope[1], endSlope[1])
    
    # Testing code
    if _debug:
        xypolynomials = []
        xzpolynomials = []
        for i, row in enumerate(np.array(XY)):
            xi = row[0]
            if i % 4 == 0:
                xypolynomials.append([])
            xypolynomials[-1].append(xi)
        for i, row in enumerate(np.array(XZ)):
            xi = row[0]
            if i % 4 == 0:
                xzpolynomials.append([])
            xzpolynomials[-1].append(xi)
        eqs = [f'x,{xycoefs[0]}+{xycoefs[1]}x+{xycoefs[2]}x^(2)+{xycoefs[3]}x^(3),' +
               f'{xzcoefs[0]}+{xzcoefs[1]}x+{xzcoefs[2]}x^(2)+{xzcoefs[3]}x^(3),x'
               for xycoefs, xzcoefs in zip(xypolynomials, xzpolynomials)]
        
        xs, ys, zs = [], [], []
        for xi, yi, zi in pts:
            xs.append(xi)
            ys.append(yi)
            zs.append(zi)
        
        bounds = []
        for i, (xi, _yi, _zi) in enumerate(pts[:-1]):
            xi_1, _yi_1, _zi_1 = pts[i + 1]
            bounds.append((xi, xi_1))
        
        graph = []
        for eq, bound in zip(eqs, bounds):
            xMin, xMax = bound
            graph.append(f'Curve({eq},{xMin},{xMax})')
        curve = ("\n".join(graph)
                 .replace('+-', '-')
                 .replace('-+', '-'))
        
        print(curve)
        print(f'({xs},{ys},{zs})'
              .replace(' ', '')
              .replace('[', '{')
              .replace(']', '}'))
    
    return XY, XZ

def parameterizedSpline2D(pts: Points2DLike,
                          startSlope: float,
                          endSlope: float,
                          _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        graphs = []
        polynomials = []
        xs, ys = [], []
        
        for xi, yi in pts:
            xs.append(xi)
            ys.append(yi)
        
        for i, row in enumerate(np.array(X)):
            xi = row[0]
            if i % 3 == 0:
                polynomials.append([])
            polynomials[-1].append(xi)
        
        for i, (xi, yi) in enumerate(pts[:-1]):
            xi_1, _yi_1 = pts[i + 1]
            dx = xi_1 - xi
            coefs = polynomials[i]
            polynomial = f'{yi}+{coefs[0]}t+{coefs[1]}t^{{2}}+{coefs[2]}t^{{3}}'
            graphs.append(f'{dx}t+{xi},{polynomial}')
        
        graph = (('\\left(' + ('\\right)\n\\left('.join(graphs)) + '\\right)')
                 .replace('+-', '-')
                 .replace('-+', '-')
                 .replace('--', '+'))
        
        print(graph)
        print(f'({xs},{ys})'.replace(' ', ''))
    
    return X

def parameterizedSpline3D(pts: Points3DLike,
                          startSlope: Point2DLike,
                          endSlope: Point2DLike,
                          _debug: bool = False) \
                                -> Tuple[Matrix, Matrix]:
    xys = []
    xzs = []
    
    for x, y, z in pts:
        xys.append((x, y))
        xzs.append((x, z))
    
    XY = parameterizedSpline2D(xys, startSlope[0], endSlope[0])
    XZ = parameterizedSpline2D(xzs, startSlope[1], endSlope[1])
    
    # Testing code (parametric form)
    if _debug:
        graphs, xypolynomials, xzpolynomials, xs, ys, zs = [], [], [], [], [], []
        
        for xi, yi, zi in pts:
            xs.append(xi)
            ys.append(yi)
            zs.append(zi)
        
        for i, row in enumerate(np.array(XY)):
            xi = row[0]
            if i % 3 == 0:
                xypolynomials.append([])
            xypolynomials[-1].append(xi)
        for i, row in enumerate(np.array(XZ)):
            xi = row[0]
            if i % 3 == 0:
                xzpolynomials.append([])
            xzpolynomials[-1].append(xi)
        
        for i, (xi, yi, zi) in enumerate(pts[:-1]):
            xi_1, _yi_1, _zi_1 = pts[i + 1]
            dx = xi_1 - xi
            xycoefs, xzcoefs = xypolynomials[i], xzpolynomials[i]
            graphs.append(f'Curve({dx}t+{xi},{yi}+{xycoefs[0]}t+{xycoefs[1]}t^(2)+{xycoefs[2]}t^(3),' +
                          f'{zi}+{xzcoefs[0]}t+{xzcoefs[1]}t^(2)+{xzcoefs[2]}t^(3),t,0,1)')
        
        print('\n'.join(graphs)
                  .replace('+-', '-')
                  .replace('-+', '-')
                  .replace('--', '+'))
        print(f'({xs},{ys},{zs})'
              .replace(' ', '')
              .replace('[', '{')
              .replace(']', '}'))
    
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
                  closed: bool = False,
                  _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        x0, y0 = xs[0], ys[0]
        xn, yn = xs[-1], ys[-1]
        graphs = []
        xpolynomials = []
        ypolynomials = []
        
        for (i, xrow), yrow in zip(enumerate(np.array(X)), np.array(Y)):
            xi, yi = xrow[0], yrow[0]
            if i % 3 == 0:
                xpolynomials.append([])
                ypolynomials.append([])
            xpolynomials[-1].append(xi)
            ypolynomials[-1].append(yi)
        
        for i, (xi, yi) in enumerate(pts[:-1]):
            xcoefs, ycoefs = xpolynomials[i], ypolynomials[i]
            xpolynomial = f'{xi:f}+{xcoefs[0]:f}t+{xcoefs[1]:f}t^{{2}}+{xcoefs[2]:f}t^{{3}}'
            ypolynomial = f'{yi:f}+{ycoefs[0]:f}t+{ycoefs[1]:f}t^{{2}}+{ycoefs[2]:f}t^{{3}}'
            graphs.append(f'{xpolynomial},{ypolynomial}')
        
        graph = '\\left(' + ('\\right)\n\\left('.join(graphs)) + '\\right)'
        points = f'\\left({xs},{ys}\\right)'.replace(' ', '')
        data = (f'{graph}\n{points}')
        
        if not closed:
            startLine = f'\\left(t+{(x0-0.5):f},{startSlope:f}t+{(y0-startSlope/2):f}\\right)'
            endLine = f'\\left(t+{(xn-0.5):f},{endSlope:f}t+{(yn-endSlope/2):f}\\right)'
            data += f'\n{startLine}\n{endLine}'
        
        data = (data.replace('+-', '-')
                    .replace('-+', '-')
                    .replace('--', '+'))
        print(data)
    
    if closed:
        # Removing appended element
        del pts[-1]
    
    return X, Y

def spaceSpline3D(pts: Points3DLike,
                  startSlope: Optional[Point2DLike] = None,
                  endSlope: Optional[Point2DLike] = None,
                  closed: bool = False,
                  _debug: bool = False) \
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
    x0, y0, z0 = xs[0], ys[0], zs[0]
    
    X = _spaceSplineAxis(xs, closed)
    Y = _spaceSplineAxis(ys, closed,
                         None if startSlope is None else startSlope[0],
                         None if endSlope is None else endSlope[0])
    Z = _spaceSplineAxis(zs, closed,
                         None if startSlope is None else startSlope[1],
                         None if endSlope is None else endSlope[1])
    
    # Testing code (parametric form)
    if _debug:
        xn, yn, zn = xs[-1], ys[-1], zs[-1]
        xpolynomials, ypolynomials, zpolynomials = [], [], []
        curve = ''
        
        for (i, xrow), yrow, zrow in zip(enumerate(np.array(X)), np.array(Y), np.array(Z)):
            xi, yi, zi = xrow[0], yrow[0], zrow[0]
            if i % 3 == 0:
                xpolynomials.append([])
                ypolynomials.append([])
                zpolynomials.append([])
            xpolynomials[-1].append(xi)
            ypolynomials[-1].append(yi)
            zpolynomials[-1].append(zi)
        
        for i, (xi, yi, zi) in enumerate(pts[:-1]):
            xcoefs, ycoefs, zcoefs = xpolynomials[i], ypolynomials[i], zpolynomials[i]
            curve += (f'Curve({xi:f}+{xcoefs[0]:f}t+{xcoefs[1]:f}t^(2)+{xcoefs[2]:f}t^(3),' +
                      f'{yi:f}+{ycoefs[0]:f}t+{ycoefs[1]:f}t^(2)+{ycoefs[2]:f}t^(3),' +
                      f'{zi:f}+{zcoefs[0]:f}t+{zcoefs[1]:f}t^(2)+{zcoefs[2]:f}t^(3),t,0,1)\n')
        
        points = (f'({xs},{ys},{zs})'
                  .replace(' ', '')
                  .replace('[', '{')
                  .replace(']', '}'))
        data = (f'{curve}{points}')
        
        if not closed:
            startLine = (f'Curve(t+{(x0-0.5):f},' +
                         f'{startSlope[0]:f}t+{(y0-startSlope[0]/2):f},' +
                         f'{startSlope[1]:f}t+{(z0-startSlope[1]/2):f},t,0,1)')
            endLine = (f'Curve(t+{(xn-0.5):f},' +
                       f'{endSlope[0]:f}t+{(yn-endSlope[0]/2):f},' +
                       f'{endSlope[1]:f}t+{(zn-endSlope[1]/2):f},t,0,1)')
            lineData = (f'\n{startLine}\n{endLine}'
                        .replace(' ', '')
                        .replace('[', '{')
                        .replace(']', '}'))
            data += lineData
        
        data = (data.replace('+-', '-')
                    .replace('-+', '-')
                    .replace('--', '+'))
        print(data)
    
    if closed:
        # Removing appended element
        del pts[-1]
    
    return X, Y, Z

def hermiteCubicSplineSegment2D(pt1: Point2DLike,
                                pt2: Point2DLike,
                                vec1: Point2DLike,
                                vec2: Point2DLike,
                                _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        vector1 = f'\\left({vx1}t+{x1},{vy1}t+{y1}\\right)'
        vector2 = f'\\left({vx2}t+{x2},{vy2}t+{y2}\\right)'
        
        pts = f'\\left(\\left[{x1},{x2}\\right],\\left[{y1},{y2}\\right]\\right)'
        
        xs = f'{X[0, 0]}+{X[1, 0]}t+{X[2, 0]}t^{{2}}+{X[3, 0]}t^{{3}}'
        ys = f'{Y[0, 0]}+{Y[1, 0]}t+{Y[2, 0]}t^{{2}}+{Y[3, 0]}t^{{3}}'
        graph = f'\\left({xs},{ys}\\right)'
        
        data = (f'{graph}\n{vector1}\n{vector2}\n{pts}'
                .replace('+-', '-')
                .replace('-+', '-'))
        print(data)
    
    return X, Y

def hermiteCubicSplineSegment3D(pt1: Point3DLike,
                                pt2: Point3DLike,
                                vec1: Point3DLike,
                                vec2: Point3DLike,
                                _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        vectors = (f'Curve({vx1}t+{x1},{vy1}t+{y1},{vz1}t+{z1},t,0,1)\n' +
                   f'Curve({vx2}t+{x2},{vy2}t+{y2},{vz2}t+{z2},t,0,1)')
        
        pts = f'({{{x1},{x2}}},{{{y1},{y2}}},{{{z1},{z2}}}])'
        
        curve = (f'Curve({X[0, 0]}+{X[1, 0]}t+{X[2, 0]}t^(2)+{X[3, 0]}t^(3)' +
                 f'{Y[0, 0]}+{Y[1, 0]}t+{Y[2, 0]}t^(2)+{Y[3, 0]}t^(3)' +
                 f'{Z[0, 0]}+{Z[1, 0]}t+{Z[2, 0]}t^(2)+{Z[3, 0]}t^(3),t,0,1)')
        
        print((f'{curve}\n{vectors}\n{pts}'
               .replace('+-', '-')
               .replace('-+', '-')))
    
    return X, Y, Z

def hermiteCubicSpline2D(pts: Points2DLike,
                         vecs: Points2DLike,
                         _debug: bool = False) \
                                -> Tuple[List[Matrix], List[Matrix]]:
    if len(pts) != len(vecs):
        raise SplineException('Must have equal quantities of control points and vectors in a hermite cubic spline')
    
    Xs, Ys = [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], vecs[:-1])):
        pt2, vec2 = pts[i + 1], vecs[i + 1]
        X, Y = hermiteCubicSplineSegment2D(pt1, pt2, vec1, vec2)
        Xs.append(X)
        Ys.append(Y)
    
    # Testing code (parametric form)
    if _debug:
        vectors = '\n'.join([f'\\left({vx}t+{x},{vy}t+{y}\\right)' for (x, y), (vx, vy) in zip(pts, vecs)])
        
        xs = [str(x) for x, _y in pts]
        ys = [str(y) for _x, y in pts]
        points = f'\\left(\\left[{",".join(xs)}\\right],\\left[{",".join(ys)}\\right]\\right)'
        
        xPolynomials = [f'{x[0, 0]}+{x[1, 0]}t+{x[2, 0]}t^{{2}}+{x[3, 0]}t^{{3}}' for x in Xs]
        yPolynomials = [f'{y[0, 0]}+{y[1, 0]}t+{y[2, 0]}t^{{2}}+{y[3, 0]}t^{{3}}' for y in Ys]
        graph = '\n'.join([f'\\left({x},{y}\\right)' for x, y in zip(xPolynomials, yPolynomials)])
        
        data = (f'{graph}\n{vectors}\n{points}'
                .replace('+-', '-')
                .replace('-+', '-'))
        print(data)
    
    return Xs, Ys

def hermiteCubicSpline3D(pts: Points3DLike,
                         vecs: Points3DLike,
                         _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        vectors = '\n'.join([f'Curve({vx}t+{x},{vy}t+{y},{vz}t+{z},t,0,1)'
                             for (x, y, z), (vx, vy, vz) in zip(pts, vecs)])
        
        xs = [str(x) for x, _y, _z in pts]
        ys = [str(y) for _x, y, _z in pts]
        zs = [str(z) for _x, _y, z in pts]
        points = f'({{{",".join(xs)}}},{{{",".join(ys)}}}),{{{",".join(zs)}}})'
        
        xPolynomials = [f'{x[0, 0]}+{x[1, 0]}t+{x[2, 0]}t^(2)+{x[3, 0]}t^(3)' for x in Xs]
        yPolynomials = [f'{y[0, 0]}+{y[1, 0]}t+{y[2, 0]}t^(2)+{y[3, 0]}t^(3)' for y in Ys]
        zPolynomials = [f'{z[0, 0]}+{z[1, 0]}t+{z[2, 0]}t^(2)+{z[3, 0]}t^(3)' for z in Zs]
        curve = '\n'.join([f'Curve({x},{y},{z},t,0,1)'
                           for x, y, z in zip(xPolynomials, yPolynomials, zPolynomials)])
        
        print(f'{curve}\n{vectors}\n{points}'
              .replace('+-', '-')
              .replace('-+', '-'))
    
    return Xs, Ys, Zs

def catmullRomSplineSegment2D(pt1: Point2DLike,
                              pt2: Point2DLike,
                              normal1: Point2DLike,
                              normal2: Point2DLike,
                              _debug: bool = False) \
                                    -> Tuple[Matrix, Matrix]:
    return cardinalSplineSegment2D(pt1, pt2, normal1, normal2, _debug)

def catmullRomSplineSegment3D(pt1: Point3DLike,
                              pt2: Point3DLike,
                              normal1: Point3DLike,
                              normal2: Point3DLike,
                              _debug: bool = False) \
                                    -> Tuple[Matrix, Matrix]:
    return cardinalSplineSegment3D(pt1, pt2, normal1, normal2, _debug)

def catmullRomSpline2D(pts: Points2DLike,
                       normals: Points2DLike,
                       _debug: bool = False) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    return cardinalSpline2D(pts, normals, 0, _debug)

def catmullRomSpline3D(pts: Points3DLike,
                       normals: Points3DLike,
                       _debug: bool = False) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    return cardinalSpline3D(pts, normals, 0, _debug)

def cardinalSplineSegment2D(pt1: Point2DLike,
                            pt2: Point2DLike,
                            normal1: Point2DLike,
                            normal2: Point2DLike,
                            t: float,
                            _debug: bool = False) \
                                    -> Tuple[Matrix, Matrix]:
    p = (1 - t) / 2
    vec1 = (p * (pt2[0] - normal1[0]),
            p * (pt2[1] - normal1[1]))
    vec2 = (p * (normal2[0] - pt1[0]),
            p * (normal2[1] - pt1[1]))
    return hermiteCubicSplineSegment2D(pt1, pt2, vec1, vec2, _debug)

def cardinalSplineSegment3D(pt1: Point3DLike,
                            pt2: Point3DLike,
                            normal1: Point3DLike,
                            normal2: Point3DLike,
                            t: float,
                            _debug: bool = False) \
                                    -> Tuple[Matrix, Matrix, Matrix]:
    p = (1 - t) / 2
    vec1 = (p * (pt2[0] - normal1[0]),
            p * (pt2[1] - normal1[1]),
            p * (pt2[2] - normal1[2]))
    vec2 = (p * (normal2[0] - pt1[0]),
            p * (normal2[1] - pt1[1]),
            p * (normal2[2] - pt1[2]))
    return hermiteCubicSplineSegment3D(pt1, pt2, vec1, vec2, _debug)

def cardinalSpline2D(pts: Points2DLike,
                     normals: Points2DLike,
                     t: float,
                     _debug: bool = False) \
                            -> Tuple[List[Matrix], List[Matrix]]:
    if len(normals) != len(pts):
        raise SplineException('Cardinal splines require equal quantities of points and normals')
    
    Xs, Ys = [], []
    for i, (pt1, vec1) in enumerate(zip(pts[:-1], normals[:-1])):
        pt2, vec2 = pts[i + 1], normals[i + 1]
        X, Y = cardinalSplineSegment2D(pt1, pt2, vec1, vec2, t)
        Xs.append(X)
        Ys.append(Y)
    
    # Testing code (parametric form)
    if _debug:
        vectors = '\n'.join([f'\\left({vx}t+{x},{vy}t+{y}\\right)' for (x, y), (vx, vy) in zip(pts, normals)])
        
        xs = [str(x) for x, _y in pts]
        ys = [str(y) for _x, y in pts]
        points = f'\\left(\\left[{",".join(xs)}\\right],\\left[{",".join(ys)}\\right]\\right)'
        
        xPolynomials = [f'{x[0, 0]}+{x[1, 0]}t+{x[2, 0]}t^{{2}}+{x[3, 0]}t^{{3}}' for x in Xs]
        yPolynomials = [f'{y[0, 0]}+{y[1, 0]}t+{y[2, 0]}t^{{2}}+{y[3, 0]}t^{{3}}' for y in Ys]
        graph = '\n'.join([f'\\left({x},{y}\\right)' for x, y in zip(xPolynomials, yPolynomials)])
        
        data = (f'{graph}\n{vectors}\n{points}'
                .replace('+-', '-')
                .replace('-+', '-'))
        print(data)
    
    return Xs, Ys

def cardinalSpline3D(pts: Points3DLike,
                     normals: Points3DLike,
                     t: float,
                     _debug: bool = False) \
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
    
    # Testing code (parametric form)
    if _debug:
        vectors = '\n'.join([f'Curve({vx}t+{x},{vy}t+{y},{vz}t+{z},t,0,1)'
                             for (x, y, z), (vx, vy, vz) in zip(pts, normals)])
        
        xs = [str(x) for x, _y, _z in pts]
        ys = [str(y) for _x, y, _z in pts]
        zs = [str(z) for x, _y, z in pts]
        points = f'({{{",".join(xs)}}},{{{",".join(ys)}}},{{{",".join(zs)}}})'
        
        xPolynomials = [f'{x[0, 0]}+{x[1, 0]}t+{x[2, 0]}t^(2)+{x[3, 0]}t^(3)' for x in Xs]
        yPolynomials = [f'{y[0, 0]}+{y[1, 0]}t+{y[2, 0]}t^(2)+{y[3, 0]}t^(3)' for y in Ys]
        zPolynomials = [f'{z[0, 0]}+{z[1, 0]}t+{z[2, 0]}t^(2)+{z[3, 0]}t^(3)' for z in Zs]
        graph = '\n'.join([f'Curve({x},{y},{z},t,0,1)'
                           for x, y, z in zip(xPolynomials, yPolynomials, zPolynomials)])
        
        print(f'{graph}\n{vectors}\n{points}'
                .replace('+-', '-')
                .replace('-+', '-'))
    
    return Xs, Ys, Zs

# TODO Make private function adapters for consistency and widespread usability
# E.g: xxToPolynomial2D, xxToParametric3D

if __name__ == '__main__':
    print('Running main...')
    
    # Fixed point testing code
    '''
    from myRenderer.Point import Point, Point2D, Point3D
    
    assert tuple(Point2D(1, 3).coords) == (1, 3)
    
    assert Point2D((1, 2)) == Point2D(1, 2)
    
    assert Point2D(4, 3) * Point2D(5, 7) == Point2D(20, 21)
    
    x = Point2D(3, 4)
    x *= [5, 7]
    assert x == Point2D(15, 28)
    
    assert Point2D(7, 9) * 2 == Point2D(14, 18)
    
    assert 2 * Point2D(7, 9) == Point2D(14, 18)
    
    
    assert tuple(Point3D(1, 2, 3).coords) == (1, 2, 3)
    
    assert Point3D((1, 2, 3)) == Point3D(1, 2, 3)
    
    assert Point3D(4, 3, 1) * Point3D(5, 7, 8) == Point3D(20, 21, 8)
    
    x = Point3D(3, 4, 9)
    x *= [5, 7, 2]
    assert x == Point3D(15, 28, 18)
    
    assert Point2D(8, 7) * 3 == Point2D(24, 21)
    
    assert 3 * Point2D(8, 7) == Point2D(24, 21)
    
    assert Point2D(1, 2) != Point3D(1, 2, 0)
    
    assert Point3D(1, 2, 0) != Point2D(1, 2)
    
    
    assert tuple(Point(1, 2, 3, 4).coords) == (1, 2, 3, 4)
    
    assert Point((1, 2, 3, 4)) == Point(1, 2, 3, 4)
    
    assert Point(4, 3, 1, 9) * Point(5, 7, 8, 6) == Point(20, 21, 8, 54)
    
    x = Point(3, 4, 9, 6)
    x *= [5, 7, 2, 9]
    assert x == Point(15, 28, 18, 54)
    
    assert Point2D(19, 5) * 4 == Point2D(76, 20)
    
    assert 4 * Point2D(19, 5) == Point2D(76, 20)
    
    assert Point2D(1, 2) != Point(1, 2)
    
    assert Point(1, 2) != Point2D(1, 2)
    
    assert Point3D(1, 2, 3) != Point(1, 2, 3)
    
    assert Point(1, 2, 3) != Point3D(1, 2, 3)
    
    print('Testing successful!')
    '''
    
    # Fixed 2D spline testing code
    '''
    pts = [
        (0, 1),
        (2, 2),
        (5, 0),
        (8, 0),
        # (8, -2),
    ]
    
    out = polynomialSpline2D(pts, _debug=True)
    '''
    
    # Randomized 2D spline testing code
    '''
    size = 50
    samples = [i for i in range(-size, size + 1)]
    pts = list(zip(
        sample(samples, size),
        sample(samples, size)))
    vecs = list(zip(
        sample(samples, size),
        sample(samples, size)))
    t = random()
    # pts.sort()
    
    cardinalSpline2D(pts, vecs, t, _debug=True)
    '''
    
    # Fixed 3D spline testing code
    '''
    pts = [
        (0, 1, 3),
        (2, 2, 4),
        (5, -2, 6),
        (5, -2, -2),
    ]
    normals = [(1, 2, 0), (-1, 0, 0), (2, 1, -1), (0, -1, 3)]
    t = -0.5
    
    out = cardinalSpline3D(pts, normals, t, _debug=True)
    '''
    