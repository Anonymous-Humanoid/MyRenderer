from numpy import mat as Matrix

# Order matters!
from .SplineException import SplineException
from .KnotVector import KnotVector
from .Point import Point2D, Point3D, Point, PointLike, Point2DLike, Point3DLike, PointsLike, Points2DLike, Points3DLike

from .splines import polynomialSpline2D, polynomialSpline3D, piecewisePolynomialSpline2D, piecewisePolynomialSpline3D, parameterizedSpline2D, parameterizedSpline3D, spaceSpline2D, spaceSpline3D, hermiteCubicSplineSegment2D, hermiteCubicSplineSegment3D, hermiteCubicSpline2D, hermiteCubicSpline3D, catmullRomSplineSegment2D, catmullRomSplineSegment3D, catmullRomSpline2D, catmullRomSpline3D, cardinalSplineSegment2D, cardinalSplineSegment3D, cardinalSpline2D, cardinalSpline3D, bSplineBasis, bSpline, bSplineSurfaceBasis, bSplineSurface, nurbsBasis, nurbs, nurbsSurfaceBasis, nurbsSurface, tSplineBasis, tSpline


# from splineTests import bezierCurve, bezierSpline, bezierSurfaceBasis, bezierSurface
