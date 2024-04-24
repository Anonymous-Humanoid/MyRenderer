'''
Created on Oct. 20, 2023

@author: Matthew
'''

from . import SplineException
from numbers import Number
from typing import Iterator, overload, List, Self, Sequence, SupportsIndex, Tuple

# TODO Integrate Point in testing code
class Point:
    @overload
    def __init__(self, pt: SupportsIndex):
        self.coords: List[Number] = [c for c in pt]
    
    @overload
    def __init__(self, *coords: Number):
        self.coords = [c for c in coords]
    
    def __init__(self, *args: Number | SupportsIndex):
        if len(args) == 1:
            self.coords: List[Number] = [c for c in args[0]]
        elif len(args) > 1:
            self.coords: List[Number] = [c for c in args]
        else:
            raise SplineException('Expected at least 1 argument')
    
    def __mul__(self, pt: Number | SupportsIndex) -> Self:
        if isinstance(pt, Number):
            return Point((c * pt for c in self.coords))
        return Point((c * p for c, p in zip(self.coords, pt)))
    
    def __rmul__(self, pt: Number | SupportsIndex) -> Self:
        return self.__mul__(pt)
    
    def __imul__(self, pt: Number | SupportsIndex) -> Self:
        if isinstance(pt, Number):
            for i in range(len(self.coords)):
                self.coords[i] *= pt
        elif len(self) != len(pt):
            raise SplineException(f'Expected an iterable of exactly {len(self)} numbers')
        else:
            for i in range(len(self.coords)):
                self.coords[i] *= pt[i]
        return self
    
    def __index__(self, i: int) -> Number:
        return self.coords[i]
    
    def __iter__(self) -> Iterator:
        return iter(self.coords)
    
    def __len__(self) -> int:
        return len(self.coords)
    
    def __eq__(self, obj) -> bool:
        if not isinstance(obj, Point):
            return False
        elif self is obj or super.__eq__(self, obj):
            return True
        for c, p in zip(self.coords, obj):
            if c != p:
                return False
        return True
    
    def __str__(self) -> str:
        return str(tuple(self))
    
    def __repr__(self) -> str:
        return f'Point({self.coords})'

class Point2D(Point):
    @property
    def coords(self):
        return self.__coords
    
    @coords.setter
    def coords(self, value: SupportsIndex):
        if len(value) != 2:
            raise SplineException('Expected an iterable containing 2 numbers')
        self.__coords = [c for c in value]
    
    @overload
    def __init__(self, pt: SupportsIndex):
        self.coords = pt
    
    @overload
    def __init__(self, x: Number, y: Number):
        self.__coords = [x, y]
    
    def __init__(self, *args: Number | SupportsIndex):
        if len(args) == 1:
            if len(args[0]) != 2:
                raise SplineException('Expected argument to contain 2 numbers')
        elif len(args) != 2:
            raise SplineException('Expected either 2 numbers or an iterable containing 2 numbers')
        self.__coords = []
        super().__init__(*args)
    
    @property
    def x(self):
        return self.__coords[0]
    
    @x.setter
    def x(self, value: Number):
        self.__coords[0] = value
    
    @property
    def y(self):
        return self.__coords[1]
    
    @y.setter
    def y(self, value: Number):
        self.__coords[1] = value
    
    def __mul__(self, pt: Number | SupportsIndex) -> Self:
        if isinstance(pt, Number):
            return Point2D(self.x * pt, self.y * pt)
        return Point2D(self.x * pt.x, self.y * pt.y)
    
    def __eq__(self, obj) -> bool:
        if not isinstance(obj, Point2D):
            return False
        elif self is obj:
            return True
        return self.x == obj.x and self.y == obj.y
    
    def __repr__(self) -> str:
        return f'Point2D({self.__coords})'

class Point3D(Point):
    @property
    def coords(self):
        return self.__coords
    
    @coords.setter
    def coords(self, value: SupportsIndex):
        if len(value) != 3:
            raise SplineException('Expected an iterable containing 3 numbers')
        self.__coords = [c for c in value]
    
    @overload
    def __init__(self, pt: SupportsIndex):
        self.coords = pt
    
    @overload
    def __init__(self, x: Number, y: Number, z: Number):
        self.__coords = [x, y, z]
    
    def __init__(self, *args: Number | SupportsIndex):
        if len(args) == 1:
            if len(args[0]) != 3:
                raise SplineException('Expected argument to contain 3 numbers')
        elif len(args) != 3:
            raise SplineException('Expected either 3 numbers or an iterable containing 3 numbers')
        self.__coords = []
        super().__init__(*args)
    
    @property
    def x(self):
        return self.__coords[0]
    
    @x.setter
    def x(self, value: Number):
        self.__coords[0] = value
    
    @property
    def y(self):
        return self.__coords[1]
    
    @y.setter
    def y(self, value: Number):
        self.__coords[1] = value
    
    @property
    def z(self):
        return self.__coords[2]
    
    @z.setter
    def z(self, value: Number):
        self.__coords[2] = value
    
    def __mul__(self, pt: Number | SupportsIndex) -> Self:
        if isinstance(pt, Number):
            return Point3D(self.x * pt, self.y * pt, self.z * pt)
        return Point3D(self.x * pt.x, self.y * pt.y, self.z * pt.z)
    
    def __eq__(self, obj) -> bool:
        if not isinstance(obj, Point3D):
            return False
        elif self is obj:
            return True
        return self.x == obj.x and self.y == obj.y and self.z == obj.z
    
    def __repr__(self) -> str:
        return f'Point3D({self.__coords})'

PointLike = Tuple[Number, ...]
Point2DLike = Tuple[Number, Number]
Point3DLike = Tuple[Number, Number, Number]
PointsLike = Sequence[Point]
Points2DLike = Sequence[Point2D]
Points3DLike = Sequence[Point3D]
