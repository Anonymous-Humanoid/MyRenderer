'''
Created on Oct. 20, 2023

@author: Matthew
'''

from myRenderer.Point import Point, Point2D, Point3D

if __name__ == '__main__':
    assert tuple(Point2D(1, 3).coords) == (1, 3)
    
    assert Point2D((1, 2)) == Point2D(1, 2)
    
    assert Point2D(4, 3) * Point2D(5, 7) == Point2D(20, 21)
    
    x = Point2D(3, 4)
    x *= [5, 7]
    assert x == Point2D(15, 28)
    
    
    assert tuple(Point3D(1, 2, 3).coords) == (1, 2, 3)
    
    assert Point3D((1, 2, 3)) == Point3D(1, 2, 3)
    
    assert Point3D(4, 3, 1) * Point3D(5, 7, 8) == Point3D(20, 21, 8)
    
    x = Point3D(3, 4, 9)
    x *= [5, 7, 2]
    assert x == Point3D(15, 28, 18)
    
    assert Point2D(1, 2) != Point3D(1, 2, 0)
    
    
    assert tuple(Point(1, 2, 3, 4).coords) == (1, 2, 3, 4)
    
    assert Point((1, 2, 3, 4)) == Point(1, 2, 3, 4)
    
    assert Point(4, 3, 1, 9) * Point(5, 7, 8, 6) == Point(20, 21, 8, 54)
    
    x = Point(3, 4, 9, 6)
    x *= [5, 7, 2, 9]
    assert x == Point(15, 28, 18, 54)
    
    assert Point2D(1, 2) != Point(1, 2)
    
    assert Point3D(1, 2, 3) != Point(1, 2, 3)
    
    print('Testing successful!')
    