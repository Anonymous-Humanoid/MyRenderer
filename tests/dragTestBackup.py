'''
Created on Nov. 18, 2023

@author: Matthew
'''

# Newer approximateSpline backup
'''
def approximateSpline(U: KnotVector,
                        control_pts: List[Tuple[float, float]],
                        lower_knot: float = 0,
                        upper_knot: float = 1) \
                            -> Optional[Tuple
                                List[List[Tuple[float, float]]],
                                KnotVector]]:
    n = len(control_pts) - 1

    if n < 2:
        return None

    # Approximating spline sample size needed
    spline_strips = [[]]
    sample = sample_rate * n
    
    # Approximating spline
    if sample > 0:
        sample = 1 / sample
        
        u = lower_knot + sample
        i = U.getKnotIntervalIndices(u)
        if i is not None:
            i = i[1]
            
        while round(u, float_precision) < upper_knot:
            pt = bSpline(control_pts, U, u)
            spline_strips[-1].append(pt)
            u += sample
            j = U.getKnotIntervalIndices(u + sample)
            if j is not None and i < j[1]:
                i += 1
                spline_strips.append([])
    
    return spline_strips, U
'''

if __name__ == '__main__':
    # Disabling pygame support message
    from dotenv import load_dotenv
    load_dotenv()

from math import dist
from typing import List, Optional, Tuple

import pygame
from myRenderer import KnotVector, SplineException, bSpline

BLACK = (0x00, 0x00, 0x00)
WHITE = (0xFF, 0xFF, 0xFF)
RED   = (0xFF, 0x00, 0x00)
GREEN = (0x00, 0xFF, 0x00)
BLUE  = (0x00, 0x00, 0xFF)
    
# No UI toggle for these options
fps = 30
controlPtRadius = 15
convexHullThickness = 4
splineThickness = 8
maxDegree = 3
closed = False # TODO Left click a point to change the closure?


# TODO Refactor to allow approximating spline segments
def approximateSpline(controlPts: List[Tuple[float, float]]) \
            -> Optional[List[List[Tuple[float, float]]]]:
    if len(rects) <= 1:
        return None
    
    # Approximating spline sample size needed
    splineStrips = [[]]
    sample = 0
    n = len(controlPts) - 1
    
    for i in range(n):
        pt1, pt2 = controlPts[i], controlPts[i + 1]
        sample += dist(pt1, pt2)
    
    if closed and n > 1:
        sample += dist(controlPts[0], controlPts[-1])
    
    # Approximating spline
    if sample > 0:
        sample = 10 / sample
        
        if closed and n > 1:
            # Bridging spline endpoints with C1 continuity
            midpoint = ((controlPts[0][0] + controlPts[-1][0]) / 2,
                        (controlPts[0][1] + controlPts[-1][1]) / 2)
            controlPts.insert(0, midpoint)
            controlPts.append(midpoint)
            n += 2

        # Intersecting spline end control points
        p = min(n, maxDegree)
        A = min(n + 1, maxDegree)
        B = 2 * A - p - 1
        multiplicities = {i: 1 for i in range(n - B + 1)}
        multiplicities[0] = multiplicities[n] = A
        
        U = KnotVector(multiplicities)
        U.normalizeAll()
        
        # TODO Bounding spline
        if closed and n > 1:
            # lowerBound, upperBound = U[p], U[n - p]
            # lowerBound, upperBound = U[0], U[n]
            lowerBound, upperBound = U[0], U[-1]
            # lowerBound, upperBound = U[1], U[-2]
        else:
            lowerBound, upperBound = U[0], U[-1]
        # lowerBound, upperBound = U[p], U[n - p]
        
        u = lowerBound + sample
        i = U.getKnotIntervalIndices(u)
        if i is not None:
            i = i[1]
            
        while u <= upperBound:
            pt = bSpline(controlPts, U, u)
            splineStrips[-1].append(pt)
            u += sample
            j = U.getKnotIntervalIndices(u + sample)
            if j is not None and i < j[1]:
                i += 1
                splineStrips.append([])

    return splineStrips

def repaint(screen: pygame.surface.Surface,
            splineStrips: List[List[Tuple[float, float]]]) \
                    -> None:
    # Repainting background
    screen.fill(BLACK)
    
    if splineStrips is not None:
        for i, strip in enumerate(splineStrips):
            for j in range(len(strip) - 1):
                pt1, pt2 = strip[j], strip[j + 1]
                pygame.draw.line(screen, BLUE, pt1, pt2, splineThickness)
            if i < len(splineStrips) - 1:
                nextStrip = splineStrips[i + 1]
                pygame.draw.line(screen, BLUE, strip[-1],
                                 nextStrip[0], splineThickness)
    
    # Filling visual gap
    if closed and len(rects) > 2:
        pygame.draw.line(screen, BLUE, splineStrips[0][0],
                         splineStrips[-1][-1], splineThickness)
    
    # Painting convex hull
    for i in range(len(rects) - 1):
        rectangle1, rectangle2 = rects[i], rects[i + 1]
        pygame.draw.line(screen, RED, rectangle1.center,
                         rectangle2.center, convexHullThickness)
    
    # Painting closed convex hull
    if closed and len(rects) > 1:
        pygame.draw.line(screen, GREEN, rects[0].center,
                        rects[-1].center, convexHullThickness)
        
    
    # Painting control points
    for rectangle in rects:
        pygame.draw.circle(screen, WHITE, rectangle.center,
                           dist(rectangle.center, rectangle.bottomright))
    
    # TODO Only draw and update as necessary
    pygame.display.update()
    clock.tick(fps)


if __name__ == '__main__':
    pygame.init()
    
    screen = pygame.display.set_mode((430, 410), pygame.RESIZABLE)
    # screen_rect = screen.get_rect()
    pygame.display.set_caption("Drag Test")
    
    rects = []
    rectangle_dragging = None
    rectangle_created = False
    clock = pygame.time.Clock()
    running = True
    
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == pygame.BUTTON_LEFT:
                    # Selecting rect for dragging
                    for rectangle in rects[::-1]:
                        if rectangle.collidepoint(event.pos):
                            rectangle_dragging = rectangle
                            mouse_x, mouse_y = event.pos
                            offset_x = rectangle.x - mouse_x
                            offset_y = rectangle.y - mouse_y
                            break
                
                # Creating draggable rect
                if event.button == pygame.BUTTON_RIGHT \
                            and not rectangle_created:
                    mouse_x, mouse_y = event.pos
                    rectangle = pygame.rect.Rect(mouse_x, mouse_y,
                                                    controlPtRadius,
                                                    controlPtRadius)
                    rects.append(rectangle)
                    rectangle_created = True
            
            elif event.type == pygame.MOUSEBUTTONUP:
                # Deselect rect for dragging
                if event.button == pygame.BUTTON_LEFT:
                    rectangle_dragging = None
                
                # Preventing spammed rect creation
                elif event.button == pygame.BUTTON_RIGHT:
                    rectangle_created = False
            
            elif event.type == pygame.MOUSEMOTION:
                # Dragging rect
                if rectangle_dragging is not None:
                    mouse_x, mouse_y = event.pos
                    rectangle_dragging.x = mouse_x + offset_x
                    rectangle_dragging.y = mouse_y + offset_y
        
        controlPts = [rect.center for rect in rects]
        splineStrips = approximateSpline(controlPts)
        
        # Repainting spline
        repaint(screen, splineStrips)
