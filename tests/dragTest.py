'''
Created on Nov. 18, 2023

@author: Matthew
'''

'''Credit to:
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-closed.html
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-mv-ctlpt.html
'''

if __name__ == '__main__':
    # Disabling pygame support message
    from dotenv import load_dotenv
    load_dotenv()

import random
from math import dist
from typing import Any, Iterator, List, Optional, Self, Sequence, Tuple, TypeVar
import pygame
from myRenderer import KnotVector, bSpline, bSplineBasis, Point2DLike

# TODO Testing imports
import numpy as np
import logging
import sys

rng = np.random.default_rng(random.seed())
logger = logging.getLogger('dragTest')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(name)s:%(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def log(msg: str = '', *args, **kwargs) -> None:
    logger.log(logging.DEBUG, msg, *args, **kwargs)

Colour = pygame.Color | int | str | Tuple[int, int, int] | Tuple[int, int, int, int] | Sequence[int]
Segment = List[Point2DLike]
T = TypeVar('T')

RED = pygame.color.THECOLORS['red']
VIOLET = pygame.color.THECOLORS['violet']
BLUE = pygame.color.THECOLORS['blue']
TURQUOISE = pygame.color.THECOLORS['turquoise']
WHITE = pygame.color.THECOLORS['white']
GREEN = pygame.color.THECOLORS['green']
BLACK = pygame.color.THECOLORS['black']

# TODO Refactor KnotVector to store non-normalized values and normalize them just-in-time

# TODO Write documentation for all methods
class InteractiveBSpline():
    def __init__(self: Self,
                 screen: pygame.surface.Surface,
                 max_degree: int = 3,
                 sample_rate: int = 10) \
                    -> None:
        self._screen = screen
        self._control_pt_rects: List[pygame.Rect] = []
        self._control_pts: List[Point2DLike] = []
        self._spline_strips: List[Segment] = []
        
        # No UI toggle for these options
        self._max_degree = max_degree
        self._U = KnotVector(first_outer_knot_mul=1, last_outer_knot_mul=1)
        self._sample_rate = sample_rate # Points per spline segment approximation
        self.control_pt_radii = 15
        self.convex_hull_thickness = 4
        self.spline_thickness = 8
        self.first_control_pt_colour = TURQUOISE
        self.middle_control_pts_colour = WHITE
        self.last_control_pt_colour = VIOLET
        self.convex_hull_colour = RED
        self.spline_colour = BLUE
        self.spline_segment_join_colour = GREEN # TODO BLUE
        self.background_colour = BLACK
    
    def _n(self: Self) -> int:
        return len(self._control_pts) - 1
    
    def _p(self: Self, n: int) -> int:
        return min(n, self._max_degree)
    
    @staticmethod
    def _reversedEnumerate(collection: List[T]) -> Iterator[Tuple[int, T]]:
        for i in range(len(collection) - 1, -1, -1):
            yield i, collection[i]
    
    def _setKnotVector(self: Self) -> None:
        n = self._n()
        old_p, p = self._p(n - 1), self._p(n)
        
        if old_p > 1:
            # Incrementing degree to clamp spline endpoints
            self._U.first_outer_knot_mul += 1
            self._U.last_outer_knot_mul += 1
        
        # Defining open uniform spline
        self._U.insertKnot(n if n - p < 1 else n - p)
    
    def _repaint(self: Self) -> None:
        # Repainting background
        self._screen.fill(self.background_colour) # TODO Repaint specified slices of spline_strips
        
        if self._spline_strips is not None:
            n = self._n()
            for i, strip in enumerate(self._spline_strips):
                # Drawing spline
                for j in range(len(strip) - 1):
                    pt1, pt2 = strip[j], strip[j + 1]
                    pygame.draw.line(self._screen, self.spline_colour,
                                     pt1, pt2, self.spline_thickness)
                
                # Joining spline segments
                if i < len(self._spline_strips) - 1:
                    next_strip = self._spline_strips[i + 1]
                    pygame.draw.line(self._screen, self.spline_segment_join_colour,
                                     strip[-1], next_strip[0], self.spline_thickness)
         
        # Painting convex hull
        for i in range(n):
            rectangle1, rectangle2 = self._control_pt_rects[i], self._control_pt_rects[i + 1]
            pygame.draw.line(self._screen, self.convex_hull_colour, rectangle1.center,
                            rectangle2.center, self.convex_hull_thickness)
        
        # Painting first control point
        rectangle = self._control_pt_rects[0]
        pygame.draw.circle(self._screen, self.first_control_pt_colour, rectangle.center,
                           dist(rectangle.center, rectangle.bottomright))
        
        if n > 0:
            # Painting middle control points
            for rectangle in self._control_pt_rects[1:-1]:
                pygame.draw.circle(self._screen, self.middle_control_pts_colour, rectangle.center,
                                dist(rectangle.center, rectangle.bottomright))
            
            # Painting last control point
            rectangle = self._control_pt_rects[-1]
            pygame.draw.circle(self._screen, self.last_control_pt_colour, rectangle.center,
                               dist(rectangle.center, rectangle.bottomright))
        
        pygame.display.flip()
    
    def _repaintControlPt(self: Self, control_pt_index: int, new_control_pt: Tuple[float, float]) -> None:
        n = self._n()
        sample = self._sample_rate * n
        
        if n < 1:
            # Updating only control point
            self._control_pts[control_pt_index] = new_control_pt
        elif sample > 0:
            # Approximating spline sample size needed
            sample = 1 / sample
            
            # Calculating control point shift
            v = (new_control_pt[0] - self._control_pts[control_pt_index][0],
                 new_control_pt[1] - self._control_pts[control_pt_index][1])
            self._control_pts[control_pt_index] = new_control_pt
            
            # Defining affected segments
            p = self._p(n)
            lower_index = max(0, control_pt_index - p - 1)
            upper_index = min(n, control_pt_index + p + 2)
            #u = self._U[lower_index + p]
            #u = max(self._U[0], self._U[lower_index + 1] - (upper_index - lower_index) / n)
            #u = max(self._U[0], self._U[lower_index + p])
            u = lower_index / n
            
            # affected_strips = []
            
            # Approximating spline
            for i in range(lower_index, upper_index):
                # affected_strips.append([])
                for j in range(self._sample_rate):
                    # affected_strips[-1].append(u)
                    basis = bSplineBasis(control_pt_index, p, self._U, u)
                    self._spline_strips[i][j] = (self._spline_strips[i][j][0] + basis * v[0],
                                                 self._spline_strips[i][j][1] + basis * v[1])
                    u += sample
            '''
            for spline_strip in self._spline_strips[lower_index:upper_index]:
                for j in range(len(spline_strip)):
                    basis = bSplineBasis(control_pt_index, p, self._U, u)
                    spline_strip[j] = (spline_strip[j][0] + basis * v[0],
                                       spline_strip[j][1] + basis * v[1])
                    u += sample
            '''
            
            # # TODO TESTING
            # DECIMAL_PRECISION = 8
            # all_strips = [
            #     [
            #         sample * (self._sample_rate * i + j)
            #         for j in range(self._sample_rate)
            #     ]
            #     for i in range(n)
            # ]
            # for i, confirmed_strip in enumerate(all_strips[lower_index: upper_index]):
            #     unconfirmed_strip = affected_strips[i]
            #     for confirmed_knot, unconfirmed_knot in zip(confirmed_strip, unconfirmed_strip):
            #         log(f'{control_pt_index} | {unconfirmed_knot} == {confirmed_knot}')
            #         assert round(unconfirmed_knot, DECIMAL_PRECISION) == round(confirmed_knot, DECIMAL_PRECISION)
            # # log(f'-1 | {all_strips[-1][-1]} == {1 - sample}')
            # assert round(all_strips[-1][-1], DECIMAL_PRECISION) == round(1 - sample, DECIMAL_PRECISION)
            
            # # TODO TESTING
            # log(f'Start u: {lower_index / n} - End u: {u} - Expected end u: {upper_index / n}')
            # log(f'Current u: {control_pt_index / n} - End u: {u} - Next u: {min(1, (control_pt_index + 1) / n)}')
            # log(f'End u: {u} - Expected: {min(1, (control_pt_index + p) / n)}')
        
        self._repaint()
    
    # TODO Fix cragginess:
    # * The segment join lines are a result of the issue
    # * Last p control points work fine
    # * Is an issue at least partially caused by _paintNewControlPt:
    #    - Setting lower_index = 0 fixes the issue, at a performance loss;
    #      note that this causes the join segments to change position
    #      as control points are added
    # 
    # - Not dependent on self._sample_rate
    # - Not dependent on self._max_degree
    # - Not a floating-point issue:
    #   - Increased self._sample_rate
    #   - Rounded u to 5 and 8 decimals
    # - v should not be recalculated for each segment
    # - Bases do need to be recalculated for each point
    # - Isn't an issue with u or sample_rate:
    #    - Ran continuity tests
    #    - Ran consistency tests
    # - Doesn't seem to be an issue with _repaintControlPt:
    #    - The first max(0, n - p - 2) join lines are wrongly stretched
    #      given n control points, which can be noticed as soon as a
    #      new control point is added (See: B-Spline Expected vs Actual.png)
    # - Knot values should be in the range [0, 1) and should not include 1:
    #    - Consistency tests fail
    #    - Lines connected to the origin may appear
    def _paintNewControlPt(self: Self, new_control_pt: Tuple[float, float]) -> None:
        self._control_pts.append(new_control_pt)
        self._setKnotVector()
        n = self._n()
        sample = self._sample_rate * n
        
        if sample > 0:
            # Approximating spline sample size needed
            sample = 1 / sample
            self._spline_strips.append([])
            
            # TODO Defining affected segments
            p = self._p(n)
            lower_index = max(0, n - p - 1)
            # lower_index = 0
            #lower_index = max(0, n - p - 2)
            #lower_index = max(0, n - 1)
            #upper_index = n
            #u = self._U[lower_index + p]
            #u = max(self._U[0], self._U[lower_index + 1] - (upper_index - lower_index) / n)
            #u = max(self._U[0], self._U[lower_index + p])
            u = lower_index / n
            
            affected_strips = []
            
            # Approximating spline
            for i in range(lower_index, n):
                affected_strips.append([])
                self._spline_strips[i].clear()
                for _ in range(self._sample_rate):
                    affected_strips[-1].append(u)
                    pt = bSpline(self._control_pts, self._U, u, self._max_degree)
                    self._spline_strips[i].append(pt)
                    u += sample
            
            # # TODO TESTING
            # DECIMAL_PRECISION = 8
            # all_strips = [
            #     [
            #         sample * (self._sample_rate * i + j)
            #         for j in range(self._sample_rate)
            #     ]
            #     for i in range(n)
            # ]
            # for i, confirmed_strip in enumerate(all_strips[lower_index:]):
            #     unconfirmed_strip = affected_strips[i]
            #     for confirmed_knot, unconfirmed_knot in zip(confirmed_strip, unconfirmed_strip):
            #         log(f'{n} | {unconfirmed_knot} == {confirmed_knot}')
            #         assert round(unconfirmed_knot, DECIMAL_PRECISION) == round(confirmed_knot, DECIMAL_PRECISION)
            # # log(f'-1 | {all_strips[-1][-1]} == {1 - sample}')
            # assert round(all_strips[-1][-1], DECIMAL_PRECISION) == round(1 - sample, DECIMAL_PRECISION)
            
            # # TODO TESTING
            # log(f'Start u: {lower_index / n} - End u: {u} - Expected end u: 1')
        
        self._repaint()
    
    @staticmethod
    def leftMouseDownHandler(renderer: 'InteractiveBSpline', pos: Tuple[int, int]) -> Optional[Tuple[int, int, int]]:
        # Selecting rect for dragging
        for i, rectangle in InteractiveBSpline._reversedEnumerate(renderer._control_pt_rects):
            if rectangle.collidepoint(pos):
                mouse_x, mouse_y = pos
                return i, rectangle.x - mouse_x, rectangle.y - mouse_y
    
    @staticmethod
    def rightMouseDownHandler(renderer: 'InteractiveBSpline', pos: Tuple[int, int]) -> None:
        # Creating draggable rect
        mouse_x, mouse_y = pos
        rectangle = pygame.rect.Rect(mouse_x, mouse_y,
                                     renderer.control_pt_radii,
                                     renderer.control_pt_radii)
        renderer._control_pt_rects.append(rectangle)
        renderer._paintNewControlPt(rectangle.center)
    
    @staticmethod
    def run(screen: pygame.surface.Surface, fps: int = 30) -> None:
        renderer = InteractiveBSpline(screen)
        clock = pygame.time.Clock()
        running = True
        dragged_rectangle_index = None
        rectangle_created = False
        offset_x, offset_y = 0, 0 # Smooths initial dragging
        
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
                elif event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == pygame.BUTTON_LEFT:
                        temp = InteractiveBSpline.leftMouseDownHandler(renderer, event.pos)
                        if temp is not None:
                            dragged_rectangle_index, offset_x, offset_y = temp
                    if event.button == pygame.BUTTON_RIGHT and not rectangle_created:
                        InteractiveBSpline.rightMouseDownHandler(renderer, event.pos)
                        rectangle_created = True
                elif event.type == pygame.MOUSEBUTTONUP:
                    if event.button == pygame.BUTTON_LEFT:
                        # Deselecting rect for dragging
                        dragged_rectangle_index = None
                    elif event.button == pygame.BUTTON_RIGHT:
                        # Preventing spammed rect creation
                        rectangle_created = False
                elif event.type == pygame.MOUSEMOTION and dragged_rectangle_index is not None:
                    # Dragging rect
                    mouse_x, mouse_y = event.pos
                    dragged_rectangle = renderer._control_pt_rects[dragged_rectangle_index]
                    dragged_rectangle.x = mouse_x + offset_x
                    dragged_rectangle.y = mouse_y + offset_y
                    renderer._repaintControlPt(dragged_rectangle_index, dragged_rectangle.center)
            clock.tick(fps)

if __name__ == '__main__':
    pygame.init()
    screen = pygame.display.set_mode((430, 410), pygame.RESIZABLE)
    pygame.display.set_caption("Interactive B-spline Test")
    InteractiveBSpline.run(screen)
    pygame.quit()