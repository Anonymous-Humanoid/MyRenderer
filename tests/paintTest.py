'''
Created on Apr. 22, 2023

@author: Matthew
'''

if __name__ == '__main__':
    # Disabling pygame support message
    from dotenv import load_dotenv
    load_dotenv()

from typing import Optional, Sequence, Tuple

import pygame


RGBAOutput = Tuple[int, int, int, int]
Colour = pygame.Color | int | str | Tuple[int, int, int] | RGBAOutput | Sequence[int]
Surface = pygame.Surface
Coordinate = Tuple[float, float] | Sequence[float] | pygame.math.Vector2
Coordinates = Sequence[Coordinate]
Rect = pygame.Rect


BLACK = (0x00, 0x00, 0x00)
RED   = (0xFF, 0x00, 0x00)
BLUE  = (0x00, 0xFF, 0x00)
GREEN = (0x00, 0x00, 0xFF)
WHITE = (0xFF, 0xFF, 0xFF)


def drawCircle(surface: Surface,
               pos: Coordinate = (0, 0),
               colour: Colour = WHITE,
               radius: float = 1,
               draw_top_right: bool = False, 
               draw_top_left: bool = False, 
               draw_bottom_left: bool = False, 
               draw_bottom_right: bool = False) \
                    -> Optional[Rect]:
    if radius > 0:
        return pygame.draw.circle(surface,
                                  colour,
                                  pos,
                                  radius,
                                  draw_top_right,
                                  draw_top_left,
                                  draw_bottom_left,
                                  draw_bottom_right)

def drawLine(surface: Surface,
             start_pos: Coordinate,
             end_pos: Coordinate = (0, 0),
             colour: Colour = WHITE,
             thickness: int = 1) \
                    -> Optional[Rect]:
    if thickness > 0:
        return pygame.draw.line(surface, colour, start_pos, end_pos, thickness)

def drawLines(surface: Surface,
              pts: Coordinates = (0, 0),
              colour: Colour = WHITE,
              closed: bool = False,
              thickness: int = 1) \
                    -> Optional[Rect]:
    if thickness > 0:
        return pygame.draw.lines(surface, colour, closed, pts, thickness)

def drawGradientLine(surface: Surface,
                     start_pos: Coordinate,
                     end_pos: Coordinate = (0, 0),
                     colour: Colour = WHITE,
                     blend: int = 1) \
                            -> Rect:
    return pygame.draw.aaline(surface, colour, start_pos, end_pos, blend)

def drawGradientLines(surface: Surface,
                      pts: Coordinates,
                      colour: Colour = WHITE,
                      closed: bool = False,
                      blend: int = 1) \
                            -> Rect:
    return pygame.draw.aalines(surface, colour, closed, pts, blend)

if __name__ == '__main__':
    pygame.init()
    
    WIDTH, HEIGHT = 500, 500
    
    windowFocused = False
    mouseOnScreen = False
    mouseIsPressed = False
    lastMousePos = None
    colour = None
    thickness = 5
    
    surface = pygame.display.set_mode((WIDTH, HEIGHT))
    surface.fill(BLACK)
    clock = pygame.time.Clock()
    
    while True:
        for event in pygame.event.get():
            # print(event)
            
            if event.type == pygame.QUIT:
                raise SystemExit
            elif event.type == pygame.WINDOWENTER:
                mouseOnScreen = True
                lastMousePos = pygame.mouse.get_pos()
            elif event.type == pygame.WINDOWLEAVE:
                mouseOnScreen = False
            elif event.type == pygame.WINDOWFOCUSLOST:
                windowFocused = False
                mouseOnScreen = False
            elif event.type == pygame.WINDOWFOCUSGAINED:
                windowFocused = True
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if mouseOnScreen:
                    mouseIsPressed = True
                    if event.button == pygame.BUTTON_LEFT:
                        colour = WHITE
                        thickness = 5
                    elif event.button == pygame.BUTTON_RIGHT:
                        colour = BLACK
                        thickness = 8
                    else:
                        colour = None
            elif event.type == pygame.MOUSEBUTTONUP:
                mouseIsPressed = False
                colour = None
            elif event.type == pygame.MOUSEMOTION:
                if mouseOnScreen:
                    if mouseIsPressed and colour != None:
                        pos = pygame.mouse.get_pos()
                        drawLine(surface, lastMousePos, pos, colour, thickness)
            
            if mouseOnScreen:
                lastMousePos = pygame.mouse.get_pos()
        pygame.display.update()
        clock.tick(60)
    