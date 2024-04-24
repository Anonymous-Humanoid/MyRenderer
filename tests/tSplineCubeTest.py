'''
import numpy as np

# Define the 3D coordinates of the vertices, edge points, face points, and center point
vertices = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1]
])

edge_points = np.array([
    [0.5, 0, 0],
    [1, 0.5, 0],
    [0.5, 1, 0],
    [0, 0.5, 0],
    [0.5, 0, 1],
    [1, 0.5, 1],
    [0.5, 1, 1],
    [0, 0.5, 1],
    [0, 0, 0.5],
    [1, 0, 0.5],
    [1, 1, 0.5],
    [0, 1, 0.5]
])

face_points = np.array([
    [0.5, 0.5, 0],
    [1, 0.5, 0.5],
    [0.5, 0.5, 1],
    [0, 0.5, 0.5],
    [0.5, 0, 0.5],
    [0.5, 1, 0.5]
])

center_point = np.array([0.5, 0.5, 0.5])

# Combine all control points
control_points = np.concatenate([vertices, edge_points, face_points, [center_point]])

# Define the knot vectors for each edge and face (example, you need to adjust based on your needs)
knots_edges = np.linspace(0, 1, len(edge_points) + 2)[1:-1]
knots_faces_u = np.linspace(0, 1, len(face_points) + 2)[1:-1]
knots_faces_v = np.linspace(0, 1, len(face_points) + 2)[1:-1]

# Example: For simplicity, assuming each edge has the same knot vector
knots_for_edges = [knots_edges] * len(edge_points)
# Example: For simplicity, assuming each face has the same knot vectors in both u and v directions
knots_for_faces_u = [knots_faces_u] * len(face_points)
knots_for_faces_v = [knots_faces_v] * len(face_points)

# Combine all knot vectors
knots = np.concatenate([knots_for_edges, knots_for_faces_u, knots_for_faces_v])

# Convert knots to a list of lists for better readability
knots = knots.tolist()

# Now, control_points and knots can be used in your T-spline representation
'''

from myRenderer import KnotVector, tSpline

if __name__ == '__main__':
    '''
    # Define the 3D coordinates of the vertices, edge points, face points, and center point
    vertices = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1]
    ]

    edge_points = [
        [0.5, 0, 0],
        [1, 0.5, 0],
        [0.5, 1, 0],
        [0, 0.5, 0],
        [0.5, 0, 1],
        [1, 0.5, 1],
        [0.5, 1, 1],
        [0, 0.5, 1],
        [0, 0, 0.5],
        [1, 0, 0.5],
        [1, 1, 0.5],
        [0, 1, 0.5]
    ]

    face_points = [
        [0.5, 0.5, 0],
        [1, 0.5, 0.5],
        [0.5, 0.5, 1],
        [0, 0.5, 0.5],
        [0.5, 0, 0.5],
        [0.5, 1, 0.5]
    ]

    center_point = [0.5, 0.5, 0.5]

    # Combine all control points
    control_points = []
    control_points.extend(vertices)
    control_points.extend(edge_points)
    control_points.extend(face_points)
    # control_points.append(center_point)

    # TODO Equal knot vectors generates a point: draw the knot vector diagrams and fix

    # Define the knot vectors for each edge and face (example, you need to adjust based on your needs)
    knots_for_vertices = [KnotVector({i: 1 for i in range(len(vertices))}) for _ in range(len(vertices))]

    # Example: For simplicity, assuming each edge has the same knot vector
    knots_for_edges_u = [KnotVector({i: 1 for i in range(len(edge_points))}) for _ in range(len(edge_points))]
    knots_for_edges_v = [KnotVector({i: 1 for i in range(len(edge_points))}) for _ in range(len(edge_points))]

    # Example: For simplicity, assuming each face has the same knot vectors in both u and v directions
    knots_for_faces_u = [KnotVector({i: 1 for i in range(len(face_points))}) for _ in range(len(face_points))]
    knots_for_faces_v = [KnotVector({i: 1 for i in range(len(face_points))}) for _ in range(len(face_points))]
    '''

    # Define the 3D coordinates of the vertices, edge points, face points, and center point
    vertices = [
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, 0),
        (0, 1, 1),
        (1, 0, 0),
        (1, 0, 1),
        (1, 1, 0),
        (1, 1, 1)
    ]

    edge_points = [
        (-1, 0, 1),
        (-1, 1, 0),
        (0, 1, -1),
        (0, -1, 1),
        (-1, -1, 0),
        (0, -1, -1),
        (-1, 0, -1),
        (1, 0, 1),
        (0, 1, 1),
        (1, 1, 0),
        (1, -1, 0),
        (1, 0, -1),
    ]

    face_points = [
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
        (0, 1, 0),
        (1, 0, 0)
    ]

    # Combine all control points
    control_points = []
    control_points.extend(vertices)
    control_points.extend(edge_points)
    control_points.extend(face_points)

    knots_for_vertices_u = [KnotVector({i: 1 for i in range(len(vertices))}) for _ in range(len(vertices))]
    knots_for_vertices_v = [KnotVector({i: 1 for i in range(len(vertices))}) for _ in range(len(vertices))]

    knots_for_edges_u = [KnotVector({i: 1 for i in range(len(edge_points))}) for _ in range(len(edge_points))]
    knots_for_edges_v = [KnotVector({i: 1 for i in range(len(edge_points))}) for _ in range(len(edge_points))]

    knots_for_faces_u = [KnotVector({i: 1 for i in range(len(face_points))}) for _ in range(len(face_points))]
    knots_for_faces_v = [KnotVector({i: 1 for i in range(len(face_points))}) for _ in range(len(face_points))]

    # Combine all knot vectors
    knots_u = []
    knots_u.extend(knots_for_vertices_u)
    knots_u.extend(knots_for_edges_u)
    knots_u.extend(knots_for_faces_u)

    knots_v = []
    knots_v.extend(knots_for_vertices_v)
    knots_v.extend(knots_for_edges_v)
    knots_v.extend(knots_for_faces_v)

    # Convert knots to a list of lists for better readability

    for T in knots_u:
        T.normalizeAll()
    for T in knots_v:
        T.normalizeAll()

    # Now, control_points and knots can be used in your T-spline representation

    stepU = 1 / 15
    stepV = 1 / 15
    
    interpolation = []
    knot0, knotn = 0, 1
    
    u = knot0 + stepU
    while u < knotn:
        print(f'{100*u:.2f}%')
        v = knot0 + stepV
        while v < knotn:
            pt = tSpline(control_points, knots_u, knots_v, u, v)
            interpolation.append(f'({pt[0] :f},{pt[1] :f},{pt[2] :f})')
            v += stepV
        u += stepU
    
    interpolationStr = '{' + (','.join(interpolation)) + '}'
    controlPtsStr = '{' + (','.join([f'({pt[0] :f},{pt[1] :f},{pt[2] :f})' for pt in control_points])) + '}'
    UsStr = '[' + (', '.join([str(U) for U in knots_u])) + ']'
    VsStr = '[' + (', '.join([str(V) for V in knots_v])) + ']'
    
    print(f'\nInterpolation:\n{interpolationStr}\n\nKnot Vectors (u):\n{UsStr}\n\nKnot Vectors (v):\n{VsStr}\n\nControl points:\n{controlPtsStr}')
