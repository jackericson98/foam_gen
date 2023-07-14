from calcs import calc_dist, box_search, get_bubbles
from numpy import sqrt, random, pi, cbrt
from pandas import DataFrame


def make_foam(mu, sd, n, density, num_boxes=None, open_cell=True):
    # Get the radii for the bubbles in the foam
    bubble_radii = [_ - 1 for _ in random.lognormal(mu, sd, n)]
    max_bub_radius = max(bubble_radii)
    total_bubble_volume = 0
    # Calculate the total bubble volume
    for bub in bubble_radii:
        total_bubble_volume += 4/3 * pi * bub ** 3
    # Calculate the retaining cube size
    cube_vol = total_bubble_volume / density
    # Calculate the cube width
    cube_width = cbrt(cube_vol)

    # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
    if num_boxes is None:
        num_boxes = int(0.5 * sqrt(n)) + 1
    else:
        num_boxes = int(cbrt(num_boxes)) + 1
    num_splits = num_boxes
    # Instantiate the grid structure of lists is locations representing a grid
    bubble_matrix = {(-1, -1, -1): [n]}
    # Get the cell size
    sub_box_size = [round(cube_width / n, 3) for i in range(3)]

    bubbles = []
    # Place the spheres
    for i, bub in enumerate(bubble_radii):
        if open_cell:
            my_loc = random.rand(3) * cube_width
            bubbles.append({'loc': my_loc, 'rad': bub, 'num': i, 'name': 'b' + str(i), 'asurfs': [], 'residue': str(i)})
        else:
            while True:
                my_loc = random.rand(3) * cube_width
                my_box = box_search(my_loc, num_splits, box_verts=[[0,0,0], [cube_width, cube_width, cube_width]])
                my_close_bubs = get_bubbles(bubble_matrix, my_box, sub_box_size, max_bub_radius)
                for bubble in my_close_bubs:
                    if calc_dist(my_loc, bubble['loc']) > bub + bubble['rad']:
                        bubbles.append({'loc': my_loc, 'rad': bub, 'num': i, 'name': 'b' + str(i), 'asurfs': [], 'residue': str(i)})
                        break
    verts = [[0, 0, 0], [cube_width, cube_width, cube_width]]
    return DataFrame(bubbles), verts
