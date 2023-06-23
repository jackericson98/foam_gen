from calcs import calc_dist
import numpy as np
from pandas import DataFrame


def make_foam(mu, sd, n, density):
    # Get the radii for the bubbles in the foam
    bubble_radii = [_ - 1 for _ in np.random.lognormal(mu, sd, n)]

    total_bubble_volume = 0
    # Calculate the total bubble volume
    for bub in bubble_radii:
        total_bubble_volume += 4/3 * np.pi * bub ** 3
    # Calculate the retaining cube size
    cube_vol = total_bubble_volume / density
    # Calculate the cube width
    cube_width = np.cbrt(cube_vol)
    bubbles = []
    # Place the spheres
    for i, bub in enumerate(bubble_radii):
        while True:
            my_loc = np.random.rand(3) * cube_width
            overlap = False
            for bubble in bubbles:
                if calc_dist(my_loc, bubble['loc']) < bub + bubble['rad']:
                    overlap = True
            if not overlap:
                bubbles.append({'loc': my_loc, 'rad': bub, 'num': i, 'name': 'b' + str(i), 'asurfs': [], 'residue': str(i)})
                break
    return DataFrame(bubbles)
