import time
from numpy import sqrt, array, random, linspace, pi, exp, cbrt
from scipy import stats
from pandas import DataFrame
from System.calcs import get_bubbles, calc_dist, calc_dist_numba, calc_tot_vol
from scipy.integrate import quad
from scipy.interpolate import interp1d
from numba import jit


def get_bubble_raddi(dist, mu, sd, n):

    # Calculate the cumulative distribution function (CDF)
    def calculate_cdf(pdf, x_values):
        cdf_values = array([quad(pdf, 0, x)[0] for x in x_values])
        cdf_values /= cdf_values[-1]  # Normalize to [0, 1]
        return cdf_values

    # Generate random samples using inverse transform sampling
    def inverse_transform_sampling(pdf, x_values, n_samples):
        cdf_values = calculate_cdf(pdf, x_values)
        inverse_cdf = interp1d(cdf_values, x_values, kind='linear', fill_value='extrapolate')
        u = random.rand(n_samples)
        return inverse_cdf(u)

    # Log normal distribution of radius sizes
    if dist == 'lognormal':
        bubble_radii = []
        while len(bubble_radii) < n:
            rad = random.lognormal(mu, sd, 1)[0] - 1
            if rad > 0:
                bubble_radii.append(rad)
    # Gamma distribution of radius sizes
    if dist == 'gamma':
        def gamma(r):
            return 0.5 * stats.gamma.pdf(r, 4 * mu, scale=1 / sd)

        x_values = linspace(0, 5, n)
        bubble_radii = inverse_transform_sampling(gamma, x_values, n)
    # Half normal distribution of radius sizes
    elif dist == 'halfnormal':
        bubble_radii = [abs(_) for _ in random.normal(0, sd, n)]
    # Normal distribution of radius sizes cut off at 0
    elif dist == 'normal':
        bubble_radii = []
        # Keep creating radii randomly
        while len(bubble_radii) < n:
            rad = random.normal(mu, sd, 1)[0]
            # Only accept radii over 0
            if rad > 0:
                bubble_radii.append(rad)
    # Geometric distribution
    elif dist == 'geometric':
        bubble_radii = [abs(_) for _ in random.geometric(1 / mu, n)]
    # devries
    elif dist == 'physical1':
        # Define the pdf for Devries
        def pdf(r):
            return 2.082 * r / (1 + 0.387 * r ** 2) ** 4

        x_values = linspace(0, 5, n)
        bubble_radii = inverse_transform_sampling(pdf, x_values, n)
    # Ranadive & Lemlich
    elif dist == 'physical2':
        # Define the pdf for Devries
        def pdf(r):
            return (32 / pi ** 2) * r ** 2 * exp(-(4 / pi) * r ** 2)

        x_values = linspace(0, 5, n)
        bubble_radii = inverse_transform_sampling(pdf, x_values, n)
        # Gal-Or & Hoelscher
    elif dist == 'physical3':
        # Define the pdf for Devries
        def pdf(r):
            return (16 / pi) * r ** 2 * exp(-sqrt(16 / pi) * r ** 2)

        x_values = linspace(0, 5, n)
        bubble_radii = inverse_transform_sampling(pdf, x_values, n)
    # Defaults to Normal with a absolute value for less than 0 applicants
    else:
        bubble_radii = [abs(_) for _ in random.normal(mu, sd, n)]
    # By sorting the bubbles they are able to be inserted more quickly into the box
    bubble_radii = sorted(bubble_radii, reverse=True)
    # Return the bubble radii
    return bubble_radii


def wall_overlap(my_loc, bub, cube_width):
    # Create the overlap tracking variable
    # Check to see if the bubble overlaps with the wall
    for k in range(3):
        if my_loc[k] - bub < 0 or my_loc[k] + bub > cube_width:
            return True


@jit
def overlap(my_loc, bub, close_locs, close_rads, open_cell):
    # Loop through the close bubbles
    for i in range(len(close_locs)):
        # In non-open cell case, check for overlap -> distance less than the sum of radii
        if not open_cell and calc_dist_numba(my_loc, close_locs[i]) < bub + close_rads[i]:
            return True
        # In open cell case, check for encapsulation -> distance less than the difference of radii
        elif open_cell and calc_dist_numba(my_loc, close_locs[i]) < 0.5 * abs(bub - close_rads[i]):
            return True
    return False


def find_bubs(bubble_radii, num_boxes, cube_width, sub_box_size, max_bub_radius, num_tries, open_cell, n, print_actions):
    # Create the bubble list
    bubbles = []
    # Instantiate the grid structure of lists is locations representing a grid
    bubble_matrix = {(-1, -1, -1): [num_boxes]}
    time_start = time.perf_counter()
    break_all = False
    # Place the spheres
    for i, bub_rad in enumerate(bubble_radii):
        # Print the loading bar
        if print_actions:
            my_time = round(time.perf_counter() - time_start)
            print("\rCreating bubbles {} - {:.2f}% - {} s".format(num_tries, 100 * (i + 1) / n, my_time), end="")
        # Keep trying to place the bubble into a spot where it won't overlap
        while True:
            # Calculate a random bubble location
            my_loc = random.rand(3) * cube_width
            # Find the box that the bubble would belong to
            my_box = [int(my_loc[j] / sub_box_size[j]) for j in range(3)]
            # Check to see if the bubble overlaps with the wall
            if wall_overlap(my_loc, bub_rad, cube_width):
                continue
            # # Find all bubbles within range of the
            bub_ints = get_bubbles(bubble_matrix, my_box, sub_box_size, max_bub_radius, bub_rad)
            close_bubs = [bubbles[_] for _ in bub_ints]
            # Sort the bubbles by size
            close_bubs.sort(key=lambda x: x['rad'], reverse=True)
            # Check for overlap with close bubs
            close_locs, close_rads = array([_['loc'] for _ in close_bubs]), array([_['rad'] for _ in close_bubs])
            if not overlap(my_loc, bub_rad, close_locs, close_rads, open_cell):
                break
        if break_all:
            break
        # Set the default residue and chain
        residue, chain = 'BUB', 'A'
        # Check the location of the bubble
        if any([my_loc[i] < bub_rad or my_loc[i] + bub_rad > cube_width for i in range(3)]):
            residue, chain = 'OUT', ''
        # Create the bubble
        bubbles.append({'chain': chain, 'loc': my_loc, 'rad': bub_rad, 'num': i, 'name': str(hex(i))[2:], 'asurfs': [],
                        'residue': residue, 'box': my_box})
        # Add the atom to the box
        try:
            bubble_matrix[my_box[0], my_box[1], my_box[2]].append(i)
        except KeyError:
            bubble_matrix[my_box[0], my_box[1], my_box[2]] = [i]

        if break_all:
            break
    return bubbles, bubble_matrix


def make_foam(sys, print_actions):
    # Get the variables
    mu, sd, n, density, open_cell, dist = sys.data['bubble size'], sys.data['bubble sd'], sys.data['bubble num'], \
        sys.data['bubble density'], sys.data['open cell'], sys.data['distribution']
    # Get the bubble radii
    bubble_radii = get_bubble_raddi(dist, mu, sd, n)
    # Get the maximum bubble radius
    max_bub_radius = max(bubble_radii)
    # Calculate the cube width by getting the cube root of total volume of the bubble radii over the density
    cube_width = cbrt(calc_tot_vol(bubble_radii) / density)
    # Create the box
    sys.box = [[0, 0, 0], [cube_width, cube_width, cube_width]]
    # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
    num_boxes = int(0.5 * sqrt(n)) + 1
    # Get the cell size
    sub_box_size = [round(cube_width / num_boxes, 3) for i in range(3)]
    num_tries = 0
    while True:
        bubbles, sys.bubble_matrix = find_bubs(bubble_radii, num_boxes, cube_width, sub_box_size, max_bub_radius, num_tries, open_cell, n, print_actions)
        num_tries += 1
        if num_tries > 3:
            return
        if len(bubble_radii) == len(bubbles):
            break

    # Create the dataframe for the bubbles
    sys.bubbles = DataFrame(bubbles)
