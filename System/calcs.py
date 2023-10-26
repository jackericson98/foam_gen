from numpy import sqrt, square, inf, array


def calc_dist(l0, l1):
    """
    Calculate distance function used to simplify code
    :param l0: Point 0 list, array, n-dimensional must match point 1
    :param l1: Point 1 list, array, n-dimensional must match point 0
    :return: float distance between the two points
    """
    # Pythagorean theorem
    return sqrt(sum(square(l0 - l1)))


def box_search(loc, num_splits, box_verts):
    # Calculate the size of the sub boxes
    sub_box_size = [round((box_verts[1][i] - box_verts[0][i]) / num_splits, 3) for i in range(3)]
    # Find the sub box for the atom
    box_ndxs = [int((loc[j] - box_verts[0][j]) / sub_box_size[j]) for j in range(3)]
    if box_ndxs[0] >= num_splits or box_ndxs[1] >= num_splits or box_ndxs[2] >= num_splits:
        return
    # Return the box indices
    return box_ndxs


def get_bubbles(bubble_matrix, cells, sub_box_size, max_atom_rad, dist=0):
    """
    Takes in the cells and the number of additional cells to search and returns an atom list
    :param cells: The initial boxes in the network to stem from
    :param dist: The number of cells out from the initial set of cells to search
    """
    # Calculate the number of boxes out needed to find any potential overlapping sphere
    reach = int((dist + max_atom_rad) / min(sub_box_size)) + 4
    n = bubble_matrix[-1, -1, -1][0]
    # If a single cell is entered
    if type(cells[0]) is int:
        cells = [cells]
    # Get the min and max of the cells
    ndx_min = [inf, inf, inf]
    ndx_max = [-inf, -inf, -inf]
    # Go through the cells and set the minimum and maximum indexes for xyz for a rectangle containing the atoms
    for cell in cells:
        # Check each xyz index to see if they are larger or smaller than the max or min
        for i in range(3):
            if cell[i] < ndx_min[i]:
                ndx_min[i] = cell[i]
            if cell[i] > ndx_max[i]:
                ndx_max[i] = cell[i]
    xs = [x for x in range(max(0, -reach + ndx_min[0] + 1), reach + ndx_max[0])]
    ys = [y for y in range(max(0, -reach + ndx_min[1] + 1), reach + ndx_max[1])]
    zs = [z for z in range(max(0, -reach + ndx_min[2] + 1), reach + ndx_max[2])]
    atoms = []
    # Get atoms
    for i in xs:
        if 0 <= i < n:
            for j in ys:
                if 0 <= j < n:
                    for k in zs:
                        if 0 <= k < n:
                            try:
                                atoms += bubble_matrix[i, j, k]
                            except KeyError:
                                pass
    return atoms


def calc_box(self, locs, rads):
    """
    Determines the dimensions of a box x times the size of the atoms
    :return: Sets the box attribute with the correct values as well as atoms_box
    """
    # Set up the minimum and maximum x, y, z coordinates
    min_vert = array([inf, inf, inf])
    max_vert = array([-inf, -inf, -inf])
    # Loop through each atom in the network
    for loc in locs:
        # Loop through x, y, z
        for i in range(3):
            # If x, y, z values are less replace the value in the mins list
            if loc[i] < min_vert[i]:
                min_vert[i] = loc[i]
            # If x, y, z values are greater replace the value in the maxes list
            if loc[i] > max_vert[i]:
                max_vert[i] = loc[i]
    # Get the vector between the minimum and maximum vertices for the defining box
    r_box = max_vert - min_vert
    # If the atoms are in the same plane adjust the atoms
    for i in range(3):
        if r_box[i] == 0 or abs(r_box[i]) == inf:
            r_box[i], min_vert[i], max_vert[i] = 4 * rads[0], locs[0][i], locs[0][i]
    # Set the atoms box value
    self.atoms_box = [min_vert.tolist(), max_vert.tolist()]
    # Set the new vertices to the x factor times the vector between them added to their complimentary vertices
    min_vert, max_vert = max_vert - r_box * self.box_size, min_vert + r_box * self.box_size
    # Return the list of array turned list vertices
    self.box = [[round(_, 3) for _ in min_vert], [round(_, 3) for _ in max_vert]]


def pdb_line(atom="ATOM", ser_num=0, name="", alt_loc=" ", res_name="", chain="A", res_seq=0, cfir="", x=0, y=0, z=0,
             occ=1, tfact=0, seg_id="", elem="h", charge=""):
    if chain == 'Z':
        chain = ' '
    return "{:<6}{:>5} {:<4}{:1}{:>3} {:^1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{}\n"\
        .format(atom, ser_num, name, alt_loc, res_name, chain, res_seq, cfir, x, y, z, occ, tfact, seg_id, elem, charge)


def calc_stats()
    import numpy as np

    import matplotlib.pyplot as plt

    from scipy import stats
    # x=np.linspace(0.1,0.1,200)
    # bubs = []
    # while len(bubs) < 1000:
    #     data = stats.lognorm(1, loc=0, scale=1).rvs(size=1)[0]
    #     if data > 0:
    #         bubs.append(data)
    #
    # plt.hist(bubs, bins=100)
    # plt.show()
    sds = [0.01, 0.05, 0.09]
    numbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900,
               1000]
    mu = 0.1
    num_its = 100
    for sd in sds:
        sd1s, sd2s, sd3s = [], [], []
        sd1gsl, sd2gsl, sd3gsl = [], [], []
        for i in numbers:
            sd1, sd2, sd3 = [], [], []
            sd1gs, sd2gs, sd3gs = 0, 0, 0
            for j in range(num_its):
                bubs = []
                while len(bubs) < i:
                    data = stats.lognorm(sd, loc=mu).rvs(size=1)[0] - 1
                    if data > 0:
                        bubs.append(data)
                sd1_count, sd2_count, sd3_count = 0, 0, 0
                for bub in bubs:
                    if abs(bub - mu) < sd:
                        sd1_count += 1
                    if abs(bub - mu) < 2 * sd:
                        sd2_count += 1
                    if abs(bub - mu) < 3 * sd:
                        sd3_count += 1
                sd1.append(sd1_count)
                sd2.append(sd2_count)
                sd3.append(sd3_count)
                sd1g, sd2g, sd3g = 0, 0, 0
                if sd1_count / i > 0.68:
                    sd1g = 1
                if sd2_count / i > 0.95:
                    sd2g = 1
                if sd3_count / i > 0.99:
                    sd3g = 1
                sd1gs += sd1g
                sd2gs += sd2g
                sd3gs += sd3g
            # print("{} Done".format(i))
            print(i, (sum(sd1) / num_its) / i, (sum(sd2) / num_its) / i, (sum(sd3) / num_its) / i, sd1gs, sd2gs, sd3gs)
import numpy as np

import matplotlib.pyplot as plt

from scipy import stats
# x=np.linspace(0.1,0.1,200)
# bubs = []
# while len(bubs) < 1000:
#     data = stats.lognorm(1, loc=0, scale=1).rvs(size=1)[0]
#     if data > 0:
#         bubs.append(data)
#
# plt.hist(bubs, bins=100)
# plt.show()
sds = [0.01, 0.05, 0.09]
numbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
mu = 0.1
num_its = 100
for sd in sds:
    sd1s, sd2s, sd3s = [], [], []
    sd1gsl, sd2gsl, sd3gsl = [], [], []
    for i in numbers:
        sd1, sd2, sd3 = [], [], []
        sd1gs, sd2gs, sd3gs = 0, 0, 0
        for j in range(num_its):
            bubs = []
            while len(bubs) < i:
                data = stats.lognorm(sd, loc=mu).rvs(size=1)[0] - 1
                if data > 0:
                    bubs.append(data)
            sd1_count, sd2_count, sd3_count = 0, 0, 0
            for bub in bubs:
                if abs(bub - mu) < sd:
                    sd1_count += 1
                if abs(bub - mu) < 2 * sd:
                    sd2_count += 1
                if abs(bub - mu) < 3 * sd:
                    sd3_count += 1
            sd1.append(sd1_count)
            sd2.append(sd2_count)
            sd3.append(sd3_count)
            sd1g, sd2g, sd3g = 0, 0, 0
            if sd1_count / i > 0.68:
                sd1g = 1
            if sd2_count / i > 0.95:
                sd2g = 1
            if sd3_count / i > 0.99:
                sd3g = 1
            sd1gs += sd1g
            sd2gs += sd2g
            sd3gs += sd3g
        # print("{} Done".format(i))
        print(i, (sum(sd1)/num_its) / i, (sum(sd2)/num_its) / i, (sum(sd3)/num_its) / i, sd1gs, sd2gs, sd3gs)
