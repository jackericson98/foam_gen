from numpy import random, pi, cbrt, sqrt
from pandas import DataFrame
from System.calcs import box_search, get_bubbles, calc_dist
import os
from System.output import set_sys_dir, output_all
from Visualize.mpl_visualize import plot_atoms, plt
import scipy as sp


class System:
    def __init__(self, args=None, bubbles=None, output_directory=None, gui=None, root_dir=None, print_actions=False):
        """
        Class used to import files of all types and return a System
        :param bubbles: List holding the atom objects
        :param output_directory: Directory for export files to be output to
        :param gui: The GUI object (tkinter) associated with loading the system and loading/creating the network
        """

        # Names
        self.name = None                    # Name                :   Name describing the system

        # Loadable objects
        self.args = args                    # Args                :   Pre loaded arguments to run faster
        self.bubbles = bubbles              # Atoms               :   List holding the atom objects
        self.bubble_matrix = None           # Bubble Matrix       :   3D matrix holding the bubbles for tracking overlap
        self.box = None                     # Bubble box          :   Vertices of the box holding the bubbles

        # Set up the file attributes
        self.data = {'open cell': False}     # Data                :   Additional data provided by the base file
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = root_dir             # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Read the arguments of the terminal
        self.read_argv()

    def read_argv(self):
        # Check if argv
        if len(self.args) == 5:
            bs, bsd, bn, bd = self.args[1:]
            self.data = {'bubble size': float(bs), 'bubble sd': float(bsd), 'bubble num': int(bn), 'bubble density': float(bd), 'open cell': self.data['open cell']}
            self.make_foam()
            output_all(self)
        elif len(self.args) == 6:
            bs, bsd, bn, bd, oc = self.args[1:]
            if oc.lower() in ['n', 'no', 'f', 'false']:
                oc = False
            self.data = {'bubble size': float(bs), 'bubble sd': float(bsd), 'bubble num': int(bn), 'bubble density': float(bd), 'open cell': oc}
            self.make_foam()
            output_all(self)
        elif len(self.args) > 1 and self.args[1].lower() == 'multi':
            bs, bsd, bn, bd = self.args[3:]
            self.data = {'bubble size': float(bs), 'bubble sd': float(bsd), 'bubble num': int(bn), 'bubble density': float(bd)}
            my_dir = set_sys_dir('Data/user_data/multi_foam')
            os.chdir(my_dir)
            for i in range(int(self.args[2])):
                self.make_foam()
                output_all(self)
                os.chdir('..')
        else:
            self.prompt()
            self.make_foam()
            output_all(self)

    def prompt(self, bubble_size=None, bubble_sd=None, bubble_num=None, bubble_density=None, open_cell=None):
        # Get the system information
        if bubble_size is None:
            self.data['bubble size'] = float(input("Enter mean bubble size - "))
        if bubble_sd is None:
            self.data['bubble sd'] = float(input("Enter bubble standard deviation - "))
        if bubble_num is None:
            self.data['bubble num'] = int(input("Enter number of bubbles - "))
        if bubble_density is None:
            self.data['bubble density'] = float(input("Enter bubble density - "))
        if open_cell is None:
            opc = input("Open cell (overlapping)? - ")
            # If user says yes, default is True so no need to catch those cases
            if opc.lower() in ['n', 'no', 'f', 'false']:
                self.data['open cell'] = False

    def make_foam(self, print_actions=True, dist='lognormal'):
        # Get the variables
        mu, sd, n, density, open_cell = self.data['bubble size'], self.data['bubble sd'], self.data['bubble num'], \
            self.data['bubble density'], self.data['open cell']
        # Log normal distribution of radius sizes
        if dist == 'lognormal':
            bubble_radii = []
            while len(bubble_radii) < n:
                rad = random.lognormal(mu, sd, 1)[0] - 1
                if rad > 0:
                    bubble_radii.append(rad)
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
            bubble_radii = [abs(_) for _ in random.geometric(1/mu, n)]
        # Defaults to Normal with a absolute value for less than 0 applicants
        else:
            bubble_radii = [abs(_) for _ in random.normal(mu, sd, n)]
        # By sorting the bubbles they are able to be inserted more quickly into the box
        bubble_radii.sort(reverse=True)
        # Get the maximum bubble radius
        max_bub_radius = max(bubble_radii)
        total_bubble_volume = 0
        # Calculate the total bubble volume
        for bub in bubble_radii:
            # In the open cell case calculate the sum of the actual volumes for the
            if open_cell:
                total_bubble_volume += 4 / 3 * pi * bub ** 3
            # In the closed cell case calculate the density off of the cube surrounding the bubble for extra cushion
            else:
                total_bubble_volume += (2 * bub) ** 3
        # Calculate the retaining cube size
        cube_vol = total_bubble_volume / density
        # Calculate the cube width
        cube_width = cbrt(cube_vol)
        # Create the box
        self.box = [[0, 0, 0], [cube_width, cube_width, cube_width]]
        # Set the number of boxes to roughly 5x the number of atoms must be a cube for the of cells per row/column/aisle
        num_boxes = int(0.5 * sqrt(n)) + 1
        # Instantiate the grid structure of lists is locations representing a grid
        self.bubble_matrix = {(-1, -1, -1): [num_boxes]}
        # Get the cell size
        sub_box_size = [round(cube_width / num_boxes, 3) for i in range(3)]
        # Create the bubble list
        bubbles = []
        # Place the spheres
        for i, bub in enumerate(bubble_radii):
            # Print the loading bar
            if print_actions:
                print("\rCreating bubbles - {:.2f} %".format(100 * (i + 1) / n), end="")
            # Keep trying to place the bubble into a spot where it won't overlap
            while True:
                # Calculate a random bubble location
                my_loc = random.rand(3) * cube_width
                # Find the box that the bubble would belong to
                my_box = [int(my_loc[j] / sub_box_size[j]) for j in range(3)]
                # Find all bubbles within range of the
                bub_ints = get_bubbles(self.bubble_matrix, my_box, sub_box_size, max_bub_radius, bub)
                close_bubs = [bubbles[_] for _ in bub_ints]
                # Create the overlap tracking variable
                overlap = False
                # Loop through the close bubbles
                for bubble in bubbles:
                    # In non-open cell case, check for overlap -> distance less than the sum of radii
                    if not open_cell and calc_dist(my_loc, bubble['loc']) < bub + bubble['rad']:
                        overlap = True
                        break
                    # In open cell case, check for encapsulation -> distance less than the difference of radii
                    elif open_cell and calc_dist(my_loc, bubble['loc']) < abs(bub - bubble['rad']):
                        overlap = True
                        break
                # Skip this location if it overlaps following the overlap criteria
                if overlap:
                    continue
                break
            # Set the default residue and chain
            residue, chain = 'BUB', 'A'
            # Check the location of the bubble
            if any([my_loc[i] < bub or my_loc[i] + bub > cube_width for i in range(3)]):
                residue, chain = 'OUT', ''
            # Create the bubble
            bubbles.append({'chain': chain, 'loc': my_loc, 'rad': bub, 'num': i, 'name': str(hex(i))[2:], 'asurfs': [],
                            'residue': residue, 'box': my_box})
            # Add the atom to the box
            try:
                self.bubble_matrix[my_box[0], my_box[1], my_box[2]].append(i)
            except KeyError:
                self.bubble_matrix[my_box[0], my_box[1], my_box[2]] = [i]
        # Create the dataframe for the bubbles
        self.bubbles = DataFrame(bubbles)
