from numpy import random, pi, cbrt, sqrt
from pandas import DataFrame
from System.calcs import box_search, get_bubbles, calc_dist
import os
from System.output import set_sys_dir, output_all


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
        self.box = None
        self.outer_spheres = None           # Outer Spheres       :   Spheres that touch the outside of the box

        # Set up the file attributes
        self.data = {}                      # Data                :   Additional data provided by the base file
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = root_dir             # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not
        
        self.read_argv()

    def read_argv(self):
        # Check if argv
        if len(self.args) == 5:
            bs, bsd, bn, bd = self.args[1:]
            self.data = {'bubble size': float(bs), 'bubble sd': float(bsd), 'bubble num': int(bn), 'bubble density': float(bd)}
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

    def prompt(self, bubble_size=None, bubble_sd=None, bubble_num=None, bubble_density=None):
        # Get the system information
        if bubble_size is None:
            self.data['bubble size'] = float(input("Enter mean bubble size - "))
        if bubble_sd is None:
            self.data['bubble sd'] = float(input("Enter bubble standard deviation - "))
        if bubble_num is None:
            self.data['bubble num'] = int(input("Enter number of bubbles - "))
        if bubble_density is None:
            self.data['bubble density'] = float(input("Enter bubble density - "))

    def make_foam(self, num_boxes=None, open_cell=True):
        # Get the variables
        mu, sd, n, density = self.data['bubble size'], self.data['bubble sd'], self.data['bubble num'], self.data['bubble density']
        # Get the radii for the bubbles in the foam
        bubble_radii = [_ - 1 for _ in random.lognormal(mu, sd, n)]
        max_bub_radius = max(bubble_radii)
        total_bubble_volume = 0
        # Calculate the total bubble volume
        for bub in bubble_radii:
            total_bubble_volume += 4 / 3 * pi * bub ** 3
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
                bubbles.append({'chain': 'A', 'loc': my_loc, 'rad': bub, 'num': i, 'name': 'b' + str(i), 'asurfs': [],
                                'residue': str(i)})
            else:
                while True:
                    my_loc = random.rand(3) * cube_width
                    my_box = box_search(my_loc, num_splits, box_verts=[[0, 0, 0], [cube_width, cube_width, cube_width]])
                    my_close_bubs = get_bubbles(bubble_matrix, my_box, sub_box_size, max_bub_radius)
                    for bubble in my_close_bubs:
                        if calc_dist(my_loc, bubble['loc']) > bub + bubble['rad']:
                            bubbles.append(
                                {'chain': 'A', 'loc': my_loc, 'rad': bub, 'num': i, 'name': 'b' + str(i), 'asurfs': [],
                                 'residue': str(i)})
                            break
        verts = [[0, 0, 0], [cube_width, cube_width, cube_width]]
        self.bubbles, self.box = DataFrame(bubbles), verts