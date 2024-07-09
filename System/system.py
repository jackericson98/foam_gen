import time

from numpy import random, pi, cbrt, sqrt, array, linspace, exp, inf
from pandas import DataFrame
from System.calcs import box_search, get_bubbles, calc_dist
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy import stats
import os
from System.output import set_sys_dir, output_all
from Visualize.GUI import settings_gui
from Visualize.mpl_visualize import plot_atoms, plt
import scipy as sp
from System.make_foam import make_foam


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
        self.vpy_dir = os.getcwd()          # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Gui
        self.gui = gui                      # GUI                 :   GUI Vorpy object that can be updated through sys
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Read the arguments of the terminal
        if bubbles is None:
            self.read_argv()
        else:
            self.set_loaded_bubs()
            output_all(self)

    def read_argv(self):
        # Set up the data dictionary
        self.data = {'bubble size': 1, 'bubble sd': 0.1, 'bubble num': 100,
                     'bubble density': 0.25, 'open cell': False, 'distribution': 'lognormal'}

        # Check to see if argv have been made
        if len(self.args) > 1:
            args = self.args[1:]
            for i, data in enumerate(self.data):
                if i >= len(args):
                    break
                self.data[data] = args[i]

        # If we want to prompt the user
        else:
            self.data = settings_gui()
        # Check the open call condition:
        try:
            if self.data['open cell'].lower() in ['true', 't', '1']:
               open_cell = True
            else:
                open_cell = False
        except AttributeError:
            open_cell = False
        self.data = {'bubble size': float(self.data['bubble size']), 'bubble sd': float(self.data['bubble sd']), 'bubble num': int(self.data['bubble num']),
                     'bubble density': float(self.data['bubble density']), 'open cell': open_cell, 'distribution': self.data['distribution']}
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

    def set_loaded_bubs(self):
        # Set the default residue and chain
        residue, chain, my_box = 'BUB', 'A', None
        # Set up the bubbles list
        bubbles = []
        box = [[inf, inf, inf], [-inf, -inf, -inf]]
        # Bubbles loaded need to be (loc, rad)
        for i, bub in enumerate(self.bubbles):
            # Pull the loc and the rad
            my_loc, bub_rad = bub
            # Set the box
            for j in range(3):
                if my_loc[j] - bub_rad < box[0][j]:
                    box[0][j] = my_loc[j] - bub_rad
                if my_loc[j] + bub_rad > box[1][j]:
                    box[1][j] = my_loc[j] + bub_rad
            # Create the bubble
            bubbles.append(
                {'chain': chain, 'loc': my_loc, 'rad': bub_rad, 'num': i, 'name': str(hex(i))[2:], 'asurfs': [],
                 'residue': residue, 'box': my_box})
        # Make the data_frame
        self.bubbles = DataFrame(bubbles)
        self.box = box

    def make_foam(self, print_actions=True):
        make_foam(self, print_actions)
