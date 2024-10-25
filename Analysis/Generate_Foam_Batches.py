import tkinter as tk
from tkinter import filedialog
import numpy as np
import os



if __name__ == '__main__':
    # Open the folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()

    # Set the settings
    cv_vals = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
    density_vals = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5}
    my_string = "python foam_gen.py cv {:.3f} den {:.3f} num 1000 olp 0.0 dist gamma pbc t\n"
    number_of_foam_files = 20
    number_of_files = 5

    # Create the dictionary
    dir_dict = {}
    # Get the quantity of each of the required generations by hand
    for rroot, directories, files in os.walk(folder):
        for directory in directories:
            print(directory)
            if '.csv' in directory:
                continue
            dir_comps = directory.split('_')
            cv, den = round(float(dir_comps[1]), 3), round(float(dir_comps[3]), 3)
            if (cv, den) in dir_dict:
                dir_dict[(cv, den)].append(directory)
            else:
                dir_dict[(cv, den)] = [directory]
    # Create the new strings list
    new_strings = []
    # Now loop through the cv_vals and density vals
    for cv in cv_vals:
        for density in density_vals:
            if (cv, density) in dir_dict:
                for _ in range(number_of_foam_files - len(dir_dict[(cv, density)])):
                    new_strings.append(my_string.format(cv, density))
            else:
                for _ in range(20):
                    new_strings.append(my_string.format(cv, density))

    for i, stank in enumerate(new_strings):
        file_number = i // (len(new_strings) / number_of_files)
        with open('../foam_gen_{}.sh'.format(file_number), 'a') as foam_writer:
            foam_writer.write(stank)