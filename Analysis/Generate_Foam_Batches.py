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
    count_dict = {}
    for rooot, folders, files in os.walk(folder):
        for folder in folders:
            folder_vals = folder.split('_')
            if len(folder_vals) < 7:
                continue
            try:
                mean, cv, number, density, olap, dist, pbc = float(folder_vals[0]), float(folder_vals[1]), int(folder_vals[2]), float(folder_vals[3]), float(folder_vals[4]), folder_vals[5], folder_vals[6]
                if (mean, cv, number, density, olap, dist, pbc) in count_dict:
                    count_dict[(mean, cv, number, density, olap, dist, pbc)] += 1
                else:
                    count_dict[(mean, cv, number, density, olap, dist, pbc)] = 1
            except ValueError:
                continue

    # Set the settings
    means = [1.0]
    cv_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    density_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    olap_vals = [0.05]
    num_vals = [1000]
    pbc_vals = ['False']
    dist_vals = ['gamma']
    my_string = "python foam_gen.py mean {:.3f} cv {:.3f} num {} den {:.3f} olp {} dist {} pbc {} \n"
    number_of_foam_files = 20
    number_of_files = 4

    # Create the new strings list
    new_strings = []
    # Now loop through the cv_vals and density vals
    for mean in means:
        for dist in dist_vals:
            for num in num_vals:
                for olap in olap_vals:
                    for cv in cv_vals:
                        for pbc in pbc_vals:
                            for density in density_vals:
                                if (mean, cv, num, density, olap, dist, pbc) in count_dict:
                                    num_files_to_make = number_of_foam_files - count_dict[(mean, cv, num, density, olap, dist, pbc)]
                                else:
                                    num_files_to_make = number_of_foam_files
                                for _ in range(num_files_to_make):

                                    new_strings.append(my_string.format(mean, cv, num, density, olap, dist, pbc))

    for i, stank in enumerate(new_strings):
        file_number = int(i // (len(new_strings) / number_of_files))
        j = file_number + 1
        with open('../foam_gen_{}.sh'.format(j), 'a') as foam_writer:
            foam_writer.write(stank)
