from make_foam import make_foam
from output import write_pdb, set_pymol_atoms


def prompt():
    # Get the system information
    bubble_size = input("Enter mean bubble size - ")
    bubble_sd = input("Enter bubble standard deviation - ")
    bubble_num = input("Enter number of bubbles - ")
    bubble_density = input("Enter bubble density - ")
    # Make the foam
    bubbles = make_foam(float(bubble_size), float(bubble_sd), int(bubble_num), float(bubble_density))
    # Write the output files
    file_name = bubble_size + '_' + bubble_sd + '_' + bubble_num + '_' + bubble_density
    write_pdb(bubbles, file_name)
    set_pymol_atoms(bubbles)
