from make_foam import make_foam
from output import write_pdb, set_pymol_atoms, write_box, set_sys_dir


def prompt(bubble_size=None, bubble_sd=None, bubble_num=None, bubble_density=None):
    # Get the system information
    if bubble_size is None:
        bubble_size = input("Enter mean bubble size - ")
    if bubble_sd is None:
        bubble_sd = input("Enter bubble standard deviation - ")
    if bubble_num is None:
        bubble_num = input("Enter number of bubbles - ")
    if bubble_density is None:
        bubble_density = input("Enter bubble density - ")
    return float(bubble_size), float(bubble_sd), int(bubble_num), float(bubble_density)

