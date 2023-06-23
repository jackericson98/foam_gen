import os


def write_pdb(atoms, file_name, directory=None):
    """
    Creates a pdb file type in the current working directory
    :param atoms: List of atom type objects for writing
    :param file_name: Name of the output file
    :param sys: System object used for writing the whole pbd file
    :param directory: Output directory for the file
    :return: Writes a pdb file for the set of atoms
    """
    # Catch empty atoms cases
    if atoms is None or len(atoms) == 0:
        return
    # Make note of the starting directory
    start_dir = os.getcwd()
    # Change to the specified directory
    if directory is not None:
        os.chdir(directory)

    # Open the file for writing
    with open(file_name + ".pdb", 'w') as pdb_file:
        # Go through each atom in the system
        for i in atoms:
            a = i
            i = a['num']
            # Get the location string
            loc = ["{:.3f}".format(_) for _ in a['loc']]
            # Get the information from the atom in writable format
            ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
            file_name = a['name'] + " " * (4 - len(a['name']))
            if 'residue' in a:
                res = " " * (3 - len(a['residue'])) + a['residue']
            else:
                res = "   "
            if 'chn' not in a or a['chn'].name.lower() == "zz" or a['chn'].name.lower() == 'mol' or a['chn'].name.lower() == 'sol':
                chain = " "
            else:
                chain = str(a.chn.name)
            if 'res_seq' in a:
                res_seq = " " * (3 - len(str(a['res_seq']))) + str(a['res_seq'])
            else:
                res_seq = "   "
            loc_strs = [" " * (7 - len(_)) + _ for _ in loc]
            occupancy = " " * 5
            t_fact = " " * 5
            if 'seg_id' in a:
                seg_id = a['seg_id'] + " " * (3 - len(a['seg_id']))
            else:
                seg_id = "   "
            if 'element' in a:
                symbol = a['element']
            else:
                symbol = 'h'
            charge = ''
            # Write the atom information
            pdb_file.write("ATOM  " + ser_num + " " + file_name + " " + res + " " + chain + res_seq + "    " +
                           " ".join(loc_strs) + occupancy + t_fact + "      " + seg_id + symbol + charge + "\n")
    # Change back to the starting directory
    os.chdir(start_dir)


def set_pymol_atoms(bubbles):

    """
    Creates a script to set the radii of the spheres in pymol
    :param sys:
    :return:
    """
    special_radii = {}
    # Check to see if the atoms in the system are all accounted for
    for bubble in bubbles.iterrows():
        special_radii[bubble['res_name']] = {bubble['name']: round(bubble['rad'], 2)}
    # Create the file
    with open('set_atoms.pml', 'w') as file:
        # Change the radii for special atoms
        for res in special_radii:
            for atom in special_radii[res]:
                res_str = "residue {} ".format(res) if res != "" else ""
                file.write("alter ({}name {}), vdw={}\n".format(res_str, atom, special_radii[res][atom]))
        # Rebuild the system
        file.write("\nrebuild")
