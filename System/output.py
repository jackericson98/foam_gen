from System.draw import draw_line
import numpy as np
import os
from System.calcs import pdb_line


def output_all(sys, dir=None):
    bubbles, verts = sys.bubbles, sys.box
    if dir is None:
        # Write the output files
        file_name = 'foam'
        my_dir = set_sys_dir('Data/user_data/' + file_name)
    else:
        file_name = 'foam'
        my_dir = set_sys_dir('foam')

    write_pdb(bubbles, file_name, directory=my_dir, box=verts)
    set_pymol_atoms(bubbles, directory=my_dir)
    write_box(verts, file_name='retaining_box', directory=my_dir)


def write_pdb(bubbles, file_name, directory=None, box=None):
    """
    Creates a pdb file type in the current working directory
    :param bubbles: List of atom type objects for writing
    :param file_name: Name of the output file
    :param sys: System object used for writing the whole pbd file
    :param directory: Output directory for the file
    :return: Writes a pdb file for the set of atoms
    """
    # Catch empty atoms cases
    if bubbles is None or len(bubbles) == 0:
        return
    # Make note of the starting directory
    start_dir = os.getcwd()
    # Change to the specified directory
    if directory is not None:
        os.chdir(directory)
    if box is None:
        min_vals = [np.inf, np.inf, np.inf]
        max_vals = [-np.inf, -np.inf, -np.inf]
        for j, bubble in bubbles.iterrows():
            for i in range(3):
                if bubble['loc'][i] < min_vals[i]:
                    min_vals[i] = bubble['loc'][i]
                if bubble['loc'][i] > max_vals[i]:
                    max_vals[i] = bubble['loc'][i]
        box = [min_vals, max_vals]
    # Open the file for writing
    with open(file_name + ".pdb", 'w') as pdb_file:
        # Write the header that lets vorpy know it is a foam pdb
        pdb_file.write('REMARK foam_gen {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(*box[0], *box[1]))
        # Go through each atom in the system
        for i, a in bubbles.iterrows():
            # Get the location string
            x, y, z = a['loc']
            occ = 1

            # Write the atom information
            pdb_file.write(pdb_line(ser_num=i, name=a['name'], res_name=a['residue'], chain=a['chain'],
                                    x=x, y=y, z=z, occ=occ, tfact=a['rad']))
    # Change back to the starting directory
    os.chdir(start_dir)


def set_pymol_atoms(bubbles, directory=None):

    """
    Creates a script to set the radii of the spheres in pymol
    :param sys:
    :return:
    """
    # Make note of the starting directory
    start_dir = os.getcwd()
    if directory is not None:
        os.chdir(directory)
    special_radii = {}
    # Check to see if the atoms in the system are all accounted for
    for i, bubble in bubbles.iterrows():
        special_radii[bubble['name']] = {bubble['name']: round(bubble['rad'], 2)}
    # Create the file
    with open('set_atoms.pml', 'w') as file:
        # Change the radii for special atoms
        for res in special_radii:
            for atom in special_radii[res]:
                res_str = "residue {} ".format(res) if res != "" else ""
                file.write("alter ({}name {}), vdw={}\n".format(res_str, atom, special_radii[res][atom]))
        # Rebuild the system
        file.write("\nrebuild")
    os.chdir(start_dir)


def write_box(verts, file_name, color=None, directory=None):
    """
    Writes an off file for the edges specified
    :param edges: Edges to be output
    :param file_name: Name for the output file
    :param color: Color for the edges
    :param directory: Output directory
    :return: None
    """
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no color is given, make the color random
    if color is None:
        color = [0.5, 0.5, 0.5]
    # Check that the edge has been drawn
    edges_draw_points, edges_draw_tris = [], []
    lines = [[[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 1]], [[1, 0, 0], [1, 1, 0]],
             [[1, 0, 0], [1, 0, 1]], [[0, 1, 0], [1, 1, 0]], [[0, 1, 0], [0, 1, 1]], [[0, 0, 1], [1, 0, 1]],
             [[0, 0, 1], [0, 1, 1]], [[1, 1, 0], [1, 1, 1]], [[1, 0, 1], [1, 1, 1]], [[0, 1, 1], [1, 1, 1]]]
    points = []
    for line in lines:
        p0, p1 = [verts[line[0][i]][i] for i in range(3)], [verts[line[1][i]][i] for i in range(3)]
        points.append([p0, p1])
        draw_points, draw_tris = draw_line([p0, p1])
        edges_draw_points.append(draw_points)
        edges_draw_tris.append(draw_tris)
    num_verts, num_tris = 72, 72
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for line in edges_draw_points:
            # Go through the points on the surface
            for point in line:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for line in edges_draw_tris:
            # Go through the triangles in the surface
            for tri in line:
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += 6


def set_sys_dir(dir_name=None):
    """
    Sets the directory for the output data. If the directory exists add 1 to the end number
    :param sys: System to assign the output directory to
    :param dir_name: Name for the directory
    :return:
    """
    if dir_name is None:
        # If no outer directory was specified use the directory outside the current one
        dir_name = os.getcwd() + 'foam'

    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    return dir_name + i_str
