import sys
from make_foam import make_foam
from cmnd import prompt
from output import set_sys_dir, write_pdb, set_pymol_atoms, write_box


if __name__ == '__main__':
    if len(sys.argv) == 5:
        bs, bsd, bn, bd = sys.argv[1:]
        bs, bsd, bn, bd = float(bs), float(bsd), int(bn), float(bd)
    else:
        bs, bsd, bn, bd = prompt()
    bubbles, verts = make_foam(bs, bsd, bn, bd)
    print(bubbles.head())
    # Write the output files
    file_name = 'foam'
    my_dir = set_sys_dir('Data/user_data/' + file_name)
    write_pdb(bubbles, file_name, directory=my_dir)
    set_pymol_atoms(bubbles, directory=my_dir)
    write_box(verts, file_name='retaining_box', directory=my_dir)
