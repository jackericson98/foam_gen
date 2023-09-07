import sys
import os
from make_foam import make_foam
from cmnd import prompt
from output import set_sys_dir, write_pdb, set_pymol_atoms, write_box, output_all


if __name__ == '__main__':
    if len(sys.argv) == 5:
        bs, bsd, bn, bd = sys.argv[1:]
        bs, bsd, bn, bd = float(bs), float(bsd), int(bn), float(bd)
        bubbles, verts = make_foam(bs, bsd, bn, bd)
        output_all(bubbles, verts)
    elif len(sys.argv) > 1 and sys.argv[1].lower() == 'multi':
        bs, bsd, bn, bd = sys.argv[3:]
        bs, bsd, bn, bd = float(bs), float(bsd), int(bn), float(bd)
        my_dir = set_sys_dir('Data/user_data/multi_foam')
        os.chdir(my_dir)
        for i in range(int(sys.argv[2])):
            bubbles, verts = make_foam(bs, bsd, bn, bd)
            output_all(bubbles, verts, dir='foam')
            os.chdir('..')
    else:
        bs, bsd, bn, bd = prompt()
        bubbles, verts = make_foam(bs, bsd, bn, bd)
        output_all(bubbles, verts)
