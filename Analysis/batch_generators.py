import numpy as np


# # Batch of physical foams
# for j in range(20):
#     with open('physical_runs_{}.bat'.format(j % 5), 'w') as write_file:
#         for k in {'gal_or'}:
#             for i in range(1, 21):
#                 write_file.write('py foam_gen.py 1 0 1000 {} False {}\n'.format(round(i*0.025, 5), k))

# def density_update(in_den):
#     return -0.42220023003338425 * in_den ** 2 + 0.9483439624428522 * in_den + 0.0036657838928034754
#
# def density_update1(in_den):
#     return 1.1231 - 1.5390 * np.sqrt(0.5362 - in_den)


# Batch of lognormal foams
for k in range(20):
    # with open('../foam_gen_{}.sh'.format(k), 'w') as peepee_dingus:
    #     peepee_dingus.write('#!/bin/sh\n')
    for density in np.linspace(0.05, 0.5, 10):
        for cv in np.linspace(0.55, 1.0, 10):
            with open('../foam_gen_{}.bat'.format(k % 10), 'a') as poopy_dingus:
                poopy_dingus.write("py foam_gen.py 1 {:.3f} 1000 {:.3f} False gamma\n".format(cv, density))
            # print("py foam_gen.py 1 {:.3f} 1000 {:.3f} True".format((i+4)*0.025, density_update1((j+1)*0.025)))
#
# Batch for Gamma foams

# with open('../foam_gene.bat', 'w') as write_file:
#     for i in range(20):
#         for j in range(17):
#             write_file.write("py foam_gen.py 10 {:.3f} 1000 {:.3f} True gamma\n".format((j+4)*0.5, (i+1)*0.025))

# for i in range(18):
#     for j in range(20):
#         for k in range(10):
# for _ in range(10):
#     for i in range(2):
#         print("python3 foam_gen.py 1 0.5 300 {:.3f} False lognormal".format((i+19)*0.025))