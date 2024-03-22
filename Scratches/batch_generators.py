
# # Batch of physical foams
# for i in range(1, 21):
#     print("python3 foam_gen.py 1 0.01 300 {} False physical1".format(round(i*0.025, 5)))

# Batch of lognormal foams
for i in range(16):
    for j in range(20):
        print("python3 foam_gen.py 1 {:.3f} 1000 {:.3f} True".format((i+4)*0.025, (j+1)*0.025))
#
# # Batch for Gamma foams
# for i in range(18):
#     for j in range(21):
#         print("python3 foam_gen.py 1 {:.3f} 1000 {:.3f} False gamma".format((j+6)*0.25, (i+1)*0.025))

# for i in range(18):
#     for j in range(20):
#         for k in range(10):
# for _ in range(10):
#     for i in range(2):
#         print("python3 foam_gen.py 1 0.5 300 {:.3f} False lognormal".format((i+19)*0.025))