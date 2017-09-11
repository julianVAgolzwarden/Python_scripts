import matplotlib.pyplot as plt
import numpy as np
import time

t0 = time.time()

with open('determinantData.txt') as f:
    content = f.readlines()

n_ci = len(content) - 1
norbs = 16
nspinorbs = 2*norbs
nelec = 12

dets = np.zeros((n_ci, nelec))

for i in range(1, n_ci + 1):
    for j in range(nelec):
        dets[i-1, j] = int(content[i].split()[j + 2])

np.savetxt('hamilton_matrix.txt', dets, fmt='%1i')
print "Total Time:", time.time() - t0
