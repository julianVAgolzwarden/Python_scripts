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

Ed = np.zeros(n_ci,)

for i in range(1, n_ci + 1):
    Ed[i-1] = float(content[i].split()[1])

np.savetxt('Determinant_Energies.txt', Ed)
print "Total Time:", time.time() - t0
exit(0)


hf_det = dets[0]
occ = np.zeros(n_ci, dtype = np.int64)
for i in range(0, n_ci):
    mask = np.in1d(hf_det, dets[i])
    occ[i] = nelec - np.sum(mask)

np.savetxt('Excitation_levels.txt', occ, fmt='%1i')
print "Total Time:", time.time() - t0
