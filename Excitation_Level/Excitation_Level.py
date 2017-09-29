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
ci = np.zeros(n_ci,)

for i in range(1, n_ci + 1):
    ci[i-1] = np.fabs(float(content[i].split()[0]))
    for j in range(nelec):
        dets[i-1, j] = int(content[i].split()[j + 2])

ci_new = ci.tolist()

hf_det = dets[0]
occ = np.zeros(n_ci,)
for i in range(0, n_ci):
    mask = np.in1d(hf_det, dets[i])
    occ[i] = int(nelec - np.sum(mask))

print mask
exit(0)
occ_new = occ.tolist()

plt.title('ci coefficients vs excitation level')
plt.xlabel('excitation level')
plt.ylabel('ci coefficients')
plt.axis([0, 10, -0.05, 0.7])
plt.plot(occ_new, ci_new, 'rs')
plt.savefig('occupation_final.png')
plt.show()

print "Total Time:", time.time() - t0
