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
ci =[]

for i in range(1, n_ci + 1):
    a = content[i].split()
    ci.append(a[0])
    for j in range(nelec):
        dets[i-1, j] = int(content[i].split()[j + 2])

ci_new = []
for numbers in ci:
    ci_new.append(np.absolute(float(numbers)))

hf_det = dets[0]
occ = []
for i in range(0, n_ci):
    mask = np.in1d(hf_det, dets[i])
    occ.append(nelec - np.sum(mask))

plt.title('ci coefficients vs excitation level')
plt.xlabel('excitation level')
plt.ylabel('ci coefficients')
plt.axis([0, 10, -0.05, 0.7])
plt.plot(occ, ci_new, 'rs')
plt.savefig('occupation_final.png')
plt.show()

print "Total Time:", time.time() - t0
