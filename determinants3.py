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
ci = []

for i in range(1, n_ci + 1):
    a = content[i].split()
    ci.append(a[0])
    for j in range(nelec):
        dets[i-1, j] = int(content[i].split()[j + 2])

ci_new = []
for numbers in ci:
    ci_new.append(np.absolute(float(numbers)))

pair = [6]*n_ci
for i in range(0, n_ci):
    for j in range(nelec - 1):
        for k in range(0 , 31, 2):
            if dets[i, j] == k and dets[i, j + 1] == k + 1:
                pair[i] = pair[i] - 1
    pair[i] = pair[i]*2

plt.title('ci coefficients vs unpaired electrons')
plt.xlabel('unpaired electrons')
plt.ylabel('ci coefficients')
plt.plot(pair, ci_new, 'gD')
plt.savefig('pair_final.png')
plt.show()

print "Total Time:", time.time() - t0
