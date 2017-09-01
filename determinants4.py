import matplotlib.pyplot as plt
import numpy as np
import time

t0 = time.time()

with open('ci_coefficients.txt') as c:
    content = c.readlines()

with open('Determinant_Energies.txt') as d:
    items = d.readlines()

with open('Excitation_levels.txt') as e:
    objects = e.readlines()

n_ci = len(content)

ci = np.zeros(n_ci,)
Ed = np.zeros(n_ci,)
occ = np.zeros(n_ci,)

for i in range(n_ci):
    ci[i] = content[i]
    Ed[i] = items[i]
    occ[i] = objects[i]

PT = np.zeros(n_ci,)
for j in range(0, n_ci):
    PT[j] = occ[j]/(Ed[0] - Ed[j])

ci_new = ci.tolist()
PT_new = PT.tolist()
ci_new.pop(0)
PT_new.pop(0)

plt.title('PT approximation')
plt.xlabel(r'excitation level/deltaEd')
plt.ylabel('ci coefficients')
plt.plot(PT_new, ci_new, 'm^')
plt.savefig('PT_approximation.png')
plt.show()

print "Total Time:", time.time() - t0
