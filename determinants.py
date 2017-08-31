import matplotlib.pyplot as plt
import numpy as np

with open('determinantData.txt') as f:
    data = f.readlines()

n_ci = len(content) - 1

ci = []
Ed = []
for line in data:
    a = line.split()
    ci.append(a[0])
    Ed.append(a[1])

ci.pop(0)
Ed.pop(0)

ci_new = []
for numbers in ci:
    ci_new.append(np.absolute(float(numbers)))

Ed_new = []
for items in Ed:
    Ed_new.append(float(items))

plt.title('ci Coefficients vs Determinant Energy')
plt.xlabel('Determinant Energy')
plt.ylabel('ci Coefficients')
plt.plot(Ed_new, ci_new, 'bo')
plt.savefig('energy_final')
plt.show()
