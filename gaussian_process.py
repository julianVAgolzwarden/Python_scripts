import matplotlib.pyplot as plt
import numpy as np
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, ConstantKernel

#accessing datas
with open('ci_coefficients.txt') as C:
    content = C.readlines()

y = np.zeros((1000-1, 1))
for i in range(1, 1000):
    y[i-1] = float(content[i])


with open('PT_approximation.txt') as P:
    data = P.readlines()

X = np.zeros((1000-1, 1))
for j in range(1000-1):
    X[j] = float(data[j])

#gaussian process
kernel = ConstantKernel() + DotProduct(sigma_0=1) + WhiteKernel(noise_level=1)

x = np.atleast_2d(np.linspace(np.amax(X), np.amin(X), 1000)).T

gp = gaussian_process.GaussianProcessRegressor(kernel=kernel)
gp.fit(X, y)
y_pred, sigma = gp.predict(x, return_std=True)

plt.title(r'Dot Product, $\sigma_{0}$ = 1')
plt.plot(X, y, 'm.', label='values')
plt.plot(x, y_pred, 'k-', label='prediction')
plt.xlabel(r'excitation level/$\Delta$Ed')
plt.ylabel('ci coefficients')
plt.legend(loc='upper right')
plt.savefig('gaussian_PT_DP2.png')
plt.show()
