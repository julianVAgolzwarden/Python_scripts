import matplotlib.pyplot as plt
import numpy as np
import time

t0 = time.time()

with open('hamiltonianData.txt') as h:
    content = h.readlines()

n_hamiltonian = len(content) - 1

#construct hamiltonian matrix
hij = np.zeros((n_hamiltonian+1, 3))

for x in range(1, n_hamiltonian + 1):
    for y in range(3):
        hij[x, y] = float(content[x].split()[y])

#access Determinants
with open('Determinants.txt') as O:
    data = O.readlines()

n_ij = int(np.amax(hij) + 1)
n_elec = 12

dets = np.zeros((n_ij, n_elec), dtype = np.int64)

for i in range(n_ij):
    for j in range(n_elec):
        dets[i, j] = int(data[i].split()[j])

#access excitation distances
with open("Excitation_levels.txt") as E:
    values = E.readlines()

excitation = []
for i in range(n_ij):
    excitation.append(int(values[i]))

#define creation and annhilation operators
def a(L, n):
    L[n] = 40
    return L

def a_dagger(L, m):
    L[np.argwhere(L == 40)] = m
    return L

#calculate excitation paths
paths = np.zeros(n_ij,)
hf_det = dets[0]
dets_new = dets.tolist()

def ha(A):
    for x in range(len(hij)):
        mask = np.in1d(A, hij[x])
        if np.sum(mask) == 2:
            return x
        else:
            return n_hamiltonian

#function finding the Determinant
def find_index(z):
    if len(z) == len(np.unique(z)):
        for x in range(len(dets)):
            mask = np.in1d(z, dets[x])
            if np.sum(mask) == 12:
                return dets_new.index(z.tolist())
            else:
                return n_ij
    else:
        return n_ij

for e in excitation:
    if e <= 1:
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[ha(np.array([0, excitation.index(e)])), 2]
    elif e == 2:
        hf_detx2 = np.tile(hf_det, (2, 1))
        mask2 = np.in1d(hf_det, dets[excitation.index(e)])
        change2 = np.where(mask2 == False)
        det201 = a_dagger(a(hf_detx2[0], change2[0][0]), dets[excitation.index(e), change2[0][0]])
        det202 = a_dagger(a(hf_detx2[1], change2[0][1]), dets[excitation.index(e), change2[0][1]])
        f_det201 = np.sort(det201)
        f_det202 = np.sort(det202)
        index201 = find_index(f_det201)
        index202 = find_index(f_det202)
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[ha(np.array([0, index201])), 2]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[ha(np.array([0, index202])), 2]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[ha(np.array([index201, excitation.index(e)])), 2]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[ha(np.array([index202, excitation.index(e)])), 2]

np.savetxt('excitation_pathways2.txt', paths, fmt='%1f')
print "Total Time:", time.time() - t0
