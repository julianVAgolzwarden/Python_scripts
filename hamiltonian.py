import matplotlib.pyplot as plt
import numpy as np
import time

t0 = time.time()

with open('hamiltonianData.txt') as h:
    content = h.readlines()

n_hamiltonian = len(content) - 1

#construct hamiltonian matrix
index_i = np.zeros(n_hamiltonian, dtype = np.int64)
index_j = np.zeros(n_hamiltonian, dtype = np.int64)

for x in range(1, n_hamiltonian + 1):
    index_i[x-1] = int(content[x].split()[0])
    index_j[x-1] = int(content[x].split()[1])

n_ij = np.amax(index_i) + 1

print index_j[14467]
exit(0)

hij = np.zeros((n_ij, n_ij))
for y in range(1, n_hamiltonian + 1):
    hij[index_i[y-1], index_j[y-1]] = float(content[y].split()[2])
    hij[index_j[y-1], index_i[y-1]] = float(content[y].split()[2])

#access Determinants
with open('occupation.txt') as O:
    data = O.readlines()

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

for e in excitation:
    if e <= 1:
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, excitation.index(e)]
    elif e == 2:
        hf_detx2 = np.tile(hf_det, (2, 1))
        mask2 = np.in1d(hf_det, dets[excitation.index(e)])
        change2 = np.where(mask2 == False)
        det201 = a_dagger(a(hf_detx2[0], change2[0][0]), dets[excitation.index(e), change2[0][0]])
        det201 = a_dagger(a(hf_detx2[1], change2[0][1]), dets[excitation.index(e), change2[0][1]])
        f_det201 = np.sort(det201).tolist()
        f_det202 = np.sort(det202).tolist()
        index201 = dets_new.index(f_det201)
        index202 = dets_new.index(f_det202)
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index201]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index202]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index201, excitation.index(e)]
        paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index202, excitation.index(e)]
    elif e == 3:
        hf_detx3 = np.tile(hf_det, (3, 1))
        mask3 = np.in1d(hf_det, dets[excitation.index(e)])
        change3 = np.where(mask3 == False)
        for i in range(3):
            det3i = a_dagger(a(hf_detx3[i], change3[0][i]), dets[excitation.index(e), change3[0][i]])
            det3ix = np.tile(det3i, (3, 1))
            f_det3i = np.sort(det3ix[0]).tolist()
            index3i = dets_new.index(f_det3i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index3i]
            for j in range(3):
                if j != i:
                    det3j = a_dagger(a(det3ix[j], change3[0][j]), dets[excitation.index(e), change3[0][j]])
                    f_det3j = np.sort(det3j).tolist()
                    index3j = dets_new.index(f_det3j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index3i, index3j]
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index3j, excitation.index(e)]
    elif e == 4:
        hf_detx4 = np.tile(hf_det, (4, 1))
        mask4 = np.in1d(hf_det, dets[excitation.index(e)])
        change4 = np.where(mask4 == False)
        for i in range(4):
            det4i = a_dagger(a(hf_detx4[i], change4[0][i]), dets[excitation.index(e), change4[0][i]])
            det4ix = np.tile(det4i, (4, 1))
            f_det4i = np.sort(det4ix[0]).tolist()
            index4i = dets_new.index(f_det4i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index4i]
            for j in range(4):
                if j != i:
                    det4j = a_dagger(a(det4ix[j], change4[0][j]), dets[excitation.index(e), change4[0][j]])
                    det4jx = np.tile(det4j, (4, 1))
                    f_det4j = np.sort(det4jx[0]).tolist()
                    index4j = dets_new.index(f_det4j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index4i, index4j]
                    for k in range(4):
                        if k != j and k != i:
                            det4k = a_dagger(a(det4jx[k], change4[0][k]), dets[excitation.index(e), change4[0][k]])
                            f_det4k = np.sort(det4k).tolist()
                            index4k = dets_new.index(f_det4k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index4j, index4k]
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index4k, excitation.index(e)]
    elif e == 5:
        hf_detx5 = np.tile(hf_det, (5, 1))
        mask5 = np.in1d(hf_det, dets[excitation.index(e)])
        change5 = np.where(mask5 == False)
        for i in range(5):
            det5i = a_dagger(a(hf_detx5[i], change5[0][i]), dets[excitation.index(e), change5[0][i]])
            det5ix = np.tile(det5i, (5, 1))
            f_det5i = np.sort(det5ix[0]).tolist()
            index5i = dets_new.index(f_det5i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index5i]
            for j in range(5):
                if j != i:
                    det5j = a_dagger(a(det5ix[j], change5[0][j]), dets[excitation.index(e), change5[0][j]])
                    det5jx = np.tile(det5j, (5, 1))
                    f_det5j = np.sort(det5jx[0]).tolist()
                    index5j = dets_new.index(f_det5j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index5i, index5j]
                    for k in range(5):
                        if k != j and k != i:
                            det5k = a_dagger(a(det5jx[k], change5[0][k]), dets[excitation.index(e), change5[0][k]])
                            det5kx = np.tile(det5k, (5, 1))
                            f_det5k = np.sort(det5kx[0]).tolist()
                            index5k = dets_new.index(f_det5k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index5j, index5k]
                            for b in range(5):
                                if b != j and b != i and b != k:
                                    det5b = a_dagger(a(det5kx[b], change5[0][b]), dets[excitation.index(e), change5[0][b]])
                                    f_det5b = np.sort(det5b).tolist()
                                    index5b = dets_new.index(f_det5b)
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index5k, index5b]
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index5b, excitation.index(e)]
    elif e == 6:
        hf_detx6 = np.tile(hf_det, (6, 1))
        mask6 = np.in1d(hf_det, dets[excitation.index(e)])
        change6 = np.where(mask6 == False)
        for i in range(6):
            det6i = a_dagger(a(hf_detx6[i], change6[0][i]), dets[excitation.index(e), change6[0][i]])
            det6ix = np.tile(det6i, (6, 1))
            f_det6i = np.sort(det6ix[0]).tolist()
            index6i = dets_new.index(f_det6i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index6i]
            for j in range(6):
                if j != i:
                    det6j = a_dagger(a(det6ix[j], change6[0][j]), dets[excitation.index(e), change6[0][j]])
                    det6jx = np.tile(det6j, (6, 1))
                    f_det6j = np.sort(det6jx[0]).tolist()
                    index6j = dets_new.index(f_det6j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index6i, index6j]
                    #excitation j -> k
                    for k in range(6):
                        if k != j and k != i:
                            det6k = a_dagger(a(det6jx[k], change6[0][k]), dets[excitation.index(e), change6[0][k]])
                            det6kx = np.tile(det6k, (6, 1))
                            f_det6k = np.sort(det6kx[0]).tolist()
                            index6k = dets_new.index(f_det6k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index6j, index6k]
                            for b in range(6):
                                if b != j and b != i and b != k:
                                    det6b = a_dagger(a(det6kx[b], change6[0][b]), dets[excitation.index(e), change6[0][b]])
                                    det6bx = np.tile(det6b, (6, 1))
                                    f_det6b = np.sort(det6bx[0]).tolist()
                                    index6b = dets_new.index(f_det6b)
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index6k, index6b]
                                    for c in range(6):
                                        if c != j and c!= i and c != k and c != b:
                                            det6c = a_dagger(a(det6bx[b], change6[0][b]), dets[excitation.index(e), change6[0][b]])
                                            f_det6c = np.sort(det6c).tolist()
                                            index6c = dets_new.index(f_det6c)
                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index6b, index6c]
                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index6c, excitation.index(e)]
    elif e == 7:
        hf_detx7 = np.tile(hf_det, (7, 1))
        mask7 = np.in1d(hf_det, dets[excitation.index(e)])
        change7 = np.where(mask7 == False)
        for i in range(7):
            det7i = a_dagger(a(hf_detx7[i], change7[0][i]), dets[excitation.index(e), change7[0][i]])
            det7ix = np.tile(det7i, (7, 1))
            f_det7i = np.sort(det7ix[0]).tolist()
            index7i = dets_new.index(f_det7i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index7i]
            for j in range(7):
                if j != i:
                    det7j = a_dagger(a(det7ix[j], change7[0][j]), dets[excitation.index(e), change7[0][j]])
                    det7jx = np.tile(det7j, (7, 1))
                    f_det7j = np.sort(det7jx[0]).tolist()
                    index7j = dets_new.index(f_det7j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7i, index7j]
                    for k in range(7):
                        if k != j and k != i:
                            det7k = a_dagger(a(det7jx[k], change7[0][k]), dets[excitation.index(e), change7[0][k]])
                            det7kx = np.tile(det7k, (7, 1))
                            f_det7k = np.sort(det7kx[0]).tolist()
                            index7k = dets_new.index(f_det7k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7j, index7k]
                            for b in range(7):
                                if b != j and b != i and b != k:
                                    det7b = a_dagger(a(det7kx[b], change7[0][b]), dets[excitation.index(e), change7[0][b]])
                                    det7bx = np.tile(det7b, (7, 1))
                                    f_det7b = np.sort(det7bx[0]).tolist()
                                    index7b = dets_new.index(f_det7b)
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7k, index7b]
                                    for c in range(7):
                                        if c != j and c != i and c != k and c != b:
                                            det7c = a_dagger(a(det7bx[c], change7[0][c]), dets[excitation.index(e), change7[0][c]])
                                            det7cx = np.tile(det7c, (7, 1))
                                            f_det7c = np.sort(det7cx[0]).tolist()
                                            index7c = dets_new.index(f_det7c)
                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7b, index7c]
                                            for d in range(7):
                                                if d != j and d != i and d != k and d != b and d != c:
                                                    det7d = a_dagger(a(det7cx[d], change7[0][d]), dets[excitation.index(e), change7[0][d]])
                                                    f_det7d = np.sort(det7d).tolist()
                                                    index7d = dets_new.index(f_det7d)
                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7c, index7d]
                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index7d, excitation.index(e)]
    elif e == 8:
        hf_detx8 = np.tile(hf_det, (8, 1))
        mask8 = np.in1d(hf_det, dets[excitation.index(e)])
        change8 = np.where(mask8 == False)
        for i in range(8):
            det8i = a_dagger(a(hf_detx8[i], change8[0][i]), dets[excitation.index(e), change8[0][i]])
            det8ix = np.tile(det8i, (8, 1))
            f_det8i = np.sort(det8ix[0]).tolist()
            index8i = dets_new.index(f_det8i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index8i]
            for j in range(8):
                if j != i:
                    det8j = a_dagger(a(det8ix[j], change8[0][j]), dets[excitation.index(e), change8[0][j]])
                    det8jx = np.tile(det8j, (8, 1))
                    f_det8j = np.sort(det8jx[0]).tolist()
                    index8j = dets_new.index(f_det8j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8i, index8j]
                    for k in range(8):
                        if k != j and k != i:
                            det8k = a_dagger(a(det8jx[k], change8[0][k]), dets[excitation.index(e), change8[0][k]])
                            det8kx = np.tile(det8k, (8, 1))
                            f_det8k = np.sort(det8kx[0]).tolist()
                            index8k = dets_new.index(f_det8k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8j, index8k]
                            for b in range(8):
                                if b != j and b != i and b != k:
                                    det8b = a_dagger(a(det8kx[b], change8[0][b]), dets[excitation.index(e), change8[0][b]])
                                    det8bx = np.tile(det8b, (8, 1))
                                    f_det8b = np.sort(det8bx[0]).tolist()
                                    index8b = dets_new.index(f_det8b)
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8k, index8b]
                                    for c in range(8):
                                        if c != j and c != i and c != k and c != b:
                                            det8c = a_dagger(a(det8bx[c], change8[0][c]), dets[excitation.index(e), change8[0][c]])
                                            det8cx = np.tile(det8c, (8, 1))
                                            f_det8c = np.sort(det8cx[0]).tolist()
                                            index8c = dets_new.index(f_det8c)
                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8b, index8c]
                                            for d in range(8):
                                                if d != j and d != i and d != k and d != b and d != c:
                                                    det8d = a_dagger(a(det8cx[d], change8[0][d]), dets[excitation.index(e), change8[0][d]])
                                                    det8dx = np.tile(det8d, (8, 1))
                                                    f_det8d = np.sort(det8dx[0]).tolist()
                                                    index8d = dets_new.index(f_det8d)
                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8c, index8d]
                                                    for s in range(8):
                                                        if s != j and s != i and s != k and s != b and s != c and s != d:
                                                            det8s = a_dagger(a(det8dx[s], change8[0][s]), dets[excitation.index(e), change8[0][s]])
                                                            f_det8s = np.sort(det8s).tolist()
                                                            index8s = dets_new.index(f_det8s)
                                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8d, index8s]
                                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8s, excitation.index(e)]
    elif e == 9:
        hf_detx9 = np.tile(hf_det, (9, 1))
        mask9 = np.in1d(hf_det, dets[excitation.index(e)])
        change9 = np.where(mask9 == False)
        for i in range(9):
            det9i = a_dagger(a(hf_detx9[i], change9[0][i]), dets[excitation.index(e), change9[0][i]])
            det9ix = np.tile(det9i, (9, 1))
            f_det9i = np.sort(det9ix[0]).tolist()
            index9i = dets_new.index(f_det9i)
            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[0, index9i]
            for j in range(9):
                if j != i:
                    det9j = a_dagger(a(det9ix[j], change9[0][j]), dets[excitation.index(e), change9[0][j]])
                    det9jx = np.tile(det9j, (9, 1))
                    f_det9j = np.sort(det9jx[0]).tolist()
                    index9j = dets_new.index(f_det9j)
                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9i, index9j]
                    for k in range(9):
                        if k != j and k != i:
                            det9k = a_dagger(a(det9jx[k], change9[0][k]), dets[excitation.index(e), change9[0][k]])
                            det9kx = np.tile(det9k, (9, 1))
                            f_det9k = np.sort(det9kx[0]).tolist()
                            index9k = dets_new.index(f_det9k)
                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[inde9j, index9k]
                            for b in range(9):
                                if b != j and b != i and b != k:
                                    det9b = a_dagger(a(det9kx[b], change9[0][b]), dets[excitation.index(e), change9[0][b]])
                                    det9bx = np.tile(det9b, (9, 1))
                                    f_det9b = np.sort(det9bx[0]).tolist()
                                    index9b = dets_new.index(f_det9b)
                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9k, index9b]
                                    for c in range(9):
                                        if c != j and c != i and c != k and c != b:
                                            det9c = a_dagger(a(det9bx[c], change9[0][c]), dets[excitation.index(e), change9[0][c]])
                                            det9cx = np.tile(det9c, (9, 1))
                                            f_det9c = np.sort(det9cx[0]).tolist()
                                            index9c = dets_new.index(f_det9c)
                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9b, index9c]
                                            for d in range(9):
                                                if d != j and d != i and d != k and d != b and d != c:
                                                    det9d = a_dagger(a(det9cx[d], change9[0][d]), dets[excitation.index(e), change9[0][d]])
                                                    det9dx = np.tile(det9d, (9, 1))
                                                    f_det9d = np.sort(det9dx[0]).tolist()
                                                    index9d = dets_new.index(f_det9d)
                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index8c, index8d]
                                                    for s in range(9):
                                                        if s != j and s != i and s != k and s != b and s != c and s != d:
                                                            det9s = a_dagger(a(det9dx[s], change9[0][s]), dets[excitation.index(e), change9[0][s]])
                                                            det9sx = np.tile(det9s, (9, 1))
                                                            f_det9s = np.sort(det9sx[0]).tolist()
                                                            index9s = dets_new.index(f_det9s)
                                                            paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9d, index9s]
                                                            for t in range(9):
                                                                if t != j and t != i and t != k and t != b and t != c and t != d and t != s:
                                                                    det9t = a_dagger(a(det9sx[t], change9[0][t]), dets[excitation.index(e), change9[0][t]])
                                                                    f_det9t = np.sort(det9t).tolist()
                                                                    index9t = dets_new.index(f_det9t)
                                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9s, index9t]
                                                                    paths[excitation.index(e)] = paths[excitation.index(e)] + hij[index9t, excitation.index(e)]

np.savetxt('excitation_pathways1.txt', paths, fmt='%1f')
print "Total Time:", time.time() - t0
