#!/bin/env python3

import math

def p_var_backbone(path_size, p, path_dist):
    if path_size == 0:
        return -math.inf
    elif path_size == 1:
        return 0

    N = 1
    Nn = path_size
    while True:
        Nn >>= 1
        if Nn == 0:
            break
        N += 1
    
    ind_offset = [0] * N
    ip = path_size
    for n in range(0, N):
        ip = (ip % 2) + (ip >> 1)
        ind_offset[n] = ip + (0 if n==0 else ind_offset[n-1])
    ind_size = ind_offset[N-1]
    for n in range(0, N):
        ind_offset[n] = ind_size - ind_offset[n]
    ind = [0.0] * ind_size
    def ind_n(j, n):
        return ind_offset[n-1] + (j >> n)
    def ind_k(j, n):
        return min(((j >> n) << n) + (1 << (n-1)), path_size - 1);

    max_p_var = 0.0
    run_p_var = [0.0] * path_size

    for j in range(0, path_size):
        for n in range(1, N + 1):
            ind[ind_n(j,n)] = max(ind[ind_n(j, n)], path_dist(ind_k(j, n), j))
        if j == 0:
            continue

        m = j - 1
        delta = 0.0
        delta_m = j
        n = N
        while True:
            skip = False
            if n > 0:
                iid = ind[ind_n(m, n)] + path_dist(ind_k(m, n), j)
                if delta >= iid:
                    skip = True
                elif m < delta_m:
                    delta = pow(max_p_var - run_p_var[m], 1. / p)
                    delta_m = m
                    if delta >= iid:
                        skip = True

            if skip:
                k = (m >> n) << n
                if k > 0:
                    m = k - 1
                    while n < N and (k >> n) % 2 == 0:
                        n += 1
                else:
                    break
            else:
                if n > 1:
                    n -= 1
                else:
                    d = path_dist(m, j)
                    if d > delta:
                        max_p_var = max(max_p_var, run_p_var[m] + pow(d, p))
                    if m > 0:
                        while n < N and (m >> n) % 2 == 0:
                            n += 1
                        m -= 1
                    else:
                        break
        run_p_var[j] = max_p_var;
    return run_p_var[-1]

path = [[0,0], [0,1], [1,1], [1,0], [0,0]]
dist = lambda a, b: math.sqrt(pow(path[b][0] - path[a][0], 2) + pow(path[b][1] - path[a][1], 2))

print("Square path: ", path)
p = 1.0
while p <= 4.0:
    pv = p_var_backbone(len(path), p, dist)
    pv_ref = 4 if p <= 2 else pow(2, p/2 + 1)
    pv_err = abs(pv - pv_ref)
    print(f'{p:5f}-variation: {pv:8f}, error {pv_err:8e}')
    p += 0.125
