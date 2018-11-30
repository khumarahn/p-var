#!/bin/env python3

import math

def p_var_backbone(path_size, p, path_dist):
    if path_size == 0:
        return -math.inf
    elif path_size == 1:
        return 0

    s = path_size - 1
    N = 1
    while s >> N != 0:
        N += 1

    ind = [0.0] * s
    def ind_n(j, n):
        return (s >> n) + (j >> n)
    def ind_k(j, n):
        return min(((j >> n) << n) + (1 << (n-1)), s);

    max_p_var = 0.0
    run_p_var = [0.0] * path_size

    for j in range(0, path_size):
        for n in range(1, N + 1):
            if not(j >> n == s >> n and (s >> (n-1)) % 2 == 0):
                ind[ind_n(j, n)] = max(ind[ind_n(j, n)], path_dist(ind_k(j, n), j))
        if j == 0:
            continue

        m = j - 1
        delta = 0.0
        delta_m = j
        n = N
        while True:
            while n > 0 and m >> n == s >> n and (s >> (n-1)) % 2 == 0:
                n -= 1;

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
