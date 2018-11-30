#!/bin/env python3

import math

# for tests
import random
import time

def p_var_backbone(path_size, p, path_dist):
    # Input:
    # * path_size >= 0 integer
    # * p >= 1 real
    # * path_dist: metric on the set {0,...,path_dist-1};
    #   I.e. path_dist(a,b) needs to be defined and nonnegative
    #   for all integer 0 <= a,b < path_dist, be symmetric and
    #   satisfy the triangle inequality:
    #   * path_dist(a,b) = path_dist(b,a)
    #   * path_dist(a,b) + path_dist(b,c) >= path_dist(a,c)
    # Output:
    # * max sum_k path_dist(a_k, a_{k-1})^p
    # over all strictly increasing subsequences a_k of 0,...,path_size-1
    # Notes:
    # * if path_size == 0, the result is -math.inf
    # * if path_size == 1, the result is 0
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

def p_var_backbone_ref(path_size, p, path_dist):
    # Reference implementation of p_var_backbone, does not need the triangle inequality
    # but may be slow; obviously correct.
    if path_size == 0:
        return -math.inf
    elif path_size == 1:
        return 0
    cum_p_var = [0.0] * path_size
    for j in range(1, path_size):
        for m in range(0, j):
            cum_p_var[j] = max(cum_p_var[j], cum_p_var[m] + pow(path_dist(m, j), p));
    return cum_p_var[-1]

def ex_sq():
    # Example: unit square
    path = [[0,0], [0,1], [1,1], [1,0], [0,0]]
    dist = lambda a, b: math.sqrt(pow(path[b][0] - path[a][0], 2) + pow(path[b][1] - path[a][1], 2))

    print(f'\nSquare path: {path}\nwith L^2 distance')
    p = 1.0
    while p <= 4.0:
        pv = p_var_backbone(len(path), p, dist)
        pv_ref = p_var_backbone_ref(len(path), p, dist)
        pv_err = abs(pv - pv_ref)
        print(f'{p:5f}-variation: {pv:8f}, error {pv_err:8e}')
        p += 0.5

def ex_bm():
    # Example: Broownian motion made of iid -1/+1 increments
    n = 2500
    print(f'\nPoor man\'s Brownian path with {n} steps:')
    path = [0.0] * (n + 1)
    sigma = 1. / math.sqrt(n)
    for k in range(1, n + 1):
        path[k] = path[k - 1] + random.choice([-1, 1]) * sigma
    dist = lambda a, b: abs(path[b] - path[a])

    for p in [1.0, math.sqrt(2), 2.0, math.exp(1)]:
        pv_start = time.time()
        pv = p_var_backbone(len(path), p, dist)
        pv_time = time.time() - pv_start
        pv_ref_start = time.time()
        pv_ref = p_var_backbone_ref(len(path), p, dist)
        pv_ref_time = time.time() - pv_ref_start
        pv_err = abs(pv - pv_ref)
        print(f'{p:5f}-variation: {pv:8e}, error {pv_err:8e}; time: {pv_time:6f}; reference time: {pv_ref_time:6f}')

def ex_bm_long():
    # Example: Broownian motion made of iid -1/+1 increments
    print('\nVery poor man\'s very long Brownian motion:')
    for n in [0,1,10,100,1000,10000,100000,1000000, 10000000]:
        path = [0.0] * (n + 1)
        sigma = 1. / math.sqrt(max(n,1))
        for k in range(1, n + 1):
            path[k] = path[k - 1] + random.choice([-1, 1]) * sigma
        dist = lambda a, b: abs(path[b] - path[a])

        p = 2.25
        pv_start = time.time()
        pv = p_var_backbone(len(path), p, dist)
        pv_time = time.time() - pv_start
        print(f'{n:10d} steps: {p:5f}-variation: {pv:8e}, time: {pv_time:6f}')

ex_sq()
ex_bm()
ex_bm_long()
