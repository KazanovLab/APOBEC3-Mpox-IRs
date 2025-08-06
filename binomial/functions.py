import os
import re
import random
from math import comb
import matplotlib.pyplot as plt
from load import *
from decimal import *
import sys


def threshold(hist, percentile, length=None):

    if length == None:
        length = sum(hist.values())
    limit = length * (percentile / 100)
    summ = 0
    for k in sorted(hist.keys()):
        summ += hist[k]
        if summ >= limit:
            return k


# function for finding groups


def find_groups(hairpins_list):

    hairpins_list.sort(key=lambda x: [x.start, x.end])
    groups = []
    i = 0
    current_group = []
    n = len(hairpins_list)

    while i < n:
        hairpin = hairpins_list[i]
        current_group.append(hairpin)
        max_right = hairpin.end

        j = i + 1
        while j < n and hairpins_list[j].start <= max_right:
            second_hairpin = hairpins_list[j]
            max_right = max(max_right, second_hairpin.end)
            current_group.append(second_hairpin)
            j += 1

        i = j
        groups.append(current_group)
        current_group = []

    return groups


def percentile(hist, value, length=None):
    if length is None:
        length = sum(hist.values())
    total_sum = 0
    for i in sorted(hist, reverse=True):
        if i < value:
            break
        total_sum += hist[i]
    return (total_sum / length) * 100


# function for hit


def hit_or_not(hairpin, mut, indexes_list=None):
    if hit_type == "spacer":
        spacer_index = hairpin.spacer_index
        l, r = spacer_index
        return l <= mut <= r
    elif hit_type == "hairpin":
        return hairpin.start <= mut <= hairpin.end
    elif hit_type == "ct_end":
        return mut in indexes_list
    elif hit_type == "c_end":
        return mut in indexes_list
    else:
        print('Structure_type not in "hairpin", "spacer", "ct_end", "c_end"')
        sys.exit(1)


# function for choosing hairpins from group, max coverage area, dp

def max_coverage(hairpins_group):

    n = len(hairpins_group)
    hairpins_group.sort(key=lambda x: x.end)
    j_indexes = [-1] * n
    for i in range(n):
        h = hairpins_group[i]
        for j in range(i - 1, -1, -1):
            h2 = hairpins_group[j]
            if h2.end < h.start:
                j_indexes[i] = j
                break
    dp = [(0, -1)] * n
    dp[0] = (hairpins_group[0].length, 1)

    for i in range(1, n):
        h = hairpins_group[i]
        not_use = dp[i - 1][0]
        use = h.length + dp[j_indexes[i]][0]
        if use >= not_use:
            dp[i] = (use, 1)
        else:
            dp[i] = (not_use,   0)

    used = []
    i = n - 1
    while i >= 0:
        u = dp[i][1]
        if u == 1:
            used.append(hairpins_group[i])
            i = j_indexes[i]
        else:
            i = i - 1
    used = used[::-1]

    return used


# spacer max coverage
def max_coverage_spacer(hairpins_group):

    n = len(hairpins_group)
    hairpins_group.sort(key=lambda x: x.end)
    j_indexes = [-1] * n
    for i in range(n):
        h = hairpins_group[i]
        for j in range(i - 1, -1, -1):
            h2 = hairpins_group[j]
            if h2.end < h.start:
                j_indexes[i] = j
                break
    dp = [(0, -1)] * n
    dp[0] = (hairpins_group[0].spacer_length, 1)

    for i in range(1, n):
        h = hairpins_group[i]
        not_use = dp[i - 1][0]
        use = h.spacer_length + dp[j_indexes[i]][0]
        if use >= not_use:
            dp[i] = (use, 1)
        else:
            dp[i] = (not_use, 0)

    used = []
    i = n - 1
    while i >= 0:
        u = dp[i][1]
        if u == 1:
            used.append(hairpins_group[i])
            i = j_indexes[i]
        else:
            i = i - 1
    used = used[::-1]

    return used


# binom

def binom(p, n=56):

    hist = {i:0 for i in range(n + 1)}
    for m in range(n + 1):
        c = comb(n, m)
        result = c * (p ** m) * ((1 - p) ** (n - m))
        hist[m] = result
    return hist

































# random mutations function
#
# def random_mutations(hairpins_list, iterations, positions=None):
#     if positions == None:
#         positions = list(range(genome.length))
#     t = []
#     hits = dict([(i, 0) for i in range(mut_cnt + 1)])
#     for i in range(iterations):
#         # if i % 2000 == 0:
#         #     print(i)
#         random_mut = random.sample(positions, mut_cnt)
#         total_hits = 0
#         for mut in random_mut:
#             hit = 0
#             for hairpin in hairpins_list:
#                 if hairpin.start <= mut <= hairpin.end:
#                     hit = 1
#                     break
#             total_hits += hit
#         t.append(total_hits)
#         hits[total_hits] += 1
#     print(sum(t) / len(t))
#     print(hits)


