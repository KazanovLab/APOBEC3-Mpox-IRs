from functions import *
import random
from load import *
from copy import deepcopy
import matplotlib.pyplot as plt
from itertools import product

# parameters
choice_type = "max_cov" # options: greedy, max_cov, min_cov, most_stable

# min len hairpin
def min_coverage(group):

    return [min(group, key=lambda x: (x.length, x.score, x.start))]


# function for choosing hairpins from group, greedy algorithm

def greedy_choose(hairpins_group, taken=None, skip_sort=False):

    if taken is None:
        taken = []
    hairpins_group_copy = deepcopy(hairpins_group)
    taken_copy = taken.copy()

    if not skip_sort:
        k = -1 if priority_is_max else 1
        hairpins_group_copy.sort(key=lambda x: x.score * k)
        n = len(hairpins_group_copy)
        i = 0
        equal_groups = []
        while i < n:
            current_group = [hairpins_group_copy[i]]
            j = i + 1
            while j < n and hairpins_group_copy[i].score == hairpins_group_copy[j].score:
                current_group.append(hairpins_group_copy[j])
                j += 1
            current_group.sort(key=lambda x: x.length * -1)
            equal_groups.append(current_group)
            i = j
        hairpins_group_copy = equal_groups

    for p in range(len(hairpins_group_copy)):

        equal_group = hairpins_group_copy[p]
        if equal_group == []:
            continue
        usable = []
        for hairpin in equal_group:
            if all([hairpin.can_exist(j) for j in taken_copy]):
                usable.append(hairpin)
        n = len(usable)
        if n == 0:
            continue
        elif n == 1:
            taken_copy.append(usable[0])
        elif n <= 11:
            variants = []
            for hairpin in usable:
                new_hairpins_group = deepcopy(hairpins_group_copy)
                new_hairpins_group[p] = [h for h in new_hairpins_group[p] if h is not hairpin]
                temp = greedy_choose(new_hairpins_group[p:], taken_copy + [hairpin], True)
                variants.append(temp)
            m = max(variants, key=lambda x: len(x))
            return m
        else:
            hairpin = usable[0]
            new_hairpins_group = deepcopy(hairpins_group_copy)
            new_hairpins_group[p] = [h for h in new_hairpins_group[p] if h is not hairpin]
            m = greedy_choose(new_hairpins_group[p:], taken_copy + [hairpin], True)
            return m

    return taken_copy


# choose most stable


def most_stable(group):

    k = -1 if priority_is_max else 1
    return [min(group, key=lambda x: (x.score * k, k * x.length, x.start))]



# fixing hairpins

groups = find_groups(all_hairpins_list)
print(len(groups))
fixed_hairpins = []
for group in groups:
    choice = []
    if hit_type == "hairpin":
        if choice_type == "greedy":
            choice = greedy_choose(group)
        if choice_type == "max_cov":
            choice = max_coverage(group)
        if choice_type == "most_stable":
            choice = most_stable(group)
        if choice_type == "min_cov":
            choice = min_coverage(group)
    else:
        if choice_type == "greedy":
            choice = greedy_choose(group)
        if choice_type == "max_cov":
            choice = max_coverage_spacer(group)
        if choice_type == "most_stable":
            choice = most_stable(group)
        if choice_type == "min_cov":
            choice = min_coverage(group)

    fixed_hairpins += choice
    # print("group")
    # for i in group:
    #     print(i.score, i.start, i.end, i.length)
    # print("choice")
    # for i in choice:
    #     print(i.score, i.start, i.end, i.length)

fixed_hairpins.sort(key=lambda x: x.start)
total_len = sum([h.length for h in fixed_hairpins])
spacer_len = sum([h.spacer_length for h in fixed_hairpins])

# tp count

end_targets = []
targets_p = 0
if hit_type == "hairpin" or hit_type == "spacer":
    j = 0
    for t in genome.targets():
        hit = 0
        for h in fixed_hairpins:
            if hit_or_not(h, t):
                hit = 1
                break
        j += hit
    targets_p = j / len(genome.targets())

if hit_type == "ct_end":
    end_targets, all_positions = genome.end_targets(fixed_hairpins)
    targets_p = len(end_targets) / len(all_positions)

if hit_type == "c_end":
    end_targets, all_positions = genome.only_c_g(fixed_hairpins)
    targets_p = len(end_targets) / len(all_positions)

print(len(end_targets))
print("tp:", targets_p * 100)
print(f"{len(fixed_hairpins)} of {len(all_hairpins_list)} hairpins are chose")
print(f"coverage percentage of fixed hairpins: {total_len * 100/genome.length}")
print(f"coverage percentage of fixed hairpins(spacers): {spacer_len * 100/genome.length}")


# real hits


hairpin_hits = []
total_real_hits = 0
for mut in mutations_list:
    hit = 0
    for hairpin in fixed_hairpins:
        if hit_or_not(hairpin, mut, end_targets):
            # print(mut, genome.sequence[mut-1:mut + 2])
            # hairpin_hits.append((mut, hairpin))
            hit = 1
            break
    total_real_hits += hit
print("total real hits", total_real_hits)

# check

# for i in hairpin_hits:
#     h = i[1]
    # print("mutation:", i[0], "  ", "hairpin:", h.length, h.start, h.end, h.score)
# file = open(r"C:\Users\fedor\Downloads\check_pins (1) (1)")
# nums = []
# for line in file:
#     line = line.replace("Start:", "")
#     line = line.replace(", End:", " ")
#     nums.append([int(i) - 1 for i in line.split()])
# nums.sort()
# fixed_hairpins.sort(key=lambda x: x.start)
# flag = True
# for i in range(3791):
#     s1, e1 = nums[i][0], nums[i][1]
#     h = fixed_hairpins[i]
#     s2, e2 = h.start, h.end
#     # print(s1, s2, e1, e2)
#     if s1 != s2 or e1 != e2:
#         print(s1, e1, s2, e2)
#         flag = False
# print(flag)


# random mutations

# hits_r = {i: 0 for i in range(mut_cnt + 1)}
# for i in range(10_000):
#     if i % 1000 == 0 and i:
#         print(i)
#     random_mut = random.sample(genome.targets(), mut_cnt)
#     total_hits = 0
#     for mut in random_mut:
#         hit = 0
#         for hairpin in fixed_hairpins:
#             if hit_or_not(hairpin, mut):
#                 hit = 1
#                 break
#         total_hits += hit
#     hits_r[total_hits] += 1
# print(hits_r)


# binom

hits_b = binom(targets_p)


# hist

keys = list(hits_b.keys())
values = list(hits_b.values())

plt.bar(keys, values)
# plt.xlabel("mut")
# plt.ylabel("experiment")
# plt.title("random")
positions = [i for i in range(0, len(keys))]
labels = [keys[i] for i in range(0, len(keys))]
plt.xticks(positions, labels, rotation=90)
plt.grid(True, axis="y", linestyle="--", alpha=0.7)


# threshold count random

# print("random")
# right = threshold(hits_r, 95)
# left = threshold(hits_r, 5)
# print(f"left threshold: {left}")
# print(f"right threshold: {right}")
# p = percentile(hits_r, total_real_hits)
# print(f"p-value: {p}")
# print()

# threshold count binom

print("binom")
right = threshold(hits_b, 95)
left = threshold(hits_b, 5)
print(f"left threshold: {left}")
print(f"right threshold: {right}")
p = percentile(hits_b, total_real_hits)
print(f"p-value: {p}")

plt.show()









