from functions import *
import re
import os
import pandas as pd

# parameters
hit_type = "hairpin"  # options: hairpin, spacer, ct_end, c_end
priority_is_max = False  # max or min hirpin energy, now always False

# paths
hairpins_file_path = r"../input/palindrom_analyzer_output.txt"
genome_file_path = r"../input/mpox_genome_seq.fna"
mutation_path = r"../input/snp_locations.xlsx"

# classes


class Hairpin:

    def __init__(self, l_s_m, position, score, palindrome):

        lsm = list(map(int, l_s_m.split("-")))
        self.spacer_length = lsm[1]
        self.stem_length = lsm[0]
        self.miss = lsm[2]
        self.length = self.spacer_length + self.stem_length * 2
        self.position = position
        self.score = score
        self.palindrome = palindrome
        self.sequence = palindrome.replace(" ", "")
        self.start = position - 1
        self.end = position - 2 + self.length
        self.stem_indexes = [(self.start, self.start + self.stem_length - 1), (self.end - self.stem_length + 1, self.end)]
        self.spacer_index = (self.start + self.stem_length, self.end - self.stem_length)

    def __eq__(self, other):
        return self.palindrome == other.palindrome and self.start == other.start and self.end == other.end

    def can_exist(self, second_hairpin):

        for i in self.stem_indexes:
            for j in second_hairpin.stem_indexes:
                if max(i[0], j[0]) <= min(i[1], j[1]):
                    return False
        return True

    def print_data(self):
        print(self.start, self.end, self.length, self.palindrome, self.spacer_index)



class Genome:

    def __init__(self, genome_file_path):

        s = ""
        if os.path.exists(genome_file_path):
            file = open(genome_file_path)
        else:
            raise FileNotFoundError(f"File not found: {genome_file_path}")
        self.name = file.readline()
        for line in file.readlines():
            line = line.rstrip()
            s += line
        self.sequence = s
        self.path = genome_file_path
        self.length = len(s)

    def targets(self):
        c_tc = [i.end() - 1 for i in re.finditer("TC", self.sequence)]
        g_ga = [i.start() for i in re.finditer("GA", self.sequence)]
        targets = sorted(list(set(c_tc + g_ga)))
        return targets

    def end_targets(self, used_hairpins):

        tc_end = []
        just_tc = [i.end() - 1 for i in re.finditer("TC", self.sequence)]
        for i in just_tc:
            for h in used_hairpins:
                if h.length == h.stem_length * 2 + 1 or h.length == h.stem_length * 2:
                    continue
                l, r = h.spacer_index
                if i == r:
                    tc_end.append(i)
                    # print(i
                    # h.print_data()
                    break

        ga_end = []
        just_ga = [i.start() for i in re.finditer("GA", self.sequence)]
        for i in just_ga:
            for h in used_hairpins:
                if h.length == h.stem_length * 2 + 1 or h.length == h.stem_length * 2:
                    continue
                l, r = h.spacer_index
                if i == l:
                    ga_end.append(i)
                    # print(i)
                    # h.print_data()
                    break
        return (sorted(tc_end + ga_end), sorted(just_tc + just_ga))

    def only_c_g(self, used_hairpins):
        c = [i.start() for i in re.finditer("C", self.sequence)]
        g = [i.start() for i in re.finditer("G", self.sequence)]
        c_end = []
        for i in c:
            for h in used_hairpins:
                if h.length == h.stem_length * 2:
                    continue
                l, r = h.spacer_index
                if i == r:
                    c_end.append(i)
                    # print(i)
                    # h.print_data()
                    break
        g_end = []
        for i in g:
            for h in used_hairpins:
                if h.length == h.stem_length * 2:
                    continue
                l, r = h.spacer_index
                if i == l:
                    g_end.append(i)
                    # print(i)
                    # h.print_data()
                    break

        return (sorted(c_end + g_end), sorted(c + g))


# function for loading hairpin to list


def load_hairpins(hairpins_file_path):

    hairpins_list = []
    if os.path.exists(hairpins_file_path):
        hairpins_file = open(hairpins_file_path)
    else:
        print("File not found")
        return 0
    hairpins_file.readline()
    for line in hairpins_file.readlines():
        line = line.split()
        h = Hairpin(line[0], int(line[1]), float(line[2]), " ".join(line[3:]))
        hairpins_list.append(h)
    hairpins_file.close()
    return hairpins_list


# load genome from file
genome = Genome(genome_file_path)

# load mutations from article
mutxls = pd.read_excel(mutation_path)
mutAPOBEC = mutxls[mutxls['isAPOBEC'] == 1]
#mutations_list = [int(i) - 1 for i in open(mutation_path).readlines()[1:]]
mutations_list = mutAPOBEC["Position"].to_list()
mut_cnt = len(mutations_list)

# load hairpins from file
all_hairpins_list = load_hairpins(hairpins_file_path)



