from load import *
from hairpin_groups import *
from functions import *
hlist = fixed_hairpins
hlist.sort(key=lambda x: x.start)
type = "max_coverage"
output_path = r"C:\Mpox_hairpins\used_hairpins"
output_file = open(output_path + rf"\{type}", "w")
output_file.write("start\tend\n")
for h in hlist:
    output_file.write(f"{h.start}\t{h.end}\n")
output_file.close()

