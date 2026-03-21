import yaml
import subprocess
from all_hairpins import *
from hairpin_groups import *

with open("config.yaml", "r") as file:
    params = yaml.safe_load(file)

if params["hairpin_selection_type"] == "all":
    #subprocess.run(["python", "all_hairpins.py"])
    pval = calculate_pval_all()
else:
    #subprocess.run(["python", "hairpin_groups.py"])
    pval = calculate_pval_groups()
    
print(pval)
