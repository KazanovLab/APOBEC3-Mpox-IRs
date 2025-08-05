import yaml
import subprocess

with open("config.yaml", "r") as file:
    params = yaml.safe_load(file)

if params["hairpin_selection_type"] == "all":
    subprocess.run(["python", "all_hairpins.py"])
else:
    subprocess.run(["python", "hairpin_groups.py"])
