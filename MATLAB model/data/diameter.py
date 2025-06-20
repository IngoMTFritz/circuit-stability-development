# Creating a Python script that reads a .hoc file and computes the mean diameter

script_content = """\
import re

def extract_diameters(filename):
    diameters = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'pt3dadd\\(.*?,.*?,.*?,(.*?)\\)', line)
            if match:
                diameters.append(float(match.group(1)))
    return diameters

def compute_mean_diameter(filename):
    diameters = extract_diameters(filename)
    if diameters:
        mean_diameter = sum(diameters) / len(diameters)
        print(f"Mean Diameter: {mean_diameter:.3f}")
    else:
        print("No diameter values found in the file.")

if __name__ == "__main__":
    filename = "SkeletonTree_ORR_aCC_48h1_NL.hoc"  # Ensure this file is in the same directory as this script
    compute_mean_diameter(filename)
"""

# Saving the script
script_path = "/mnt/data/compute_mean_diameter.py"
with open(script_path, "w") as script_file:
    script_file.write(script_content)

# Providing the download link
script_path
