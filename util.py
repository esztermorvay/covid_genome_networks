import os
import random
import zipfile
import json

util_dir = "util_files"
graphs_dir = "graphs/"
gml_dir = "gml_files"

computer_using = "lenovo_laptop"
groups_to_use = [1,2,3]
if computer_using == "lenovo_laptop":
    groups_to_use = [7,8,9]
elif computer_using == "dell_laptop":
    groups_to_use = [4,5,6]

def get_all_files_in_dir_as_list(dir):
    # get all the files in the directory
    files = os.listdir(dir)
    # get all the files in the directory as a list
    files = list(files)
    return files

def extract_file_from_zip(zip_file_path, dir_to_save, file_to_extract = 'ncbi_dataset/data/genomic.fna'):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        # zip_ref.extractall('temp')
        zip_ref.extract(file_to_extract, dir_to_save)

def take_random_sample(x, n):
    # get a random sample of size n from list x
    return random.sample(x, n)

def divide_names_into_groups(file_names, n):
    # return a list of lists dividing file_names into size n separate lists
    size = n
    start_slice = 0
    end_slice = size
    groups = []
    while (start_slice < len(file_names)):
        if end_slice > len(file_names):
            end_slice = len(file_names)
        groups.append(file_names[start_slice:end_slice])
        start_slice += size
        end_slice += size
    i = 1
    # print("[")
    # for group in groups:
    #     print("\t", group, ",")
    #     i += 1
    # print("]")
    with open("util_files/multiprocessing_groups.json", 'w') as f:
        # json_str = json.dumps(groups)
        # json.dump(json_str, f)
        json.dump(groups, f, indent=4)
    return groups

def create_combos(groups):
    combos = {}
    i = 1
    print("[")
    for group in groups:
        j = 1
        for group2 in groups:
            comb = f"group{i}_group{j}"
            combos[comb] = False
            j += 1
        i += 1
    with open(util_dir + "/combinations.json", "w") as f:
        json.dump(combos, f, indent=4)

    # get all individual group combos
    i = 1
    for group in groups:
        print("getting combintaions for group ", i)
        individual_combos = get_combos_btwn_groups(group)
        with open(f"util_files/combinations/group{i}.json", "w") as f:
            json.dump(individual_combos, f, indent=4)
        i+=1
    # get all combos between groups
    for i in range(0,len(groups)):
        for j in range(i+1, len(groups)):
            print("getting combintaions for group ", i+1, j+1)

            group_combos = get_combos_btwn_groups(groups[i], groups[j])
            with open(f"util_files/combinations/group{i+1}group{j+1}.json", "w") as f:
                json.dump(group_combos, f, indent=4)
    return True

def get_combos_btwn_groups(group1, group2=None):
    # return a list w all combos
    output = []
    if group2 is None:
        for i in range(0, len(group1)):
            for j in range(i+1, len(group1)):
                combo = [group1[i], group1[j]]
                output.append(combo)
    else:
        for group in group1:
            for other in group2:
                combo = [group,other]
                output.append(combo)
    return output

