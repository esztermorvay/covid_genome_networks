import os
import random
import zipfile
import json

util_dir = "util_files"
graphs_dir = "graphs/"
gml_dir = "gml_files"
# scores_dir = "scores/"
# combinations_dir = util_dir + "/combinations/"
scores_dir = "garbage_collection/scores/"
combinations_dir = "garbage_collection/combinations/"
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

def get_incomplete_files():
    files = get_all_files_in_dir_as_list(combinations_dir)
    done = {}
    score_files = get_all_files_in_dir_as_list(scores_dir)
    for file in files:
        with open(combinations_dir + file, 'r') as f:
            data = json.load(f)
            for pair in data:
                node1 = pair[0][11:-4]
                node2 = pair[1][11:-4]
                key = node1 + "_" + node2
                done[key] = False
         # get the corresponding group
        file = file[:-5]
        group1_name = file[:6]
        if len(file) > 6:
            group2_name = file[6:]
        else:
            group2_name = group1_name
        group = group1_name + "_" + group2_name
        # get the corresponding files from scores
        for score_file in score_files:
            if group in score_file:
                try:
                    with open(scores_dir + score_file, 'r') as f:
                        data = json.load(f)
                        for key in data:
                            done[key] = True
                except:
                    print("error with file: ", score_file)
    # get the keys that are false
    incomplete = []
    for key in done:
        if not done[key]:
            incomplete.append(key)
    incomplete_reversed = []
    for key in incomplete:
        node1 = key.split("_")[0]
        node2 = key.split("_")[1]
        incomplete_reversed.append(node2 + "_" + node1)
    # check if any of the ones in incomplete reversed are in done
    # for key in incomplete_reversed:
    #     if key in done:
    #         if done[key]:
    #             incomplete.remove(key)

    return incomplete
    # "SARS-CoV-2-" + name + ".zip"

def create_incomplete_combos(incomplete):
    incomplete_combos = []
    for key in incomplete:
        node1 = key.split("_")[0]
        node2 = key.split("_")[1]
        name1 = "SARS-CoV-2-" + node1 + ".zip"
        name2 = "SARS-CoV-2-" + node2 + ".zip"
        incomplete_combos.append([name1, name2])
    return incomplete_combos

def create_files_from_combos(num_files, combos):
    # divide combos into num_fies lists
    # create a file for each list
    # return a list of the file names
    groups = []

    size = num_files
    start_slice = 0
    end_slice = size
    groups = []
    while (start_slice < len(combos)):
        if end_slice > len(combos):
            end_slice = len(combos)
        groups.append(combos[start_slice:end_slice])
        start_slice += size
        end_slice += size
    # check the size of groups
    total_size = 0
    i = 1
    for group in groups:
        total_size += len(group)
        # save it in the garbage collection dir
        name = "group" + str(i) + ".json"
        with open(f"garbage_collection2/combinations/{name}", "w") as f:
            json.dump(group, f, indent=4)
        i += 1
    if total_size != len(combos):
        print("error in create_files_from_combos")
    return groups

def main():
    incomplete = get_incomplete_files()
    incomplete_combos = create_incomplete_combos(incomplete)
    # create_files_from_combos(len(incomplete_combos)//8, incomplete_combos)
    create_files_from_combos(len(incomplete_combos), incomplete_combos)
    print("hi")

if __name__ == "__main__":
    main()

