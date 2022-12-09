import json
import os

util_dir = "util_files"
web_scraping_dir = "genomes"
lineages_file = os.path.join(os.getcwd(), util_dir, "lineage_data.full.json")




def run_cmd_script(directory, cmd_to_run):

    # Run the script in the directory
    # print the current director
    # print(os.getcwd())
    os.system(cmd_to_run)


def read_json_as_list(json_file_path):
    # read the json file into a dictionary
    with open(json_file_path, 'r') as f:
        data = json.load(f)

    # get all the keys in the dictionary as a list
    keys = list(data.keys())
    return keys



def install_commands():
    os.chdir(web_scraping_dir)
    command1 = 'curl -o datasets.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe"'
    command2 = 'curl -o dataformat.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/dataformat.exe"'
    # run the commands
    os.system(command1)
    os.system(command2)


def main():
    # read in the list of pango lineages we want to scrape
    # install_commands()
    lineages = read_json_as_list(lineages_file)
    os.chdir(web_scraping_dir)

    length = len(lineages)
    iteration = 0
    #create a file and use it to keep track of progress
    with open("progress.txt", "w") as f:
        for lineage in lineages:
            f.write("iteration {iteration} of {length}, {lineage} \n".format(iteration=iteration, length=length, lineage=lineage))
            f.write(f"{iteration/length*100}% done \n")
            print("iteration {iteration} of {length}, {lineage}".format(iteration=iteration, length=length, lineage=lineage))
            command = f"datasets download virus genome taxon SARS2 --complete-only --lineage {lineage} --filename SARS-CoV-2-{lineage}.zip"
            run_cmd_script(web_scraping_dir, command)
            iteration+=1


if __name__ == "__main__":
    main()