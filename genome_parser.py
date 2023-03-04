import os
import sys
import time
import zipfile
from Bio import SeqIO
import util
from Bio import pairwise2
import networkx as nx
from multiprocessing import Pool
import threading
import json
import traceback

"""
https://biopython.org/docs/1.75/api/Bio.pairwise2.html
creates raw gml file without normalization
"""

zips_dir = "genomes"
fna_file = "/ncbi_dataset/data/genomic.fna"


def get_longest_sequence_from_fasta(fasta_file_path):
    fasta_sequences = SeqIO.parse(open(fasta_file_path), 'fasta')
    longest_length = 0
    longest_sequence = ""
    count = 0
    for fasta in fasta_sequences:
        count += 1
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) > longest_length:
            longest_length = len(sequence)
            longest_sequence = sequence
    return longest_sequence, count


def get_similarity_score(sequence1, sequence2):
    # get the similarity score between two sequences
    # https://biopython.org/docs/1.75/api/Bio.pairwise2.html
    score = pairwise2.align.globalxx(sequence1, sequence2, score_only=True)
    return score
    # print(alignments)
    # print(pairwise2.format_alignment(*alignments[0]))


def remove_empty_files(file_names, dir):
    for file_name in file_names:
        file_path = dir + "/" + file_name
        try:
            util.extract_file_from_zip(file_path, "temp")
        except:

            os.remove(file_path)
            print("removed " + file_name)

def get_similarity_scores_subprocess(combinations_info, debugging=False):
    # tuple of thread_num, combinations
    thread_num, group1name, group2name, combinations = combinations_info
    print("thread " + str(thread_num) + " started")
    log_file_name = group1name + "_" + group2name + "_thread" +  str(thread_num) + ".json"
    # inititialize a list of empty strings with size length
    counts = {}
    results = {}
    benchmark = 5
    counter = 0
    # logging = thread_num == 1
    logging = debugging
    if logging:
        print("logging for thread " + str(thread_num), len(combinations))
    for combination in combinations:
        counter += 1
        if counter % benchmark == 0 or counter == len(combinations):
            if not debugging:
                with open("garbage_collection/counts/" + log_file_name, "w") as counts_file:
                    json.dump(counts, counts_file, indent=4)
                with open("garbage_collection/scores/" + log_file_name, "w") as counts_file:
                    json.dump(results, counts_file, indent=4)
            # print current timestamp as str
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            print("\tcombination " + str(thread_num) + " " + str(counter) + "/" + str(len(combinations)), 100*counter/len(combinations), "% done")
        file_name1 = combination[0][11:-4]
        file_name2 = combination[1][11:-4]
        file_path1 = zips_dir + "/" + combination[0]
        file_path2 = zips_dir + "/" + combination[1]
        try:
            temp1 = "temp" + group1name + group2name +  str(thread_num)
            temp2 = "temp" + group1name + group2name + str(thread_num*100)
            util.extract_file_from_zip(file_path1, temp1)
            util.extract_file_from_zip(file_path2, temp2)
            # time how long this takes
            start_time = time.time()
            if logging:
                print("attempting to get longest seq")
            sequence1, counts[file_name1] = get_longest_sequence_from_fasta(temp1 + fna_file)

            sequence2, counts[file_name2] = get_longest_sequence_from_fasta(temp2 + fna_file)
            if logging:
                print("got longest seq")
            score = get_similarity_score(sequence1, sequence2)
            end_time = time.time()

            # print("time to get similarity score: " + str(end_time - start_time))
            # key = (file_name1, file_name2)
            key = file_name1 + "_" + file_name2
            # print(key, ":", score)
            results[key] = score
            if logging:
                print("iteration completed: " + str(counter) + "/" + str(len(combinations)))
            # FOR DEBUGGING ONLY
            # return
        except Exception as e:
            print("Exception occured in thread " + str(thread_num), combination)
            print(traceback.format_exc())
            continue
    print("done with thread " + str(thread_num))
    if not debugging:
        with open("garbage_collection/counts/" + log_file_name, "w") as counts_file:
            json.dump(counts, counts_file, indent=4)
        with open("garbage_collection/scores/" + log_file_name, "w") as counts_file:
            json.dump(results, counts_file, indent=4)
    print("done writing to file for thread " + str(thread_num))
    return 0


def run_multithreading(group1name, group2name, num_threads=4):
    to_process = get_subgroups(group1name, group2name, num_threads=num_threads)
    with Pool(num_threads) as p:
        p.map(get_similarity_scores_subprocess, to_process)
    print("threads joined")
    # with Pool(2) as p:
    #     p.map(get_similarity_scores_subprocess, to_process)


def get_subgroups(group1name, group2name, num_threads=4, directory="garbage_collection/combinations/"):
    """ return thread name"""
    # open the file corresponding to the combinations we need to do
    file_name = directory + group1name
    if group2name != group1name:
        file_name += group2name
    file_name += ".json"
    # open the combinations file
    combinations = []
    with open(file_name, "r") as f:
        combinations = json.load(f)
    # divide the total amount of combinations into num_threads
    interval_size = len(combinations)//num_threads
    to_process = []
    # last_interval = 0
    for i in range(num_threads):
        to_process.append((i+1, group1name, group2name, combinations[i*interval_size:(i+1)*interval_size]))
        # last_interval = (i+1)*interval_size
    # if last_interval < len(combinations) - 1:
    #     to_process.append((num_threads+1, group1name, group2name, combinations[last_interval:]))
    # print(len(to_process))
    # interval_size = len(combinations)//4
    # to_process2 = [(1,group1name, group2name, combinations[0:interval_size]), (2,group1name, group2name, combinations[interval_size:2*interval_size]),
    #               (3,group1name, group2name, combinations[2*interval_size:3*interval_size]), (4,group1name, group2name, combinations[3*interval_size:])]
    return to_process
def main():
    """
    USAGE
    python3 genome_parser.py [num_threads] [optional: combination to do]


    :return:
    """
    # test()
    # get commad line arg
    numthreads = int(sys.argv[1])
    combination = None
    if len(sys.argv) > 2:
        combination = sys.argv[2]
    # testing
    # util.extract_file_from_zip("genomes/SARS-CoV-2-BA.1.13.1.zip", "temp1")
    # util.extract_file_from_zip("genomes/SARS-CoV-2-BA.1.1.15.zip", "temp2")
    #
    # se1 = get_longest_sequence_from_fasta("temp1/ncbi_dataset/data/genomic.fna")
    # se2 = get_longest_sequence_from_fasta("temp2/ncbi_dataset/data/genomic.fna")
    # load the combinations file
    with open("garbage_collection/combinations.json", "r") as f:
        combinations_done = json.load(f)

    group1name = ""
    group2name = ""
    current = ""
    if combination is not None:
        current = combination
        groups = combination.split("_")
        group1name = groups[0]
        group2name = groups[1]

    else:
        # get the next combination which is false
        for combination in combinations_done:
            if not combinations_done[combination]:
                current = combination
                groups = combination.split("_")
                group1name = groups[0]
                group2name = groups[1]
                break
    if group1name == "":
        print("all combinations done")
        return
    print("starting with " + current)

    combinations_done[current] = True
    with open("garbage_collection/combinations.json", "w") as f:
        json.dump(combinations_done, f, indent=4)
    run_multithreading(group1name, group2name, numthreads)

    print("Done with multithreading for " + current + "")
    # update the combinations file
    combinations_done[current] = True
    with open("garbage_collection/combinations.json", "w") as f:
        json.dump(combinations_done, f, indent=4)


def test():
    to_test = [['SARS-CoV-2-AY.117.zip', 'SARS-CoV-2-AY.4.7.zip']]
    com_info = (1, "SARS-CoV-2", "SARS-CoV-2", to_test)
    get_similarity_scores_subprocess(com_info, debugging=True)


def orig_function():
    G = nx.Graph()

    file_names = util.get_all_files_in_dir_as_list("genomes")
    file_names = file_names[2:]
    length = len(file_names)
    # inititialize a list of empty strings with size length
    sequences = [""] * length
    counts = [0] * length
    sum_ = 0
    # from file names check how many are  more than 5 digits
    for file in file_names:
        # find the first occurance of a period
        index = file.find(".")
        # get the substring from the period to the end
        substr = file[index + 1:]
        # get the count of numerical characters in substr
        # count = sum(c.isdigit() for c in substr)
        count = len(substr)
        # if the count is more than 5 digits
        if count > 10:
            sum_ += 1
    print(sum_)

    temp1 = "temp1"
    temp2 = "temp2"
    # nested for loop to iterate through all pairs of file names
    # each calculation of similarity score takes ~13 seconds
    for i in range(length):
        print()
        print(f"{i / length * 100}% done")
        print()
        # extract the first file
        file_name1 = file_names[i]
        file_path1 = zips_dir + "/" + file_name1
        try:
            util.extract_file_from_zip(file_path1, temp1)
            if sequences[i] == "":
                sequences[i], counts[i] = get_longest_sequence_from_fasta(temp1 + fna_file)
            sequence1 = sequences[i]
            # sequence1 = get_longest_sequence_from_fasta(temp1 + fna_file)
            print("iteration i {i} of {length}, {file_name1}".format(i=i, length=length, file_name1=file_name1))
            # set the count of the current node
            if not G.has_node(file_name1[11:-4]):
                G.add_node(file_name1[11:-4], count=counts[i])
            for j in range(i + 1, length):
                file_name2 = file_names[j]
                file_path2 = zips_dir + "/" + file_name2
                try:
                    util.extract_file_from_zip(file_path2, temp2)
                    # get the longest sequence from each file
                    if sequences[j] == "":
                        sequences[j], counts[j] = get_longest_sequence_from_fasta(temp2 + fna_file)
                    sequence2 = sequences[j]
                    # sequence2 = get_longest_sequence_from_fasta(temp2 + fna_file)
                    print(
                        "\titeration j {j} of {length}, {file_name2}".format(j=j, length=length, file_name2=file_name2))
                    # get the similarity score between the two sequences, and the time taken to get it
                    # get the current time
                    start_time = time.time()
                    score = get_similarity_score(sequence1, sequence2)
                    end_time = time.time()
                    time_taken = end_time - start_time

                    print("\ttime taken", int(time_taken), "seconds")
                    if not G.has_node(file_name2[11:-4]):
                        G.add_node(file_name2[11:-4], count=counts[j])
                    G.add_edge(file_name1[11:-4], file_name2[11:-4], weight=score)
                    print("\tscore between " + file_name1[:-4] + " and " + file_name2[:-4] + " is " + str(score))
                except Exception as e:
                    print("\terror with " + file_name2)
                    continue
        except Exception as e:
            print("error with " + file_name1)
            continue

    nx.write_gml(G, "graph_raw.gml")
    # t1 = get_longest_sequence_from_fasta("temp/ncbi_dataset/data/genomic.fna")
    # t2 = get_longest_sequence_from_fasta("temp2/ncbi_dataset/data/genomic.fna")
    # get_similarity_score(t1, t2)

    # print(file_names)



if __name__ == "__main__":
    main()
