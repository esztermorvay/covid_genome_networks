import os
import time
import zipfile
from Bio import SeqIO
import util
from Bio import pairwise2
import networkx as nx

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


def main():
    # testing
    util.extract_file_from_zip("genomes/SARS-CoV-2-BA.1.13.1.zip", "temp1")
    util.extract_file_from_zip("genomes/SARS-CoV-2-BA.1.1.15.zip", "temp2")

    se1 = get_longest_sequence_from_fasta("temp1/ncbi_dataset/data/genomic.fna")
    se2 = get_longest_sequence_from_fasta("temp2/ncbi_dataset/data/genomic.fna")
    G = nx.Graph()
    # this is a list of all the different genome zip files
    file_names = util.get_all_files_in_dir_as_list("genomes")[3:]
    file_names = util.take_random_sample(file_names, 50)
    # remove_empty_files(file_names, zips_dir)
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
