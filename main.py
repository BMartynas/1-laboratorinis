import numpy as np
from itertools import product
from Bio import SeqIO
import os
import glob

folder_path = './data'
all_codons = [''.join(i) for i in product(['T', 'C', 'A', 'G'], repeat=3)]
all_dicodons = [''.join(i) for i in product(all_codons, repeat=2)]
all_codons_dict = {codon: 0 for codon in all_codons}
all_dicodons_dict = {dicodon: 0 for dicodon in all_dicodons}


def read_fasta_file(file):
    for seq_record in SeqIO.parse(file, "fasta"):
        return seq_record


def find_codons(sequence):
    i = 0
    codons = []
    while i < len(sequence):
        if sequence[i] == 'ATG':
            start = i
            j = i
            while j < len(sequence):
                if sequence[j] == 'TAA' or sequence[j] == 'TAG' or sequence[j] == 'TGA':
                    end = j
                    codons.append(''.join(str(e)
                                          for e in sequence[start:end + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codons


def divide_into_triplets(sequence):
    frame1 = [sequence.seq[i:i + 3]
              for i in range(0, len(sequence), 3)]
    frame2 = [sequence.seq[i:i + 3]
              for i in range(1, len(sequence), 3)]
    frame3 = [sequence.seq[i:i + 3]
              for i in range(2, len(sequence), 3)]
    frame4 = [sequence.seq.reverse_complement()[i:i + 3]
              for i in range(0, len(sequence), 3)]
    frame5 = [sequence.seq.reverse_complement()[i:i + 3]
              for i in range(1, len(sequence), 3)]
    frame6 = [sequence.seq.reverse_complement()[i:i + 3]
              for i in range(2, len(sequence), 3)]
    return frame1 + frame2 + frame3 + frame4 + frame5 + frame6


def filter_out(codons):
    filteredList = list(filter(lambda x: len(x) >= 100, codons))
    return filteredList


def get_frequencies(fragments):
    for fragment in fragments:
        codon = [fragment[i:i + 3]
                 for i in range(0, len(fragment), 3)]
        dicodon = [fragment[i:i + 6]
                   for i in range(0, len(fragment), 3)]
        for key in all_codons_dict:
            all_codons_dict[key] += codon.count(key)
        for key in all_dicodons_dict:
            all_dicodons_dict[key] += dicodon.count(key)

    codons_sum = 0
    dicodons_sum = 0
    for key in all_codons_dict:
        codons_sum += all_codons_dict[key]
    for key in all_dicodons_dict:
        dicodons_sum += all_dicodons_dict[key]

    for key in all_codons_dict:
        all_codons_dict[key] = (all_codons_dict[key]/codons_sum) * 100
    for key in all_dicodons_dict:
        all_dicodons_dict[key] = (all_dicodons_dict[key]/dicodons_sum) * 100

    return all_codons_dict.values(), all_dicodons_dict.values()


def calculate_distance(freq1, freq2):
    return np.sum(np.absolute(freq1 - freq2))


codon_frequencies = []
dicodon_frequencies = []

for filename in glob.glob(os.path.join(folder_path, '*.fasta')):
    sequence = read_fasta_file(filename)
    triplets = divide_into_triplets(sequence)
    codons = find_codons(triplets)
    filtered_fragments = filter_out(codons)
    cod_freq, dicod_freq = get_frequencies(filtered_fragments)
    codon_frequencies.append(np.fromiter(cod_freq, dtype=float))
    dicodon_frequencies.append(np.fromiter(dicod_freq, dtype=float))

codons_distance_matrix = [[0 for x in range(8)] for y in range(8)]
dicodons_distance_matrix = [[0 for x in range(8)] for y in range(8)]
for z in range(8):
    for v in range(8):
        codons_distance_matrix[z][v] = calculate_distance(
            codon_frequencies[z], codon_frequencies[v])
        dicodons_distance_matrix[z][v] = calculate_distance(
            dicodon_frequencies[z], dicodon_frequencies[v])

print('\nKodonų atstumo matrica:')
print(codons_distance_matrix)
print('\nDikodonų atstumo matrica:')
print(dicodons_distance_matrix)
