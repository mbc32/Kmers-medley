import itertools
import configparser
import argparse

dna_translation_table = str.maketrans("ACTG", "TGAC")


def reverse_complement(dna):
    """
    create the reverse complement DNA string in uppercase
    :param dna: str
    :return: str
    """
    return dna.upper().translate(dna_translation_table)[::-1]


def create_kmers_generator(alphabet, repeat_length):
    """
    Creates the Cartesian product i.e alphabet 'ab' and repeat_length 2 would give 'aa','ab','ba','bb'.
    :param alphabet: str
    :param repeat_length: int
    :return generator
    """
    kmers_generator = (''.join(i) for i in itertools.product(alphabet, repeat=repeat_length))
    return kmers_generator


def select_most_different_kmers(search_sequences,  kmers_generator, kmer_length, duplex=False):
    """
    select a set of kmers that are most different to the search sequences, based on count of match/mis-match positions
    :param search_sequences: list
    :param kmers_generator: generator object
    :param kmer_length: int
    :param duplex: boolean
    :return: set
    """
    global_min_count = kmer_length
    most_different_set = set()
    for kmer in kmers_generator:
        max_count = 0
        if duplex:
            kmer_revcomp = reverse_complement(kmer)
        for sequence in search_sequences:
            if len(kmer) == len(sequence):
                count = sum(1 for base1, base2 in zip(kmer, sequence) if base1 == base2)
                if count > max_count:
                    max_count = count
                if duplex:
                    revcomp_count = sum(1 for base1, base2 in zip(kmer_revcomp, sequence) if base1 == base2)
                    if revcomp_count > max_count:
                        max_count = revcomp_count
                if max_count > global_min_count:
                    # Break early if the kmer is more similar to the sequence than the present most different kmers.
                    break
            else:
                raise ValueError("The Kmer {} and Sequence {} are not the same length".format(kmer, sequence))
        if max_count < global_min_count:
            most_different_set = {kmer}
            global_min_count = max_count
        elif max_count == global_min_count:
            most_different_set.add(kmer)
    return most_different_set


def compare_kmer_sets(sets):
    """
    find the intersection of sets
    :param sets: list
    :return: set
    """
    return sets.pop().intersection(*sets)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', nargs=1)
    args = parser.parse_args()

    config_file = args.config
    config = configparser.ConfigParser()
    config.read(config_file)

    dna_bases = config['SETUP']['alphabet']
    kmer_repeat = int(config['SETUP']['kmer_length'])
    dna_duplex = config['SETUP']['duplex']
    output_file = config['SETUP']['output_file']
    file_paths = config['INPUT_FILES']['file_path']
    result_kmer_sets = list()

    file_list = file_paths.split(',')

    for file_path in file_list:
        with open(file_path, 'r') as file:
            sequences = file.read().rstrip().splitlines()
        kmers = create_kmers_generator(dna_bases, kmer_repeat)
        most_different_kmer_set = select_most_different_kmers(sequences, kmers, kmer_repeat, dna_duplex)
        result_kmer_sets.append(most_different_kmer_set)

    intersection_set = compare_kmer_sets(result_kmer_sets)

    with open(output_file, 'w') as file:
        for element in intersection_set:
            file.write(element + '\n')
