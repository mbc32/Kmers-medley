import itertools


def create_kmers_generator(alphabet, kmer_length):
    """
    Creates the Cartesian product i.e alphabet 'ab' and kmer_length 2 would give 'aa','ab','ba','bb'.
    :param alphabet: str
    :param kmer_length: int
    :return generator
    """
    kmers_generator = (''.join(i) for i in itertools.product(alphabet, repeat=kmer_length))
    return kmers_generator


def select_most_different_kmers(search_sequences,  kmers_generator, kmer_length):
    """
    select a set of kmers that are most different to the search sequences, based on count of match/mis-match positions
    :param search_sequences: list
    :param kmers_generator: generator object
    :param kmer_length: int
    :return: set
    """
    global_min_count = kmer_length
    base_match_counter = {'0': set(), '1': set(), '2': set(), '3': set(), '4': set(), '5': set(), '6': set()}
    for kmer in kmers_generator:
        max_count = 0
        for sequence in search_sequences:
            if len(kmer) == len(sequence):
                count = sum(1 for base1, base2 in zip(kmer, sequence) if base1 == base2)
                if count > max_count:
                    max_count = count
                if count > global_min_count:
                    break
            else:
                print('Error')
        if max_count < global_min_count:
            base_match_counter[str(max_count)].add(kmer)
            global_min_count = max_count


if __name__ == '__main__':
    pass
