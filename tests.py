import unittest
import kmers_medley as medley


class KmersTestCase(unittest.TestCase):

    def test_Generator(self):
        expected_list = ['aa', 'ab', 'ba', 'bb']
        observed_list = list(medley.create_kmers_generator('ab', 2))
        self.assertEqual(expected_list, observed_list, msg='generator return correct values')

    def test_Select_kmers(self):
        expected_set = {'bb'}
        sequences = ['aa', 'ab', 'ba']
        generator = medley.create_kmers_generator('ab', 2)
        observed_set = medley.select_most_different_kmers(sequences, generator, 2)
        self.assertEqual(expected_set, observed_set,  msg='the most different kmer was selected')

    def test_Reverse_complement(self):
        expected_dna = 'CGTA'
        observed_dna = medley.reverse_complement('TaCg')
        self.assertEqual(expected_dna, observed_dna, msg='reverse complement DNA is correct')

    def test_Duplex_kmers(self):
        expected_set = {'CCC', 'CCG', 'CGC', 'CGG', 'GCC', 'GCG', 'GGC', 'GGG'}

        sequences = ['AAA']

        generator = medley.create_kmers_generator('ATCG', 3)
        observed_set = medley.select_most_different_kmers(sequences, generator, 3, True)
        self.assertEqual(expected_set, observed_set, msg='the most different duplex kmer was selected')

        expected_set2 = {'AAA', 'TTT'}
        sequences2 = list(medley.create_kmers_generator('ATCG', 3))
        sequences2.remove('AAA')
        sequences2.remove('TTT')
        sequences2.remove('CCC')

        generator2 = medley.create_kmers_generator('ATCG', 3)
        observed_set2 = medley.select_most_different_kmers(sequences2, generator2, 3, True)
        self.assertEqual(expected_set2, observed_set2, msg='the most different duplex kmer was selected')

    def test_Compare_sets(self):
        expected_set = {'AA'}
        set1 = {'AA', 'TT'}
        set2 = {'AA', 'CC'}
        set3 = {'AA', 'GG'}
        observed_set = medley.compare_kmer_sets([set1, set2, set3])
        self.assertEqual(expected_set, observed_set, msg='union is correct')


if __name__ == '__main__':
    unittest.main()
