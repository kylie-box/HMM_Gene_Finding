from __future__ import division

import math
import os

from parser import Parser
from hmm import HMM
from gene_stock import GeneStock

local = True
if local:
    root_path = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3/HMM_Gene_Finding/asm3_source"

else:
    root_path = "/network/home/hejingyi/Documents/561/asm3_source"

gff_file = os.path.join(root_path, "Vibrio_cholerae.GFC_11.37.gff3")
seq_file = os.path.join(root_path, "Vibrio_cholerae.GFC_11.dna.toplevel.fa")

# TEST FILES
# test_anno_0 = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_0.annot"
# test_seq_0 = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_0.fa"
# test_anno_1 = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_1.annot"
# test_seq_1 = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_1.fa"
# test_anno_2 ="/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_2.annot"
# test_seq_2 = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/fake_level_2.fa"

test_seq = os.path.join(root_path, "fake_level_{}.fa".format(2))
test_anno = os.path.join(root_path, "fake_level_{}.annot".format(2))
# examiner = GFFExaminer()
# in_handle = open(in_file)
# pprint.pprint(examiner.available_limits(in_handle))
# in_handle.close()

# in_handle = open(in_file)
# for rec in GFF.parse(in_handle):
#     print(rec)
# in_handle.close()

################
#    Q1.1      #
################

NUCLEOTIDES = ['A', 'T', 'C', 'G']
STATES = ['I', 'S', 'M', 'P']




def q1(parser, verbose=False):
    parser.cal_len()
    parser.sequence_loading(test_seq)
    parser.validation()
    for gs in parser.gene_stock.values():
        gs.state_sequence_from_loc()
    parser.cal_bp_freq()
    parser.cal_codon_freq()

    if verbose:
        print("average length gene sequence is: {}".format(
            parser.gene_len))
        print("average length intergenic sequence is: {}".format(
            parser.intergenic_len))

        print("intergenic nucleotide frequency:")
        print(parser.intergenic_nucleotide_freq)

        print("gene codons frequency:")
        print(parser.gene_codon_freq)
        print("start codons frequency: ")
        print(parser.start_codons_freq)
        print("stop codons frequency: ")
        print(parser.stop_codons_freq)


def q2(parser):
    hmm = HMM(parser, test_seq)
    hmm.build_transition_prob()
    hmm.build_emission_prob()
    print(hmm.transition_prob)
    print(hmm.emission_prob)
    for i in hmm.input_gene_stock.keys():
        # print(i)
        hmm.viterbi(i)
        # print(len(hmm.input_gene_stock[i].prediction[:998]))
        GeneStock.validate_state_seq(hmm.input_gene_stock[i].prediction)
        path = os.path.join(os.getcwd(), "test_gff.txt")
        hmm.gen_gff_file(path)


def main():
    parser = Parser(test_anno)
    # parser.sequence_loading(test_seq_0)
    for gn, gs in parser.gene_stock.items():
        gs.build_intergenic_loc()
    q1(parser, verbose=False)
    q2(parser)


if __name__ == '__main__':
    main()
    # print(-math.inf * 0)
