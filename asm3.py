from __future__ import division
import sys
import os

from gene_stock import GeneStock
from hmm import HMM
from parser import Parser

local = True
if local:
    root_path = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3/HMM_Gene_Finding/asm3_source"

else:
    root_path = "/network/home/hejingyi/Documents/561/asm3_source"

vibrio_cholerae_gff = os.path.join(root_path, "Vibrio_cholerae.GFC_11.37.gff3")
vibrio_cholerae_seq = os.path.join(root_path,
                                   "Vibrio_cholerae.GFC_11.dna.toplevel.fa")

vibrio_vulnificus_gff = os.path.join(root_path,
                                     "Vibrio_vulnificus.ASM74310v1.37.gff3")
vibrio_vulnificus_seq = os.path.join(root_path,
                                     "Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa")
# TEST FILES
test_seq = os.path.join(root_path, "fake_level_{}.fa".format(2))
test_anno = os.path.join(root_path, "fake_level_{}.annot".format(2))

test_seq_0 = os.path.join(root_path, "fake_level_{}.fa".format(0))
test_anno_0 = os.path.join(root_path, "fake_level_{}.annot".format(0))

ANNOTATION_FILE = vibrio_cholerae_gff
TRAINING_SEQ_FILE = vibrio_cholerae_seq
EVALUATE_FILE = vibrio_vulnificus_seq

TEST_ANNOTATION_FILE = test_anno
TEST_TRAINING_SEQ_FILE = test_seq
TEST_EVALUATE_FILE = test_seq
# OUTPUT_FILE = os.path.join(os.getcwd(), "test.gff")
# OUTPUT_FILE = os.path.join(os.getcwd(), "Vibrio_vulnificus_prediction.gff")


################
#    Q1.1      #
################

NUCLEOTIDES = ['A', 'T', 'C', 'G']
STATES = ['I', 'S', 'M', 'P']


def q1(parser, verbose=False):
    parser.cal_len()
    # parser.sequence_loading(TRAINING_SEQ_FILE)
    parser.validation()
    for gs in parser.gene_stock.values():
        gs.state_sequence_from_loc()
    parser.cal_bp_freq()
    parser.cal_codon_freq()

    if verbose:
        print("annotation file name: {}".format(os.path.split(ANNOTATION_FILE)[-1]))
        print("----------------- Q1 -----------------")

        print("(i)  average length intergenic sequence is: \n{}".format(
            parser.intergenic_len))
        print("(ii)   average length gene sequence is: \n{}".format(
            parser.gene_len))
        print("(iii)   intergenic nucleotide frequency:")
        to_print = ""
        for k, i in parser.intergenic_nucleotide_freq.items():
            to_print += "\t{} --- {:2.4f}\n".format(k, i)
        print(to_print)

        print("(iv)   gene codons frequency:")
        to_print = ""
        for k, i in parser.gene_codon_freq.items():
            if i > 0:
                to_print += "\t{} --- {:2.4f}\n".format(k, i)
        print(to_print)

        print("\t start codons frequency:")

        to_print = ""
        for k, i in parser.start_codons_freq.items():
            if i > 0:
                to_print += "\t{} --- {:2.4f}\n".format(k, i)
        print(to_print)

        print("\tstop codons frequency:")
        to_print = ""
        for k, i in parser.stop_codons_freq.items():
            if i > 0:
                to_print += "\t{} --- {:2.4f}\n".format(k, i)

        print(to_print)

def q2(parser, evaluate_file, verbose=True):
    if verbose: print("----------------- Q2 -----------------")

    for gn, gs in parser.gene_stock.items():
        gs.build_intergenic_loc()
    q1(parser, verbose=False)

    hmm = HMM(parser, evaluate_file)
    hmm.build_transition_prob()
    hmm.build_emission_prob()
    if verbose:
        print("HMM transition probability: ")
        to_print = ""
        for k, i in hmm.transition_prob.items():
            if i > 0:
                to_print += "\t{} --- {:2.4f}\n".format(k, i)
        print(to_print)
        print("HMM emission probability: ")
        to_print = ""
        for k, i in hmm.emission_prob.items():
            if i > 0:
                to_print += "\t{} --- {:2.4f}\n".format(k, i)
        print(to_print)

    for i in hmm.input_gene_stock.keys():
        if verbose: print("predicting genes in sequence... {}".format(i))
        hmm.viterbi(i)
        # print(len(hmm.input_gene_stock[i].prediction[:998]))
        GeneStock.validate_state_seq(hmm.input_gene_stock[i].prediction)
    write_file = evaluate_file.split('.')[0] + "_predicted.gff"
    hmm.generate_gff_file(write_file)
    print("Generated prediction file in ...{}\n\n".format(
        os.path.split(write_file)[-1]))

    return hmm

def comparison(true_loc, pred_loc):
    all_true_start = [i[0] for i in true_loc]
    all_true_end = [i[1] for i in true_loc]

    perfect_match_cnt = 0
    perfect_stop_cnt = 0
    perfect_start_cnt = 0
    all_miss = 0

    for i, j in pred_loc:
        if i in all_true_start: perfect_start_cnt += 1
        if j in all_true_end: perfect_stop_cnt += 1
        if i in all_true_start and j in all_true_end: perfect_match_cnt += 1
        if i not in all_true_start and j not in all_true_end: all_miss += 1

    return perfect_match_cnt, perfect_stop_cnt, perfect_start_cnt, all_miss


def q4(true_gene_stock, predicted_gene_stock, sequence_name, my_true=False,verbose=True):
    perfect_match_cnt = 0
    perfect_stop_cnt = 0
    perfect_start_cnt = 0
    all_miss_cnt = 0
    gene_cnt = 0

    if verbose:
        print("----------------- Q4 -----------------")
        print("Statistics for gene found in .. {}\n".format(sequence_name))

    if my_true:
        print(" Annotated match my prediction: (Using my prediction as gold standard)")
    else:
        print(
            " My prediction match annotated: (Using annotated as gold standard)")

    for key in true_gene_stock.keys():
        true_gene_loc = true_gene_stock[key].gene_loc
        gene_cnt += len(true_gene_loc)
        predicted_gene_loc = predicted_gene_stock[key].gene_loc
        perfect_match_cnt_, perfect_stop_cnt_, perfect_start_cnt_, all_miss_cnt_ = comparison(
            true_gene_loc, predicted_gene_loc)
        perfect_match_cnt += perfect_match_cnt_
        perfect_start_cnt += perfect_start_cnt_
        perfect_stop_cnt += perfect_stop_cnt_
        all_miss_cnt += all_miss_cnt_

    print("perfect match: {:10.2f}%".format(perfect_match_cnt / gene_cnt  * 100))
    print("perfect start:  {:10.2f}%".format(perfect_start_cnt / gene_cnt * 100))
    print("perfect stop:  {:10.2f}%".format(perfect_stop_cnt / gene_cnt * 100))
    print("all miss:  {:10.2f}%".format(all_miss_cnt / gene_cnt * 100))

def q3(vv_seq, annotated_parser):
    print("----------------- Q3 -----------------")
    hmm = q2(annotated_parser, vv_seq, verbose=False)
    return hmm


def main():

    parser = Parser(ANNOTATION_FILE, TRAINING_SEQ_FILE)
    # parser.sequence_loading(test_seq_0)
    for gn, gs in parser.gene_stock.items():
        gs.build_intergenic_loc()
    q1(parser, verbose=True)


    test_parser = Parser(TEST_ANNOTATION_FILE, TEST_TRAINING_SEQ_FILE)
    q2(test_parser, TEST_EVALUATE_FILE, verbose=True)

    hmm = q3(EVALUATE_FILE, parser)

    # compute the true location from vv gff
    parser_vv = Parser(vibrio_vulnificus_gff, vibrio_vulnificus_seq)
    for gn, gs in parser_vv.gene_stock.items():
        gs.build_intergenic_loc()
    q1(parser_vv, verbose=True)
    true_gene_stock = parser_vv.gene_stock
    predicted_gene_stock = hmm.input_gene_stock
    evaluate_seq_name = os.path.split(EVALUATE_FILE)[-1]
    q4(true_gene_stock, predicted_gene_stock, evaluate_seq_name, my_true=False)
    q4(predicted_gene_stock, true_gene_stock, evaluate_seq_name, my_true=True,verbose=False)


if __name__ == '__main__':

    old_stdout = sys.stdout
    log_file = open("message.log", "w")
    sys.stdout = log_file
    main()

    log_file.close()
    sys.stdout = old_stdout


    # print(-math.inf * 0)
