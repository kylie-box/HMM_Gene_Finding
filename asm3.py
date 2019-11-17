from __future__ import division

from collections import Counter
import math
import os

local = True
if local:
    root_path = "/Users/kyliehe/Documents/words/study/class/comp/561/asm3_source/"

else:
    root_path = "/network/home/hejingyi/Documents/561/asm3_source"

gff_file = os.path.join(root_path, "Vibrio_cholerae.GFC_11.37.gff3")
seq_file = os.path.join(root_path,"Vibrio_cholerae.GFC_11.dna.toplevel.fa")

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


class Parser:
    def __init__(self, file_name):
        self.file_name = file_name
        self.gene_stock = {}
        self.file_parser()
        self.intergenic_nucleotide_freq = None
        self.gene_codon_freq = None
        self.start_codons_cnt = None
        self.stop_codons_cnt = None
        self.intergenic_len = None
        self.gene_len = None
        self.gene_nucleotide_freq = None

    def validation(self):
        for gs in self.gene_stock.values():
            if gs.seq_len == 0 and gs.seq is not None:
                gs.seq_len = len(gs.seq)
            if len(gs.gene_loc)>0 and gs.gene_loc[-1][1] != gs.seq_len:
                gs.intergenic_loc.append((gs.gene_loc[-1][1] + 1, gs.seq_len))
                # print(gs.gene_len)


    def file_parser(self):
        print(self.file_name)
        with open(self.file_name, 'r') as w:
            cnt = 0
            for l in w:
                # init buckets
                if l.startswith("##sequence-region  "):
                    gene_name = \
                    l.strip("##sequence-region  ").strip('\n').split(' ')[0]
                    gs = GeneStock()
                    gs.gene_name = gene_name
                    self.gene_stock[gene_name] = gs
                    gs.seq_len = int(
                        l.strip("##sequence-region  ").strip('\n').split(' ')[
                            2])
                    continue

                # genetic region
                elif l.startswith("#"):
                    continue
                else:
                # if l.startswith("DN38"):

                    gene_lst = l.strip('\n').split('\t')
                    gene_name = gene_lst[0]
                    gene_starts = int(gene_lst[3])
                    gene_ends = int(gene_lst[4])
                    try:
                        crt_gene_stock = self.gene_stock[gene_name]
                    except KeyError:
                        gs = GeneStock()
                        gs.gene_name = gene_name
                        self.gene_stock[gene_name] = gs
                        crt_gene_stock = gs

                    # process gene regions
                    if gene_lst[6] == "+" \
                            and gene_lst[2] == "CDS":
                        loc = (gene_starts, gene_ends)
                        # overlapping gene
                        if len(crt_gene_stock.gene_loc) > 0 \
                                and gene_starts < crt_gene_stock.gene_loc[-1][1]:
                            # print("\n\n",gene_starts, crt_gene_stock.loc[-1][1])
                            cnt += 1
                        crt_gene_stock.gene_loc.append(loc)
                        crt_gene_stock.gene_count += 1

                    # elif gene_lst[6] == "." and crt_gene_stock.loc == []:
                    #     span = (gene_starts, gene_ends)
                    #     crt_gene_stock.span = span
                    # else:
                    #     continue

        # print("number of overlapping gene regions: ", cnt)

    def sequence_loading(self, fasta_file):
        with open(fasta_file, 'r') as f:
            gene_name = None
            gene_seq = ""
            for l in f:
                if l.startswith(">"):
                    if gene_name is not None:
                        self.gene_stock[gene_name].seq = gene_seq
                        if self.gene_stock[gene_name].seq_len == 0 \
                                and self.gene_stock[gene_name].seq is not None:
                            self.gene_stock[gene_name].seq_len = len(gene_seq)
                    gene_name = l.strip('\n').strip('>').split(" ")[0]
                    gene_seq = ""
                    continue
                else:
                    gene_seq += l.strip('\n')
            # last seqeuence in gene_seq needs to get into the bucket!

            self.gene_stock[gene_name].seq = gene_seq
            if self.gene_stock[gene_name].seq_len == 0 \
                    and self.gene_stock[gene_name].seq is not None:
                self.gene_stock[gene_name].seq_len = len(gene_seq)

    def cal_len(self):
        gene_len = 0
        intergenic_len = 0
        gene_cnts = 0
        intergenic_cnts = 0

        for gn, gs in self.gene_stock.items():
            gs.build_intergenic_loc()
            gs.gene_len = gs.cal_len_from_loc(gs.gene_loc)
            gene_len += gs.gene_len
            gene_cnts += gs.gene_count
            # print(gs.gene_name, gs.loc)
            gs.intergenic_len = gs.cal_len_from_loc(gs.intergenic_loc)
            intergenic_len += gs.intergenic_len
            intergenic_cnts += gs.intergenic_count

        self.intergenic_len = intergenic_len // intergenic_cnts
        self.gene_len = gene_len // gene_cnts

    def cal_bp_freq(self):
        inter_cnt = Counter()
        gene_cnt = Counter()
        # res = {'A':0, 'T':0,'C':0,'G':0}
        for gn, gs in self.gene_stock.items():
            gs.cal_bp_freq(inter_cnt, gene=False)
            gs.cal_bp_freq(gene_cnt, gene=True)
        tot = sum(inter_cnt.values())
        self.intergenic_nucleotide_freq = {n: cnt / tot for n, cnt in
                                           inter_cnt.items()}

        if len(self.intergenic_nucleotide_freq) < 4:
            self.intergenic_nucleotide_freq = self.suppliment_nucleotide_usage(
                        self.intergenic_nucleotide_freq)

        tot = sum(gene_cnt.values())
        self.gene_nucleotide_freq = {n: cnt / tot for n, cnt in
                                     gene_cnt.items()}

        if len(self.gene_nucleotide_freq) < 4:
            self.gene_nucleotide_freq = self.suppliment_nucleotide_usage(
                        self.gene_nucleotide_freq)

    def cal_codon_freq(self):
        res = Counter()
        start_codons = Counter()
        stop_codons = Counter()
        for gn, gs in self.gene_stock.items():
            res, start_codons, stop_codons = gs.cal_gene_codon_freq(res,
                                                                    start_codons,
                                                                    stop_codons)

        tot = sum(res.values())
        self.gene_codon_freq = {n: cnt / tot for n, cnt in res.items()}
        if len(self.gene_codon_freq) < 64:
            self.gene_codon_freq = self.suppliment_codon_usage(
                self.gene_codon_freq)

        # tot = sum(start_codons.values())
        # self.start_codons_cnt = {n: cnt / tot for n, cnt in
        #                           start_codons.items()}
        self.start_codons_cnt = start_codons

        # if len(self.start_codons_freq) < 64:
        #     self.start_codons_freq = self.suppliment_codon_usage(
        #         self.start_codons_freq)

        # tot = sum(stop_codons.values())
        # self.stop_codons_cnt = {n: cnt / tot for n, cnt in stop_codons.items()}
        self.stop_codons_cnt = stop_codons
        # if len(self.stop_codons_freq) < 64:
        #     self.stop_codons_freq = self.suppliment_codon_usage(
        #         self.stop_codons_freq)

    @staticmethod
    def suppliment_codon_usage(freq_table):
        for i in NUCLEOTIDES:
            for j in NUCLEOTIDES:
                for k in NUCLEOTIDES:
                    try:
                        _ = freq_table[str(i + j + k)]
                    except KeyError:
                        freq_table[str(i + j + k)] = 0
                        continue
        return freq_table

    @staticmethod
    def suppliment_nucleotide_usage(freq_table):
        for i in NUCLEOTIDES:
            if freq_table.get(i) is None:
                freq_table[i] = 0

        return freq_table



class GeneStock:
    def __init__(self):
        self.gene_loc = []
        self.intergenic_loc = []
        self.gene_count = 0
        self.intergenic_count = 0
        self.gene_name = None
        self.intergenic_len = 0
        self.gene_len = 0
        self.seq_len = 0
        self.seq = None
        self.prediction = []

    def build_intergenic_loc(self):
        next_start = 1
        intergenic_len = 0
        if len(self.gene_loc) == 0:
            self.intergenic_loc = [(1, self.seq_len)]
        else:
            for start, end in self.gene_loc:
                if start > next_start:
                    anchor = next_start
                    next_start = end + 1
                assert anchor < start

                # print(anchor, start, end)
                self.intergenic_loc.append((anchor, start - 1))
                intergenic_len += start - anchor
                self.intergenic_count += 1

            # interval of last end, seq len
            if self.seq_len != 0:
                self.intergenic_loc.append((end + 1, self.seq_len))
                self.intergenic_count += 1

    def cal_len_from_loc(self, loc):
        length = 0
        for start, end in loc:
            length += end - start + 1
        return length


    def cal_bp_freq(self, counter, gene=True):
        if gene:
            loc = self.gene_loc
        else:
            loc = self.intergenic_loc
        for start, end in loc:
            # inclusive end
            region = self.seq[slice(start - 1, end)]
            counter.update(region)
        return counter

    def cal_gene_codon_freq(self, codon_counter,
                            start_codon_counter,
                            stop_codon_coutner):
        def count_codons(seq):
            codon_counter.update([seq[i: i + 3] for i in range(0, len(seq), 3)])

        for start, end in self.gene_loc:
            gene_region = self.seq[slice(start - 1, end)]
            start_codon_counter.update([self.seq[start - 1:start + 2]])
            stop_codon_coutner.update([self.seq[end - 3: end]])
            count_codons(gene_region)

        return codon_counter, start_codon_counter, stop_codon_coutner

    def state_sequence_from_loc(self):
        # do it in reverse direction
        is_gene = False
        state_seq = []
        gene_loc = self.gene_loc
        intergenic_loc = self.intergenic_loc
        anchor = 0
        i = 0
        for start, end in gene_loc:

            while anchor < start:
                state_seq.append('I')
                anchor += 1
                i += 1
            while anchor <= end:
                state_seq.append('M')
                anchor += 1
                i += 1
            # anchor = end
        while anchor < self.seq_len:
            state_seq.append('I')
            anchor += 1

        # print(gene_loc[-1])
        # print(intergenic_loc[-1])
        # if self.seq_len  == gene_loc[-1][1]:
        #     is_gene = True
        #
        # gene_loc_idx = -1
        # intergenic_loc_idx = -1
        # while len(state_seq) < self.seq_len:
        #     if is_gene:
        #         for i in range(gene_loc[gene_loc_idx][0]-1, gene_loc[gene_loc_idx][1]):
        #             state_seq.append('M')
        #         gene_loc_idx -= 1
        #         is_gene = False
        #     else:
        #         for j in range(intergenic_loc[intergenic_loc_idx][0]-1, intergenic_loc[intergenic_loc_idx][1]):
        #             state_seq.append('I')
        #         intergenic_loc_idx -= 1
        #         is_gene = True
        # state_seq.reverse()

        assert len(state_seq) == self.seq_len
        # self.validate_state_seq(state_seq)
        self.prediction = state_seq

    @staticmethod
    def validate_state_seq(state_seq):
        gene_loc = []
        seq = ""
        for i in range(len(state_seq) -1):
            seq += state_seq[i]
            if state_seq[i] == 'I' and state_seq[i + 1] == 'S':
                gene_loc.append(i + 1)
            elif state_seq[i] == 'P' and state_seq[i + 1] == 'I':
                gene_loc.append(i + 1)

        print(gene_loc)

    def prediction_sequence_to_loc(self):
        gene_turing_pt = []
        gene_loc = []
        gene_start_idx = 0 if self.prediction[0] == 'S' else None
        gene_end_idx = self.seq_len if self.prediction[-1] == 'P' else None

        if gene_start_idx is not None:
            gene_turing_pt.append(gene_start_idx)

        for i in range(len(self.prediction) -1):
            if self.prediction[i] == 'I' and self.prediction[i+1] == 'S':
                gene_turing_pt.append(i+1)
            elif self.prediction[i] == 'P' and self.prediction[i+1] == 'I':
                gene_turing_pt.append(i+1)
        if gene_end_idx is not None:
            gene_turing_pt.append(gene_end_idx)

        for j in range(0, len(gene_turing_pt)-1, 2):
            gene_loc.append((gene_turing_pt[j], gene_turing_pt[j+1]))

        print(gene_loc)

        print(self.prediction[999:1021])



class Cell:
    def __init__(self):
        """
                o0   o1   o2   o3
        --------------------------
        I   |
        S   |
        M   |
        P   |
        """
        self.coor = None  # (State, Emission)
        self.value = None   # store the log value of the prob
        self.pointer = None  # previous cell

    def __str__(self):
        output_msg = "coor: {}\nvalue: {}\npointer: {}\n".format(self.coor,
                                                                self.value,
                                                                self.pointer)
        return output_msg

    def __abs__(self):
        return self.value


class HMM:
    def __init__(self, parser, input_fasta):
        self.states = ({s:i for s, i in zip(STATES, range(4))})
        self.observations = NUCLEOTIDES
        # I: intergenic; S: Start; M: Middle; P: Stop
        self.observations = set()
        self.emission_prob = {}     # s:obs given state s observe obs
        self.transition_prob = {}       # s0:s1 given s0 transition to s1
        self.initial_prob = {'I': 1, 'S': 0, 'M': 0, 'P': 0}
        self.parser = parser
        self.input_gene_stock = {}
        self.sequence_loading(fasta_file=input_fasta)



    def sequence_loading(self, fasta_file):
        with open(fasta_file, 'r') as f:
            gene_name = None
            gene_seq = ""
            i=0
            for l in f:
                if l.startswith(">"):
                    if gene_name is not None:
                        self.input_gene_stock[gene_name] = GeneStock()
                        self.input_gene_stock[gene_name].seq = gene_seq
                        i += 1
                    gene_name = l.strip('\n').strip('>').split(" ")[0]
                    gene_seq = ""
                    continue
                else:
                    gene_seq += l.strip('\n')
            # last seqeuence in gene_seq needs to get into the bucket!
            self.input_gene_stock[gene_name] = GeneStock()
            self.input_gene_stock[gene_name].seq = gene_seq

    def build_transition_prob(self):
        # for s0 in STATES:
        #     for s1 in STATES:
        #         prob = 0
        #         if s0 == 'I':
        #             if s1 == 'I':
        #                 prob = 1 - 1 / self.parser.intergenic_len
        #             elif s1 == 'S':
        #                 prob = 1 / self.parser.intergenic_len
        #             else:
        #                 prob = 0
        #         elif s0 == 'S':
        #             if s1 == 'S':
        #                 prob = 2/3
        #             elif s1 == 'M':
        #                 prob = 1/3
        #             else:
        #                 prob = 0
        #         elif s0 == 'M':
        #             if s1 == 'M':
        #                 prob = 1 - 1 / self.parser.gene_len
        #             elif s1 == 'P':
        #                 prob = 1 / self.parser.gene_len
        #             else:
        #                 prob = 0
        #         elif s0 == 'P':
        #             if s1 == 'I':
        #                 prob = 1/3
        #             elif s1 == 'P':
        #                 prob = 2/3
        #             else:
        #                 prob = 0
        #         else:
        #             raise ValueError("unknown state")
        #
        #         self.transition_prob[(s0, s1)] = prob
        for s0 in STATES:
            for s1 in STATES:
                prob = 0
                if s0 == 'I':
                    if s1 == 'I':
                        prob = 1 - 1 / self.parser.intergenic_len
                    elif s1 == 'S':
                        prob = 1 / self.parser.intergenic_len
                    else:
                        prob = 0
                elif s0 == 'S':
                    if s1 == 'M':
                        prob = 1
                    else:
                        prob = 0
                elif s0 == 'M':
                    if s1 == 'M':
                        prob = 1 - 1 / self.parser.gene_len
                    elif s1 == 'P':
                        prob = 1 / self.parser.gene_len
                    else:
                        prob = 0
                elif s0 == 'P':
                    if s1 == 'I':
                        prob = 1
                    else:
                        prob = 0
                else:
                    raise ValueError("unknown state")

                self.transition_prob[(s0, s1)] = prob

        assert len(self.transition_prob) == 16

    def build_emission_prob(self):
        """
        Given state, probability of observation

        :param state: choice I S M P
        :param observation: choice A T C G
        :return:
        """
        for s in STATES:
            for n in NUCLEOTIDES:
                prob = 0
                if s == 'I':
                    prob = self.parser.intergenic_nucleotide_freq[n]
                elif s == 'S':
                    freq = self.codon_freq_to_true_freq(
                        self.parser.start_codons_cnt
                    )
                    prob = freq[n]

                elif s == 'M':
                    prob = self.parser.gene_nucleotide_freq[n]

                elif s == 'P':
                    freq = self.codon_freq_to_true_freq(
                        self.parser.stop_codons_cnt
                    )
                    prob = freq[n]

                else:
                    raise ValueError("Unknown key! Got {}".format(s))
                self.emission_prob[(s, n)] = prob

        assert len(self.emission_prob) == 16


    def viterbi(self, seq_name):

        table = [[Cell() for _ in range(len(self.input_gene_stock[seq_name].seq))]
                 for _ in range(len(self.states))]
        seq = self.input_gene_stock[seq_name].seq

        def get_previous_cell(i, this_state):
            """
            compute V(this_state, i) = max_{s in S}(V[s, i-1] +
            Transition(s, this state)

            :param i: index of observations
            :param this_state:
            :return:
            """
            #
            max_value = -math.inf
            prev_coor = (0, 0)

            for _s, _sid in self.states.items():
                # get the max previous state in the table

                if table[_sid][i-1].value + \
                    self.get_transition_prob(_s, this_state)  >= max_value:
                        prev_coor = (_sid, i-1)
                        max_value = table[_sid][i-1].value

            return prev_coor, max_value

        def backtracking():
            max_value = -math.inf
            max_coor = (None, None)
            cur_cell = None
            out_seq = []
            for _s, _sid in self.states.items():
                if table[_sid][len(seq) -1].value > max_value:
                    max_coor = (_s, len(seq) -1)
                    max_value = table[_sid][len(seq) -1].value
                    cur_cell = table[_sid][len(seq) -1]
            # print(max_coor)
            assert cur_cell is not None

            while cur_cell.pointer is not None:
                out_seq.append(cur_cell.coor[0])
                prev_cell = table[cur_cell.pointer[0]][cur_cell.pointer[1]]
                cur_cell = prev_cell

            # last cell!
            out_seq.append(cur_cell.coor[0])
            out_seq.reverse()
            print(len(out_seq))
            print(len(seq))
            assert len(out_seq) == len(seq)
            return out_seq

        # first column
        for s, sid in self.states.items():
            this_cell = table[sid][0]
            this_cell.coor = (s, 0)
            this_cell.value = self.get_initial_prob(s) + self.get_emission_prob(s, seq[0])

        # filling the table column by column
        for obs_idx in range(1, len(seq)):
            for s, sid in self.states.items():
                # print(s)
                this_cell = table[sid][obs_idx]
                this_cell.coor = (s, obs_idx)
                this_cell.pointer, prev_max = get_previous_cell(obs_idx, s)
                this_cell.value = prev_max + self.get_emission_prob(s, seq[obs_idx])

        self.input_gene_stock[seq_name].prediction = backtracking()


    def gen_gff_file(self, write_file):
        with open(write_file, 'w') as f:
            for gn, gs in self.input_gene_stock.items():
                gs.prediction_sequence_to_loc()



    @staticmethod
    def codon_freq_to_true_freq(codon_freq):

        nu_freq = {}
        for codon, cnt in codon_freq.items():
            for nu in codon:
                try:
                    nu_freq[nu] += cnt
                except KeyError:
                    nu_freq[nu] = cnt

        if len(nu_freq) < 4:
            for i in NUCLEOTIDES:
                if nu_freq.get(i) is None:
                    nu_freq[i] = 0

        return {i : f / sum(nu_freq.values()) for i, f in nu_freq.items()}


    def get_emission_prob(self, state, observation):
        """
        Given state, probability of observation

        :param state: choice I S M P
        :param observation: choice A T C G
        :return:
        """
        if self.emission_prob[(state, observation)] == 0:
            return -math.inf
        return math.log(self.emission_prob[(state, observation)])

    def get_initial_prob(self, state):
        if self.initial_prob[state] == 0:
            return -math.inf
        return math.log(self.initial_prob[state])

    def get_transition_prob(self, state_0, state_1):
        if self.transition_prob[(state_0, state_1)] == 0:
            return -math.inf
        return math.log(self.transition_prob[(state_0, state_1)])



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
        print(parser.start_codons_cnt)
        print("stop codons frequency: ")
        print(parser.stop_codons_cnt)


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
    q1(parser, verbose=True)
    q2(parser)


if __name__ == '__main__':
    main()
    # print(-math.inf * 0)
