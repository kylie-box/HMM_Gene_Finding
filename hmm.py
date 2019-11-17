import math

from gene_stock import GeneStock

NUCLEOTIDES = ['A', 'T', 'C', 'G']
STATES = ['I', 'S', 'M', 'P']


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
        self.value = None  # store the log value of the prob
        self.pointer = None  # previous cell
        self.end_codon = False

    def __str__(self):
        output_msg = "coor: {}\nvalue: {}\npointer: {}\n".format(self.coor,
                                                                 self.value,
                                                                 self.pointer)
        return output_msg

    def __abs__(self):
        return self.value


class HMM:
    def __init__(self, parser, input_fasta):
        self.states = ({s: i for s, i in zip(STATES, range(4))})
        self.observations = NUCLEOTIDES
        # I: intergenic; S: Start; M: Middle; P: Stop
        self.observations = set()
        self.emission_prob = {}  # s:obs given state s observe obs
        self.transition_prob = {}  # s0:s1 given s0 transition to s1
        self.initial_prob = {'I': 1, 'S': 0, 'M': 0, 'P': 0}
        self.parser = parser
        self.input_gene_stock = {}
        self.sequence_loading(fasta_file=input_fasta)

    def sequence_loading(self, fasta_file):
        with open(fasta_file, 'r') as f:
            gene_name = None
            gene_seq = ""
            i = 0
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
                        prob = 1 - 3 / self.parser.gene_len
                    elif s1 == 'P':
                        prob = 3 / self.parser.gene_len
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
        # for s in STATES:
        #     for n in NUCLEOTIDES:
        #         prob = 0
        #         if s == 'I':
        #             prob = self.parser.intergenic_nucleotide_freq[n]
        #         elif s == 'S':
        #             freq = self.codon_freq_to_true_freq(
        #                 self.parser.start_codons_cnt
        #             )
        #             prob = freq[n]
        #
        #         elif s == 'M':
        #             prob = self.parser.gene_nucleotide_freq[n]
        #
        #         elif s == 'P':
        #             freq = self.codon_freq_to_true_freq(
        #                 self.parser.stop_codons_cnt
        #             )
        #             prob = freq[n]
        #
        #         else:
        #             raise ValueError("Unknown key! Got {}".format(s))
        #         self.emission_prob[(s, n)] = prob

        s = 'I'
        for n in NUCLEOTIDES:
            prob = self.parser.intergenic_nucleotide_freq[n]
            self.emission_prob[(s, n)] = prob

        s = 'S'
        print(len(self.parser.start_codons_freq))
        assert len(self.parser.start_codons_freq) == 64
        for codon in self.parser.start_codons_freq.keys():
            prob = self.parser.start_codons_freq[codon]
            self.emission_prob[(s, codon)] = prob

        s = 'M'
        print(len(self.parser.gene_codon_freq))
        for codon in self.parser.gene_codon_freq.keys():
            prob = self.parser.gene_codon_freq[codon]
            self.emission_prob[(s, codon)] = prob

        s = 'P'
        print(len(self.parser.stop_codons_freq))
        for codon in self.parser.stop_codons_freq.keys():
            prob = self.parser.stop_codons_freq[codon]
            self.emission_prob[(s, codon)] = prob

        print(len(self.emission_prob))  # 196

    def viterbi(self, seq_name):

        table = [
            [Cell() for _ in range(len(self.input_gene_stock[seq_name].seq))]
            for _ in range(len(self.states))]
        seq = self.input_gene_stock[seq_name].seq

        def get_previous_cell(i, this_state):
            """
            compute V(this_state, i) = max_{s in S}(V[s, i-1] +
            Transition(s, this state)

            :param i: index of observations; or start codon position
            :param this_state:
            :return:
            """
            # TODO: Change the multiple of 3... in S/M/P states
            max_value = -math.inf
            prev_coor = (0, 0)

            if i >= 3:
                # we have to consider the codon position when this state is not
                # Intergenic
                if this_state == 'I':
                    # can only come from I or P, if from P, we only store the
                    # first codon emission prob
                    _sid = self.states['I']
                    score_from_I = table[_sid][i - 1].value + \
                                self.get_transition_prob('I',this_state)

                    _sid = self.states['P']
                    score_from_P = table[_sid][i - 3].value + \
                                self.get_transition_prob('P', this_state)

                    if score_from_P > score_from_I:
                        prev_coor = (self.states['P'], i - 3)
                        max_value = score_from_P
                    else:
                        prev_coor = (self.states['I'], i - 1)
                        max_value = score_from_I

                elif this_state == 'S':
                    _sid = self.states['I']
                    score_from_I = table[_sid][i - 1].value + \
                                    self.get_transition_prob('I', this_state)

                    prev_coor = (self.states['I'], i - 1)
                    max_value = score_from_I

                elif this_state == 'M':
                    _sid = self.states['S']
                    score_from_S = table[_sid][i - 3].value + \
                                    self.get_transition_prob('S', this_state)

                    _sid = self.states['M']
                    score_from_M = table[_sid][i - 3].value + \
                                    self.get_transition_prob('M', this_state)

                    if score_from_S > score_from_M:
                        prev_coor = (self.states['S'], i - 3)
                        max_value = score_from_S
                    else:
                        prev_coor = (self.states['M'], i - 3)
                        max_value = score_from_M

                elif this_state == 'P':
                    _sid = self.states['M']
                    score_from_M = table[_sid][i - 3].value +\
                                    self.get_transition_prob('M', this_state)

                    prev_coor = (self.states['M'], i - 3)
                    max_value = score_from_M


                    # for _s, _sid in self.states.items():
                    #     # get the max previous state in the table
                    #
                    #     end_codon = True
                    #     if table[_sid][i - 1].value + \
                    #             self.get_transition_prob(_s,
                    #                                      this_state) >= max_value:
                    #         prev_coor = (_sid, i - 1)
                    #         max_value = table[_sid][i - 1].value


            # first 3 columns
            else:
                for _s, _sid in self.states.items():

                    # start filling normal viterbi
                    if table[_sid][i - 1].value + \
                            self.get_transition_prob(_s,
                                                     this_state) >= max_value:
                        prev_coor = (_sid, i - 1)
                        max_value = table[_sid][i - 1].value

            return prev_coor, max_value
            # for _s, _sid in self.states.items():
            #     # get the max previous state in the table
            #
            #     if table[_sid][i - 1].value + \
            #             self.get_transition_prob(_s, this_state) >= max_value:
            #         prev_coor = (_sid, i - 1)
            #         max_value = table[_sid][i - 1].value
            #
            # return prev_coor, max_value

        def backtracking():

            # TODO: codon traceback broken
            max_value = -math.inf
            max_coor = (None, None)
            cur_cell = None
            out_seq = []
            for _s, _sid in self.states.items():
                if table[_sid][len(seq) - 1].value > max_value:
                    max_coor = (_s, len(seq) - 1)
                    max_value = table[_sid][len(seq) - 1].value
                    cur_cell = table[_sid][len(seq) - 1]
            # print(max_coor)
            assert cur_cell is not None

            while cur_cell.pointer is not None:
                out_seq.append(cur_cell.coor[0])
                prev_cell = table[cur_cell.pointer[0]][cur_cell.pointer[1]]
                cur_cell = prev_cell

            # last cell!
            out_seq.append(cur_cell.coor[0])
            out_seq.reverse()
            # print(len(out_seq))
            # print(len(seq))
            # assert len(out_seq) == len(seq)
            return out_seq

        # first column
        for s, sid in self.states.items():
            this_cell = table[sid][0]
            this_cell.coor = (s, 0)
            if s == 'I':
                this_cell.value = self.get_initial_prob(s) + self.get_emission_prob(
                    s, seq[0])
            else:
                this_cell.value = -math.inf
            this_cell.end_codon = True


        # filling the table column by column
        for obs_idx in range(1, len(seq)):
            for s, sid in self.states.items():
                # print(s)
                #                 this_cell = table[sid][obs_idx]
                # this_cell.coor = (s, obs_idx)
                # this_cell.pointer, prev_max, this_cell.end_codon = \
                #     get_previous_cell(obs_idx, s)
                # if s == 'I':
                #
                #     this_cell.value = prev_max + self.get_emission_prob(s,seq[obs_idx])
                # elif s != 'I' :
                #     if this_cell.end_codon:
                #     # reading the frame
                #         codon = "".join(seq[obs_idx:obs_idx + 3])
                #         if len(codon) < 3:
                #             this_cell.pointer, prev_max, this_cell.end_codon = \
                #                 get_previous_cell(obs_idx, s)
                #             this_cell.value = prev_max - math.inf
                #         else:
                #             this_cell.pointer, prev_max, this_cell.end_codon = \
                #                 get_previous_cell(obs_idx, s)
                #             this_cell.value = prev_max + self.get_emission_prob(s,codon)
                #
                #     else:
                #         this_cell.value = prev_max
                #         this_cell.pointer = (sid, obs_idx - 1)
                #         if get_previous_cell()
                this_cell = table[sid][obs_idx]
                this_cell.coor = (s, obs_idx)

                if s == 'I':
                    prev_max_pointer, prev_max = get_previous_cell(obs_idx, s)
                    this_cell.value = prev_max + self.get_emission_prob(s, seq[
                        obs_idx])
                    this_cell.pointer = prev_max_pointer
                else:
                    # codon = "".join(seq[obs_idx:obs_idx + 3])
                    # if len(codon) < 3:
                    #     this_cell.pointer, prev_max, this_cell.end_codon = \
                    #         get_previous_cell(obs_idx, s)
                    #     this_cell.value = prev_max - math.inf
                    #     continue
                    # if this_cell.end_codon:
                    #
                    #     this_cell.value = prev_max + self.get_emission_prob(s, codon)
                    #     this_cell.pointer = prev_max_pointer
                    #
                    # else:
                    #     if prev_max < table[sid][obs_idx - 1].value:
                    #         this_cell.value = table[sid][obs_idx - 1].value
                    #         this_cell.pointer = (sid, obs_idx - 1)
                    #         this_cell.end_codon = True
                    #
                    #     else:
                    #         this_cell.value = prev_max + self.get_emission_prob(s, codon)
                    #         this_cell.pointer = prev_max_pointer
                    this_cell.pointer, prev_max = get_previous_cell(obs_idx, s)
                    codon = "".join(seq[obs_idx:obs_idx + 3])
                    if len(codon) < 3:
                        this_cell.value = prev_max - math.inf
                        continue
                    this_cell.value = prev_max + self.get_emission_prob(s,codon)

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

        return {i: f / sum(nu_freq.values()) for i, f in nu_freq.items()}

    def get_emission_prob(self, state, observation):
        """
        Given state, probability of observation

        :param state: choice I S M P
        :param observation: choice A T C G for intergenic or codon usage table
        for S/M/P
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
