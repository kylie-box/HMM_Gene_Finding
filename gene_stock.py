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
        for i in range(len(state_seq) - 1):
            seq += state_seq[i]
            if state_seq[i] == 'I' and state_seq[i + 1] == 'S':
                gene_loc.append(i + 1)
            elif state_seq[i] == 'P' and state_seq[i + 1] == 'I':
                gene_loc.append(i + 1)

        # print(gene_loc)

    def prediction_sequence_to_loc(self):
        gene_turing_pt = []
        gene_loc = []
        gene_start_idx = 0 if self.prediction[0] == 'S' else None
        gene_end_idx = self.seq_len if self.prediction[-1] == 'P' else None

        if gene_start_idx is not None:
            gene_turing_pt.append(gene_start_idx)

        for i in range(len(self.prediction) - 1):
            if self.prediction[i] == 'I' and self.prediction[i + 1] == 'S':
                gene_turing_pt.append(i + 1)
            elif self.prediction[i] == 'P' and self.prediction[i + 1] == 'I':
                gene_turing_pt.append(i + 1)
        if gene_end_idx is not None:
            gene_turing_pt.append(gene_end_idx)

        for j in range(0, len(gene_turing_pt) - 1, 2):
            gene_loc.append((gene_turing_pt[j], gene_turing_pt[j + 1]))

        print(gene_turing_pt)

        print("\n\n\n", gene_turing_pt[997:1010],gene_turing_pt[1150:1056])