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
        anchor = 1
        if len(self.gene_loc) == 0:
            self.intergenic_loc = [(1, self.seq_len)]
        else:
            for start, end in self.gene_loc:
                if start >= next_start:
                    anchor = next_start
                    next_start = end + 1
                if anchor == 1:
                    continue
                assert anchor <= start
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

        state_seq = []
        gene_loc = self.gene_loc

        anchor = 1
        i = 0
        for start, end in gene_loc:

            while anchor < start -1:
                state_seq.append('I')
                anchor += 1
                i += 1
            while anchor <= end :
                state_seq.append('M')
                anchor += 1
                i += 1
            # anchor = end
        while anchor < self.seq_len:
            state_seq.append('I')
            anchor += 1

        # print(len(state_seq), self.seq_len)

        # assert len(state_seq) == self.seq_len
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
    # TODO: in SME states multiple of 3
    def prediction_sequence_to_loc(self):
        gene_turning_pt = []
        gene_loc = []
        gene_start_idx = 0 if self.prediction[0] == 'S' else None
        gene_end_idx = self.seq_len if self.prediction[-1] == 'P' else None

        if gene_start_idx is not None:
            gene_turning_pt.append(gene_start_idx)

        # start recovering index
        for i in range(len(self.prediction) - 1):
            if self.prediction[i] == 'I' and self.prediction[i + 1] == 'S':
                gene_turning_pt.append(i + 2)
            elif self.prediction[i] == 'P' and self.prediction[i + 1] == 'I':
                gene_turning_pt.append(i + 1)

        if gene_end_idx is not None:
            gene_turning_pt.append(gene_end_idx)

        for j in range(0, len(gene_turning_pt) - 1, 2):
            gene_loc.append((gene_turning_pt[j], gene_turning_pt[j + 1]))
        self.gene_loc = gene_loc
        return gene_loc
