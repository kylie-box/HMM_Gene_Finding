from gene_stock import GeneStock
from collections import Counter


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
        self.start_codons_freq = None
        self.stop_codons_cnt = None
        self.stop_codons_freq = None
        self.intergenic_len = None
        self.gene_len = None
        self.gene_nucleotide_freq = None

    def validation(self):
        for gs in self.gene_stock.values():
            if gs.seq_len == 0 and gs.seq is not None:
                gs.seq_len = len(gs.seq)
            if len(gs.gene_loc) > 0 and gs.gene_loc[-1][1] != gs.seq_len:
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
                                and gene_starts < crt_gene_stock.gene_loc[-1][
                            1]:
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

        tot = sum(start_codons.values())
        self.start_codons_freq = {n: cnt / tot for n, cnt in
                                  start_codons.items()}
        self.start_codons_cnt = start_codons

        if len(self.start_codons_freq) < 64:
            self.start_codons_freq = self.suppliment_codon_usage(
                self.start_codons_freq)

        tot = sum(stop_codons.values())
        self.stop_codons_freq = {n: cnt / tot for n, cnt in stop_codons.items()}
        self.stop_codons_cnt = stop_codons
        if len(self.stop_codons_freq) < 64:
            self.stop_codons_freq = self.suppliment_codon_usage(
                self.stop_codons_freq)

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