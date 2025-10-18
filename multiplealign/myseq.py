class MySeq:
    """ Class for Biological sequences. """
    def __init__(self, seq, seq_type = "DNA"):
        self.seq = seq
        self.seq_type = seq_type
        
    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.seq
    
    def __getitem__(self, n):
        return self.seq[n]
    
    def __getslice__(self, i, j):
        return self.seq[i:j]
    
    def print_sequence(self):
        print("Sequence: " + self.seq)

    def get_seq_biotype(self):
        return self.seq_type

    def show_info_seq(self):
        print("Sequence: " + self.seq + " biotype: " + self.seq_type)

    def count_occurrences(self, seq_search):
        return self.seq.count(seq_search)

    # class method to validate biotype update
    def set_seq_biotype(self, bt):
        biotype = bt.upper()
        if biotype == "DNA" or biotype == "RNA" or biotype == "PROTEIN":
            self.seq_type = biotype
        else:
            print("Non biological sequence type!")

    def alphabet(self):
        if (self.seq_type=="DNA"): return "ACGT"
        elif (self.seq_type=="RNA"): return "ACGU"
        elif (self.seq_type=="PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
        else: return None
        
    def validate(self):
        alp = self.alphabet()
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp: res = False
            else: i += 1
        return res

    def transcription(self):
        if (self.seq_type == "DNA"):
            return MySeq(self.seq.replace("T","U"), "RNA")
        else:
            return None


