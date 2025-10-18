
# version with PWA as dist

from nummatrix import NumMatrix
from hierarchicalclustering import HierarchicalClustering
from myseq import MySeq
from pairwisealignment import PairwiseAlignment
from substmatrix import SubstMatrix
from myalign import MyAlign
from multiplealign import MultipleAlignment

class UPGMA:
    
    def __init__(self, seqs, alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.create_mat_dist()
        
    def create_mat_dist(self):
        self.matdist = NumMatrix(len(self.seqs), len(self.seqs))
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needleman_Wunsch(s1, s2)
                alin = self.alseq.recover_align()
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)
                    if (col[0] != col[1]): ncd += 1
                self.matdist.set_value(i, j, ncd)

    def run(self):
        ch = HierarchicalClustering(self.matdist)
        t = ch.execute_clustering()
        return t
    
    def create_mat_dist_align(self, align):
        self.matdist = NumMatrix(len(align.listseqs), len(align.listseqs))
        for i in range(len(align.listseqs)):
            for j in range(i, len(align.listseqs)):
                s1 = align.listseqs[i]
                s2 = align.listseqs[j]
                alin = MyAlign([s1,s2], align.al_type)
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)
                    if (col[0] != col[1]): ncd += 1
                self.matdist.set_value(i, j, ncd/len(alin))

    def run_align(self, align):
        self.create_mat_dist_align(align)
        ch = HierarchicalClustering(self.matdist)
        t = ch.execute_clustering()
        return t
    
def test():
    seq1 = MySeq("ATAGCGAT")
    seq2 = MySeq("ATAGGCCT")
    seq3 = MySeq("CTAGGCCC")
    seq4 = MySeq("CTAGGCCT")
    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    alseq = PairwiseAlignment(sm, -2)
    up = UPGMA([seq1, seq2, seq3, seq4], alseq)
    arv = up.run()
    arv.print_tree()
    
def test_align():
    seq1 = MySeq("ACATATCAT")
    seq2 = MySeq("AGATATTAG")
    seq3 = MySeq("AACAGATCT")
    seq4 = MySeq("GCATCGATT")
    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    alseq = PairwiseAlignment(sm, -2)
    ma = MultipleAlignment([seq1, seq2, seq3, seq4], alseq)
    al = ma.align_consensus()
    print(al)
    up = UPGMA([seq1, seq2, seq3, seq4], alseq)
    arv = up.run_align(al)
    arv.print_tree()
 