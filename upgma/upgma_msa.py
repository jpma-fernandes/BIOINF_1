# version with distance = num of diff pos from MSA

from nummatrix import NumMatrix
from hierarchicalclustering import HierarchicalClustering
from myseq import MySeq

class UPGMA:
    
    def __init__(self, seqs):
        self.seqs = seqs
        self.create_mat_dist()
        
    def create_mat_dist(self):
        self.matdist = NumMatrix(len(self.seqs), len(self.seqs))
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                ncd = 0
                for k in range(len(s1)):
                    if (s1[k] != s2[k]): ncd += 1
                self.matdist.set_value(i, j, ncd)

    def run(self):
        ch = HierarchicalClustering(self.matdist)
        t = ch.execute_clustering()
        return t
    
    
    
def test():
    seq1 = MySeq("ATAGCGAT")
    seq2 = MySeq("ATAGGCCT")
    seq3 = MySeq("CTAGGCCC")
    seq4 = MySeq("CTAGGCCT")
    up = UPGMA([seq1, seq2, seq3, seq4])
    arv = up.run()
    arv.print_tree()

test()

 