from .pairwisealignment import PairwiseAlignment
from .myalign import MyAlign
from .myseq import MySeq
from .substmatrix import SubstMatrix

class MultipleAlignment():
    def __init__(self, seqs, alignseq):
        self.seqs = seqs
        self.alignpars = alignseq
        
    def add_seq_alignment(self, alignment , seq):
        res = []
        for i in range(len(alignment.listseqs)+1):
            res.append("")
        cons = MySeq(alignment.consensus(),alignment.al_type)
        #print("cons",cons)
        self.alignpars.needleman_Wunsch(cons, seq)
        align2 = self.alignpars.recover_align()
        #print("align2",align2)
        orig = 0
        for i in range(len(align2)):
            if align2[0,i]== "-":
                for k in range(len(alignment.listseqs)):
                    res[k] += "-"
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k,orig]
                orig+=1
        res[len(alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res, alignment.al_type)

    def align_consensus(self):
        self.alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recover_align()
        #print(res)
        for i in range(2, len(self.seqs)):
            res = self.add_seq_alignment(res, self.seqs[i])
            #print(res)
        return res
    
def testMSA():
    s1 = MySeq("ATATCCG")
    s2 = MySeq("TCCG")
    s3 = MySeq("ATGTACTG")
    s4 = MySeq("ATGTCTG")
    sm = SubstMatrix()
    sm.create_submat(0,-1,"ACGT")
    aseq = PairwiseAlignment(sm,-1)
    ma = MultipleAlignment([s1,s2,s3,s4], aseq)
    al = ma.align_consensus()
    print(al)
    print(ma)

testMSA()
    
