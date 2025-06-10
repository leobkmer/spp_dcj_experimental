import unittest as ut
import data_utils as du

class TestUnimogUtils(ut.TestCase):
    def test_invert(self):
        chrm = [(du.ORIENT_NEGATIVE,'A'),(du.ORIENT_POSITIVE,'B'),(du.ORIENT_POSITIVE,'C')]
        chrm_ = du.invert_chromosome(chrm)
        self.assertEqual(len(chrm_),len(chrm))
        self.assertEqual([(du.ORIENT_NEGATIVE,'C'),(du.ORIENT_NEGATIVE,'B'),(du.ORIENT_POSITIVE,'A')],chrm_)
    def test_canonicize_circular(self):
        chrm = [(du.ORIENT_NEGATIVE,'m3'),(du.ORIENT_NEGATIVE,'m1'),(du.ORIENT_POSITIVE,'m2')]
        chrm_ = du.canonicize_circular(chrm)
        self.assertEqual(len(chrm_),len(chrm))
        self.assertEqual(chrm_,[(du.ORIENT_POSITIVE,'m1'),(du.ORIENT_POSITIVE,'m3'),(du.ORIENT_NEGATIVE,'m2')])
    def test_canonicize_linear(self):
        chrm = [(du.ORIENT_NEGATIVE,'A'),(du.ORIENT_POSITIVE,'B'),(du.ORIENT_NEGATIVE,'C')]
        chrm_ = du.invert_chromosome(chrm)
        chrm__ = du.canonicize_linear(chrm)
        self.assertEqual(chrm_,chrm__)
    def test_canonicize_circular_nothing(self):
        chrm = [(du.ORIENT_POSITIVE,'A'),(du.ORIENT_POSITIVE,'B')]
        chrm_ = du.canonicize_circular(chrm)
        self.assertEqual(chrm,chrm_)
    def test_canonicize_linear_nothing(self):
        chrm = [(du.ORIENT_POSITIVE,'A')]
        chrm_ = du.canonicize_linear(chrm)
        self.assertEqual(chrm,chrm_)
    def test_canonicize_single_marker(self):
        chrm = [(du.ORIENT_NEGATIVE,'x')]
        chrm_l = du.canonicize_linear(chrm)
        chrm_c = du.canonicize_circular(chrm)
        self.assertEqual([(du.ORIENT_POSITIVE,'x')],chrm_l)
        self.assertEqual([(du.ORIENT_POSITIVE,'x')],chrm_c)
    def test_canonicize_circular_strict(self):
        chrm = [(du.ORIENT_NEGATIVE,'y'),(du.ORIENT_POSITIVE,'y'),(du.ORIENT_POSITIVE,'z'),(du.ORIENT_NEGATIVE,'x'),(du.ORIENT_POSITIVE,'x'),(du.ORIENT_NEGATIVE,'y'),(du.ORIENT_POSITIVE,'z'),(du.ORIENT_POSITIVE,'x')]
        canonizations = []
        for i in range(len(chrm)):
            rotation = chrm[i::]+chrm[0:i]
            canonizations.append(du.canonicize_circular(rotation,strict=True))
            canonizations.append(du.canonicize_circular(du.invert_chromosome(rotation),strict=True))
        first = canonizations[0]
        for c in canonizations:
            self.assertEqual(first,c)
    def test_canonicize_unimog(self):
        chr1 = [(du.ORIENT_NEGATIVE,'y'),(du.ORIENT_POSITIVE,'x'),(du.ORIENT_POSITIVE,'z')]
        chr1_lc = chr1
        chr1_cc = chr1[1::]+chr1[0:1]
        chr2 = [(du.ORIENT_NEGATIVE,'a')]
        chr2_c = du.invert_chromosome(chr2)
        chr3 = [(du.ORIENT_POSITIVE,'c'),(du.ORIENT_NEGATIVE,'b')]
        chr3_c = du.invert_chromosome(chr3)

        
        unimog = [('C',[(du.CHR_LINEAR,chr2),(du.CHR_CIRCULAR,chr1)]),('B',[(du.CHR_CIRCULAR,chr3)])]
        man_canon = [('B',[(du.CHR_CIRCULAR, chr3_c)]),('C',[(du.CHR_CIRCULAR,chr1_cc),(du.CHR_LINEAR,chr2_c)])]
        canon_unimog = du.canonicize_unimog(unimog)
        self.assertEqual(man_canon,canon_unimog)

class TestTrees(ut.TestCase):
    def test_cp_to_pc(self):
        tree = dict([('A','X'),('B','Y'),('C','Y'),('Y','X')])
        root,pc_tree = du.cp_to_pc(tree)
        self.assertEqual(root,'X')
        self.assertEqual(set(pc_tree.keys()),set(['X','Y']))
        self.assertEqual(set(pc_tree['X']),{'A','Y'})
        self.assertEqual(set(pc_tree['Y']),{'B','C'})
    def test_lca_trace(self):
        tree = dict([('A','X'),('B','Y'),('C','Y'),('Y','X'),('X','Z'),('D','Z')])
        traces = {('A','B') : {'A','B','Y','X'},
                  ('B','D'):{'B','X','Y','Z','D'},
                  ('A','D'):{'A','X','Z','D'},
                  ('B','Y'):{'B','Y'},
                  ('X','X'):{'X'}}
        for x in "ABCDXYZ":
            for y in "ABCDXYZ":
                trace = du.lca_trace_cp_tree(tree,x,y)
                trace_r = du.lca_trace_cp_tree(tree,y,x)
                if (x,y) in traces:
                    self.assertEqual(set(trace),traces[(x,y)])
                    self.assertEqual(set(trace_r),traces[(x,y)])
                self.assertEqual(set(trace_r),set(trace))
            #no duplicate nodes
            self.assertEqual(len(set(trace)),len(trace))
            self.assertEqual(len(set(trace_r)),len(trace_r))
        


if __name__ == '__main__':
    ut.main()