import unittest
from BITS.seq.align import EdlibRunner

from BITS.seq.cigar import Cigar
from BITS.seq.util import revcomp_seq


class TestEdlibRunner(unittest.TestCase):
    def setUp(self):
        self.x = "acagttaccgt"
        self.y = "cagatacc"
        self.z = "accgacag"

        self.cigar_x2y_global = Cigar("1D3=1X4=2D")

        self.cigar_x2y_glocal = Cigar("3=1X4=")
        self.x_glocal_start = 1
        self.x_glocal_end = len(self.x) - 2
        self.x_glocal = self.x[self.x_glocal_start:self.x_glocal_end]

        self.cigar_x2y_prefix = Cigar("1D3=1X4=")
        self.x_prefix_end = len(self.x) - 2
        self.x_prefix = self.x[:self.x_prefix_end]

        self.cigar_xx2z = (Cigar("4=1D4=2D"), Cigar("2D4=1D4="))
        self.x_cyclic_start = (6, 4)
        self.x_cyclic = ("accgtacagtt", "ttaccgtacag")
        self.cigar_zz2x = (Cigar("4=2I4=1I"), Cigar("1I4=2I4="))
        self.z_cyclic_start = (4, 3)
        self.z_cyclic = ("acagaccg",)

        """
        x vs y
        acagttaccgt
        *|||*||||**
        -cagatacc--

        xx vs z
        acagttaccgtacagttacgt
              ||||*||||**
              accg-acag--
        or
        acagttaccgtacagttacgt
            **||||*||||
            --accg-acag

        zz vs x
        accgacag--accg-acag
            ||||**||||*
            acagttaccgt
        or
        accg-acag--accgacag
            *||||**||||
            tacagttaccg
        """

    def test_global_forward(self):
        er = EdlibRunner("global", revcomp=False)

        aln = er.align(self.y, self.x)
        self.assertEqual(aln.cigar, self.cigar_x2y_global)
        self.assertEqual(aln.a_aligned_seq, self.y)
        self.assertEqual(aln.b_aligned_seq, self.x)

        aln = er.align(self.x, self.y)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.swap())
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertEqual(aln.b_aligned_seq, self.y)

    def test_global_revcomp(self):
        er = EdlibRunner("global", revcomp=True)

        aln = er.align(revcomp_seq(self.y), self.x)
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.revcomp())
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, 0)
        self.assertEqual(aln.a_end, len(self.y))
        self.assertEqual(aln.b_end, len(self.x))
        self.assertEqual(aln.a_aligned_seq, revcomp_seq(self.y))
        self.assertEqual(aln.b_aligned_seq, revcomp_seq(self.x))

        aln = er.align(self.y, revcomp_seq(self.x))
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_global)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, 0)
        self.assertEqual(aln.a_end, len(self.y))
        self.assertEqual(aln.b_end, len(self.x))
        self.assertEqual(aln.a_aligned_seq, self.y)
        self.assertEqual(aln.b_aligned_seq, self.x)

        aln = er.align(self.x, revcomp_seq(self.y))
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.swap())
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, 0)
        self.assertEqual(aln.a_end, len(self.x))
        self.assertEqual(aln.b_end, len(self.y))
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertEqual(aln.b_aligned_seq, self.y)

    def test_glocal_forward(self):
        er = EdlibRunner("glocal", revcomp=False)

        aln = er.align(self.y, self.x)
        self.assertEqual(aln.cigar, self.cigar_x2y_glocal)
        self.assertEqual(aln.a_aligned_seq, self.y)
        self.assertEqual(aln.b_aligned_seq, self.x_glocal)

        # inappropriate use-case (mapped seq is shorter)
        aln = er.align(self.x, self.y)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.swap())
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertEqual(aln.b_aligned_seq, self.y)

    def test_glocal_revcomp(self):
        er = EdlibRunner("glocal", revcomp=True)

        aln = er.align(revcomp_seq(self.y), self.x)
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_glocal.revcomp())
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, self.x_glocal_start)
        self.assertEqual(aln.a_end, len(self.y))
        self.assertEqual(aln.b_end, self.x_glocal_end)   # position is forward
        self.assertEqual(aln.a_aligned_seq, revcomp_seq(self.y))
        self.assertEqual(aln.b_aligned_seq, revcomp_seq(self.x_glocal))   # seq is revcomp

        aln = er.align(self.y, revcomp_seq(self.x))
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_glocal)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, len(self.x) - self.x_glocal_end)
        self.assertEqual(aln.a_end, len(self.y))
        self.assertEqual(aln.b_end, len(self.x) - self.x_glocal_start)   # rc(x) is forward
        self.assertEqual(aln.a_aligned_seq, self.y)
        self.assertEqual(aln.b_aligned_seq, self.x_glocal)   # seq is forward

        # inappropriate use-case (mapped seq is shorter)
        aln = er.align(self.x, revcomp_seq(self.y))
        self.assertEqual(aln.strand, 1)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.swap())
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.b_start, 0)
        self.assertEqual(aln.a_end, len(self.x))
        self.assertEqual(aln.b_end, len(self.y))
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertEqual(aln.b_aligned_seq, self.y)

    def test_prefix(self):
        er = EdlibRunner("prefix", revcomp=False)

        aln = er.align(self.y, self.x)
        self.assertEqual(aln.cigar, self.cigar_x2y_prefix)
        self.assertEqual(aln.a_aligned_seq, self.y)
        self.assertEqual(aln.b_aligned_seq, self.x_prefix)

        # inappropriate use-case (mapped seq is shorter)
        aln = er.align(self.x, self.y)
        self.assertEqual(aln.cigar, self.cigar_x2y_global.swap())
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertEqual(aln.b_aligned_seq, self.y)

    def test_global_forward_cyclic(self):
        er = EdlibRunner("global", revcomp=False, cyclic=True)

        aln = er.align(self.z, self.x)   # z vs xx
        self.assertIn(aln.cigar, self.cigar_xx2z)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.z))
        self.assertIn(aln.b_start, self.x_cyclic_start)
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, self.z)
        self.assertIn(aln.b_aligned_seq, self.x_cyclic)

        aln = er.align(self.x, self.z)   # x vs zz
        self.assertIn(aln.cigar, self.cigar_zz2x)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.x))
        self.assertIn(aln.b_start, self.z_cyclic_start)
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertIn(aln.b_aligned_seq, self.z_cyclic)

    def test_global_revcomp_cyclic(self):
        er = EdlibRunner("global", revcomp=True, cyclic=True)

        aln = er.align(revcomp_seq(self.z), self.x)   # rc(z) vs xx
        self.assertEqual(aln.strand, 1)
        self.assertIn(aln.cigar, [x.revcomp() for x in self.cigar_xx2z])
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.z))
        self.assertIn(aln.b_start, self.x_cyclic_start)
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, revcomp_seq(self.z))
        self.assertIn(aln.b_aligned_seq, [revcomp_seq(x) for x in self.x_cyclic])

        aln = er.align(self.z, revcomp_seq(self.x))   # z vs rc(xx)
        self.assertEqual(aln.strand, 1)
        self.assertIn(aln.cigar, self.cigar_xx2z)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.z))
        self.assertIn(aln.b_start, [len(self.x) - x for x in self.x_cyclic_start])
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, self.z)
        self.assertIn(aln.b_aligned_seq, self.x_cyclic)

        aln = er.align(revcomp_seq(self.x), self.z)   # rc(x) vs zz
        self.assertEqual(aln.strand, 1)
        self.assertIn(aln.cigar, [x.revcomp() for x in self.cigar_zz2x])
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.x))
        self.assertIn(aln.b_start, self.z_cyclic_start)
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, revcomp_seq(self.x))
        self.assertIn(aln.b_aligned_seq, [revcomp_seq(x) for x in self.z_cyclic])

        aln = er.align(self.x, revcomp_seq(self.z))   # x vs rc(zz)
        self.assertEqual(aln.strand, 1)
        self.assertIn(aln.cigar, self.cigar_zz2x)
        self.assertEqual(aln.a_start, 0)
        self.assertEqual(aln.a_end, len(self.x))
        self.assertIn(aln.b_start, [len(self.z) - x for x in self.z_cyclic_start])
        self.assertEqual(aln.b_start, aln.b_end)
        self.assertEqual(aln.a_aligned_seq, self.x)
        self.assertIn(aln.b_aligned_seq, self.z_cyclic)


if __name__ == "__main__":
    unittest.main()
