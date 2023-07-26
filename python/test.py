import unittest
import numpy as np
from solver import Solver, Determinant, trace
from configuration import Configuration

from segment import is_segment_proper, Segment, segments

from segment import AntiSegment

S = Solver()

class TestTinycthyb(unittest.TestCase):

    def test_det(self):
        c = Configuration([1.0], [3.0]) 
        d = Determinant(c, S.e) 
        np.testing.assert_allclose(d.mat, np.array([[0.06230103101396599]]))
        np.testing.assert_allclose(d.value, 0.06230103101396599)

        c = Configuration([19.0], [1.0]) 
        d = Determinant(c, S.e) 
        np.testing.assert_allclose(d.mat, np.array([[-0.06230103101396599]]))
        np.testing.assert_allclose(d.value, -0.06230103101396599)

        c = Configuration([1.0, 5.0], [3.0, 7.0])
        d = Determinant(c, S.e)
        np.testing.assert_allclose(d.mat, np.array([[0.06230103101396599,-0.06230103101396599],
                                                    [0.03003466159282729, 0.06230103101396599]])                                   )
        np.testing.assert_allclose(d.value, 0.005752608848791858)

        c = Configuration([2.0, 19.0], [5.0, 1.0])
        d = Determinant(c, S.e)
        np.testing.assert_allclose(d.mat, np.array([[0.048396641468896745,-0.03003466159282729],
                                                    [-0.08525324452511202,-0.06230103101396599]])                                   )
        np.testing.assert_allclose(d.value, -0.005575713010127824)


    def test_trace(self):
        c = Configuration([1.0], [3.0]) 
        np.testing.assert_allclose(trace(c,S.e), 1.0)
        c = Configuration([19.0], [1.0]) 
        np.testing.assert_allclose(trace(c,S.e), -1.0)
        c = Configuration([1.0, 5.0], [3.0, 7.0])
        np.testing.assert_allclose(trace(c,S.e), 1.0)
        c = Configuration([2.0, 19.0], [5.0, 1.0])
        np.testing.assert_allclose(trace(c,S.e), -1.0)

    def test_is_segment_proper(self):
        c = Configuration([1.0], [3.0]) 
        self.assertEqual(is_segment_proper(c), True)
        c = Configuration([19.0], [1.0]) 
        self.assertEqual(is_segment_proper(c), True)
        c = Configuration([1.0, 5.0], [3.0, 7.0])
        self.assertEqual(is_segment_proper(c), True)
        c = Configuration([2.0, 19.0], [5.0, 1.0])
        self.assertEqual(is_segment_proper(c), True)
        c = Configuration([1.0, 3.0], [2.0, 4.0])
        self.assertEqual(is_segment_proper(c), True)
        c = Configuration([1.0, 3.0, 5.0], [0.5, 2.0, 4.0])
        self.assertEqual(is_segment_proper(c), True)

    def test_segment_iterator(self):
        c = Configuration([1.0, 3.0], [2.0, 4.0])

        segs  = [Segment(1.0, 2.0), Segment(3.0, 4.0)]

        for s1, s2 in zip(segments(c), segs):
            self.assertEqual(s1.t_i, s2.t_i)
            self.assertEqual(s1.t_f, s2.t_f)

        c = Configuration([1.0, 3.0, 5.0], [0.5, 2.0, 4.0])

    def test_segment_length(self):
        self.assertEqual(Segment(1.0, 2.0).len(10.0), 1.0)
        self.assertEqual(Segment(9.5, 0.5).len(10.0), 1.0)

    def test_antisegment_length(self):
        self.assertEqual(AntiSegment(1.0, 2.0).len(10.0), 1.0)
        self.assertEqual(AntiSegment(9.5, 0.5).len(10.0), 1.0)


if __name__ == "__main__":
    unittest.main()





#print("Anti-Segments")
#for s in antisegments(c):
#    print(s)

#
#print("Segments")
#for s in segments(c):
#    print(s)

#print("Anti-Segments")
#for s in antisegments(c):
#    print(s)
##
#print(onsegment(1.5, c))
#print(onsegment(3.5, c))
#print(onsegment(6.0, c))
#print(onsegment(0.2, c))
#print(onsegment(0.6, c))
#print(onsegment(2.5, c))
#print(onsegment(4.5, c))
#
#print(onantisegment(1.5, c))
#print(onantisegment(3.5, c))
#print(onantisegment(6.0, c))
#print(onantisegment(0.2, c))
#print(onantisegment(0.6, c))
#print(onantisegment(2.5, c))
#print(onantisegment(4.5, c))
#
#   
#print(segments(c).indices(0))
#print(segments(c).indices(1))
#print(segments(c).indices(2))
#
#for idx in range(len(c)):
#    c_tmp = deepcopy(c)
#    remove_segment(c_tmp, idx)
#    for s in segments(c_tmp): print(s)
#
#for idx in range(len(c)):
#    c_tmp = deepcopy(c)
#    remove_antisegment(c_tmp, idx)
#    for s in antisegments(c_tmp): print(s)
#
#c = Configuration([1.0, 2.0], [3.0, 4.0])

#print(c)
#print(is_segment_proper(c))

#assert not is_segment_proper(c)
