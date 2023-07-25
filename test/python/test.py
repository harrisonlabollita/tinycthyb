import unittest

from tinycthyb.solver import Determinant, trace
from tinycthyb.configuration import Configuration

class TestTinycthyb(unittest.TestCase):

    def test_config(self):
        c = Configuration([1.0], [3.0]) 
        d = Determinant(c, e) 
        np.testing.assert_allclose(d.mat, np.array([0.06230103101396599]))
        np.testing.assert_allclose(d.mat, 0.06230103101396599)
        np.testing.assert_allclose(trace(c,e), 1.0)



c = Configuration([19.0], [1.0]) 
d = Determinant(c, e) 
t = trace(c,e)

c = Configuration([1.0, 5.0], [3.0, 7.0])
d = Determinant(c, e)
t = trace(c,e)

c = Configuration([2.0, 19.0], [5.0, 1.0])
d = Determinant(c, e)
t = trace(c,e)
assert is_segment_proper(c)


c = Configuration([1.0, 3.0], [2.0, 4.0])

print(c)
assert is_segment_proper(c)

print("Segments")
for s in segments(c):
    print(s)

print("Anti-Segments")
for s in antisegments(c):
    print(s)

c = Configuration([1.0, 3.0, 5.0], [0.5, 2.0, 4.0])
print(c)
assert is_segment_proper(c)

print("Segments")
for s in segments(c):
    print(s)

print("Anti-Segments")
for s in antisegments(c):
    print(s)

print(onsegment(1.5, c))
print(onsegment(3.5, c))
print(onsegment(6.0, c))
print(onsegment(0.2, c))
print(onsegment(0.6, c))
print(onsegment(2.5, c))
print(onsegment(4.5, c))

print(onantisegment(1.5, c))
print(onantisegment(3.5, c))
print(onantisegment(6.0, c))
print(onantisegment(0.2, c))
print(onantisegment(0.6, c))
print(onantisegment(2.5, c))
print(onantisegment(4.5, c))

   
print(segments(c).indices(0))
print(segments(c).indices(1))
print(segments(c).indices(2))

for idx in range(len(c)):
    c_tmp = deepcopy(c)
    remove_segment(c_tmp, idx)
    for s in segments(c_tmp): print(s)

for idx in range(len(c)):
    c_tmp = deepcopy(c)
    remove_antisegment(c_tmp, idx)
    for s in antisegments(c_tmp): print(s)

c = Configuration([1.0, 2.0], [3.0, 4.0])

#print(c)
#print(is_segment_proper(c))

assert not is_segment_proper(c)
assert Segment(1.0, 2.0).len(10.0) == 1.0
assert Segment(9.5, 0.5).len(10.0) == 1.0
assert AntiSegment(1.0, 2.0).len(10.0) == 1.0
assert AntiSegment(9.5, 0.5).len(10.0) == 1.0
