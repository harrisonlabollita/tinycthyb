class Segment:
    def __new__(cls, t_i, t_f):
        obj = super().__new__(cls)
        obj.t_i, obj.t_f = t_i, t_f;
        return obj

    def __str__(self):
        return f"Segment({self.t_i}, {self.t_f})"

    def len(self, beta):
        if self.t_i < self.t_f:
            return self.t_f - self.t_i
        else:
            return beta - self.t_i + self.t_f

    def onsegment(self, t):
        if self.t_i < self.t_f:
            return (t > self.t_i) and (t <  self.t_f)
        else:
            return (t < self.t_f) or (t >  self.t_i)

class SegmentIterator:
    def __init__(self, configuration):
        self.c = configuration
        self._state = 0

    def __len__(self): return len(self.c)

    def indices(self, idx):
        if self.c.t_f[0] > self.c.t_i[0]:
            return (idx, idx)
        else:
            return (idx, (idx+1) % len(self.c))

    def getindex(self, idx):
        i_idx, f_idx = self.indices(idx)
        return Segment(self.c.t_i[i_idx], self.c.t_f[f_idx])

    def __iter__(self): return self

    def __next__(self):
        l = len(self.c)
        i_idx, f_idx = self.indices(self._state)
        if i_idx < l:
            self._state += 1
            return Segment(self.c.t_i[i_idx], self.c.t_f[f_idx])
        else:
            raise StopIteration

def segments(c): return SegmentIterator(c)

def onsegment(t, c):
    for s in segments(c):
        if s.onsegment(t):
            return s
    return None


def remove_segment(c, segment_idx):
    assert 0 <= segment_idx < len(c)
    i_idx, f_idx = segments(c).indices(segment_idx)
    c.t_i = np.delete(c.t_i, i_idx)
    c.t_f = np.delete(c.t_f, f_idx)

class AntiSegment:
    def __new__(cls, t_i, t_f):
        obj = super().__new__(cls)
        obj.t_i, obj.t_f = t_i, t_f;
        return obj

    def len(self, beta):
        if self.t_i < self.t_f:
            return self.t_f - self.t_i
        else:
            return beta - self.t_i + self.t_f

    def __str__(self):
        return f"AntiSegment({self.t_i}, {self.t_f})"

    def onantisegment(self, t):
        if self.t_i < self.t_f:
            return (t > self.t_i) and (t <  self.t_f)
        else:
            return (t < self.t_f) or (t >  self.t_i)

class AntiSegmentIterator:
    def __init__(self, configuration):
        self.c = configuration
        self._state = 0

    def __len__(self): return len(self.c)

    def indices(self, idx):
        if self.c.t_f[0] < self.c.t_i[0]:
            return (idx, idx)
        else:
            return (idx, (idx+1) % len(self.c))

    def getindex(self, idx):
        f_idx, i_idx = self.indices(idx)
        return AntiSegment(self.c.t_f[f_idx], self.c.t_i[i_idx])

    def __iter__(self): return self

    def __next__(self):
        l = len(self.c)
        f_idx, i_idx = self.indices(self._state)
        if f_idx < l:
            s = AntiSegment(self.c.t_f[f_idx], self.c.t_i[i_idx])
            self._state += 1
            return s
        else:
            raise StopIteration

def antisegments(c): return AntiSegmentIterator(c)

def onantisegment(t, c):
    for s in antisegments(c):
        if s.onantisegment(t):
            return s
    return None

def remove_antisegment(c, segment_idx):
    assert 0 <= segment_idx < len(c)
    f_idx, i_idx = antisegments(c).indices(segment_idx)
    c.t_i = np.delete(c.t_i, i_idx)
    c.t_f = np.delete(c.t_f, f_idx)

def is_segment_proper(c):

    if c.t_i[0] < c.t_f[0]:
        for idx in range(1, len(c)):
            if not (c.t_f[idx-1] < c.t_i[idx]) or not (c.t_i[idx] < c.t_f[idx]):
                return False
    else:
        for idx in range(1, len(c)):
            if not (c.t_i[idx-1] < c.t_f[idx]) or not (c.t_f[idx] < c.t_i[idx]):
                return False
    return True
