

import os
from copy import deepcopy, copy
import numpy as np
from pydlr import kernel
from scipy.integrate import quad

import matplotlib.pyplot as plt

seed = 1234
np.random.seed(seed)

getenv = lambda x, default : int(os.getenv(x)) if os.getenv(x) is not None else default

DEBUG = getenv('DEBUG', 0)

class Configuration:
    def __init__(self, t_i, t_f):
        assert len(t_i) == len(t_f)
        if len(t_i) == 0:
            self.t_i, self.t_f = [], []
        else:
            self.t_i, self.t_f = sorted(t_i), sorted(t_f)

    def __len__(self): return len(self.t_i)

    def __add__(self, move):

        if move.type == 'Insert':
            t_i = np.hstack((self.t_i, move.t_i))
            t_f = np.hstack((self.t_f, move.t_f))
            return Configuration(t_i, t_f)

        elif move.type == 'Remove':
            if move.i_idx > 0:
                t_i = np.delete(copy(self.t_i), move.i_idx)
                t_f = np.delete(copy(self.t_f), move.f_idx)
                return Configuration(t_i, t_f)
            else:
                return self
    
    def __str__(self):
        return f"Configuration({self.t_i}, {self.t_f})"


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
            # TODO: check?
            return (idx, (idx+1) % len(self.c))

    def getindex(self, idx):
        i_idx, f_idx = self.indices(idx)
        return Segment(self.c.t_i[i_idx], self.c.t_f[f_idx])

    def __iter__(self): return self

    def __next__(self):
        l = len(self.c)-1
        i_idx, f_idx = self.indices(self._state)
        if i_idx <= l:
            s = Segment(self.c.t_i[i_idx], self.c.t_f[f_idx])
            self._state += 1
            return s
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


class Hybridization:
    def __init__(self, times, values, beta):
        self.times = times
        self.values = values
        self.beta = beta

    def __call__(self, time):
        def call(t):
            s = 1.0
            if t < 0.0:
                s = -1.0
                t += self.beta
            idx = np.searchsorted(self.times, t)
            idx = idx-1 if idx == len(self.times) else idx
            ti, tf = self.times[idx-1], self.times[idx]
            vi, vf = self.values[idx-1], self.values[idx]
            return s * (vi + (t-ti) * (vf-vi) / (tf - ti))

        if hasattr(time, '__iter__'): return list(map(call, time))
        else: return call(time)


class Expansion:

    def __init__(self, beta, h, Delta, moves):
        self.beta = beta
        self.h = h
        self.Delta = Delta
        self.moves = moves
        n = len(self.moves)
        self.move_prop = np.zeros(n)
        self.move_acc =  np.zeros(n)


class Determinant:

    def __init__(self, c, e):

        self.t_f = np.copy(c.t_f)
        self.t_i = np.copy(c.t_i)

        if len(self.t_f) > 0 and self.t_f[0] < self.t_i[0]:
            self.t_f = np.roll(self.t_f, -1)

        self.mat = np.zeros((len(self.t_i), len(self.t_f)), dtype=float)
        for itf, tf in enumerate(self.t_f):
            for iti, ti in enumerate(self.t_i):
                self.mat[itf,iti] = e.Delta(tf-ti)
        self.value = np.linalg.det(self.mat)


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


def trace(c, e):
    if len(c) == 0:
        return 2.0
    if not is_segment_proper(c):
        return 0.0

    value = -1.0 if c.t_f[0] < c.t_i[0] else +1.0
    for s in segments(c):
        value *= np.exp(-e.h * s.len(e.beta))
    return value

def eval(c, e): return trace(c,e) * Determinant(c,e).value

class InsertMove:
    def __init__(self, t_i, t_f, l):
        self.t_i = t_i
        self.t_f = t_f
        self.l = l
        self.type = 'Insert'

class RemovalMove:
    def __init__(self, i_idx, f_idx, l):
        self.i_idx= i_idx
        self.f_idx= f_idx
        self.l = l
        self.type = 'Remove'

def new_insertion_move(c, e):
    t_i = e.beta * np.random.rand()
    t_f = e.beta * np.random.rand()
    return InsertMove(t_i, t_f, 0.0)

def new_segment_insertion_move(c, e):
    t_i = e.beta * np.random.rand()

    if len(c) == 0:
        t_f = e.beta * np.random.rand()
        l = e.beta
    else:
        s = onantisegment(t_i, c)
        l = 0.0 if s is None else Segment(t_i, s.t_f).len(e.beta) 
        t_f = t_i + l * np.random.rand() % e.beta

    assert l >= 0 and l <= e.beta
    return InsertMove(t_i, t_f, l)

def new_antisegment_insertion_move(c, e):
    t_i = e.beta * np.random.rand()

    if len(c) == 0:
        t_f = e.beta * np.random.rand()
        l = e.beta
    else:
        s = onsegment(t_i, c)
        l = 0.0 if s is None else AntiSegment(t_i, s.t_f).len(e.beta) 
        t_f = t_i + l * np.random.rand() % e.beta

    assert l >= 0 and l <= e.beta
    return InsertMove(t_f, t_i, l)


def new_removal_move(c, e):
    l = len(c)
    if l > 0:
        i_idx = np.random.randint(0,l)
        f_idx = np.random.randint(0,l)
        return RemovalMove(i_idx, f_idx, l)
    else:
        return RemovalMove(0, 0, 0.0)


def new_segment_removal_move(c, e):
    if len(c) > 0:
        idx = np.random.randint(0, len(c))
        i_idx, f_idx = segments(c).indices(idx)

        s = Segment(c.t_i[i_idx], c.t_i[(idx+1) % len(c)])
        l = s.len(e.beta)
        return RemovalMove(i_idx, f_idx, l)
    else:
        return RemovalMove(0, 0, 0.0)

def new_antisegment_removal_move(c, e):
    if len(c) > 0:
        idx = np.random.randint(0, len(c))
        f_idx, i_idx = antisegments(c).indices(idx)
        s = AntiSegment(c.t_f[f_idx], c.t_f[(f_idx+1) % len(c)])
        l = s.len(e.beta)
        return RemovalMove(i_idx, f_idx, l)
    else:
        return RemovalMove(0, 0, 0.0)

def propose(move, c_old, e):
    c_new = c_old + move

    if move.type == 'Insert':
        R = move.l * e.beta / len(c_new) * (np.abs(eval(c_new, e) / eval(c_old, e) ))
    elif move.type == 'Remove':
        if move.l == 0: return np.nan
        R = len(c_old) / e.beta / move.l * (np.abs(eval(c_new, e) / eval(c_old, e) ))
    return R

def finalize(move, c):
    c_new = c + move
    return c_new

def metropolis_hastings_update(c, e):
    move_idx = np.random.randint(0, len(moves))
    move = e.moves[move_idx](c,e)
    e.move_prop[move_idx] += 1
    R = propose(move, c, e)
    if R > np.random.rand():
        c = finalize(move, c)
        e.move_acc[move_idx] += 1
    return c


class GreensFunction:

    def __init__(self, beta, data=None, sign=None, N=None):
        self.beta = beta
        if data is not None:
            self.data = data
            self.sign = 0.0
        if N is not None:
            self.data = np.zeros(N)
            self.sign = 0.0

    def __len__(self): return len(self.data)

    def __mul__(self, scalar):
        return GreensFunction(self.beta, 
                              self.data * scalar, 
                              self.sign)

    def __truediv__(self, scalar):
        return GreensFunction(self.beta, 
                              self.data / scalar, 
                              self.sign)

def accumulate(g, time, value):

    if time < 0.0:
        value *= -1
        time += g.beta
    idx = int(np.floor(len(g) * time / g.beta))
    g.data[idx] += value

def sample_greens_function(g, c, e):

    d = Determinant(c, e)
    M = np.linalg.inv(d.mat)

    w = trace(c,e) * d.value

    g.sign += np.sign(w)

    nt = len(c)
    for i in range(nt):
        for j in range(nt):
            accumulate(g, d.t_f[i]-d.t_i[j], M[j,i])


def eval_semi_circular_g_tau(times, t, h, beta):
    I = lambda x : -2 / np.pi / t**2 * kernel(np.array([times])/beta, 
                                              beta*np.array([x]))[0,0]
    g, res = quad(I, -t+h, t+h, weight='alg', wvar=(0.5, 0.5))
    return g

eval_semi_circular_g_tau = np.vectorize(eval_semi_circular_g_tau)

def check(c, d, t):
    print(c)
    print('d.mat = ', d.mat)
    print('d.value = ', d.value)
    print('t = ', t)
    print('is_segment_proper(c) = ', is_segment_proper(c))
	
if __name__ == "__main__":
    beta = 20
    h = 0.0
    t = 1.0
    nt = 200
    times = np.linspace(0, beta, nt)
    g_ref = eval_semi_circular_g_tau(times, t, h, beta)

    Δ = Hybridization(times, -0.25*t**2 * g_ref, beta)

    moves = [ 
                new_segment_insertion_move,
                new_segment_removal_move, 
                new_antisegment_insertion_move,
                new_antisegment_removal_move
            ]

    e = Expansion(beta, h, Δ, moves)

    if DEBUG == 1:

        c = Configuration([1.0], [3.0])
        d = Determinant(c, e)
        t = trace(c,e)
        check(c,d,t)

        c = Configuration([19.0], [1.0])
        d = Determinant(c, e)
        t = trace(c,e)
        check(c,d,t)

        c = Configuration([1.0, 5.0], [3.0, 7.0])
        d = Determinant(c, e)
        t = trace(c,e)
        check(c,d,t)

        c = Configuration([2.0, 19.0], [5.0, 1.0])
        d = Determinant(c, e)
        t = trace(c,e)
        check(c,d,t)
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

        print(c)
        print(is_segment_proper(c))

        assert not is_segment_proper(c)
        assert Segment(1.0, 2.0).len(10.0) == 1.0
        assert Segment(9.5, 0.5).len(10.0) == 1.0
        assert AntiSegment(1.0, 2.0).len(10.0) == 1.0
        assert AntiSegment(9.5, 0.5).len(10.0) == 1.0

    print("Starting CT-HYB QMC")
    epoch_steps = 10
    warmup_epochs = 1000
    sampling_epochs = int(1e4)

    g = GreensFunction(beta=beta, N=nt)
    c = Configuration([], [])

    print(f"Warmup epochs {warmup_epochs} with {epoch_steps} steps.")
    for epoch in range(warmup_epochs):
        for step in range(epoch_steps):
            c = metropolis_hastings_update(c, e)

    print(f"Sampling epochs {sampling_epochs} with {epoch_steps} steps.")
    for epoch in range(sampling_epochs):
        for step in range(epoch_steps):
            c = metropolis_hastings_update(c, e)
        sample_greens_function(g, c, e)

    print('<sign> = ', g.sign )

    dt = beta / len(g)
    g  = g / (-g.sign * beta * dt)
    print('move prob = ', e.move_prop )
    print('move acc  = ', e.move_acc  )

    dt = beta/nt
    t = np.linspace(0.5*dt, beta-0.5*dt, nt)
    plt.figure()
    plt.plot(t, g.data, ".", label = "G")
    plt.plot(times, g_ref, "-", label = "G (ref)")
    times = np.linspace(0, beta, 1001)
    plt.plot(times, Δ(times), '-', label='Delta (ref)')
    plt.plot(Δ.times, Δ.values, '.', label='Delta')
    plt.legend()
    plt.show()


