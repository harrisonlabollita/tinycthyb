import numpy as np
from copy import deepcopy, copy

from segment import Segment, SegmentIterator,segments, onsegment, remove_segment, is_segment_proper
from segment import AntiSegment, AntiSegmentIterator, antisegments, onantisegment, remove_antisegment

from configuration import Configuration

from green import GreensFunction, eval_semi_circular_g_tau

from hybridization import Hybridization


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

        self.t_f = copy(c.t_f)
        self.t_i = copy(c.t_i)

        if len(self.t_f) > 0 and self.t_f[0] < self.t_i[0]:
            self.t_f = np.roll(self.t_f, -1)

        self.mat = np.zeros((len(self.t_i), len(self.t_f)), dtype=float)
        for (i, ti) in enumerate(self.t_i):
            for (j, tf) in enumerate(self.t_f):
                self.mat[j,i] = e.Delta(tf-ti)

        self.value = np.linalg.det(self.mat)

def eval(c, e): return trace(c,e) * Determinant(c,e).value


def trace(c, e):
    if len(c) == 0:
        return 2.0
    if not is_segment_proper(c):
        return 0.0

    value = -1.0 if c.t_f[0] < c.t_i[0] else +1.0
    for s in segments(c):
        value *= np.exp(-e.h * s.len(e.beta))
    return value

# Moves

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

## Insertion Moves
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
        t_f = (t_i + l * np.random.rand()) % e.beta

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
        t_f = (t_i + l * np.random.rand()) % e.beta

    assert l >= 0 and l <= e.beta
    return InsertMove(t_f, t_i, l)


## Removal Moves
def new_removal_move(c, e):
    l = len(c)
    if l >= 0:
        i_idx = np.random.randint(0,l)
        f_idx = np.random.randint(0,l)
        return RemovalMove(i_idx, f_idx, l)
    else:
        return RemovalMove(0, 0, 0.0)


def new_segment_removal_move(c, e):
    if len(c) > 0:
        idx = np.random.randint(0, len(c))
        i_idx, f_idx = segments(c).indices(idx)
        s = Segment(c.t_i[i_idx], c.t_i[(i_idx+1) % len(c)])
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


class Solver:

    def __init__(self, beta=20, h=0.0, t=1.0, nt=200):

        self.beta = beta
        self.h = h
        self.t = t
        self.nt = nt
        self.times = np.linspace(0, self.beta, self.nt)

        self.g_ref = eval_semi_circular_g_tau(self.times, self.t, self.h, self.beta)

        self.Δ = Hybridization(self.times, -0.25*t**2 * self.g_ref, self.beta)
	
        self.moves = [ 
                    new_segment_insertion_move,
                    new_segment_removal_move, 
                    new_antisegment_insertion_move,
                    new_antisegment_removal_move
                ]

        self.e = Expansion(self.beta, self.h, self.Δ, self.moves)


    def metropolis_hastings_update(self):

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

        move_idx = np.random.randint(0, len(self.moves))
        move = self.e.moves[move_idx](self.c, self.e)
        self.e.move_prop[move_idx] += 1
        R = propose(move, self.c, self.e)
        if R > np.random.rand():
            c = finalize(move, self.c)
            self.e.move_acc[move_idx] += 1
            return c
        else:
            return self.c

    def sample_greens_function(self):


        d = Determinant(self.c, self.e)
        M = np.linalg.inv(d.mat)

        w = trace(self.c,self.e) * d.value

        self.g.sign += np.sign(w)

        nt = len(self.c)
        for i in range(nt):
            for j in range(nt):
                self.g.accumulate(d.t_f[i]-d.t_i[j], M[j,i])

    def solve(self):

        print("Starting CT-HYB QMC")
        epoch_steps = 10
        warmup_epochs = 1000
        sampling_epochs = int(1e4)

        self.g = GreensFunction(self.beta, N=self.nt)
        self.c = Configuration([], [])

        print(f"Warmup epochs {warmup_epochs} with {epoch_steps} steps.")
        for epoch in range(warmup_epochs):
            for step in range(epoch_steps):
                self.c = self.metropolis_hastings_update() #self.c, self.e)

        print(f"Sampling epochs {sampling_epochs} with {epoch_steps} steps.")
        for epoch in range(sampling_epochs):
            for step in range(epoch_steps):
                self.c = self.metropolis_hastings_update() #self.c, self.e)
            self.sample_greens_function() #g, c, e)

        print('<sign> = ', self.g.sign/(sampling_epochs))

        dt = self.beta / len(self.g)
        self.g  = self.g / (-self.g.sign * self.beta * dt)
        print('move prob = ', self.e.move_prop/(epoch_steps*(sampling_epochs+warmup_epochs)))
        print('move acc  = ', self.e.move_acc/(epoch_steps*(sampling_epochs+warmup_epochs)))
