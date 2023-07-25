import numpy as np
from copy import copy

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
            if move.i_idx >= 0 and len(self) > 0:
                t_i = np.delete(copy(self.t_i), move.i_idx)
                t_f = np.delete(copy(self.t_f), move.f_idx)
                return Configuration(t_i, t_f)
            else:
                return self #Configuration(self.t_i, self.t_f)
    
    def __str__(self):
        return f"Configuration({self.t_i}, {self.t_f})"

