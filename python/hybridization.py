import numpy as np

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
            idx = 1 if idx == 0 else idx
            ti, tf = self.times[idx-1], self.times[idx]
            vi, vf = self.values[idx-1], self.values[idx]
            return s * (vi + (t-ti) * (vf-vi) / (tf - ti))

        if hasattr(time, '__iter__'): return list(map(call, time))
        else: return call(time)
