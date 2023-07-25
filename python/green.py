import numpy as np
from scipy.integrate import quad
from pydlr import kernel

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
                              data=self.data * scalar, 
                              sign=self.sign)

    def __truediv__(self, scalar):
        return GreensFunction(self.beta, 
                              data=self.data / scalar, 
                              sign=self.sign)



    def accumulate(self, time, value):
        if time < 0.0:
            value *= -1
            time += self.beta

        idx = int(np.floor(len(self) * time / self.beta))
        self.data[idx] += value



def eval_semi_circular_g_tau(tau, t, h, beta):
    I = lambda x : -(2.0/np.pi)*kernel(np.array([tau])/beta, beta*np.array([x]))[0,0]*np.sqrt(1-x**2)
    g, res = quad(I, -1, 1)
    return g

eval_semi_circular_g_tau = np.vectorize(eval_semi_circular_g_tau)
