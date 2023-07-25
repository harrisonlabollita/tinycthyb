import os, sys
import matplotlib.pyplot as plt
import numpy as np

from .solver import Solver

getenv = lambda x, default : int(os.getenv(x)) if os.getenv(x) is not None else default

DEBUG = getenv('DEBUG', 0)

def main():

    S = Solver();
    S.solve()


    dt = S.beta/S.nt
    t = np.linspace(0.5*dt, S.beta-0.5*dt, len(S.g))
    plt.figure()
    plt.plot(t, S.g.data, ".", label = "G")
    plt.plot(S.times, S.g_ref, "-", label = "G (ref)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()





