from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import hypercube as hc

def main():
    nps = np.array([10, 10**2, 10**3, 10**4, 10**5])
    for N in range(2, 7):
        vars = []
        for n_p in nps:
            sample_arr = []
            for _ in range(20):
                sample_arr.append(hc.hc_monte_carlo(N, n_p))
            vars.append(np.var(sample_arr)/n_p)
        y = np.log(np.array(vars))
        x = np.log(nps)
        z = np.polyfit(x, y, 1)
        plt.scatter(x, y)
        plt.plot(x, z[1]+z[0]*x, ls='--', 
                 label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)}\,\,\, N={N}$')
    plt.xlabel(r'$\log(n_p)$')
    plt.ylabel(r'$\log(\tilde{\sigma}^2)$')
    plt.title('Monte Carlo Approximation Sample variance')
    plt.legend()
    plt.show()

    return None

if __name__ == '__main__':
    main()