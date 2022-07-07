import numpy as np
import matplotlib.pyplot as plt
import hypercube as hc

def main():
    y1 = np.zeros(5)
    y2 = np.zeros(5)  
    "For N=2"
    nps = np.array([3**2, 10**2, 31**2, 100**2, 330**2])
    for i, n_p in enumerate(nps):
        y1[i] = hc.hc_analytic(2)
        y2[i] = hc.hc_rect(2, n_p)
    y = np.log10(abs(y1-y2))
    x = np.log10(nps)
    plt.scatter(x, y)
    z = np.polyfit(x, y, 1)
    plt.plot(x, z[1]+z[0]*x, ls='--',
             label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)} \,\, N=2$')
    "For N=3"
    nps = np.array([2**3, 5**3, 10**3, 24**3, 50**3])
    for i, n_p in enumerate(nps):
        y1[i] = hc.hc_analytic(3)
        y2[i] = hc.hc_rect(3, n_p)
    y = np.log10(abs(y1-y2))
    x = np.log10(nps)
    plt.scatter(x, y)
    z = np.polyfit(x, y, 1)
    plt.plot(x, z[1]+z[0]*x, ls='--',
             label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)}\,\, N=3$')
    "For N=4"
    nps = np.array([2**4, 3**4, 6**4, 10**4, 18**4])
    for i, n_p in enumerate(nps):
        y1[i] = hc.hc_analytic(4)
        y2[i] = hc.hc_rect(4, n_p)
    y = np.log10(abs(y1-y2))
    x = np.log10(nps)
    plt.scatter(x, y)
    z = np.polyfit(x, y, 1)
    plt.plot(x, z[1]+z[0]*x, ls='--',
             label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)}\,\, N=4$')
    "For N=5"
    nps = np.array([2**5, 3**5, 4**5, 6**5, 10**5])
    for i, n_p in enumerate(nps):
        y1[i] = hc.hc_analytic(5)
        y2[i] = hc.hc_rect(5, n_p)
    y = np.log10(abs(y1-y2))
    x = np.log10(nps)
    plt.scatter(x, y)
    z = np.polyfit(x, y, 1)
    plt.plot(x, z[1]+z[0]*x, ls='--',
             label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)}\,\, N=5$')
    "For N=6"
    nps = np.array([4**6, 2**6, 3**6, 5**6, 7**6])
    for i, n_p in enumerate(nps):
        y1[i] = hc.hc_analytic(6)
        y2[i] = hc.hc_rect(6, n_p)
    y = np.log10(abs(y1-y2))
    x = np.log10(nps)
    plt.scatter(x, y)
    z = np.polyfit(x, y, 1)
    plt.plot(x, z[1]+z[0]*x, ls='--',
             label=f'$f(x)={round(z[0], 2)}x+{round(z[1], 2)}\,\, N=6$')
    
    plt.xlabel(r'$\log(n_p)$')
    plt.xlim(2, 12)
    plt.title('Rectangular Approximation Error')
    plt.ylabel(r'$\log|y_{rect}-y_{analyt}|$')
    plt.legend()
    plt.show()
    return None

if __name__ == '__main__':
    main()