import numpy as np
from tqdm import tqdm
from random import random


def hc_rect(N: int, n_p: int):
    num = round(n_p**(1/N))
    count = 0
    for point in range(n_p):
        coord = np.zeros(N)
        for k in range(N):
            coord[k] = int(point/(num**k)) % num
        r = 0
        for x in coord:
            r += (x-(num-1)/2)**2
        count += (r <= ((num-1)/2)**2)
    return(count/n_p)*2**N

def hc_monte_carlo(N: int, n_p: int):
    count = 0
    for _ in range(n_p):
        coord_sqrd = [(random()*2 - 1)**2 for _ in range(N)]
        r = sum(coord_sqrd)
        count += (r<=1)
    return (count/n_p)*2**N

def hc_analytic(N:int):
    if N == 0:
        val = 1
    elif N == 1:
        val = 2
    else:
        val = 2*np.pi*hc_analytic(N-2)/N
    return val
    

# print(hc_monte_carlo(3, int(100e3)))
# print(hc_analytic(3))