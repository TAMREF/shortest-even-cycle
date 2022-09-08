from posixpath import split
import numpy as np
import random as rnd
from typing import List

'''
for k <= 4, even walk with length k == even cycle with length k
'''
def capture_girth_leq_4(adj: np.ndarray) -> int:
    n = adj.shape[0]
    for i in range(n):
        adj[i, i] = 0
    
    adj = adj @ adj
    for i in range(n):
        if adj[i, i]:
            return 2
    
    adj = adj @ adj
    for i in range(n):
        if adj[i, i]:
            return 4
    
    return -1


def generate_random(n: int) -> np.ndarray:
    return np.random.randint(2, size=(n, n))

def no_two_cycles(n: int) -> np.ndarray:
    a = np.zeros((n, n)).astype(np.int32)
    for i in range(n):
        for j in range(i+1, n):
            x = rnd.randint(0, 2)
            a[i][j] = x >> 1
            a[j][i] = x & 1
        
    return a


def no_four_cycles(n: int) -> np.ndarray:
    while True:
        a = no_two_cycles(n)
        if capture_girth_leq_4(a) == -1:
            return a


def sample_odd(a: int, b: int):
    st = a // 2
    ed = (b - 1) // 2
    assert st <= ed
    return rnd.randint(st, ed) * 2 + 1

def split_by_odds(n: int) -> List[int]:
    assert n != 4

    odds = []
    while n > 0:
        if n % 2 == 1 and n < 15:
            odds.append(n)
            break
    
        x = sample_odd(5, n-5)
        assert x <= n - 5
        odds.append(x)
        n -= x

    return odds

def bowties(n: int) -> np.ndarray:
    if n % 2 == 0:
        assert n >= 10
    ver = list(range(n))
    rnd.shuffle(ver)
    adj = np.zeros((n, n), dtype=np.int32)

    ties = split_by_odds(n)
    
    # print(ver)
    # print(ties)
    cur = 0
    for tie in ties:
        l1 = sample_odd(3, tie - 2)
        # print(l1, tie - l1)
        anchor = cur + l1 - 1
        for j in range(anchor, cur, -1):
            adj[ver[j], ver[j-1]] = 1
        adj[ver[cur], ver[anchor]] = 1

        for j in range(anchor, cur + tie - 1):
            adj[ver[j], ver[j+1]] = 1

        adj[ver[cur + tie - 1], ver[anchor]] = 1
        cur += tie

    return adj

def even_hamiltonian(n: int):
    ver = list(range(n))
    rnd.shuffle(ver)
    adj = np.zeros((n, n)).astype(np.int64)
    for i in range(n):
        adj[ver[i], ver[(i+1)%n]] = 1
    return adj

if __name__ == '__main__':
    n = 12
    print(no_four_cycles(n))