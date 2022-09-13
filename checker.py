from textwrap import fill
import numpy as np
from typing import List
from random import randint
from bjorklund import SELCSolver

from generator import bowties, even_hamiltonian, no_four_cycles, no_two_cycles


def is_even_cycle(adj: np.ndarray, path: List[int]):
    if len(path) % 2 or len(path) == 0:
        return False
    
    if len(set(path)) < len(path):
        return False

    for i in range(len(path)):
        x = path[i]
        y = path[(i + 1) % len(path)]
        if adj[x][y] == 0:
            return False

    return True


'''
Naive algorithm to find shortest even cycle
'''
def shortest_even_cycle(adj: np.ndarray) -> List[int]:
    n = adj.shape[0]
    
    vis = np.zeros((n, n, 1 << n), dtype=bool)
    par = np.full((n, n, 1 << n), fill_value=-1, dtype=np.int32)

    for i in range(n):
        vis[i, i, 1 << i] = True

    shortest = (n + 1, -1, -1, -1)

    mask = (1 << n) - 1
    for bit in range(1, 1 << n):
        
        for nxt in range(n):
            if bit >> nxt & 1:
                continue

            for st in range(n):
                for ed in range(n):
                    if vis[st, ed, bit] and not vis[st, nxt, bit | 1 << nxt] and adj[ed, nxt]:
                        vis[st, nxt, bit | 1 << nxt] = 1
                        par[st, nxt, bit | 1 << nxt] = ed
        
        cnt: int = bin(bit).count("1")
        if cnt % 2 == 0:
            for st in range(n):
                for ed in range(n):
                    if vis[st, ed, bit] and adj[ed, st]:
                        shortest = min(shortest, (cnt, st, ed, bit))
        
    if shortest[0] == n + 1:
        return []
    
    path = []
    
    st, ed, bit = shortest[1:]
    for _ in range(shortest[0]):
        path.append(ed)
        bit ^= (1 << ed)
        ed = par[st, ed, bit | 1 << ed]
    
    path.reverse()
    # print(path)
    assert len(path) % 2 == 0
    assert is_even_cycle(adj, path)
    return path

if __name__ == '__main__':
    maxn = 13
    TC = [
        no_four_cycles(12),
        bowties(12),
        bowties(10),
        even_hamiltonian(12),
        no_four_cycles(10),
        no_two_cycles(12),
        no_two_cycles(6),
        no_two_cycles(12)
    ] + [no_two_cycles(12)] * 10 + [no_four_cycles(13)] * 10

    solver = SELCSolver(maxn, iter=1)

    from time import perf_counter as pf
    from tqdm import tqdm

    pbar = tqdm(TC)
    for adj in pbar:
        t = pf()
        path = shortest_even_cycle(adj)
        t2= pf()
        # print(f"naive took {pf() - t:.3f} sec")
        path2 = solver.get_selc(adj)
        t3 = pf()
        # print(f"bjork took {pf() - t:.3f} sec")

        if len(path2) == 0:
            assert len(path) == 0
        else:
            # print(path, path2)
            assert len(path) == len(path2)
            assert is_even_cycle(adj, path2)
                
        pbar.set_description(f"naive {t2 - t:.3f}, bjork {t3 - t2:.3f}, len {len(path)}")
        # print(f"bjor: {path2}")