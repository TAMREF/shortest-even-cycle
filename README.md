# Shortest even cycle

Implementation of the paper [Shortest Even Cycle Problem is Tractible](https://arxiv.org/pdf/2111.02992.pdf) by A. Björklund et al.

This implementation is **highly inefficient**. Any kind of contribution is welcomed.

## Dependencies
- Python 3.9
- SageMath 9.6
- numpy 1.23

## Supports
- $O(n^3 \cdot 2^{n})$ exponential algorithm, finding the length and an instance of the shortest even cycle
- (Probably) $\tilde{O}(n^7)$ implementation of Björklund et al, finding the length of the shortest even cycle
  - Due to suboptimal determinant oracle, which needs about $\tilde{O}(n^4)$ field operations

## TODOs
- [ ] Providing an instance of the shortest even cycle
- [ ] Employ faster determinant oracle for polynomial matrix

