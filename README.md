# 🔐 Discrete Logarithm Solver

Provides performant Python implementations of several **Discrete Logarithm Problem (DLP)** algorithms for prime-order multiplicative subgroups of prime fields ($F_p$).


Given the following parameters:
- $p$ the prime modulus. Specifies the field $F_p$ we are working in
- $g$ a generator element of the multiplicative subgroup $\langle g \rangle$ of $F_p$
- $h$ the target element. We must have $h \in \langle g \rangle$
- $n$ the subgroup order $|\langle g \rangle|$

Output:
$x \in \{0,1,\ldots,n-1\} \text{ s.t. } g^x \equiv h \bmod p$

