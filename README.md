# PQM4_SABER
Some fast operation algorithm implementations on pqm4 SABER

# Reference
All codes are modified from [pqm4](https://github.com/mupq/pqm4) and [NIST Round3 Submission SABER](https://csrc.nist.gov/CSRC/media/Projects/post-quantum-cryptography/documents/round-3/submissions/SABER-Round3.zip).
I collaborated this project with [Gting6](https://github.com/gting6). Much of the work is completed with him.

# Works
Saber uses polynomial multiplication in $Z_{8192}[x]/(x^{256}+1)$. Reference codes uses Toom-Cook 4way + Karatsuba trick. In this project we use NTT trick.
Multiplications are used in two places:
* `MatrixVectorMul`: Multiply a matrix with a vector
* `InnerProd`: Compute inner product

Each entry in the matrix and vector is a polynomial in $Z_{8192}[x]/(x^{256}+1)$.

## NTT
$Z_{8192}[x]/(x^{256}+1)$ is an unfriendly ring for NTT, so we change modulo ring to $Z_Q[x]/(x^{256}+1)$, where Q = 25166081.
There is no primitive root of -1 in $Z_Q[x]/(x^{256}+1)$, but there is a root $\xi=1708789$ such that $\xi^{128}=-1$.
Therefore we can perform incomplete NTT with layer <= 7; we provide incomplete NTT with layer 5, 6, 7 in this package.

### Karatsuba
For base multiplication after NTT transform, we use Karatsuba multiplication.

## Reduction
Montgomery Reduction is applied in multiplication after NTT transform.

# Usage
To use this package, clone the repository to local.
```bash=
git clone https://github.com/Ergodica10002/PQM4_SABER.git
cd PQM4_SABER/
```
Build all codes by `make`, and clean all outputs by `make clean`. If some unknown error occurs, make sure doing `make clean` before every `make`.
```bash=
make clean
make
```
The result program is in `test/` folder, one can see the result by
```bash=
./test/test_kex
```
If all go well, cycles for each stage(keypair, encaps, and decaps) are shown.
If there is logical error in the code, error message will occur.
```bash=
----- ERR CCA KEM ------
```

## Parameter

### NTT

To enable/disable NTT algorithm, open `Makefile` and set the parameter `NTT` to TRUE:

```make
# Set this value to TRUE to use NTT multiplication
NTT = TRUE
```

Default is enabled.

### Karatsuba

Choose schoolbook or Karatsuba algorithm to be used in NTT base multiplication by modifying `base_multiplication.c`:

```c
// Define whether to use Karatsuba algorithm
#define Karatsuba
```

### Reduction

To enable/disable Montgomery Reduction algorithm used in NTT algorithm, open `Makefile` and set the parameter `Mont_reduce` to TRUE:

```ma
# Set this value to TRUE to use Montgomery Reduction in NTT
Mont_reduce = TRUE
```



