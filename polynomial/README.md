# Polynomial Library for Competitive Programming

A header-only C++20 library for polynomial arithmetic, designed for competitive programming and discrete mathematics. Built around algebraic abstractions using C++ concepts — operations are defined over any type satisfying `ring` or `field`, including modular arithmetic.

> **Work in progress.** FFT/NTT multiplication and new algebraic structures (matrices, quotient rings) are planned.

---

## Motivation

Polynomials are far more than a mathematical abstraction — they are a computational tool that unlocks elegant solutions to problems that would otherwise seem intractable.

### A surprising application: the Circle Flipping Problem

Consider **n points arranged in a circle**, each labeled **+1** or **−1**. You are allowed to pick any fixed subset of **k positions relative to a starting point** and flip all their signs simultaneously, then repeat this operation starting from any position on the circle. The question: *for a given flip pattern, can you always reach the all-ones configuration from any starting configuration?*

This problem looks combinatorial and hard. It is. But it reduces completely to a single algebraic check.

**The key insight:** encode each configuration as a polynomial in **Z₂[X]/(Xⁿ − 1)** — the ring of polynomials with binary coefficients, reduced modulo Xⁿ − 1. A point labeled −1 at position i contributes Xⁱ; a point labeled +1 contributes 0. Each flip operation corresponds to adding a fixed polynomial **M(X)** (shifted to the appropriate position on the circle).

Reaching the all-ones state from any configuration is possible if and only if:

> **GCD(M(X), Xⁿ − 1) = 1** in Z₂[X]

If this GCD is 1, the extended Euclidean algorithm gives an explicit sequence of moves. If not, some configurations are permanently unreachable.

This result is fully general: **M(X) can be any polynomial**, encoding any fixed flip pattern — not just consecutive blocks. The library lets you verify reachability and compute move sequences for any n, any M(X), in a few lines:

```cpp
using Z2 = modular_element<2>;
using Poly = polynomial<Z2>;

// Circle of size n=8, flip move: 3 consecutive points (1 + X + X^2)
Poly move_poly  = {Z2(1), Z2(1), Z2(1)};
Poly circle_mod = {Z2(1), Z2(0), Z2(0), Z2(0),
                   Z2(0), Z2(0), Z2(0), Z2(0), Z2(1)}; // X^8 - 1 in Z2

// Check if any configuration can be solved
Poly g = gcd(move_poly, circle_mod); // coming soon
bool always_solvable = (g == Poly{Z2(1)});
```

This is the kind of problem this library is built for — where the algebraic structure of polynomials over finite fields is not just convenient, but *necessary*.

> Full derivation and generalisation available in the LaTeX documentation (coming soon).

---

## Requirements

- **C++20 or later** (uses `<concepts>` and `requires` clauses)
- A modern compiler: GCC 10+, Clang 13+, or MSVC 19.28+

Compile with:
```bash
g++ -std=c++20 -O2 your_file.cpp -o your_program
```

---

## Quick Start

```cpp
#include "polynomial.h"

int main() {
    // Polynomials over doubles: 1 + 2X + 3X^2
    polynomial<double> P = {1.0, 2.0, 3.0};
    polynomial<double> Q = {1.0, 1.0};

    auto sum  = P + Q;
    auto prod = P * Q;
    auto quot = P / Q;
    auto rem  = P % Q;

    std::cout << prod << "\n"; // prints in descending degree order

    // Evaluate P at x = 2.0
    std::cout << P(2.0) << "\n"; // 1 + 4 + 12 = 17

    // Modular arithmetic: Z/998244353Z
    using Zp = modular_element<998244353>;
    polynomial<Zp> A = {Zp(1), Zp(2), Zp(3)};
    polynomial<Zp> B = {Zp(1), Zp(1)};
    std::cout << A * B << "\n";
}
```

---

## Algebraic Concepts

The library enforces algebraic structure at compile time via C++20 concepts.

| Concept | Requires | Supports |
|---------|----------|----------|
| `ring<T>` | `+`, `-`, `*`, `T(0)` | All polynomial operations except division |
| `field<T>` | `ring<T>` + `/`, `T(1)` | Full polynomial arithmetic including division |

Built-in types satisfying these concepts:
- `double`, `float` — satisfy `field`
- `int`, `long long` — satisfy `ring` *(integer division truncates; use with care)*
- `modular_element<mod>` — satisfies `field` when `mod` is prime

---

## `modular_element<mod>`

Represents an integer modulo a prime `mod`, with full field arithmetic. Division is implemented via Fermat's little theorem (modular inverse = `a^(mod-2) mod mod`), which requires `mod` to be **prime**.

```cpp
using Zp = modular_element<998244353>; // common CP prime
Zp a = 3, b = 7;
Zp c = a / b; // computes 3 * 7^(mod-2) mod mod
```

Common primes used in CP: `998244353`, `1e9+7`.

---

## Polynomial Operations

Coefficients are stored in **ascending degree order**: `a[i]` is the coefficient of `X^i`.

### Arithmetic operators

| Operator | Description | Requirement |
|----------|-------------|-------------|
| `P + Q` | Addition | `ring` |
| `P - Q` | Subtraction | `ring` |
| `P * Q` | Multiplication (auto-selects algorithm) | `ring` |
| `P / Q` | Quotient of Euclidean division | `field` |
| `P % Q` | Remainder of Euclidean division | `field` |
| `c * P`, `P * c` | Scalar multiplication | `ring` |
| `+=`, `-=`, `*=` | In-place variants | `ring` |

### Multiplication algorithms

`operator*` automatically selects the best algorithm based on degree:

| Degree threshold | Algorithm | Complexity |
|-----------------|-----------|------------|
| max degree < 64 | Brute force | O(n²) |
| max degree ≥ 64 | Karatsuba | O(n^1.585) |

FFT/NTT (O(n log n)) — *coming soon*.

### Other operations

| Expression | Description | Requirement |
|------------|-------------|-------------|
| `--P` | Derivative | `ring` |
| `++P` | Antiderivative (constant = 0) | `field` |
| `P(x)` | Evaluate at point `x` (Horner's method) | `ring` |
| `P == Q`, `P != Q` | Equality (leading zeros ignored) | `ring` |
| `std::cout << P` | Print in descending degree order | `ring` |

> **Note:** `--` and `++` are used as shorthand for derivative and antiderivative. This is a non-standard semantic choice that may change in a future version.

### Utility functions

| Function | Description |
|----------|-------------|
| `degree(P)` | Degree of P |
| `coefficient(i, P)` | Coefficient of X^i (0 if out of range) |
| `trim(P)` | Remove leading zero coefficients |
| `shift_left(k, P)` | Multiply by X^k |
| `shift_right(k, P)` | Divide by X^k (integer floor) |
| `subpolynomial(l, r, P)` | Extract coefficients from index l to r |

---

## Constructors

```cpp
polynomial<T> P;            // zero polynomial
polynomial<T> P(n);         // zero polynomial of degree n
polynomial<T> P = {1,2,3};  // from initializer list: 1 + 2X + 3X^2
```

---

## Roadmap

- [ ] FFT multiplication — O(n log n), floating-point coefficients
- [ ] NTT multiplication — O(n log n), exact modular arithmetic
- [ ] Polynomial composition P(Q(x))
- [ ] Higher-order derivatives, nth antiderivative
- [ ] GCD of polynomials (Euclidean algorithm over fields)
- [ ] Extended GCD — for move sequence computation in the Circle Problem
- [ ] Formal power series mod X^n (inverse, exp, log)
- [ ] Square matrix type satisfying `ring` (for linear recurrences)
- [ ] Full LaTeX reference documentation

---

## Project structure

```
polynomial.h   — full library (header-only, single file)
README.md      — this file
```

---

## License

MIT
