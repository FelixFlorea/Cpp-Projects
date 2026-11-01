# Polynomial Library for C++20

A **header-only C++20 library** for polynomial arithmetic, modular arithmetic and square matrices, designed around **algebraic structures**: rings, commutative rings and fields.

The core idea is that **every operation asks for exactly the algebraic structure it needs**. Polynomial division requires coefficients from a field, so `polynomial<int> / polynomial<int>` does **not compile**: the mistake is caught by the compiler instead of producing a wrong result at run time.

```cpp
#include "polynomial_refined.h"
#include <iostream>

int main() {
    polynomial<int> P = {1, 2, 3};      // 1 + 2X + 3X^2
    polynomial<int> Q = {4, -1};        // 4 - X
    std::cout << P * Q << "\n";         // -3X^3+10X^2+7X+4
}
```

---

## Table of contents

1. [Features](#features)
2. [Getting started](#getting-started)
3. [Quick tour](#quick-tour)
4. [Design philosophy](#design-philosophy)
5. [Algebraic concepts](#algebraic-concepts)
6. [`modular_element<mod>`](#modular_elementmod)
7. [`matrix<T, n>`](#matrixt-n)
8. [`polynomial<T>`](#polynomialt)
9. [Error handling](#error-handling)
10. [Compile-time errors by design](#compile-time-errors-by-design)
11. [Complexity](#complexity)
12. [Known limitations](#known-limitations)
13. [Roadmap](#roadmap)
14. [FAQ](#faq)
15. [Author](#author)

---

## Features

| Component | What it provides |
| --- | --- |
| `ring`, `commutative_ring`, `field` | C++20 concepts describing algebraic structures |
| `is_commutative<T>` | an explicit list of types whose multiplication is commutative |
| `modular_element<mod>` | integers modulo `mod` (assumed prime), with Fermat inversion and binary exponentiation |
| `matrix<T, n>` | n×n matrices over any ring, with identity, transpose, trace and fast exponentiation |
| `polynomial<T>` | polynomials over any ring: addition, subtraction, scalar and polynomial multiplication (schoolbook and Karatsuba), division with remainder, evaluation, derivative, antiderivative, input/output |

The library depends **only on the standard library**, and only on five headers:

| Header | Why it is included |
| --- | --- |
| `<concepts>` | `std::regular` and `std::same_as`, the building blocks of the algebraic concepts |
| `<ostream>` | `operator<<` for every type (the library never uses `std::cout` itself, so `<iostream>` is not needed) |
| `<istream>` | `operator>>` for every type |
| `<vector>` | storage for polynomial coefficients, whose number changes constantly |
| `<initializer_list>` | brace initialization: `polynomial<int> P = {1, 2, 3};` |

---

## Getting started

### Requirements

- A compiler with **C++20** support, including concepts. Tested with GCC 13.
- Nothing else: no build system, no external dependencies.

### Usage

Copy `polynomial_refined.h` next to your code and include it:

```cpp
#include "polynomial_refined.h"
```

Compile with C++20 enabled:

```bash
g++ -std=c++20 main.cpp -o main
./main
```

---

## Quick tour

### Integer polynomials

```cpp
polynomial<int> P = {1, 2, 3};      // coefficients from X^0 upwards: 1 + 2X + 3X^2
polynomial<int> Q = {4, -1};        // 4 - X

P + Q;              // 3X^2+X+5
P * Q;              // -3X^3+10X^2+7X+4
2 * P;              // 6X^2+4X+2
P(2);               // 17       (Horner evaluation)
derivative(P);      // 6X+2
degree(P);          // 2
P[1];               // 2        (coefficient of X)
P[10];              // 0        (above the degree, every coefficient is 0)
P - P;              // 0        (normalized: degree 0, no trailing zeros)
```

### Modular arithmetic

```cpp
using Z7 = modular_element<7>;

Z7 a(3), b(5);
a + b;              // 1
a - b;              // 5
a * b;              // 1
a / b;              // 2
-a;                 // 4
a.inverse();        // 5        (3 · 5 = 15 ≡ 1 mod 7)
Z7::pow_log(a, 6);  // 1        (Fermat's little theorem)
Z7(-10);            // 4        (negative values are reduced into [0, 6])
```

### Polynomials over a field: division with remainder

```cpp
polynomial<Z7> A = {Z7(1), Z7(0), Z7(2), Z7(3)};   // 1 + 2X^2 + 3X^3
polynomial<Z7> B = {Z7(2), Z7(1)};                 // 2 + X

A / B;              // 3X^2+3X+1
A % B;              // 6
A(Z7(-2));          // 6        (remainder theorem: A mod (X + 2) = A(-2))
antiderivative(A);  // 6X^4+3X^3+X
```

### Real polynomials

```cpp
polynomial<double> D1 = {-1.0, 0.0, 1.0};   // X^2 - 1
polynomial<double> D2 = {-1.0, 1.0};        // X - 1

D1 / D2;            // X+1
D1 % D2;            // 0
```

### Matrices: Fibonacci in logarithmic time

```cpp
//  [1 1]^n   [F(n+1)  F(n)  ]
//  [1 0]   = [F(n)    F(n-1)]
matrix<long long, 2> F;             // starts as the zero matrix
F.m[0][0] = 1; F.m[0][1] = 1;
F.m[1][0] = 1;

matrix<long long, 2>::pow_log(F, 10).m[0][1];   // 55
```

With `long long`, `n` can go up to 91. With `matrix<modular_element<p>, 2>`, `n` can be arbitrarily large (for example 10¹⁸) and the result is computed modulo `p`.

---

## Design philosophy

1. **Algebra first.** The code states the same rules as the mathematics. Operations are constrained by concepts, so an operation without mathematical meaning is rejected at compile time.
2. **Invariants protected by the type.** Each class that has a rule to keep makes its data private and enforces the rule in its constructors:
   - `modular_element`: the stored value is always in `[0, mod − 1]`;
   - `polynomial`: the coefficient vector is never empty and never ends in a zero (except for the zero polynomial).

   `matrix` has no such rule (any n×n values form a valid matrix), so it is a `struct` with public data.
3. **No hidden conversions.** Converting constructors are `explicit`: an `int` never silently becomes a `modular_element`, and a `std::vector` never silently becomes a `polynomial`.
4. **Minimal dependencies.** Only the five standard headers listed above.
5. **Zero runtime cost for the checks.** Concepts are evaluated entirely at compile time.

---

## Algebraic concepts

```
ring  →  commutative_ring  →  field
```

Each level contains the previous one and adds one requirement.

### `ring<T>`

```cpp
template <typename T>
concept ring =
    std::regular<T> and
    requires(T a, T b) {
        {-a}     -> std::same_as<T>;
        {a + b}  -> std::same_as<T>;
        {a - b}  -> std::same_as<T>;
        {a * b}  -> std::same_as<T>;
        {a += b} -> std::same_as<T&>;
        {a -= b} -> std::same_as<T&>;
        {a *= b} -> std::same_as<T&>;
        {T(0)}   -> std::same_as<T>;
        {T(1)}   -> std::same_as<T>;
    };
```

| Requirement | Meaning |
| --- | --- |
| `std::regular<T>` | the type behaves like a value: default-constructible, copyable, movable, equality-comparable |
| `-a` | every element has an additive inverse |
| `a + b`, `a - b`, `a * b` | closed operations that return a `T` |
| `+=`, `-=`, `*=` returning `T&` | in-place operations that return the modified object itself |
| `T(0)`, `T(1)` | the additive and multiplicative identities, built from the integers 0 and 1 |

Because `T(1)` is required, the concept describes a **ring with unity**.

The parameters `a` and `b` in `requires(T a, T b)` never exist at run time; they are only names used to describe which expressions must compile.

### `commutative_ring<T>`

```cpp
template <typename T>
concept commutative_ring = ring<T> and is_commutative<T>;
```

Commutativity cannot be checked by the compiler: for matrices, both `a * b` and `b * a` compile, but they give different results. Commutative types are therefore **declared explicitly**:

```cpp
template <typename T>
inline constexpr bool is_commutative = false;     // safe default

template <> inline constexpr bool is_commutative<int>       = true;
template <> inline constexpr bool is_commutative<long long> = true;
template <> inline constexpr bool is_commutative<float>     = true;
template <> inline constexpr bool is_commutative<double>    = true;

template <unsigned int mod>
inline constexpr bool is_commutative<modular_element<mod>> = true;
```

- The list lives **outside** the classes, so built-in types such as `double` can be added to it.
- The default is `false`. Forgetting to register a commutative type only costs speed (the schoolbook multiplication is used instead of Karatsuba); the result stays correct.

### `field<T>`

```cpp
template <typename T>
concept field = commutative_ring<T> and
    requires(T a, T b) {
        {inverse(a)} -> std::same_as<T>;
        {a / b}      -> std::same_as<T>;
        {a /= b}     -> std::same_as<T&>;
    };
```

The free function `inverse(a)` is what separates a real field from a type that merely has a `/` operator. `int` has `a / b` (integer division) but no `inverse`, so `int` is correctly **not** a field.

All fields are assumed to be commutative (quaternion-like structures are out of scope).

`inverse` overloads for `double` and `float` are declared **before** the concepts, because for built-in types the compiler only looks for functions declared earlier in the file. The overload for `modular_element` can appear later, because for user-defined types the compiler also searches next to the type (argument-dependent lookup).

### Which types satisfy which concept

| Type | `ring` | `commutative_ring` | `field` |
| --- | :---: | :---: | :---: |
| `int`, `long long` | yes | yes | no |
| `float`, `double` | yes | yes | yes |
| `modular_element<p>` | yes | yes | yes |
| `matrix<T, n>` | no | no | no |
| `polynomial<T>` | no | no | no |

### What each operation requires

| Operation | Requires |
| --- | --- |
| `polynomial<T>`, `matrix<T, n>` | `ring<T>` |
| polynomial `+`, `-`, `*`, `derivative`, evaluation | `ring<T>` |
| Karatsuba multiplication | `commutative_ring<T>` |
| polynomial `/`, `%`, `antiderivative` | `field<T>` |

---

## `modular_element<mod>`

An element of ℤ<sub>mod</sub>. The modulus is a **template parameter**, so it is part of the type:

- `modular_element<7>` and `modular_element<11>` are different types and cannot be mixed by accident;
- the compiler knows the modulus in advance and can optimize the `%` operation.

`mod` is assumed to be **prime**, which makes ℤ<sub>mod</sub> a field.

### Invariant

The stored value `v` is always in `[0, mod − 1]`. It is private, so it cannot be set to an invalid value from outside.

### Constructors

| Constructor | Access | Behaviour |
| --- | --- | --- |
| `modular_element()` | public | the element 0 |
| `explicit modular_element(long long x)` | public | the **safe** path: reduces any value, including negative ones, into `[0, mod − 1]` |
| `modular_element(int val, direct_tag)` | private | the **fast** path: stores a value that is already reduced, skipping `%` |

- In C++, `-10 % 7` is `-3`, not `4`, so the public constructor adds `mod` when the remainder is negative.
- The private constructor is selected with an empty tag type, `direct_tag`, because both constructors take a number. It is used only inside the operators, where the result is known to be in range; external code cannot bypass the reduction.
- Data coming from outside (the public constructor, `operator>>`) always takes the safe path.

### Operations

| Operation | Implementation |
| --- | --- |
| `a + b` | one addition and at most one subtraction of `mod` (no `%`) |
| `a - b` | one subtraction and at most one addition of `mod` (no `%`) |
| `a * b` | product computed in 64 bits (`1LL * ...`), then `% mod` |
| `-a` | `mod - v`, with 0 handled separately so the result stays in range |
| `a / b` | `a * b.inverse()` |
| `+=`, `-=`, `*=`, `/=` | modify the object in place |
| `==` | compares values (`!=` is generated automatically in C++20) |
| `val()` | read-only access to the stored value |
| `<<`, `>>` | print the value; read a `long long` and reduce it |

All binary operators are `friend` functions: they need access to the private value and to the fast constructor, and they treat both operands symmetrically.

### Exponentiation and inverse

```cpp
static modular_element pow_log(modular_element base, long long exp);
modular_element inverse() const;
```

- `pow_log` uses **binary exponentiation**: the exponent is halved at every step, so the cost is O(log exp), about 30 steps for an exponent around 10⁹. A negative exponent uses the inverse: a⁻ᵉ = (a⁻¹)ᵉ.
- `pow_log` is `static` because it does not work on a particular element: the base is passed as a parameter. It is called as `Z7::pow_log(a, 5)`.
- `inverse()` uses **Fermat's little theorem**: a⁻¹ = a<sup>mod − 2</sup>. It is `const` because it does not modify the element (this also allows it to be called on constant objects, as in `operator/`).
- `inverse()` **throws** for 0, which has no inverse.
- A free function `inverse(a)` is also provided, as required by the `field` concept.

---

## `matrix<T, n>`

```cpp
template <typename T, unsigned int n>
requires ring<T>
struct matrix {
    T m[n][n];
    ...
};
```

A square n×n matrix whose elements come from any ring. The size is a template parameter, so matrices of different sizes are different types. `n = 0` does not compile.

Elements are stored in a fixed-size array (no dynamic allocation) and accessed directly as `M.m[i][j]`.

### Members and functions

| Member / function | Description |
| --- | --- |
| `matrix()` | the zero matrix (every element is `T(0)`, so built-in types never contain garbage) |
| `static identity()` | the identity matrix |
| `-M`, `M1 + M2`, `M1 - M2`, `M1 * M2` | matrix operations |
| `+=`, `-=`, `*=` | in-place versions |
| `==` | element-wise comparison |
| `transpose(M)` | the transpose |
| `trace(M)` | the sum of the diagonal |
| `static pow_log(base, exp)` | binary exponentiation, starting from the identity |
| `<<` | prints the matrix row by row, as `[ a b ]` |
| `>>` | reads the n² elements row by row |

### Implementation notes

- **Loop order `i, k, j`** in multiplication: matrices are stored row by row, so with `j` in the innermost loop the memory is read sequentially, which is cache-friendly.
- **Factor order is preserved** (`m1[i][k] * m2[k][j]`), so multiplication is also correct for non-commutative element types.
- **`*=` is not in place** (`m1 = m1 * m2`): every new element needs several old ones, so writing in place would overwrite values that are still needed.
- **`>>` is all-or-nothing**: input is read into a temporary matrix, and the target is only changed if every element was read successfully.
- `transpose` and `trace` are `friend` functions not for access (the data is public) but to define free functions inside the class template without repeating its template header, and to call them in mathematical style: `transpose(M)`.

### Initialization

Because `matrix` has a user-written constructor, it cannot be filled with braces:

```cpp
matrix<long long, 2> M = {{{1, 1}, {1, 0}}};    // does not compile
```

Create it and set the elements instead:

```cpp
matrix<long long, 2> M;
M.m[0][0] = 1; M.m[0][1] = 1;
M.m[1][0] = 1;
```

A brace-list constructor is planned (see the roadmap).

---

## `polynomial<T>`

```cpp
template <typename T>
requires ring<T>
class polynomial;
```

Coefficients are stored in **increasing order**: index `i` holds the coefficient of Xⁱ. This makes the formulas natural: in a product, Xⁱ · Xʲ lands at index `i + j`.

### Invariant: always normalized

1. the coefficient vector is **never empty**; the zero polynomial is `{0}`;
2. the last coefficient is **never zero**, except for the zero polynomial.

The private function `normalize()` is called by every constructor and by `set_coef`. Consequences:

- `degree()` is always simply `size − 1`;
- `operator==` compares coefficient vectors directly; without the invariant, `{1, 2}` and `{1, 2, 0}` would be the same polynomial but different vectors.

By convention, the degree of the zero polynomial is 0.

### Constructors

| Constructor | `explicit` | Use |
| --- | :---: | --- |
| `polynomial()` | — | the zero polynomial |
| `polynomial(std::initializer_list<T>)` | no | `polynomial<int> P = {1, 2, 3};` |
| `polynomial(const std::vector<T>&)` | yes | conversion from a vector, used internally |

- Thanks to the initializer-list constructor, the type can be deduced: `polynomial P = {1, 2, 3};` is a `polynomial<int>`.
- For `modular_element` coefficients, each value must be written explicitly, because that constructor is `explicit`:

  ```cpp
  polynomial<Z7> A = {Z7(1), Z7(2)};      // correct
  polynomial<Z7> A = {1, 2};              // does not compile
  ```

### Member functions

| Member | Description |
| --- | --- |
| `degree()` | the degree, as an `int` |
| `operator[](int i)` | the coefficient of Xⁱ; returns `T(0)` above the degree; throws for a negative index |
| `set_coef(int i, T val)` | sets a coefficient, extending the polynomial if needed, then renormalizes |
| `operator()(T x)` | evaluation with **Horner's scheme**: P(x) = a₀ + x(a₁ + x(a₂ + …)) |
| `get_coefficients()` | `const` reference to the coefficient vector |

`operator[]` returns a **copy**, so coefficients can only be read through it. Writing `P[2] = 0` could break the normalization invariant; `set_coef` exists for writing and renormalizes.

### Free functions and operators

| Function / operator | Description |
| --- | --- |
| `degree(P)` | the degree, in mathematical notation |
| `P1 + P2`, `P1 - P2` | coefficient by coefficient, then normalized |
| `c * P`, `P * c` | multiplication by a scalar (factor order preserved) |
| `P1 * P2` | polynomial multiplication, with automatic algorithm choice |
| `P1 / P2`, `P1 % P2` | quotient and remainder (requires `field<T>`) |
| `+=`, `-=`, `*=` | compound assignment |
| `==` | equality |
| `derivative(P)` | the derivative |
| `antiderivative(P)` | the antiderivative with integration constant 0 (requires `field<T>`) |
| `<<` | prints from the highest degree down, skipping zero terms and coefficients equal to 1 |
| `>>` | reads the degree, then the coefficients from X⁰ upwards |

Helper functions used by Karatsuba are also available: `shift_left(k, P)` (multiply by Xᵏ), `shift_right(k, P)` (divide by Xᵏ, dropping low terms), `subpolynomial(start, stop, P)`, `add`, `subtract`, `scalar_multiply`, `multiply_brute`, `multiply_karatsuba`, `divide_quotient`, `divide_remainder`.

### Multiplication: algorithm selection by concepts

```cpp
template <typename T> requires ring<T> && (!commutative_ring<T>)
polynomial<T> operator*(const polynomial<T>&, const polynomial<T>&);   // schoolbook

template <typename T> requires commutative_ring<T>
polynomial<T> operator*(const polynomial<T>&, const polynomial<T>&);   // below degree 64: schoolbook
                                                                       // from degree 64: Karatsuba
```

The two constraints exclude each other, so the compiler always picks exactly one overload.

**Schoolbook multiplication** multiplies every coefficient with every other one: O(n · m).

**Karatsuba multiplication:**

1. base cases: a constant polynomial reduces to scalar multiplication; below degree 32 on both sides, the schoolbook method is used;
2. both polynomials are split **at the same point** `m` into a low part and a high part;
3. three recursive products: `A = high · high`, `B = low · low`, `C = (low + high) · (low + high)`;
4. the middle term is `C − A − B`;
5. the result is `A · X²ᵐ + (C − A − B) · Xᵐ + B`.

Three half-size products instead of four give roughly **O(n^1.58)** instead of O(n²).

On degree 1, with (a + bx)(c + dx): the three products are `ac`, `bd` and `(a + b)(c + d)`; the middle coefficient is `(a + b)(c + d) − ac − bd = ad + bc`.

Karatsuba is restricted to commutative rings because one of its base cases (`scalar_multiply(P2[0], P1)`) swaps the order of the factors.

### Division with remainder

`divide_quotient` and `divide_remainder` implement **long division**: at every step, the leading coefficient of the remainder is divided by the leading coefficient of the divisor, and the corresponding multiple of the divisor is subtracted. Cost: O(n · m).

- Division by the zero polynomial throws.
- The divisor's leading coefficient is never zero (normalization), so there is never a division by zero inside the loop.
- After each step, the top coefficient is removed **without checking** that it is zero: it is zero mathematically, but with `double` a value like 1e-16 may remain from rounding.

### Derivative and antiderivative

- `derivative` only needs a ring: the coefficient of Xⁱ is multiplied by `T(i)`. In ℤ<sub>p</sub>, `T(p) = 0`, so the derivative of Xᵖ is 0, which is mathematically correct in characteristic p.
- `antiderivative` needs a field: the coefficient of Xⁱ is divided by `T(i + 1)`. Before dividing, the function checks whether `T(i + 1)` became zero (in ℤ<sub>p</sub>, when `i + 1` is a multiple of `p`) and throws. Zero coefficients are skipped, so the error only appears when there is actually a term to integrate. Example: in ℤ₇, the antiderivative of X⁶ would require dividing by 7 ≡ 0.
- The function is named `antiderivative` rather than `integral` to avoid a name clash with the C++20 concept `std::integral` when `using namespace std;` is used.

### Input format

```text
2  1 2 3
```

is read as degree 2 with coefficients 1, 2, 3, that is 1 + 2X + 3X². A negative degree marks the stream as failed.

---

## Error handling

Runtime errors are reported by throwing a string literal (`const char*`).

| Situation | Function | Message |
| --- | --- | --- |
| inverse of 0 | `modular_element::inverse`, `/`, `/=`, `pow_log` with a negative exponent | `"0 has no inverse"` |
| negative index | `polynomial::operator[]`, `set_coef` | `"negative index"` |
| division by the zero polynomial | `/`, `%` on polynomials | `"Division by zero polynomial"` |
| antiderivative in characteristic p | `antiderivative` | `"Integral undefined in this characteristic"` |

The code compiles normally; the error only appears at run time, when the faulty line is executed. If the exception is not caught, the program stops:

```cpp
try {
    Z7(0).inverse();
} catch (const char* message) {
    std::cout << message;       // 0 has no inverse
}
```

---

## Compile-time errors by design

These lines are rejected by the compiler, which is the intended behaviour:

```cpp
polynomial<int> R = P / Q;                          // int is not a field
modular_element<7> c = 3;                           // explicit constructor
auto d = modular_element<7>(1) + modular_element<11>(1);   // different types
polynomial<modular_element<7>> A = {1, 2};          // ints are not converted implicitly
matrix<long long, 2> M = {{{1, 1}, {1, 0}}};        // no brace-list constructor yet
```

---

## Complexity

n and m are polynomial degrees; k is the matrix size.

| Operation | Complexity |
| --- | --- |
| addition, subtraction, scalar multiplication | O(n) |
| derivative, antiderivative | O(n) |
| evaluation (Horner) | O(n) |
| schoolbook multiplication | O(n · m) |
| Karatsuba multiplication | about O(n^1.58) |
| division with remainder | O(n · m) |
| `modular_element::pow_log` | O(log exp) |
| `modular_element::inverse` | O(log mod) |
| matrix multiplication | O(k³) |
| `matrix::pow_log` | O(k³ · log exp) |

---

## Known limitations

**Correctness**

- The modulus of `modular_element` is assumed to be prime but is not checked. For a composite modulus, addition and multiplication are correct, but `inverse` returns wrong results, and the type still satisfies `field`.
- The modulus must be at most about 2³⁰, so that the sum of two values fits in an `int`.
- `matrix::pow_log` with a negative exponent returns the identity matrix.
- `int` and `long long` coefficients can overflow silently.
- With `double`, rounding can leave tiny coefficients (around 1e-16) instead of exact zeros, which affects the degree.

**Design**

- `matrix<T, n>` does not satisfy `ring` (no constructor from a number), so polynomials with matrix coefficients are not possible yet, and the non-commutative branch of `operator*` is currently unused.
- `polynomial<T>` does not satisfy `ring` either (no unary minus, no construction from a number), so multivariate polynomials cannot be built by nesting.
- `std::complex<double>` satisfies `ring` but is not registered as commutative and has no `inverse`, so it gets neither Karatsuba nor division.
- Quotient and remainder are computed by two separate functions, so computing both does the work twice.
- Errors are thrown as `const char*` rather than `std::exception` types, so an uncaught error does not show its message.
- Everything lives in the global namespace.
- The modulus is fixed at compile time.

**Performance**

- The Karatsuba threshold (64) is not tuned: on one test machine, Karatsuba only became faster than the schoolbook method from about degree 256–512. Two different thresholds are used (64 and 32).
- Karatsuba allocates many temporary vectors at each recursion level.
- Compound assignments (`+=`, `-=`, `*=`) create a new polynomial instead of working in place.

**Output**

- Negative coefficients are printed after a `+` (for example `X+-1`), and −1 is printed as `-1X`.

---

## Roadmap

### Short term: fixes and structure

- compile-time check that the modulus is prime and small enough (`static_assert` with a `constexpr` primality test);
- a dedicated `namespace`, with helpers in a `detail` namespace;
- exceptions from `<stdexcept>` (`std::domain_error`, `std::out_of_range`);
- `divmod(A, B)` returning quotient and remainder in one pass;
- a constructor `explicit matrix(T x)` (x on the diagonal), so that `matrix` satisfies `ring`, plus brace-list initialization;
- `modular_element` storing a templated integral type instead of `int`;
- Karatsuba threshold tuned by benchmarks, and factor order preserved in its base case;
- correct printing of negative coefficients;
- more types in `is_commutative` (`unsigned`, `short`, `long`).

### Planned layout

```text
polynomial_lib/
├── core/
│   ├── concepts.hpp              ring, commutative_ring, field
│   ├── modular_element.hpp       modular_element<mod>
│   ├── matrix.hpp                matrix<T, n>
│   └── polynomial.hpp            polynomial<T>, operators, Karatsuba, GCD, derivative
├── transforms/
│   ├── fft.hpp                   Fast Fourier Transform on complex numbers
│   └── ntt.hpp                   Number Theoretic Transform on modular_element<mod>
├── series/
│   ├── newton_inversion.hpp      inverse of a power series mod X^N in O(N log N)
│   ├── formal_series.hpp         logarithm, exponential, real/complex powers
│   └── taylor_laurent.hpp        Taylor and Laurent expansions
├── roots/
│   ├── durand_kerner.hpp         all complex roots at once (Durand–Kerner / Weierstrass)
│   └── cantor_zassenhaus.hpp     factorization and roots over Z_p
├── evaluation/
│   ├── fast_interpolation.hpp    Lagrange interpolation in O(N log² N)
│   └── multipoint_eval.hpp       evaluation at M points in O(N log² N) with subproduct trees
└── polynomial_lib.hpp            umbrella header including everything
```

The `series/` and `evaluation/` modules reach their target complexity only on top of fast multiplication, so `transforms/` comes first.

### Polynomial algebra

- polynomial GCD (Euclid's algorithm) and extended GCD (Bézout coefficients);
- composition P(Q(X));
- irreducibility tests, resultant and discriminant;
- Sturm sequences (number of real roots in an interval);
- parsing from text (`"3X^2 + 2X - 1"`) and LaTeX output.

### Fast algorithms

- NTT over ℤ<sub>p</sub>; the modulus 998244353 = 119 · 2²³ + 1 is suitable;
- FFT over complex numbers for real coefficients;
- Newton iteration for series inversion, then fast division;
- logarithm, exponential and square root of formal power series;
- Taylor shift: computing P(X + a);
- allocation-free Karatsuba;
- Montgomery / Barrett reduction for a modulus chosen at run time.

### New number types

- `fraction<T>`: exact rational numbers, a true field without rounding errors;
- `complex<T>` over any ring: Gaussian integers ℤ[i], ℚ(i), and GF(p²) when p ≡ 3 (mod 4);
- finite fields GF(pⁿ) as polynomials modulo an irreducible polynomial;
- bit-packed polynomials over GF(2);
- arbitrary-precision integers.

### Linear algebra

- determinant (Gaussian elimination over fields, Bareiss over rings), inverse, rank, linear systems;
- characteristic polynomial, connecting matrices and polynomials (with Cayley–Hamilton as a test);
- Berlekamp–Massey (shortest linear recurrence) and Kitamasa's method;
- rectangular and runtime-sized matrices.

### Applications

- **Signal processing:** convolution, FIR/IIR filters as ratios of polynomials in z⁻¹, Z-transform and frequency response, poles and zeros, Chebyshev polynomials, fixed-point arithmetic.
- **Coding theory and telecommunications:** CRC (polynomial division over GF(2)), Reed–Solomon over GF(2⁸), BCH and cyclic codes, LFSRs.
- **Hardware:** using the library as a reference model for FPGA implementations (polynomial multipliers, CRC, GF(2⁸) arithmetic).

### Engineering

- property-based tests on random data (`A == (A / B) * B + A % B`, Karatsuba equals schoolbook, `(P · Q)' = P' · Q + P · Q'`);
- compile-time tests of the concepts with `static_assert`;
- sanitizers (`-fsanitize=address,undefined`) and warning-free builds with `-Wall -Wextra -Wpedantic`;
- CMake target, continuous integration, Doxygen documentation, published benchmarks;
- modern C++: move semantics, `constexpr` everywhere possible, `std::span`, ranges, `std::expected`, C++20 modules.

---

## FAQ

**Why header-only?**
Template code must be visible wherever it is used, so a single header is the natural form for a template library. The trade-off is longer compile times as the library grows.

**Why are there no destructors or copy constructors?**
They are not needed. Memory is managed by `std::vector`, and the compiler-generated copy and move operations are correct (the "rule of zero").

**Why does `T(0)` appear in the concept?**
Algebraically, `T(0)` and `T(1)` are the identities for addition and multiplication. In code, they mean building an element from the integers 0 and 1. The concept only checks that they can be built; that they really are identities is a convention respected by the constructors.

**Why is `matrix` a `struct` and the others `class`?**
`modular_element` and `polynomial` have invariants to protect, so their data is private. Any n×n values form a valid matrix, so there is nothing to protect.

**Do the concepts cost anything at run time?**
No. They are checked entirely at compile time.

---

## Author

**Felix Florea**
