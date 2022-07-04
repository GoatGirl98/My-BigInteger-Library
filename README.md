# My-BigInteger-Library

A mini lib for arbitrary-precision arithmetic of BigInteger with high performance in C++11 (without STL).

[README 中文版本](https://github.com/GoatGirl98/My-BigInteger-Library/blob/main/README_cn.md)

## `10-Radix Big-Integer.cpp`

Library for 10-radix **signed** BigInteger, including all kinds of arithmetic operation. Including:

- **add ( $A+B$ )**
- **subtract ( $A-B$ )**
- **multiply ( $A*B$ )**
- **divide (abs round down) ( $\lfloor A/B \rfloor$ )**
- **modulo ( $A\mod B$ )**
- **power ( $A^B$ )**
- **log (abs round down) ($\lfloor \log_A B\rfloor$ )**
- **root (abs round down) ($ \lfloor \sqrt[A]{B} \rfloor $)**
- **greatest common divisor ($ \gcd (A, B) $)**
- **factorial ($ A! $)**
- **combination ( $C(n, m)$ )**

## `Arbitrary-Radix Big-Integer.cpp`

Library for arbitrary-radix **unsigned** BigInteger, including Base-Conversion between radix 2 to 62, with time complexity $O(n \log ^2 n)$.

## Remark

README for source code is not finished yet. To be continued...
