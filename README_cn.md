# My-BigInteger-Library

一个可以在任意精度下完成算数运算的大整数库，采用标准为C++11(不使用STL)。

[README English Version](https://github.com/GoatGirl98/My-BigInteger-Library/blob/main/README.md)

## `10-Radix Big-Integer.cpp`

该文件是针对**10进制有符号大整数**的各种运算的大整数库。包括如下运算：

- **大数加法 ( $A+B$ )**
- **大数减法 ( $A-B$ )**
- **大数乘法 ( $A*B$ )**
- **大数除法 (结果按绝对值向下取整) ( $\lfloor A/B \rfloor$ )**
- **大数取模 ( $A\%B$ )**
- **大数次幂 ( $A^B$ )**
- **大数取对数 (结果按绝对值向下取整) ($\lfloor \log_A B\rfloor$ )**
- **大数开根 (结果按绝对值向下取整) ($ \lfloor \sqrt[A]{B} \rfloor $)**
- **最大公约数 ($ \gcd (A, B) $)**
- **阶乘 ($A!$)**
- **组合数 ( $C(n, m)$ )**

## `Arbitrary-Radix Big-Integer.cpp`

该文件是针对**2到62进制之间的无符号大整数**的进制转换运算的大整数库。其时间复杂度为 $O(n \log ^2 n)$ 。

## 备注

README针对源文件的详细解释尚未完成，有待补充...
