# dictos

Symbolic discretization tools for the finite difference method on a regular and staggered grid.

## Installation
There are three ways to install dictos.

### install using pip
The simplest way to install dictos is using pip.
After cloning the dictos repository, install dictos with the commands below:

```console
$ git clone https://github.com/degawa/dictos.git
$ cd dictos
$ pip install .
```

<details><summary>other installation ways to avoid using pip</summary><div>

### copy dictos directory to your projects
The second way to install dictos is simply copying the `dictos` directory in the cloned repository to other projects.

The central part of dictos is in `dictos/dictos` directory.

```console
$ git clone https://github.com/degawa/dictos.git
$ tree dictos -d
dictos
├── dictos
│   └── error
└── test
```

If importing dictos from `main.py`, copy `dictos/dictos` to the same location and then import in `main.py` like below:

```console
exampleproject/
└── src/
    ├── dictos/
    └── main.py
```

```python
from dictos import finite_difference as fd
```

### add path to `sys.path`
The third way to install dictos is to set the path to tell the dictos directory to Python.

When dictos is cloned to a directory `/home/hoge/dictos`, add the dictos path to `sys.path` as shown below before importing dictos:

```console
$ git clone https://github.com/degawa/dictos.git
$ cd dictos
$ pwd
/home/hoge/dictos
```

```python
import sys

sys.path.insert(1, "/home/hoge/dictos")

from dictos import finite_difference as fd
```

</div>
</details>

## A Story Telling Features

### Finite Difference
Tom, a university student assigned to a Computational Fluid Dynamics Laboratory, asks dictos.

Tom: "As part of my graduation research theme, I need to write a program simulating incompressible fluid flow using the finite difference method."
"I can find finite difference coefficients on the net. There are very thankful, but those are derived on the regular grid."
"I need the coefficients on the staggered grid. Dictos, can you provide those?"

Dictos: "Yes, I can."

Tom: "How should I do? I'm not good at computer operation and programming."

Dictos: "Import the `finite_difference` module and pass a **relative** stencil to `equation`. I derive a finite difference equation at `0`."
```Python
import dictos.finite_difference as fd

stencil = [-2, -1, 0, 1, 2]
# stencil [i-2, i-1, i, i+1, i+2] relative to i.
# o---o---x---o---o
# -2  -1  0   1   2
print(fd.equation(stencil, deriv=1))
# (f_0 - 8*f_1 + 8*f_3 - f_4)/(12*h)

stencil = [-1.5, -0.5, 0.5, 1.5]
#    +---+-x-+---+
# -1.5-0.5 0 0.5 1.5
print(fd.equation(stencil, deriv=1))
# (f_0 - 27*f_1 + 27*f_3 - f_4)/(24*h)
```

Tom: "I guess `deriv` is the order of derivate. Isn't the subscripts in the resulting equation the same as the stencil?"

Dictos: "use `same_subscripts_as_stencil` option."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0.5, 1.5]
print(fd.equation(stencil, deriv=1, same_subscripts_as_stencil=True))
# (-27*f_{-0.5} + f_{-1.5} + 27*f_{0.5} - f_{1.5})/(24*h)
```

"set `evaluate=False`, get the equation that you may want."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0, 0.5, 1.5]
print(fd.equation(stencil, deriv=1, same_subscripts_as_stencil=True, evaluate=False))
# (f_{-1.5} - 27*f_{-0.5} + 27*f_{0.5} - f_{1.5})/(24*h)
```

Tom: "Well... can I extract the coefficients?"

Dictos: "`coefficients` may be your best friends."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0.5, 1.5]
print(fd.coefficients(stencil, deriv=1))
# [1/24, -9/8, 9/8, -1/24]

print(fd.coefficients(stencil, deriv=1, as_numer_denom=True))
# ([1, -27, 27, -1], 24)
```

Tom: "Thanks!"

### Interpolation
...A few days later

Tom: "I found that interpolation is also necessary for numerical simulations on the staggered grid. Dictos, can you provide the interpolation equation?"

Dictos: "Import the `interpolation` module. The usability is almost the same as the `finite_difference`."
"Do not contain `0` in the stencil because the I derive the interpolation equation at `0`."
```Python
import dictos.interpolation as intp

stencil = [-1.5, -0.5, 0.5, 1.5]
#    +---+-x-+---+
# -1.5-0.5 0 0.5 1.5
print(intp.equation(stencil))
# -f_0/16 + 9*f_1/16 + 9*f_2/16 - f_3/16

print(intp.equation(stencil, same_subscripts_as_stencil=True))
# 9*f_{-0.5}/16 - f_{-1.5}/16 + 9*f_{0.5}/16 - f_{1.5}/16
# `evalulate` option is not provided yet.

print(intp.coefficients(stencil))
# [-1/16, 9/16, 9/16, -1/16]

print(intp.coefficients(stencil, as_numer_denom=True))
# ([-1, 9, 9, -1], 16)
```

Tom: "Thanks!"

### Extrapolation
...A few days later

Tom: "I need to extrapolate the velocity outside the boundaries when computing differences at grid points near boundaries. How complicated the staggered grid is!"
"Dictos, can you provide the extrapolation equation?"

Dictos: "It is possible to pass one-sided stencil to `interpolation`"
```Python
import dictos.interpolation as intp

stencil = [1, 2]
# x---o---o
# 0   1   2
print(intp.equation(stencil,same_subscripts_as_stencil=True))
# 2*f_{1} - f_{2}

#  o---o---o---x
# -3  -2  -1   0
stencil = [-3, -2, -1]
print(intp.equation(stencil,same_subscripts_as_stencil=True))
# 3*f_{-1} - 3*f_{-2} + f_{-3}
```

Tom: "Thanks!"

### Truncation Error
...A few months later

Tom: "My boss has requested me to derive formal truncation errors."
"Di..."

Dictos: "`truncation_error`. `finite_differene` and `interplation` have `truncation_error`"
```Python
import dictos.finite_difference as fd
import dictos.interpolation as intp

# finite difference
print(fd.truncation_error([-0.5, 0.5], deriv=1))
# -f^(3)*h**2/24
print(fd.truncation_error([-1.5, -0.5, 0.5, 1.5], deriv=1))
# 3*f^(5)*h**4/640
print(fd.truncation_error([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], deriv=1))
# -5*f^(7)*h**6/7168

# interpolation
print(intp.truncation_error([-0.5, 0.5]))
# -f^(2)*h**2/8
print(intp.truncation_error([-1.5, -0.5, 0.5, 1.5]))
# 3*f^(4)*h**4/128
print(intp.truncation_error([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]))
# -5*f^(6)*h**6/1024

# extrapolation
print(intp.truncation_error([1, 2]))
# f^(2)*h**2
print(intp.truncation_error([1, 2, 3]))
# -f^(3)*h**3
```

Tom: "Thanks!"
"It would be nice if the truncation error were provided when I asked for an equation or coefficients.

Dictos: ☺

### Features to be implemented
Tom successfully  graduated from the university and will go on to a graduate school in the new academic year.

Tom: "I'm going to research numerical simulations of aeroacoustics using the compact finite difference in the graduate school."
"Dictos, do you support the compact finite difference?"

Dictos: "Not supported yet. I will support it."

Tom: "Since the leading error term of the compact finite difference is dispersive, a low-pass filter is needed to suppress numerical oscillation."
"Dictos, do you support the low-pass filter?"

Dictos: "Not supported yet. I will support explicit and compact filters."

## feature
- derive a finite difference equation on regular or staggered grid based on a given stencil.
- derive an interpolation formula on regular or staggered grid based on a given stencil.
- derive an extrapolation formula on regular or staggered grid based on a given stencil.
- calculate a formal truncation error of a finite difference equation or an interpolation formula.

## todo
- [ ] add documents
- [ ] add examples
- [ ] update tests

## Finite Difference Coefficients
### Central Finite Difference on the Regular Grid
#### Notation

|||
|:--|:--|
|![S=\{-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5\}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+S%3D%5C%7B-5%2C+-4%2C+-3%2C+-2%2C+-1%2C+0%2C+1%2C+2%2C+3%2C+4%2C+5%5C%7D)|stencil|
|![n_s](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+n_s%0A)|numerator|
|![d](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+d)|denominator|
|![h](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+h)|interval between grid points|

#### 1st-order derivative
On the regular grid, the 1st-order derivative is calculated by the following equation:

![f_i^{(1)} = \frac{1}{d\cdot h}\sum_{s\in S} n_s f_{i+s}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f_i%5E%7B%281%29%7D+%3D+%5Cfrac%7B1%7D%7Bd%5Ccdot+h%7D%5Csum_%7Bs%5Cin+S%7D+n_s+f_%7Bi%2Bs%7D)


| Order of Accuracy |  -5   |  -4   |  -3   |  -2   |  -1   |   0   |   1   |   2   |   3   |   4   |   5   | Denominator |  Trunctaion Error  |
| :---------------: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---------: | :----------------- |
|         2         |   0   |   0   |   0   |   0   |  -1   |   0   |   1   |   0   |   0   |   0   |   0   |      2      | -f^(3)*h**2/6      |
|         4         |   0   |   0   |   0   |   1   |  -8   |   0   |   8   |  -1   |   0   |   0   |   0   |     12      | f^(5)*h**4/30      |
|         6         |   0   |   0   |  -1   |   9   |  -45  |   0   |  45   |  -9   |   1   |   0   |   0   |     60      | -f^(7)*h**6/140    |
|         8         |   0   |   3   |  -32  |  168  | -672  |   0   |  672  | -168  |  32   |  -3   |   0   |     840     | f^(9)*h**8/630     |
|        10         |  -2   |  25   | -150  |  600  | -2100 |   0   | 2100  | -600  |  150  |  -25  |   2   |    2520     | -f^(11)*h**10/2772 |


Note that these coefficients are NOT for a usual equation:

![f_i^{(1)} = a\frac{f_{i+1}-f_{i-1}}{2h} + b\frac{f_{i+2}-f_{i-2}}{4h} + c\frac{f_{i+3}-f_{i-3}}{6h}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f_i%5E%7B%281%29%7D+%3D+a%5Cfrac%7Bf_%7Bi%2B1%7D-f_%7Bi-1%7D%7D%7B2h%7D+%2B+b%5Cfrac%7Bf_%7Bi%2B2%7D-f_%7Bi-2%7D%7D%7B4h%7D+%2B+c%5Cfrac%7Bf_%7Bi%2B3%7D-f_%7Bi-3%7D%7D%7B6h%7D)

#### 2nd-order derivative
The 2nd-order derivative is calculated by the following equation:

![f_i^{(2)} = \frac{1}{d\cdot h^2}\sum_{s\in S} n_s f_{i+s}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f_i%5E%7B%282%29%7D+%3D+%5Cfrac%7B1%7D%7Bd%5Ccdot+h%5E2%7D%5Csum_%7Bs%5Cin+S%7D+n_s+f_%7Bi%2Bs%7D%0A)

| Order of Accuracy |  -5   |  -4   |  -3   |  -2   |  -1   |   0    |   1   |   2   |   3   |   4   |   5   | Denominator |  Trunctaion Error   |
| :---------------: | :---: | :---: | :---: | :---: | :---: | :----: | :---: | :---: | :---: | :---: | :---: | :---------: | :------------------ |
|         2         |   0   |   0   |   0   |   0   |   1   |   -2   |   1   |   0   |   0   |   0   |   0   |      1      | -f^(4)*h**2/12      |
|         4         |   0   |   0   |   0   |  -1   |  16   |  -30   |  16   |  -1   |   0   |   0   |   0   |     12      | f^(6)*h**4/90       |
|         6         |   0   |   0   |   2   |  -27  |  270  |  -490  |  270  |  -27  |   2   |   0   |   0   |     180     | -f^(8)*h**6/560     |
|         8         |   0   |  -9   |  128  | -1008 | 8064  | -14350 | 8064  | -1008 |  128  |  -9   |   0   |    5040     | f^(10)*h**8/3150    |
|        10         |   8   | -125  | 1000  | -6000 | 42000 | -73766 | 42000 | -6000 | 1000  | -125  |   8   |    25200    | -f^(12)*h**10/16632 |

Note that these coefficients are NOT for a usual equation:

![f_i^{(2)} = a\frac{f_{i+1}-2f_i+f_{i-1}}{h^2} + b\frac{f_{i+2}-2f_i+f_{i-2}}{(2h)^2} + c\frac{f_{i+3}-2f_i+f_{i-3}}{(3h)^2}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f_i%5E%7B%282%29%7D+%3D+a%5Cfrac%7Bf_%7Bi%2B1%7D-2f_i%2Bf_%7Bi-1%7D%7D%7Bh%5E2%7D+%2B+b%5Cfrac%7Bf_%7Bi%2B2%7D-2f_i%2Bf_%7Bi-2%7D%7D%7B%282h%29%5E2%7D+%2B+c%5Cfrac%7Bf_%7Bi%2B3%7D-2f_i%2Bf_%7Bi-3%7D%7D%7B%283h%29%5E2%7D)
### Central Finite Difference on the Staggered Grid
#### Notation

|||
|:--|:--|
|![S=\{-4.5, -3.5, -2.5, -1.5,-0.5, 0.5, 1.5, 2.5, 3.5, 4.5\}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+S%3D%5C%7B-4.5%2C+-3.5%2C+-2.5%2C+-1.5%2C-0.5%2C+0.5%2C+1.5%2C+2.5%2C+3.5%2C+4.5%5C%7D)|stencil|
|![n_s](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+n_s%0A)|numerator|
|![d](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+d)|denominator|
|![h](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+h)|interval between grid points|

#### 1st-order derivative
On the staggered grid, the 1st-order derivative is calculated by the equation with the same form as it on the regular grid:

![f_i^{(1)} = \frac{1}{d\cdot h}\sum_{s\in S} n_s f_{i+s}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+f_i%5E%7B%281%29%7D+%3D+%5Cfrac%7B1%7D%7Bd%5Ccdot+h%7D%5Csum_%7Bs%5Cin+S%7D+n_s+f_%7Bi%2Bs%7D)

| Order of Accuracy | -4.5  | -3.5  |  -2.5   |  -1.5  |   -0.5    |   0.5    |   1.5   |  2.5   |  3.5   |  4.5  | Denominator |     Trunctaion Error     |
| :---------------: | :---: | :---: | :-----: | :----: | :-------: | :------: | :-----: | :----: | :----: | :---: | :---------: | :----------------------- |
|         2         |   0   |   0   |    0    |   0    |    -1     |    1     |    0    |   0    |   0    |   0   |      1      | -f^(3)*h**2/24           |
|         4         |   0   |   0   |    0    |   1    |    -27    |    27    |   -1    |   0    |   0    |   0   |     24      | 3*f^(5)*h**4/640         |
|         6         |   0   |   0   |   -9    |  125   |   -2250   |   2250   |  -125   |   9    |   0    |   0   |    1920     | -5*f^(7)*h**6/7168       |
|         8         |   0   |  75   |  -1029  |  8575  |  -128625  |  128625  |  -8575  |  1029  |  -75   |   0   |   107520    | 35*f^(9)*h**8/294912     |
|        10         | -1225 | 18225 | -142884 | 926100 | -12502350 | 12502350 | -926100 | 142884 | -18225 | 1225  |  10321920   | -63*f^(11)*h**10/2883584 |

#### higher-order derivative
On the staggered grid, the odd-order derivative is calculated at each cell center, and the even-order derivative is calculated at each point.
This means we can reuse the equations for even-order derivative on the regular grid. However, to maintain consistency between 2nd-order derivative and 1st-order derivative of 1st-order derivate, the 2nd-order derivative is calculated as the 1st-order derivative of the 1st-order derivative. The higher-order derivatives than 2nd-order must be calculated in the same way.

The stencil width will be wider, but the consistency is more important.

## Interpolation Coefficients
### Notation

|||
|:--|:--|
|![S=\{-4.5, -3.5, -2.5, -1.5,-0.5, 0.5, 1.5, 2.5, 3.5, 4.5\}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+S%3D%5C%7B-4.5%2C+-3.5%2C+-2.5%2C+-1.5%2C-0.5%2C+0.5%2C+1.5%2C+2.5%2C+3.5%2C+4.5%5C%7D)|stencil|
|![n_s](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+n_s%0A)|numerator|
|![d](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+d)|denominator|
|![h](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+h)|interval between grid points|

### Central Interpolation on the Staggered Grid
The central interpolation is calculated by the following equation:

![\overline{f}_i = \frac{1}{d}\sum_{s\in S} n_s f_{i+s}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Coverline%7Bf%7D_i+%3D+%5Cfrac%7B1%7D%7Bd%7D%5Csum_%7Bs%5Cin+S%7D+n_s+f_%7Bi%2Bs%7D%0A)

| Order of Accuracy | -4.5  | -3.5  | -2.5  | -1.5  | -0.5  |  0.5  |  1.5  |  2.5  |  3.5  |  4.5  | Denominator |    Trunctaion Error     |
| :---------------: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---------: | :---------------------- |
|         2         |   0   |   0   |   0   |   0   |   1   |   1   |   0   |   0   |   0   |   0   |      2      | -f^(2)*h**2/8           |
|         4         |   0   |   0   |   0   |  -1   |   9   |   9   |  -1   |   0   |   0   |   0   |     16      | 3*f^(4)*h**4/128        |
|         6         |   0   |   0   |   3   |  -25  |  150  |  150  |  -25  |   3   |   0   |   0   |     256     | -5*f^(6)*h**6/1024      |
|         8         |   0   |  -5   |  49   | -245  | 1225  | 1225  | -245  |  49   |  -5   |   0   |    2048     | 35*f^(8)*h**8/32768     |
|        10         |  35   | -405  | 2268  | -8820 | 39690 | 39690 | -8820 | 2268  | -405  |  35   |    65536    | -63*f^(10)*h**10/262144 |

The above coefficients can not be used to the more general interpolation equation expressed as follows:

![\overline{f}_i = a\frac{f_{i+1}+f_{i-1}}{2} +b\frac{f_{i+2}+f_{i-2}}{2}+c\frac{f_{i+3}+f_{i-3}}{2}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Coverline%7Bf%7D_i+%3D+a%5Cfrac%7Bf_%7Bi%2B1%7D%2Bf_%7Bi-1%7D%7D%7B2%7D+%2Bb%5Cfrac%7Bf_%7Bi%2B2%7D%2Bf_%7Bi-2%7D%7D%7B2%7D%2Bc%5Cfrac%7Bf_%7Bi%2B3%7D%2Bf_%7Bi-3%7D%7D%7B2%7D)
