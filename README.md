# dictos

Symbolic discretization tools for the finite difference method on a regular and staggered grid.

## A Story Telling Features

### Finite Difference
Tom, a university student assigned to a Computational Fluid Dynamics Laboratory, asks dictos.

Tom: "As part of my graduation research theme, I need to write a program simulating incompressible fluid flow in the finite difference method."
"I can find finite difference coefficients on the net. There are very thankful, but there are derived by discretization on the regular grid."
"I need the coefficients on the staggered grid. Dictos, can you provide those?"

Dictos: "Yes, I can."

Tom: "How should I do? I'm not good at computer operation and programming."

Dictos: "Import the `finite_difference` module and pass a **relative** stencil to `equation`."
```Python
import dictos.finite_difference as fd

stencil = [-2, -1, 0, 1, 2]
# stencil [i-2, i-1, i, i+1, i+2] relative to i.
print(fd.equation(stencil, deriv=1))
# (f_0 - 8*f_1 + 8*f_3 - f_4)/(12*h)

stencil = [-1.5, -0.5, 0, 0.5, 1.5]
print(fd.equation(stencil, deriv=1))
# (f_0 - 27*f_1 + 27*f_3 - f_4)/(24*h)
```

Tom: "I guess `deriv` is the order of derivate. Isn't the subscripts in the resulting equation the same as the stencil?"

Dictos: "use `same_subscripts_as_stencil` option."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0, 0.5, 1.5]
print(fd.equation(stencil, deriv=1, same_subscripts_as_stencil=True))
# (-27*f_{-0.5} + f_{-1.5} + 27*f_{0.5} - f_{1.5})/(24*h)
```

"set `evaluate=False`, get the equation that you may want."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0, 0.5, 1.5]
print(fd.equation(stencil, deriv=1, same_subscripts_as_stencil=True, evaluate=False))
# (f_{-1.5} - 27*f_{-0.5} + 0*f_{0} + 27*f_{0.5} - f_{1.5})/(24*h)
```

Tom: "Well... can I extract the coefficients?"

Dictos: "`coefficients` may be your best friends."
```Python
import dictos.finite_difference as fd

stencil = [-1.5, -0.5, 0, 0.5, 1.5]
print(fd.coefficients(stencil, deriv=1))
# [1/24, -9/8, 0, 9/8, -1/24]

print(fd.coefficients(stencil, deriv=1, as_numr_denom=True))
# ([1, -27, 0, 27, -1], 24)
```

Tom: "Thanks!"

### Interpolation
...A few days later

Tom: "I found that interpolation is also necessary for numerical simulations on the staggered grid. Dictos, can you provide the interpolation equation?"

Dictos: "Import the `interpolation` module. The usability is almost the same as the `finite_difference`."
"Do not contain `0` in the stencil because the module derives the interpolation equation at `0`."
```Python
import dictos.interpolation as intp

stencil = [-1.5, -0.5, 0.5, 1.5]
print(intp.equation(stencil))
# -f_0/16 + 9*f_1/16 + 9*f_2/16 - f_3/16

print(intp.equation(stencil, same_subscripts_as_stencil=True))
# 9*f_{-0.5}/16 - f_{-1.5}/16 + 9*f_{0.5}/16 - f_{1.5}/16
# `evalulate` option is not provided yet.

print(intp.coefficients(stencil))
# [-1/16, 9/16, 9/16, -1/16]

print(intp.coefficients(stencil, as_numr_denom=True))
# ([-1, 9, 9, -1], 16)
```

Tom: "Thanks!"

### Extrapolation
...A few days later

Tom: "I need to extrapolate the velocity outside the boundaries when computing grid points near boundaries. How complicated the staggered grid is!"
"Dictos, can you provide the extrapolation equation?"

Dictos: "It is possible to pass one-sided stencil to `interpolation`"
```Python
import dictos.interpolation as intp

stencil = [1, 2]
print(intp.equation(stencil,same_subscripts_as_stencil=True))
# 2*f_{1} - f_{2}

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

print(fd.truncation_error([-0.5, 0, 0.5], deriv=1))
# -f^(3)*h**2/24
print(fd.truncation_error([-1.5, -0.5, 0, 0.5, 1.5], deriv=1))
# 3*f^(5)*h**4/640
print(fd.truncation_error([-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5], deriv=1))
# -5*f^(7)*h**6/7168

print(intp.truncation_error([-0.5, 0.5]))
# -f^(2)*h**2/8
print(intp.truncation_error([-1.5, -0.5, 0.5, 1.5]))
# 3*f^(4)*h**4/128
print(intp.truncation_error([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]))
# -5*f^(6)*h**6/1024

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
- [ ] error handling
- [ ] update tests
