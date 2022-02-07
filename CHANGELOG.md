# Change Log
All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Breaking Changes
- replaced the term "function" that to be differentiated with "differentiand", and changed related the parameter, functions, and arguments name. #89
- replased the argument `same_subscripts_as_stencil` to `sort`. #74
- updated the `equation` in `interpolation` module to sort the results. #74

### Changes
- replaced % format and str.format() with fstring. #71
- moved functions related to stencil to stencil module. #81
- moved default parameter to spec module. #81

### New features
- added new module for handling stencil. #81
- added funciton to convert `int` or `float` to string for subscript. #76
- added function to sort expression by subscript. #74
- added function to extract each term in numerator. #74

### Repository updates
- update README to specify install version during cloning #83
- update README to hide install methods #83

## [0.3.0] - 2022-01-25
### Breaking changes
- changed default of `dot_product`'s `evaluate` argument to `True`. #62

### Changes
- changed algorithm computing dot product in `dot_product` when `evaluate=True`. #63
- replaced list comprehension to compute dot products of 2 lists with `dot_product`

### New features
- added new module for commonly used in tests to generate random input. #46
- added new module to describe the specifitions of methods used in dictos. #42, #50, #52, #53, #54
- added CHANGELOG.md.
- added linalg module. #51
- added version number in `__version__.py` and import in `__init__.py`. #45
- added `dictos/error/__init__.py` for packaging. #14
- added installation settings using pip. #14

### Repository updates
- added installation guide. #14

## [0.2.1] - 2022-01-21
### Fixes
- fixed a bug that returns an equation in which terms are not arranged in the order of the stencil when unsorted stencil are passed. #43, #4

## [0.2.0] - 2022-01-19
### Breaking changes
- changed argument's name `as_numr_denom` to `as_numer_denom` #31
- changed return value of `derivative symbol` with `deriv=0` from `f^(0)` to `f` #41

### Changes
- refactored coefficient extraction #35
- refactored `taylor_series` #41
- refactored rationalization operation #15

### New features
- added custom error modules #26
  - introducing custom errors #2, #3, #5, #6, #7, #8, #9, #10, #11, #21, #22, #32, #34

### Fixes
- fixed test class name #28

### Repository updates
- delete 0 from staggered grid table and examples in README #29

## [0.1.1] - 2022-01-15
### Fixes
- fixed algorithm of determining limit_denominator in simplify_coefficients. #1

### Repository updates
- updated README
- added coefficients tables

## 0.1.0 - 2022-01-15
### New feature
- lagrangian_polynomial
  - construct Lagrangian basis polynomial based on given sets of coordinate and function.
  - construct Lagrangian interpolation polynomial.
- taylor_expansion
  - calculate Taylor series.
  - create symbol representing derivative like f^(2).
- finite_difference
  - based on given stencil like [-1, 0, 1], derive finite difference equation at 0.
  - derive finite difference coefficients
  - derive the leading-order of truncation error.
- interpolation
  - based on given stencil, not including 0 like [-1, 1], derive interpolation equation at 0.
  - derive interpolation coefficients
  - derive the leading-order of truncation error.

[Unreleased]: https://github.com/degawa/dictos/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/degawa/dictos/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/degawa/dictos/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/degawa/dictos/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/degawa/dictos/compare/v0.1.0...v0.1.1
