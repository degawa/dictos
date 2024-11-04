# Change Log
All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Fixes
- fixed missing the 2nd argument of `Expr._latex`. #124

## [0.5.0] - 2024-11-04
### Breaking Changes
- changed project directory structure. #108
- moved default symbols to dictos/default.py. #112

### Changes
- refactored rationalization algorithm from my own implementation to sympy `nsimplify`. #93
- added comment to a parameter `_MINIMUM_STENCIL_WIDTH`. #95
- added `is_positive_integer` and `is_not_positive_integer` methods to be used instead of `is_natural_number` and `is_not_natural_number`. #110

### New features
- added error handling when passing out-of-range values to argument `point_at` of `lagrangian_basis`. #48
- added `Expr` class that extends SymPy's `Expr` to preserve expression display order. #101
- added `generate` in `finite_difference` module for generating the finite difference equation or coefficient from deriv, acc, and grid type. #70
- added an enumerator `GridType` for spceifiying the grid type on `generate` method. #122

### Fixes
- fixed a typo (natrual to natural). #94
- fixed sympy class references from full-path to alias. #104
- fixed tests along with changes in directory structure. #107
- fixed `finite_difference.equation` to remove the coefficient from terms with a coefficient of 1. #118

### Repository updates
- fixed equations in README. #117
- update examples in README along with the current behavior. #117, #118

## [0.4.0] - 2022-02-08
### Breaking Changes
- replaced the term "function" that to be differentiated with "differentiand", and changed related the parameter, functions, and arguments name. #89
- replased the argument `same_subscripts_as_stencil` to `sort`. #74
- updated the `equation` in `interpolation` module to sort the results. #74
- replased the argument `evaluate` in the `equation` in `finite_difference` module with `keep_zero`. #80
- changed the default value of `keep_zero` to `False`. #80

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
- update README to specify install version during cloning. #83
- update README to hide install methods. #83
- update README to reflect breaking changes. #74, #80

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

[Unreleased]: https://github.com/degawa/dictos/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/degawa/dictos/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/degawa/dictos/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/degawa/dictos/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/degawa/dictos/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/degawa/dictos/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/degawa/dictos/compare/v0.1.0...v0.1.1
