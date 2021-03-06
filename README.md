[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.46597.svg)](http://dx.doi.org/10.5281/zenodo.46597)

## Introduction

This repository contains the supplementary information for the paper 

Dmitri Chernyshenko, and Hans Fangohr <br>
*Computation of the magnetostatic interaction between linearly magnetized polyhedrons*<br>
Journal of Magnetism and Magnetic Materials **412**, 132-137 (2016)<br>

The paper is available online 

- Journal's URL: http://dx.doi.org/10.1016/j.jmmm.2016.03.083
- On the arXiv: http://arxiv.org/abs/1601.04185 

This software repository is available on GitHub [[1]](#github-link)</a> and has the DOI [10.5281/zenodo.46597](http://dx.doi.org/10.5281/zenodo.46597).

## Source code for the surface integral derivation

The primary source file is `surface-int-derivation/MagnetostaticTetTetInteractions.nb`. The derivation proceeds by repeatedly applying sets of transformation rules to transform the volume integrand to the surface integrand. The source code for the transformation rules and functions are in `surface-int-derivation/PotentialIntegration.nb`.

## Numerical integration weights

The weights for the numerical integration formulas are in `numint-weights/`. 

The integration point coordinates _(x<sub>i</sub>, y<sub>i</sub>, z<sub>i</sub>)_ are barycentric, i.e. _x<sub>i</sub> + y<sub>i</sub> + z<sub>i</sub> = 1_.

For an arbitrary triangle _T = (r<sub>1</sub>, r<sub>2</sub>,r<sub>3</sub>)_, the numerical integration formula over T is defined by 

_&int;<sub>T</sub> f(r) dr_ &approx; det(T) &middot; _&sum;<sub>i</sub>  w<sub>i</sub> f(x<sub>i</sub> r<sub>1</sub> + y<sub>i</sub> r<sub>2</sub> + z<sub>i</sub> r<sub>3</sub>)_

where 

det(T) = |(r<sub>2</sub> - r<sub>1</sub>) &times; (r<sub>3</sub> - r<sub>1</sub>)|

## Reference matrix _E<sub>ref</sub>_

The 12x12 reference test matrix obtained via finite difference simulation is in `Eref-test-data/E_ref.txt` - the construction of this matrix is described in the paper.

## Acknowledgements

This work was supported by an EPSRC Doctoral Training Centre grant (EP/G03690X/1).


## References

<a name="github-link">1</a>. [https://github.com/dc0/magnetostatic-polyhedrons](https://github.com/dc0/magnetostatic-polyhedrons)
