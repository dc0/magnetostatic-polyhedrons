## Introduction

This repository contains the supplementary information for the paper "Computation of the magnetostatic interaction between linearly magnetized polyhedrons" by D.Chernyshenko and 
H. Fangohr [[1]](#ref1), currently being prepared for publication.

The latest version of this repository is available on GitHub [[2]](#ref2)</a>.


## Source code for the surface integral derivation


## Numerical integration weights

The integration point coordinates _(x<sub>i</sub>, y<sub>i</sub>, z<sub>i</sub>)_ are barycentric, i.e. _x<sub>i</sub> + y<sub>i</sub> + z<sub>i</sub> = 1_.

For an arbitrary triangle _T = (r<sub>1</sub>, r<sub>2</sub>,r<sub>3</sub>)_, the numerical integration formula over T is defined by 

_&int;<sub>T</sub> f(r) dr_ &approx; det(T) &middot; _&sum;<sub>i</sub>  w<sub>i</sub> f(x<sub>i</sub> r<sub>1</sub> + y<sub>i</sub> r<sub>2</sub> + z<sub>i</sub> r<sub>3</sub>)_

where 

det(T) = |(r<sub>2</sub> - r<sub>1</sub>) &times; (r<sub>3</sub> - r<sub>1</sub>)|

## Reference matrix _E<sub>ref</sub>_


## References

<a name="ref1">1.</a> "Computation of the magnetostatic interaction between linearly magnetized polyhedrons", D.Chernyshenko and H. Fangohr

<a name="ref1">2.</a> [https://github.com/dc0/magnetostatic-polyhedrons](https://github.com/dc0/magnetostatic-polyhedrons)

## License

A short snippet describing the license (MIT, Apache, etc.)

