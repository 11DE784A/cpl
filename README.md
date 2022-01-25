# Code for P452 Computational Physics

A numerical analysis library written in C. Contents reflect the structure of the
[Computational Physics](https://www.niser.ac.in/sps/course/p452-computational-physics) elective as taught in Spring 2022 at NISER.

The code is likely very bad and certainly very slow; please don't use it anywhere.

## Features

- A rudimentary `cpl_tensor` struct that can be used to store vectors, matrices
  and tensors
- Indices start from 1 in all interfaces to `cpl_tensor` objects
- Matrix reduction methods: Gauss-Jordan, ~~LU decomposition, Cholesky, QR
  decomposition, Jacobi iterations, Gauss-Siedel, Conjugate gradient, Minimal
  residue (MRES), Generalized minimal residue (GMRES)~~
- Literally no dependencies
- Written by me

## TODO

- [ ] Documentation
- [ ] More robust error handling
- [ ] Unit testing
