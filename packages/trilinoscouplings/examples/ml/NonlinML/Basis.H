// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
                                                                                
// 1D Linear basis function for finite element method

#ifndef _NOX_EXAMPLE_EPETRA_LINEAR_BASIS_H
#define _NOX_EXAMPLE_EPETRA_LINEAR_BASIS_H

class Basis {

 public:

  // Constructor
  Basis();

  // Destructor
  ~Basis();

  // Calculates the values of u and x at the specified gauss point
  void getBasis(int gp, double *x, double *u);

 public:
  // Variables that are calculated at the gauss point
  double *phi, *dphide; 
  double uu, xx, duu, eta, wt;
  double dx;
};

#endif
