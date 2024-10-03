// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Linear basis function for finite element method

#ifndef _NOX_EXAMPLE_EPETRA_LINEAR_BASIS_H
#define _NOX_EXAMPLE_EPETRA_LINEAR_BASIS_H

class Basis {

 public:

  // Constructor
  Basis(int numSpec);

  // Destructor
  ~Basis();

  // Calculates the values of u and x at the specified gauss point
  void getBasis(int gp, double *x, double *u, double *udot);

 private:
  // Private to prohibit copying
  Basis(const Basis&);
  Basis& operator=(const Basis&);

 public:
  int NumSpecies;
  // Variables that are calculated at the gauss point
  double *phi, *dphide;
  double xx, eta, wt;
  double *uu, *duu;
  double *uudot, *duudot;
  double dx;
};

#endif
