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
  void getBasis(int gp, double *x, double *u, double *oldu);

 public:
  // Variables that are calculated at the gauss point
  double *phi, *dphide; 
  double uu, uuold, xx, duu, duuold, eta, wt;
  double dx;
};

#endif
