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
  void getBasis(int gp, double *x, double *u, double *uold, double *aux);

 public:
  int NumSpecies;
  // Variables that are calculated at the gauss point
  double *phi, *dphide; 
  double xx, eta, wt;
  double uu, duu;
  double uuold, duuold;
  double aaux;
  double dx;
};

#endif
