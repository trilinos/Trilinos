#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>

class Basis {

 public:
  double *phi, *dphide; 
  double uu, xx, duu, eta, wt;
  Basis();
  ~Basis();
  void getBasis(int gp, double *x, double *u);
  double dx;
};


