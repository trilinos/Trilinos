// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <math.h>
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
//#include "NOX_Common.H"
#include "Basis.H"

// Constructor
Basis::Basis() {
  phi = new double[2];
  dphide = new double[2];
}

// Destructor
Basis::~Basis() {
  delete [] phi;
  delete [] dphide;
}

// Calculates a linear 1D basis
void Basis::getBasis(int gp, double *x, double *u) {
  int N = 2;
  if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;
  
  // Caculate basis function and derivative at GP.
  dx=0.5*(x[1]-x[0]);
  xx=0.0;
  uu=0.0;
  duu=0.0;
  for (int i=0; i < N; i++) {
    xx += x[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
  }

  return;
}
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO)
