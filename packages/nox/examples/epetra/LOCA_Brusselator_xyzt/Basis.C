// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "Basis.H"
#include <cstddef>
#include <cmath>

// Constructor
Basis::Basis(int numSpec) :
  NumSpecies(numSpec),
  phi(NULL),
  dphide(NULL),
  xx(0.0),
  eta(0.0),
  wt(0.0),
  uu(NULL),
  duu(NULL),
  uudot(NULL),
  duudot(NULL),
  dx(0.0)
{
  uu = new double[numSpec];
  duu = new double[numSpec];
  uudot = new double[numSpec];
  duudot = new double[numSpec];
  phi = new double[2];
  dphide = new double[2];
}

// Destructor
Basis::~Basis() {
  delete [] uu;
  delete [] duu;
  delete [] uudot;
  delete [] duudot;
  delete [] phi;
  delete [] dphide;
}

// Calculates a linear 1D basis
void Basis::getBasis(int gp, double *x, double *u, double *udot) {

  int N = 2;
  if (gp==0) {eta=-1.0/std::sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/std::sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;

  // Caculate basis function and derivative at GP.
  dx=0.5*(x[1]-x[0]);
  xx=0.0;
  for (int k=0; k<NumSpecies; k++)
    uu[k] = duu[k] = uudot[k] = duudot[k] = 0.0;

  for (int i=0; i < N; i++) {
    xx += x[i] * phi[i];
    for (int k=0; k<NumSpecies; k++) {
      uu[k] += u[NumSpecies * i + k] * phi[i];
      duu[k] += u[NumSpecies * i + k] * dphide[i];
      uudot[k] += udot[NumSpecies * i + k] * phi[i];
      duudot[k] += udot[NumSpecies * i + k] * dphide[i];
    }
  }

  return;
}
