//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
