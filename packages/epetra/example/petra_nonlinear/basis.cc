/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include<iostream.h>
#include<math.h>
#include "basis.h"

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
  if (gp==0) {eta=-sqrt(3)/sqrt(5); wt=5.0/9.0;}
  if (gp==1) {eta=0.0; wt=8.0/9.0;}
  if (gp==2) {eta=sqrt(3)/sqrt(5); wt=5.0/9.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;
  
  // Caculate basis function and derivative at GP.
  dx=x[1]-x[0];
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
