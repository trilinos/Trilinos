//@HEADER
// ************************************************************************
//
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef PRECONDITIONER_CRD_HPP
#define PRECONDITIONER_CRD_HPP

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "myzero.hpp"

using namespace std;

class preconditioner_crd 
{
public: // functions
  preconditioner_crd() { };
  virtual ~preconditioner_crd() { };
  virtual double norm2(double a[], int n)=0;
  virtual double dotprod(double a[], double b[], int n)=0;
  virtual void sum_vectors(double a[], int n, double a_sum[])=0;
  virtual int initialize_solve(double u[], double r[])=0;
  virtual void apply_preconditioner(const double r[], double z[])=0;
  virtual void A_times_x(double x[] , double Ax[])=0;
private:
protected:
  
};
#endif // PRECONDITIONER_CRD_HPP
