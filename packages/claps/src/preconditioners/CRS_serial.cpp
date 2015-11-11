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

#include <stdio.h>
#include <iostream>
#include "CRS_serial.hpp"

CRS_serial::CRS_serial(
     double*   A_,        // nonzero entries of matrix
     int*      rowbeg_,   // beginning of rows array
     int*      colidx_,   // column indices
     int       nrow_,     // number of rows in matrix
     int       ncol_,     // number of columns in matrix
     int       nnode_,    // number of nodes
     int*      nodebeg_,  // starting locations for nodal dofs
     int*      localdof_, // local dof numbers (1-7)
     double*   x_,        // x-coordinates of nodes
     double*   y_,        // y-coordinates of nodes
     double*   z_)        // z-coordinates of nodes
  : A(A_), rowbeg(rowbeg_), colidx(colidx_), nrow(nrow_), ncol(ncol_), 
    nnode(nnode_), nodebeg(nodebeg_), localdof(localdof_),
    x(x_), y(y_), z(z_)
{
}

CRS_serial::~CRS_serial()
{
}

int CRS_serial::get_nrow()
{
  return nrow;
}

int CRS_serial::get_ncol()
{
  return ncol;
}

int CRS_serial::get_nnode()
{
  return nnode;
}

int CRS_serial::get_nnz()
{
  return rowbeg[nrow];
}

int* CRS_serial::get_nodebeg()
{
  return nodebeg;
}

int* CRS_serial::get_localdof()
{
  return localdof;
}

double* CRS_serial::get_xcoord()
{
  return x;
}

double* CRS_serial::get_ycoord()
{
  return y;
}

double* CRS_serial::get_zcoord()
{
  return z;
}

int* CRS_serial::get_rowbeg()
{
  return rowbeg;
}

int* CRS_serial::get_colidx()
{
  return colidx;
}

double* CRS_serial::get_vals()
{
  return A;
}

void CRS_serial::multiply(double* xx, double* yy)
{
  int i, j;
  double sum;
  for (i=0; i<nrow; i++) {
    sum = 0;
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
      sum += A[j]*xx[colidx[j]];
    yy[i] = sum;
  }
}


