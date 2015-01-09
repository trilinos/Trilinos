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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
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


