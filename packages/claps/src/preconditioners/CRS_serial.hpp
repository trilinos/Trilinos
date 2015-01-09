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

#ifndef CRS_SERIAL_HPP
#define CRS_SERIAL_HPP

class CRS_serial 
{
 public: // functions
  CRS_serial(
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
     double*   z_);       // z-coordinates of nodes
  ~CRS_serial();
  int get_nrow();
  int get_ncol();
  int get_nnode();
  int get_nnz();
  int* get_nodebeg();
  int* get_localdof();
  double* get_xcoord();
  double* get_ycoord();
  double* get_zcoord();
  int* get_rowbeg();
  int* get_colidx();
  double* get_vals();
  void multiply(double* x, double* y);
 private: // functions

 private: // variables
  int *rowbeg, *colidx, nrow, ncol, nnode, *nodebeg, *localdof;
  double *A, *x, *y, *z;
};
#endif // CRS_SERIAL_HPP

