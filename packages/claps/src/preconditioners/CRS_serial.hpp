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

