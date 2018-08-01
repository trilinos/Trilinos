// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_ELEMENT_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_ELEMENT_HPP

#include <iomanip>
#include "Kokkos_View.hpp"

namespace TpetraExamples 
{

//
// This function does a very rough approximation of 2D finite-difference
// stencil.  This is hardwired for QUAD4.
//
// Layout:
//
//    Node connectivity in the quad is 0 - 1 - 2 - 3 - 0 
//    0 ------ 1 
//    |        |
//    |        |
//    3 ------ 2 
//
//  - Each 'node' has a self edge and is assigned 2.
//  - Neighbors on cardinal directions are assigned -1.
//  - Diagonal neighbors aren't set (i.e., 0).
//
//        0   1   2   3
//      +--------------
//    0 | 2  -1      -1
//    1 |-1   2  -1
//    2 |    -1   2  -1
//    3 |-1      -1   2
// 
void ReferenceQuad4(scalar_2d_array_t & elementMatrix) {
  size_t lr[4] = {1,0,3,2};
  size_t ud[4] = {3,2,1,0};

  for (size_t i=0; i<4; i++)  {

    // Zero everything
    for (size_t j=0; j<4; j++) {
      elementMatrix(i,j) = 0.0;
    }

    // Diagonals
    elementMatrix(i,i) = 2.0;

    // Off-diagonals
    elementMatrix(i,lr[i]) = -1.0;
    elementMatrix(i,ud[i]) = -1.0;
  }
}

// RHS vector for the reference quad.
// This can be thought of as a unit source being equally distributed to the
// 4 nodes of the quad.
void ReferenceQuad4RHS(Teuchos::Array<Scalar>& rhs) {
  rhs.resize(4);
  std::fill(rhs.begin(), rhs.end(), static_cast<Scalar>(.25));
}

//
// This function prints out the quad4 array in a nice way.
//  rows x cols?
//
void PrettyPrintQuad4(scalar_2d_array_t & elementMatrix)
{
  size_t nr = elementMatrix.extent(0);
  size_t nc = elementMatrix.extent(1);
  for(size_t row_idx=0; row_idx<nr; row_idx++)
  { 
    std::cout << "[ ";
    for(size_t col_idx=0; col_idx<nc; col_idx++) {
      std::cout << std::setw(2) << elementMatrix(row_idx, col_idx) << " ";
    }
    std::cout << "]" << std::endl;
  }
}

} // end of namespace TpetraExamples


#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_ELEMENT_HPP
