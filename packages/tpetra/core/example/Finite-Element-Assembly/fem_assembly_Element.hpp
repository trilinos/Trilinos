// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_ELEMENT_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_ELEMENT_HPP

#include <iomanip>

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
template <class ViewType>
KOKKOS_INLINE_FUNCTION void ReferenceQuad4(ViewType & elementMatrix) {
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
template<class ViewType>
KOKKOS_INLINE_FUNCTION void ReferenceQuad4RHS(ViewType& rhs) {
  for(size_t i=0; i<rhs.extent(0); i++)
    rhs[i] = static_cast<Scalar>(.25);
}

void ReferenceQuad4RHS(Teuchos::Array<Scalar>& rhs) {
  for(int i=0; (int)i<rhs.size(); i++)
    rhs[i] = static_cast<Scalar>(.25);
}

//
// This function prints out the quad4 array in a nice way.
//  rows x cols?
//
void PrettyPrintQuad4(scalar_2d_array_type & elementMatrix)
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
