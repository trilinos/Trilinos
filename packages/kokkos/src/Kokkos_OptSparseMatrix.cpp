
//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include "Kokkos_OptMatrix.hpp"
using namespace Kokkos;
//==============================================================================
template<typename OrdinalType, typename ScalarType>
 OptMatrix<OrdinalType, ScalarType>::CrsMatrix() 
  : Kokkos_CompObject(),
    allocated_(false),
    numRows_(0),
    numCols_(0),
    numEntries_(0),
    values_(0),
    allValues_(0),
    indices_(0),
    allIndices_(0),
    pntr_(0),
    profile_(0)
{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
OptMatrix<OrdinalType, ScalarType>::OptMatrix(const OptMatrix<OrdinalType, ScalarType> &matrix) 
  : Kokkos_CompObject(matrix),
    allocated_(matrix.allocated_),
    numRows_(matrix.numRows_),
    numCols_(matrix.numCols_),
    numEntries_(matrix.numEntries_),
    values_(matrix.values_),
    allValues_(matrix.allValues_),
    indices_(matrix.indices_),
    allIndices_(matrix.allIndices_),
    pntr_(matrix.pntr_),
    profile_(matrix.profile_)

{
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
void OptMatrix<OrdinalType, ScalarType>::initializeDefaults() { // Initialize all attributes that have trivial default values
  return;
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
void OptMatrix<OrdinalType, ScalarType>::allocate() { // Initialize all attributes that have trivial default values
  return;
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
void OptMatrix<OrdinalType, ScalarType>::~Kokkos_OptMatrix(){
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OptMatrix<OrdinalType, ScalarType>::apply(OrdinalType xLength, ScalarType * x, 
					      OrdinalType yLength, ScalarType * y) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OptMatrix<OrdinalType, ScalarType>::apply(OrdinalType numVectors, 
					      OrdinalType xLength, ScalarType ** x, 
					      OrdinalType yLength, ScalarType ** y) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OptMatrix<OrdinalType, ScalarType>::applyInverse(OrdinalType xLength, ScalarType * x, 
						     OrdinalType yLength, ScalarType * y) const {
}

//==============================================================================
template<typename OrdinalType, typename ScalarType>
int OptMatrix<OrdinalType, ScalarType>::applyInverse(OrdinalType numVectors, 
						     OrdinalType xLength, ScalarType ** x, 
						     OrdinalType yLength, ScalarType ** y) const {
}
