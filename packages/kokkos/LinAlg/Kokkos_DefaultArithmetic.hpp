//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTARITHMETIC_H
#define KOKKOS_DEFAULTARITHMETIC_H

namespace Kokkos {

  // FINISH: insert structs from NewAPI

  template <class MultiVector, class Vector>
  class DefaultArithmetic {
    public:
      typedef typename MultiVector::NodeType    NodeType;
      typedef typename MultiVector::OrdinalType OrdinalType;
      typedef typename MultiVector::ScalarType  ScalarType;

      //! Multiply one MultiVector by another, element-wise: B *= A
      void Multiply(const MultiVector &A, MultiVector &B) {
        // FINISH: if stride is one, one kernel invocation. otherwise, one for each column.
      };

      //! Divide one MultiVector by another, element-wise: B /= A
      void Divide(const MultiVector &A, MultiVector &B) {
        // FINISH: if stride is one, one kernel invocation. otherwise, one for each column.
      };

      // FINISH: add vector routines as well, if Vector isn't a MultiVector (i.e., doesn't have numCols())
  };

}

#endif
