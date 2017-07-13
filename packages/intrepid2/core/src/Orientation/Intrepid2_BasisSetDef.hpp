// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::BasisSet class.
    \author Created by Kyungjoo Kim
*/
#ifndef INTREPID2_BASISSETDEF_HPP
#define INTREPID2_BASISSETDEF_HPP

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
namespace Intrepid2 {

  template<class Scalar, class ArrayType>
  void
  Basis_Dummy_FEM<Scalar,ArrayType>::initializeTags() {}

  template<class Scalar, class ArrayType>
  void
  Basis_Dummy_FEM<Scalar,ArrayType>::getValues(ArrayType &       outputValues,
                                               const ArrayType & inputPoints,
                                               const EOperator   operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (Basis_Dummy_FEM): Dummy FEM basis is called getValues");
  }

  template<class Scalar, class ArrayType>
  void
  Basis_Dummy_FEM<Scalar,ArrayType>::getValues(ArrayType &       outputValues,
                                               const ArrayType & inputPoints,
                                               const ArrayType & cellVertices,
                                               const EOperator   operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (Basis_Dummy_FEM): Dummy FEM basis is called getValues");
  }

  template<class Scalar, class ArrayType>
  BasisSet<Scalar,ArrayType>::BasisSet(const EFunctionSpace space)
    : _space(space) {}

  template<class Scalar, class ArrayType>
  EFunctionSpace
  BasisSet<Scalar,ArrayType>::getFunctionSpace() const {
    return _space;
  }

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet<Scalar,ArrayType>::getCellBasis() const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (BasisSet): Dummy FEM basis is requested");
    return _dummy;
  }
  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet<Scalar,ArrayType>::getTriangleBasis() const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (BasisSet): Dummy FEM basis is requested");
    return _dummy;
  }
  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet<Scalar,ArrayType>::getQuadrilateralBasis() const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (BasisSet): Dummy FEM basis is requested");
    return _dummy;
  }
  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet<Scalar,ArrayType>::getLineBasis() const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                                ">>> ERROR (BasisSet): Dummy FEM basis is requested");
    return _dummy;
  }

  template<class Scalar, class ArrayType>
  BasisSet_HGRAD_TRI_Cn_FEM<Scalar,ArrayType>::BasisSet_HGRAD_TRI_Cn_FEM(const int n,
                                                                         const EPointType pointType)
    : BasisSet<Scalar,ArrayType>(FUNCTION_SPACE_HGRAD),
      _cellBasis(n, pointType),
      _lineBasis(n, pointType) { }

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet_HGRAD_TRI_Cn_FEM<Scalar,ArrayType>::getCellBasis() const {
    return _cellBasis;
  }

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet_HGRAD_TRI_Cn_FEM<Scalar,ArrayType>::getLineBasis() const {
    return _lineBasis;
  }

  template<class Scalar, class ArrayType>
  BasisSet_HGRAD_TET_Cn_FEM<Scalar,ArrayType>::BasisSet_HGRAD_TET_Cn_FEM(const int n,
                                                                         const EPointType pointType)
    : BasisSet<Scalar,ArrayType>(FUNCTION_SPACE_HGRAD),
      _cellBasis(n, pointType),
      _trigBasis(n, pointType),
      _lineBasis(n, pointType) {}

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet_HGRAD_TET_Cn_FEM<Scalar,ArrayType>::getCellBasis() const {
    return _cellBasis;
  }

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet_HGRAD_TET_Cn_FEM<Scalar,ArrayType>::getTriangleBasis() const {
    return _trigBasis;
  }

  template<class Scalar, class ArrayType>
  const Basis<Scalar,ArrayType>&
  BasisSet_HGRAD_TET_Cn_FEM<Scalar,ArrayType>::getLineBasis() const {
    return _lineBasis;
  }

}
#endif // INTREPID_USING_EXPERIMENTAL_HIGH_ORDER

#endif // INTREPID2_BASISSETDEF_HPP
