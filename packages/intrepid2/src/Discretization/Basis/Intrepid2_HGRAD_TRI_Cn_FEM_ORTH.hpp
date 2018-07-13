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

/** \file   Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH class.
    \author Created by Robert Kirby
            Kokkorized by Kyungjoo Kim and Mauro Perego
 */

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_HPP__

#include "Intrepid2_Basis.hpp"

namespace Intrepid2 {

/** \class  Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
      \brief  Implementation of the default H(grad)-compatible orthogonal basis (Dubiner) of
      arbitrary degree on triangle.

      \remarks

      \li   All degrees of freedom are considered to be internal (ie not assembled)
 */

namespace Impl {

/**
   \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
*/
template<typename outputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv, ordinal_type n>
struct OrthPolynomialTri {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(       outputViewType output,
            const inputViewType input,
                  workViewType work,
            const ordinal_type p );
};

/**
   \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
*/
template<typename outputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
struct OrthPolynomialTri<outputViewType,inputViewType,workViewType,hasDeriv,0> {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(       outputViewType output,
            const inputViewType input,
                  workViewType work,
            const ordinal_type p );
};

/**
   \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
*/
template<typename outputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
struct OrthPolynomialTri<outputViewType,inputViewType,workViewType,hasDeriv,1> {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(   outputViewType output,
      const inputViewType input,
            workViewType work,
      const ordinal_type p );
};

/**
   \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
*/
class Basis_HGRAD_TRI_Cn_FEM_ORTH {
public:

  /**
     \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
  */
  template<EOperator opType>
  struct Serial {
    template<typename outputViewType,
    typename inputViewType,
    typename workViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getValues(       outputViewType output,
               const inputViewType  input,
                     workViewType work,
               const ordinal_type   order);
  };

  template<typename ExecSpaceType, ordinal_type numPtsPerEval,
  typename outputValueValueType, class ...outputValueProperties,
  typename inputPointValueType,  class ...inputPointProperties>
  static void
  getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const ordinal_type order,
             const EOperator operatorType );

  /**
     \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
  */
  template<typename outputValueViewType,
  typename inputPointViewType,
  typename workViewType,
  EOperator opType,
  ordinal_type numPtsEval>
  struct Functor {
          outputValueViewType _outputValues;
    const inputPointViewType  _inputPoints;
          workViewType        _work;
    const ordinal_type        _order;

    KOKKOS_INLINE_FUNCTION
    Functor(       outputValueViewType outputValues_,
                   inputPointViewType  inputPoints_,
                   workViewType        work_,
             const ordinal_type        order_ )
    : _outputValues(outputValues_), _inputPoints(inputPoints_), _work(work_),_order(order_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type iter) const {
      const auto ptBegin = Util<ordinal_type>::min(iter*numPtsEval,    _inputPoints.extent(0));
      const auto ptEnd   = Util<ordinal_type>::min(ptBegin+numPtsEval, _inputPoints.extent(0));

      const auto ptRange = Kokkos::pair<ordinal_type,ordinal_type>(ptBegin, ptEnd);
      const auto input   = Kokkos::subview( _inputPoints, ptRange, Kokkos::ALL() );

      switch (opType) {
      case OPERATOR_VALUE : {
        auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange );
        Serial<opType>::getValues( output, input, _work, _order );  //here work is not used
        break;
      }
      case OPERATOR_GRAD :
      case OPERATOR_D1 :
      {
        const auto work = Kokkos::subview( _work, Kokkos::ALL(), ptRange, Kokkos::ALL() );
        auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
        Serial<opType>::getValues( output, input, work, _order);
        break;
      }
      case OPERATOR_D2 :
      case OPERATOR_D3 :
      case OPERATOR_D4 :
      case OPERATOR_D5 :
      case OPERATOR_D6 :
      case OPERATOR_D7 :
      case OPERATOR_D8 :
      case OPERATOR_D9 :
      case OPERATOR_D10 :
      {
        auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
        Serial<opType>::getValues( output, input, _work, _order); //here work is not used
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
            ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::Functor) operator is not supported");

      }
      }
    }
  };

};

}

template<typename ExecSpaceType = void,
    typename outputValueType = double,
    typename pointValueType = double>
class Basis_HGRAD_TRI_Cn_FEM_ORTH
    : public Basis<ExecSpaceType,outputValueType,pointValueType> {
    public:
  typedef double value_type;
  typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_1d_host ordinal_type_array_1d_host;
  typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_2d_host ordinal_type_array_2d_host;
  typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_3d_host ordinal_type_array_3d_host;

  /** \brief  Constructor.
   */
  Basis_HGRAD_TRI_Cn_FEM_ORTH( const ordinal_type order );

  typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::outputViewType outputViewType;
  typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::pointViewType  pointViewType;

  using Basis<ExecSpaceType,outputValueType,pointValueType>::getValues;

  virtual
  void
  getValues(       outputViewType outputValues,
             const pointViewType  inputPoints,
             const EOperator operatorType = OPERATOR_VALUE ) const {
    #ifdef HAVE_INTREPID2_DEBUG
          Intrepid2::getValues_HGRAD_Args(outputValues,
                                          inputPoints,
                                          operatorType,
                                          this->getBaseCellTopology(),
                                          this->getCardinality() );
    #endif
    constexpr ordinal_type numPtsPerEval = Parameters::MaxNumPtsPerBasisEval;
    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    getValues<ExecSpaceType,numPtsPerEval>( outputValues,
        inputPoints,
        this->getDegree(),
        operatorType );
  }
};

}// namespace Intrepid2

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTHDef.hpp"

#endif

