// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TET_Cn_FEM_ORTH.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH class.
    \author Created by Robert C. Kirby
    
    Kokkorized by Kyungjoo Kim and Mauro Perego
*/

#ifndef __INTREPID2_HGRAD_TET_Cn_FEM_ORTH_HPP__
#define __INTREPID2_HGRAD_TET_Cn_FEM_ORTH_HPP__

#include "Intrepid2_Basis.hpp"

namespace Intrepid2 {
  
/** \class  Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
    \brief  Implementation of the default H(grad)-compatible orthogonal basis of
            arbitrary degree on tetrahedron.
  
    \remarks
    \li   All degrees of freedom are considered to be internal (ie not assembled)
  */
  
namespace Impl {

/** \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
 */
template<typename OutputViewType,
typename inputViewType, 
typename workViewType,
bool hasDeriv, ordinal_type n>
struct OrthPolynomialTet {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(       OutputViewType output,
            const inputViewType input,
                  workViewType work,
            const ordinal_type p );
};

/** \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
 */
template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
struct OrthPolynomialTet<OutputViewType,inputViewType,workViewType,hasDeriv,0> {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(       OutputViewType output,
            const inputViewType input,
                  workViewType work,
            const ordinal_type p );
};

/** \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
 */
template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
struct OrthPolynomialTet<OutputViewType,inputViewType,workViewType,hasDeriv,1> {
  KOKKOS_INLINE_FUNCTION
  static void
  generate(   OutputViewType output,
      const inputViewType input,
            workViewType work,
      const ordinal_type p );
};

/**
  \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
*/
class Basis_HGRAD_TET_Cn_FEM_ORTH {
public:

  /**
    \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
  */
  template<EOperator opType>
  struct Serial {
    template<typename OutputViewType,
    typename inputViewType,
    typename workViewType>
    KOKKOS_INLINE_FUNCTION
    static void
    getValues(       OutputViewType output,
               const inputViewType  input,
                     workViewType work,
               const ordinal_type   order);
  };

  template<typename DeviceType, ordinal_type numPtsPerEval,
  typename outputValueValueType, class ...outputValueProperties,
  typename inputPointValueType,  class ...inputPointProperties>
  static void
  getValues( const typename DeviceType::execution_space& space,
                   Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const ordinal_type order,
             const EOperator operatorType );

  /**
    \brief See Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH
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
            ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH::Functor) operator is not supported");

      }
      }
    }
  };

};

}

template<typename DeviceType = void,
    typename outputValueType = double,
    typename pointValueType = double>
class Basis_HGRAD_TET_Cn_FEM_ORTH
    : public Basis<DeviceType,outputValueType,pointValueType> {
    public:
  typedef double value_type;

  using BasisBase = Basis<DeviceType,outputValueType,pointValueType>;
  using typename BasisBase::ExecutionSpace;

  using typename BasisBase::OrdinalTypeArray1DHost;
  using typename BasisBase::OrdinalTypeArray2DHost;
  using typename BasisBase::OrdinalTypeArray3DHost;

  using typename BasisBase::OutputViewType;
  using typename BasisBase::PointViewType ;
  using typename BasisBase::ScalarViewType;

  /** \brief  Constructor.
   */
  Basis_HGRAD_TET_Cn_FEM_ORTH( const ordinal_type order );
  
  using BasisBase::getValues;

  virtual
  void
  getValues( const ExecutionSpace& space,
                   OutputViewType  outputValues,
             const PointViewType   inputPoints,
             const EOperator       operatorType = OPERATOR_VALUE ) const override {
    #ifdef HAVE_INTREPID2_DEBUG
          Intrepid2::getValues_HGRAD_Args(outputValues,
                                          inputPoints,
                                          operatorType,
                                          this->getBaseCellTopology(),
                                          this->getCardinality() );
    #endif
    constexpr ordinal_type numPtsPerEval = Parameters::MaxNumPtsPerBasisEval;
    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
    getValues<DeviceType,numPtsPerEval>(space, 
                                        outputValues,
                                        inputPoints,
                                        this->getDegree(),
                                        operatorType);
  }
};

}// namespace Intrepid2

#include "Intrepid2_HGRAD_TET_Cn_FEM_ORTHDef.hpp"

#endif

