// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_LINE_Cn_FEM_JACOBI.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_HPP__
#define __INTREPID2_HGRAD_LINE_CN_FEM_JACOBI_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Polylib.hpp"

namespace Intrepid2 {
  
  /** \class  Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI
      \brief  Implementation of the locally H(grad)-compatible FEM basis of variable order
      on the [-1,1] reference line cell, using Jacobi polynomials. 
    
      Implements Jacobi basis of variable order \f$n\f$ on
      the reference [-1,1] line cell. Jacobi polynomials depend on three parameters
      \f$ \alpha \f$, \f$ \beta \f$, and \f$ n \f$ and are defined via the so-called
      Gamma function by
      \f[
      P_n^{(\alpha,\beta)} (z) = 
      \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)}
      \sum_{m=0}^n {n\choose m}
      \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m
      \f]
      The basis has cardinality \f$n+1\f$ and spans a COMPLETE linear polynomial space.
      Basis functions are dual to a unisolvent set of degrees of freedom (DoF) enumerated as follows:
    
      <table>
      <tr>
      <th rowspan="2"> Basis order <th colspan="4"> DoF tag table <th rowspan="2"> DoF definition
      </tr>
      <tr>
      <th> subc dim <th> subc ordinal <th> subc DoF tag <th> subc num DoFs
      </tr>

      <tr align="center"> <td> 0 <td> 1 <td> 0 <td> 0   <td> 1
      <td align="left"> \f$ P_0^{(\alpha,\beta)} \f$ </tr>
      <tr align="center"> <td> 1 <td> 1 <td> 0 <td> 0-1 <td> 2
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)} \f$ </tr>
      <tr align="center"> <td> 2 <td> 1 <td> 0 <td> 0-2 <td> 3
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, P_2^{(\alpha,\beta)} \f$ </tr>
      <tr align="center"> <td> 3 <td> 1 <td> 0 <td> 0-3 <td> 4
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, ..., P_3^{(\alpha,\beta)} \f$ </tr>
      <tr align="center"> <td> ... <td> 1 <td> 0 <td> ... <td> ...
      <td align="left"> ... </tr>
      <tr align="center"> <td> n <td> 1 <td> 0 <td> 0-n <td> n+1
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, ..., P_n^{(\alpha,\beta)} \f$ </tr>
      </table>
  
      For example, for Legendre polynomials (\f$\alpha=\beta=0\f$), the first 11 bases are given by

      <table>
      <tr>
      <th rowspan="2"> Basis order <th colspan="4"> DoF tag table <th rowspan="2"> DoF definition
      </tr>
      <tr>
      <th> subc dim <th> subc ordinal <th> subc DoF tag <th> subc num DoFs
      </tr>

      <tr align="center"> <td> 0 <td> 1 <td> 0 <td> 0   <td> 1
      <td align="left"> \f$ 1 \f$ </tr>
      <tr align="center"> <td> 1 <td> 1 <td> 0 <td> 0-1 <td> 2
      <td align="char" char=":"> and: \f$ x \f$ </tr>
      <tr align="center"> <td> 2 <td> 1 <td> 0 <td> 0-2 <td> 3
      <td align="char" char=":"> and: \f$ \frac{1}{2} (3x^2-1) \f$ </tr>
      <tr align="center"> <td> 3 <td> 1 <td> 0 <td> 0-3 <td> 4
      <td align="char" char=":"> and: \f$ \frac{1}{2} (5x^3-3x) \f$ </tr>
      <tr align="center"> <td> 4 <td> 1 <td> 0 <td> 0-4 <td> 5
      <td align="char" char=":"> and: \f$ \frac{1}{8} (35x^4-30x^2+3) \f$ </tr>
      <tr align="center"> <td> 5 <td> 1 <td> 0 <td> 0-5 <td> 6
      <td align="char" char=":"> and: \f$ \frac{1}{8} (63x^5-70x^3+15x) \f$ </tr>
      <tr align="center"> <td> 6 <td> 1 <td> 0 <td> 0-6 <td> 7
      <td align="char" char=":"> and: \f$ \frac{1}{16} (231x^6-315x^4+105x^2-5) \f$ </tr>
      <tr align="center"> <td> 7 <td> 1 <td> 0 <td> 0-7 <td> 8
      <td align="char" char=":"> and: \f$ \frac{1}{16} (429x^7-693x^5+315x^3-35x) \f$ </tr>
      <tr align="center"> <td> 8 <td> 1 <td> 0 <td> 0-8 <td> 9
      <td align="char" char=":"> and: \f$ \frac{1}{128} (6435x^8-12012x^6+6930x^4-1260x^2+35) \f$ </tr>
      <tr align="center"> <td> 9 <td> 1 <td> 0 <td> 0-9 <td> 10
      <td align="char" char=":"> and: \f$ \frac{1}{128} (12155x^9-25740x^7+18018x^5-4620x^3+315x) \f$ </tr>
      <tr align="center"> <td>10 <td> 1 <td> 0 <td> 0-10<td> 11
      <td align="char" char=":"> and: \f$ \frac{1}{128} (46189x^{10}-109395x^8+90090x^6-30030x^4+3465x^2-63) \f$ </tr>
      </table>
  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI
    */
    class Basis_HGRAD_LINE_Cn_FEM_JACOBI {
    public:

      /**
        \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI
      */
      template<EOperator opType>
      struct Serial {
        template<typename OutputViewType,
                 typename inputViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(       OutputViewType output,
                   const inputViewType  input,
                   const ordinal_type   order,
                   const double         alpha,
                   const double         beta,
                   const ordinal_type   opDn = 0 );
      };
      
      template<typename DeviceType, ordinal_type numPtsPerEval,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      static void
      getValues(  const typename DeviceType::execution_space& space,
                        Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                  const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                  const ordinal_type order,
                  const double alpha,
                  const double beta,
                  const EOperator operatorType );
      
      /**
        \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI
      */
      template<typename outputValueViewType,
               typename inputPointViewType,
               EOperator opType,
               ordinal_type numPtsEval>
      struct Functor {
              outputValueViewType _outputValues;
        const inputPointViewType  _inputPoints;
        const ordinal_type        _order;
        const double              _alpha, _beta;
        const ordinal_type        _opDn;
        
        KOKKOS_INLINE_FUNCTION
        Functor(       outputValueViewType outputValues_,
                       inputPointViewType  inputPoints_,
                 const ordinal_type        order_,
                 const double              alpha_,
                 const double              beta_,
                 const ordinal_type        opDn_ = 0 )
          : _outputValues(outputValues_), _inputPoints(inputPoints_), 
            _order(order_), _alpha(alpha_), _beta(beta_), _opDn(opDn_) {}
        
        KOKKOS_INLINE_FUNCTION
        void operator()(const size_type iter) const {
          const auto ptBegin = Util<ordinal_type>::min(iter*numPtsEval,    _inputPoints.extent(0));
          const auto ptEnd   = Util<ordinal_type>::min(ptBegin+numPtsEval, _inputPoints.extent(0));

          const auto ptRange = Kokkos::pair<ordinal_type,ordinal_type>(ptBegin, ptEnd);
          const auto input   = Kokkos::subview( _inputPoints, ptRange, Kokkos::ALL() );

          if (input.extent(0)) {
            switch (opType) {
            case OPERATOR_VALUE : {
              auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange );
              Serial<opType>::getValues( output, input, _order, _alpha, _beta );
              break;
            }
            case OPERATOR_GRAD : {
              auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
              Serial<opType>::getValues( output, input, _order, _alpha, _beta );
              break;
            }
            case OPERATOR_Dn: {
              auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
              Serial<opType>::getValues( output, input, _order, _alpha, _beta, _opDn );
              break;
            }
            default: {
              INTREPID2_TEST_FOR_ABORT( true, 
                                        ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM_JACOBI::Functor) operator is not supported");
              
            }
            }
          }
        }
      };
    };
  }

  template<typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_LINE_Cn_FEM_JACOBI
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
    Basis_HGRAD_LINE_Cn_FEM_JACOBI( const ordinal_type order, 
                                    const double alpha = 0, 
                                    const double beta = 0 );
   
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
      constexpr ordinal_type numPtsPerEval = 1;
      Impl::Basis_HGRAD_LINE_Cn_FEM_JACOBI::
        getValues<DeviceType,numPtsPerEval>(space,
                                            outputValues, 
                                            inputPoints, 
                                            this->getDegree(),
                                            this->alpha_, 
                                            this->beta_,
                                            operatorType);
    }
    
  private:
    double alpha_, beta_;

  };

}

#include "Intrepid2_HGRAD_LINE_Cn_FEM_JACOBIDef.hpp"

#endif
