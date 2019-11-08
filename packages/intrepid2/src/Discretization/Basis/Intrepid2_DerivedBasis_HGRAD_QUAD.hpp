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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HGRAD_QUAD.hpp
    \brief  Implementation of H(grad) basis on the quadrilateral that is templated on H(grad) on the line.
    \author Created by N.V. Roberts.
 
 H(grad) on the quadrilateral is defined (and implemented) as H(grad) x H(grad) on the line.
 */

#ifndef Intrepid2_DerivedBasis_HGRAD_QUAD_h
#define Intrepid2_DerivedBasis_HGRAD_QUAD_h

#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE>
  class Basis_Derived_HGRAD_QUAD
  : public Basis_TensorBasis<HGRAD_LINE, HGRAD_LINE>
  {
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;
    
    using LineBasis = HGRAD_LINE;
    using TensorBasis = Basis_TensorBasis<LineBasis,LineBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
     */
    Basis_Derived_HGRAD_QUAD(int polyOrder_x, int polyOrder_y)
    :
    TensorBasis(LineBasis(polyOrder_x),
                LineBasis(polyOrder_y))
    {
      this->functionSpace_ = FUNCTION_SPACE_HGRAD;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in both dimensions.
     */
    Basis_Derived_HGRAD_QUAD(int polyOrder) : Basis_Derived_HGRAD_QUAD(polyOrder,polyOrder) {}
    
    using Basis<ExecutionSpace,OutputValueType,PointValueType>::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType  inputPoints2,
                           bool tensorPoints) const override
    {
      Intrepid2::EOperator op1, op2;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        this->TensorBasis::getValues(outputValues,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
      }
      else if (operatorType == Intrepid2::OPERATOR_GRAD)
      {
        // to evaluate gradient, we actually need both OP_VALUE and OP_GRAD (thanks to product rule)
        // for 1D line x line, we will put derivative * value in first component, and value * derivative in second
        
        // outputValues1 and outputValues2 are computed by basis1 and basis2 -- these are tensorial components
        // outputValuesComponent1 is a slice of the final output container (similarly, outputValuesComponent2)
        // when the component basis is 1D, it expects not to have a "dimension" component in the output container
        // the int argument in the dimension component creates a subview that skips the dimension component; the std::pair argument retains it
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        
        // compute first component -- derivative happens in x, and value taken in y
        op1 = Intrepid2::OPERATOR_GRAD;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        this->TensorBasis::getValues(outputValuesComponent1,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        
        // second component -- value in x, derivative in y
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD;
        
        this->TensorBasis::getValues(outputValuesComponent2,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_QUAD_h */
