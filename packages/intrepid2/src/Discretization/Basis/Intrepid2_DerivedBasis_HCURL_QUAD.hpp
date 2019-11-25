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

/** \file   Intrepid2_DerivedBasis_HCURL_QUAD.hpp
    \brief  Implementation of H(curl) basis on the quadrilateral that is templated on H(vol) and H(grad) on the line.
    \author Created by N.V. Roberts.
 
 This class constructs the H(curl) space as the direct sum of two families of tensor-product bases on the quad:
 - family 1: H(vol)  x  H(grad), placed in the x component of vector output
 - family 2: H(grad) x  H(vol),  placed in the y component of vector output
 */

#ifndef Intrepid2_DerivedBasis_HCURL_QUAD_h
#define Intrepid2_DerivedBasis_HCURL_QUAD_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family1_QUAD
  : public Basis_TensorBasis<HVOL_LINE, HGRAD_LINE>
  {
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis = Basis_TensorBasis<LineHVolBasis,LineGradBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
     */
    Basis_Derived_HCURL_Family1_QUAD(int polyOrder_x, int polyOrder_y)
    :
    TensorBasis(LineHVolBasis(polyOrder_x-1),
                LineGradBasis(polyOrder_y))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
    }
    
    using TensorBasis::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const
    {
      Intrepid2::EOperator op1, op2;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        // family 1 goes in the x component; 0 in the y component
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        
        this->TensorBasis::getValues(outputValuesComponent1,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        // place 0 in the y component
        Kokkos::deep_copy(outputValuesComponent2,0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_CURL)
      {
        // family 1 gets a -d/dy applied to the first (nonzero) vector component
        // since this is H(VOL)(x) * H(GRAD)(y), this amounts to taking the derivative in the second tensorial component
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD;
        
        double weight = -1.0; // the minus sign in front of d/dy
        this->TensorBasis::getValues(outputValues,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints, weight);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family2_QUAD
  : public Basis_TensorBasis<HGRAD_LINE, HVOL_LINE>
  {
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis = Basis_TensorBasis<LineGradBasis,LineHVolBasis>;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
     */
    Basis_Derived_HCURL_Family2_QUAD(int polyOrder_x, int polyOrder_y)
    :
    TensorBasis(LineGradBasis(polyOrder_x),
                LineHVolBasis(polyOrder_y-1))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
    }
    
    using TensorBasis::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const
    {
      Intrepid2::EOperator op1, op2;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        // family 2 goes in the y component; 0 in the x component
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        
        // place 0 in the x component
        Kokkos::deep_copy(outputValuesComponent1, 0.0);
        this->TensorBasis::getValues(outputValuesComponent2,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        
      }
      else if (operatorType == Intrepid2::OPERATOR_CURL)
      {
        // family 2 gets a d/dx applied to the second (nonzero) vector component
        // since this is H(GRAD)(x) * H(VOL)(y), this amounts to taking the derivative in the first tensorial component
        op1 = Intrepid2::OPERATOR_GRAD;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        this->TensorBasis::getValues(outputValues,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_QUAD
  : public Basis_DirectSumBasis<Basis_Derived_HCURL_Family1_QUAD<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HCURL_Family2_QUAD<HGRAD_LINE, HVOL_LINE> >
  {
    using Family1 = Basis_Derived_HCURL_Family1_QUAD<HGRAD_LINE, HVOL_LINE>;
    using Family2 = Basis_Derived_HCURL_Family2_QUAD<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family1,Family2>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
     */
    Basis_Derived_HCURL_QUAD(int polyOrder_x, int polyOrder_y)
    :
    DirectSumBasis(Family1(polyOrder_x, polyOrder_y),
                   Family2(polyOrder_x, polyOrder_y)) {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
     */
    Basis_Derived_HCURL_QUAD(int polyOrder) : Basis_Derived_HCURL_QUAD(polyOrder, polyOrder) {}
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HCURL_QUAD_h */
