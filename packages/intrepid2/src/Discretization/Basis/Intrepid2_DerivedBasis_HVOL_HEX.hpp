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

/** \file   Intrepid2_DerivedBasis_HVOL_HEX.hpp
    \brief  Implementation of H(vol) basis on the hexahedron that is templated on H(vol) on the line.
    \author Created by N.V. Roberts.
 
 H(vol) on the hexahedron is defined (and implemented) as H(vol) on the quad times H(vol) on the line.
 
 For consistency with the H(div) and H(curl) implementations on the hexahedron, we may in the future define this as
 
 H(vol) x H(vol) x H(vol)
 
 where each component is defined on the line.
 */

#ifndef Intrepid2_DerivedBasis_HVOL_HEX_h
#define Intrepid2_DerivedBasis_HVOL_HEX_h

#include "Intrepid2_TensorBasis.hpp"
#include "Intrepid2_DerivedBasis_HVOL_QUAD.hpp"

namespace Intrepid2
{
  template<class HVOL_LINE>
  class Basis_Derived_HVOL_HEX
  :
  public Basis_TensorBasis<Intrepid2::Basis_Derived_HVOL_QUAD<HVOL_LINE>,
                           HVOL_LINE>
  // TODO: make this a subclass of TensorBasis3 instead, following what we've done for H(curl) and H(div)
  {
  public:
    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using LineBasis = HVOL_LINE;
    using QuadBasis = Intrepid2::Basis_Derived_HVOL_QUAD<HVOL_LINE>;
    using TensorBasis = Basis_TensorBasis<QuadBasis,LineBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HVOL_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis(QuadBasis(polyOrder_x,polyOrder_y),
                LineBasis(polyOrder_z))
    {
      this->functionSpace_ = FUNCTION_SPACE_HVOL;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
     */
    Basis_Derived_HVOL_HEX(int polyOrder) : Basis_Derived_HVOL_HEX(polyOrder, polyOrder, polyOrder) {}
    
    using TensorBasis::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x,y dimensions
        \param [in] inputPoints2 - input points in the z dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType inputPoints1, const PointViewType inputPoints2,
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
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HVOL_HEX_h */
