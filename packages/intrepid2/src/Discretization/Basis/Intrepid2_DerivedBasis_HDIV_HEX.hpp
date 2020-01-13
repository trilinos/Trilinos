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

/** \file   Intrepid2_DerivedBasis_HDIV_HEX.hpp
    \brief  Implementation of H(div) basis on the hexahedron that is templated on H(vol) and H(grad) on the line.
    \author Created by N.V. Roberts.
 
  This class constructs the H(div) space as the direct sum of three families of tensor-product bases on the hexahedron:
 - family 1: H(grad) x H(vol)  x H(vol),  placed in the x component of vector output
 - family 2: H(vol)  x H(grad) x H(vol),  placed in the y component of vector output
 - family 3: H(vol)  x H(vol)  x H(grad), placed in the z component of vector output
 
 Note that for compatibility with the ordering in the ESEAS libray, the direct sum is Family3 + Family1 + Family2;
 i.e., the Family3 members come first in the direct sum basis enumeration, then Family1, then Family2.
 */

#ifndef Intrepid2_DerivedBasis_HDIV_HEX_h
#define Intrepid2_DerivedBasis_HDIV_HEX_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family1_HEX
  :
  public Basis_TensorBasis3<HGRAD_LINE, HVOL_LINE, HVOL_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineGradBasis, LineHVolBasis, LineHVolBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HDIV_Family1_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineGradBasis(polyOrder_x  ),
                 LineHVolBasis(polyOrder_y-1),
                 LineHVolBasis(polyOrder_z-1))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    using TensorBasis3::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis3)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] inputPoints3 - input points in the z dimension
        \param [in] tensorPoints - if true, inputPoints1, inputPoints2, and inputPoints3 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of the three inputPoints containers).  If false, each inputPoints container should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2, const PointViewType inputPoints3,
                           bool tensorPoints) const
    {
      Intrepid2::EOperator op1, op2, op3;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        
        // family 1 goes in the x component; 0 in the y and z components
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent23 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::make_pair(1,3));
        
        this->TensorBasis3::getValues(outputValuesComponent1,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints);
        // place 0 in the y and z components
        Kokkos::deep_copy(outputValuesComponent23,0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_DIV)
      {
        // family 1 is nonzero in the x component, so the div is d/dx of the first component
        // outputValues is scalar, so no need to take subviews
        
        op1 = Intrepid2::OPERATOR_GRAD;  // d/dx
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        
        double weight = 1.0; // the plus sign in front of d/dx
        this->TensorBasis3::getValues(outputValues,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family2_HEX
  :
  public Basis_TensorBasis3<HVOL_LINE, HGRAD_LINE, HVOL_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineHVolBasis, LineGradBasis, LineHVolBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HDIV_Family2_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineHVolBasis(polyOrder_x-1),
                 LineGradBasis(polyOrder_y  ),
                 LineHVolBasis(polyOrder_z-1))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    using TensorBasis3::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis3)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] inputPoints3 - input points in the z dimension
        \param [in] tensorPoints - if true, inputPoints1, inputPoints2, and inputPoints3 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of the three inputPoints containers).  If false, each inputPoints container should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2, const PointViewType inputPoints3,
                           bool tensorPoints) const
    {
      Intrepid2::EOperator op1, op2, op3;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        
        // family 2 goes in the y component; 0 in the x and z components
        auto outputValuesComponent_x = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent_y = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        auto outputValuesComponent_z = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
       
        // 0 in x component
        Kokkos::deep_copy(outputValuesComponent_x,0.0);
        
        double weight = 1.0;
        this->TensorBasis3::getValues(outputValuesComponent_y,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
        
        // 0 in z component
        Kokkos::deep_copy(outputValuesComponent_z,0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_DIV)
      {
        // family 2 is nonzero in the y component, so the div is d/dy of the second component
        // ESEAS has -1 weight here; we follow that, to simplify verification
        // outputValues is scalar, so no need to take subviews
        
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD;  // d/dy
        op3 = Intrepid2::OPERATOR_VALUE;
        
        double weight = 1.0;
        this->TensorBasis3::getValues(outputValues,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family3_HEX
  : public Basis_TensorBasis3<HVOL_LINE, HVOL_LINE, HGRAD_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineHVolBasis, LineHVolBasis, LineGradBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HDIV_Family3_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineHVolBasis(polyOrder_x-1),
                 LineHVolBasis(polyOrder_y-1),
                 LineGradBasis(polyOrder_z  ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    using TensorBasis3::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis3)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] inputPoints3 - input points in the z dimension
        \param [in] tensorPoints - if true, inputPoints1, inputPoints2, and inputPoints3 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of the three inputPoints containers).  If false, each inputPoints container should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2, const PointViewType inputPoints3,
                           bool tensorPoints) const
    {
      Intrepid2::EOperator op1, op2, op3;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        
        // family 3 goes in the z component; 0 in the x and y components
        auto outputValuesComponent_xy = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::make_pair(0,1));
        auto outputValuesComponent_z = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
       
        // 0 in x and y components
        Kokkos::deep_copy(outputValuesComponent_xy,0.0);
        
        // z component
        this->TensorBasis3::getValues(outputValuesComponent_z,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints);
      }
      else if (operatorType == Intrepid2::OPERATOR_DIV)
      {
        // family 3 is nonzero in the z component, so the div is d/dz of the third component
        // outputValues is scalar, so no need to take subviews
        
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_GRAD;  // d/dz
        
        double weight = 1.0; // the plus sign in front of d/dz
        this->TensorBasis3::getValues(outputValues,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
  
  // ESEAS numbers its H(div) families differently, with the nonzero going in z, x, y for I,II,III.
  // To allow our interior orderings to match that of ESEAS, we put the direct sum in the same order as ESEAS,
  // which is to say that we go 3,1,2.
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family3_Family1_HEX
  : public Basis_DirectSumBasis<Basis_Derived_HDIV_Family3_HEX<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HDIV_Family1_HEX<HGRAD_LINE, HVOL_LINE>>
  {
    using Family3 = Basis_Derived_HDIV_Family3_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family1 = Basis_Derived_HDIV_Family1_HEX<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family3,Family1>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HDIV_Family3_Family1_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    DirectSumBasis(Family3(polyOrder_x, polyOrder_y, polyOrder_z),
                   Family1(polyOrder_x, polyOrder_y, polyOrder_z)) {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_HEX
  : public Basis_DirectSumBasis<Basis_Derived_HDIV_Family3_Family1_HEX<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HDIV_Family2_HEX<HGRAD_LINE, HVOL_LINE> >
  {
    using Family31 = Basis_Derived_HDIV_Family3_Family1_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family2  = Basis_Derived_HDIV_Family2_HEX        <HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family31,Family2>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HDIV_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    DirectSumBasis(Family31(polyOrder_x, polyOrder_y, polyOrder_z),
                   Family2 (polyOrder_x, polyOrder_y, polyOrder_z)) {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
     */
    Basis_Derived_HDIV_HEX(int polyOrder) : Basis_Derived_HDIV_HEX(polyOrder, polyOrder, polyOrder) {}
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HDIV_HEX_h */
