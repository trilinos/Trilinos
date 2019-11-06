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

/** \file   Intrepid2_DerivedBasis_HCURL_HEX.hpp
    \brief  Implementation of H(curl) basis on the hexahedron that is templated on H(vol) and H(grad) on the line.
    \author Created by N.V. Roberts.
 
 This class constructs the H(curl) space as the direct sum of three families of tensor-product bases on the hexahedron:
 - family 1: H(vol)  x  H(grad) x  H(grad), placed in the x component of vector output
 - family 2: H(grad) x  H(vol)  x  H(grad), placed in the y component of vector output
 - family 3: H(grad) x  H(grad) x  H(vol),  placed in the z component of vector output
 */

#ifndef Intrepid2_DerivedBasis_HCURL_HEX_h
#define Intrepid2_DerivedBasis_HCURL_HEX_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family1_HEX
  : public Basis_TensorBasis3<HVOL_LINE,
                              HGRAD_LINE,
                              HGRAD_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;

    using LineGradBasis = HGRAD_LINE;
    using LineVolBasis  = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineVolBasis, LineGradBasis, LineGradBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HCURL_Family1_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineVolBasis (polyOrder_x-1),
                 LineGradBasis(polyOrder_y  ),
                 LineGradBasis(polyOrder_z  ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
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
      else if (operatorType == Intrepid2::OPERATOR_CURL)
      {
        // family 1 is nonzero in the x component, so the curl is d/dz placed in the y component, and -d/dy placed in the z component.
        auto outputValuesComponent_x = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent_y = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        auto outputValuesComponent_z = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        // x component is zero
        Kokkos::deep_copy(outputValuesComponent_x, 0.0);
        
        // y component is d/dz
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_GRAD; // d/dz
        
        double weight = 1.0; // the plus sign in front of d/dz
        this->TensorBasis3::getValues(outputValuesComponent_y,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
        
        // z component is -d/dy
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD; // d/dy
        op3 = Intrepid2::OPERATOR_VALUE;
        weight = -1.0; // the -1 weight on d/dy
        this->TensorBasis3::getValues(outputValuesComponent_z,
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
  class Basis_Derived_HCURL_Family2_HEX
  : public Basis_TensorBasis3<HGRAD_LINE,
                              HVOL_LINE,
                              HGRAD_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineVolBasis  = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineGradBasis, LineVolBasis, LineGradBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HCURL_Family2_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineGradBasis(polyOrder_x),
                 LineVolBasis (polyOrder_y-1),
                 LineGradBasis(polyOrder_z))
    {}
    
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
        
        // place 0 in the x component
        Kokkos::deep_copy(outputValuesComponent_x,0.0);
        // evaluate y component
        this->TensorBasis3::getValues(outputValuesComponent_y,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints);
        // place 0 in the z component
        Kokkos::deep_copy(outputValuesComponent_z,0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_CURL)
      {
        // family 2 is nonzero in the y component, so the curl is -d/dz placed in the x component, and d/dx placed in the z component.
        auto outputValuesComponent_x = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent_y = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        auto outputValuesComponent_z = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        // x component is -d/dz
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_GRAD; // d/dz
        
        double weight = -1.0; // the minus sign in front of d/dz
        this->TensorBasis3::getValues(outputValuesComponent_x,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
        
        // y component is zero
        Kokkos::deep_copy(outputValuesComponent_y, 0.0);
        
        // z component is d/dx
        op1 = Intrepid2::OPERATOR_GRAD;  // d/dx
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        weight = 1.0; // the weight on d/dx
        this->TensorBasis3::getValues(outputValuesComponent_z,
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
  class Basis_Derived_HCURL_Family3_HEX
  : public Basis_TensorBasis3<HGRAD_LINE,
                              HGRAD_LINE,
                              HVOL_LINE>
  {
    using OutputViewType = typename HGRAD_LINE::outputViewType;
    using PointViewType  = typename HGRAD_LINE::pointViewType ;
    using ScalarViewType = typename HGRAD_LINE::scalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineVolBasis  = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<LineGradBasis, LineGradBasis, LineVolBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HCURL_Family3_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    TensorBasis3(LineGradBasis(polyOrder_x  ),
                 LineGradBasis(polyOrder_y  ),
                 LineVolBasis (polyOrder_z-1))
    {}
    
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
        auto outputValuesComponent_xy = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::make_pair(0,2));
        auto outputValuesComponent_z  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        // place 0 in the x and y components
        Kokkos::deep_copy(outputValuesComponent_xy,0.0);
        // evaluate z component
        this->TensorBasis3::getValues(outputValuesComponent_z,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints);
      }
      else if (operatorType == Intrepid2::OPERATOR_CURL)
      {
        // family 3 is nonzero in the z component, so the curl is d/dy placed in the x component, and -d/dx placed in the y component.
        auto outputValuesComponent_x = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent_y = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        auto outputValuesComponent_z = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        // x component is d/dy
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD; // d/dy
        op3 = Intrepid2::OPERATOR_VALUE;
        
        double weight = 1.0; // the sign in front of d/dy
        this->TensorBasis3::getValues(outputValuesComponent_x,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
        // y component is -d/dx
        op1 = Intrepid2::OPERATOR_GRAD; // d/dx
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        weight = -1.0; // the weight on d/dx
        this->TensorBasis3::getValues(outputValuesComponent_y,
                                      inputPoints1, op1,
                                      inputPoints2, op2,
                                      inputPoints3, op3, tensorPoints, weight);
        
        // z component is zero
        Kokkos::deep_copy(outputValuesComponent_z, 0.0);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family1_Family2_HEX
  : public Basis_DirectSumBasis<Basis_Derived_HCURL_Family1_HEX<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HCURL_Family2_HEX<HGRAD_LINE, HVOL_LINE> >
  {
    using Family1 = Basis_Derived_HCURL_Family1_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family2 = Basis_Derived_HCURL_Family2_HEX<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family1,Family2>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HCURL_Family1_Family2_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    DirectSumBasis(Family1(polyOrder_x, polyOrder_y, polyOrder_z),
                   Family2(polyOrder_x, polyOrder_y, polyOrder_z)) {}
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_HEX
  : public Basis_DirectSumBasis<Basis_Derived_HCURL_Family1_Family2_HEX<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HCURL_Family3_HEX<HGRAD_LINE, HVOL_LINE> >
  {
    using Family12 = Basis_Derived_HCURL_Family1_Family2_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family3  = Basis_Derived_HCURL_Family3_HEX        <HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family12,Family3>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
     */
    Basis_Derived_HCURL_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z)
    :
    DirectSumBasis(Family12(polyOrder_x, polyOrder_y, polyOrder_z),
                   Family3 (polyOrder_x, polyOrder_y, polyOrder_z)) {}
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
     */
    Basis_Derived_HCURL_HEX(int polyOrder) : Basis_Derived_HCURL_HEX(polyOrder, polyOrder, polyOrder) {}
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HCURL_HEX_h */
