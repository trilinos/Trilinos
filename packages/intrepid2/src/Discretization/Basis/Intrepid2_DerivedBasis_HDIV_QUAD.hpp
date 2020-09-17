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

/** \file   Intrepid2_DerivedBasis_HIV_QUAD.hpp
    \brief  Implementation of H(div) basis on the quadrilateral that is templated on H(vol) and H(grad) on the line.
    \author Created by N.V. Roberts.
 
 This class constructs the H(div) space as the direct sum of two families of tensor-product bases on the quad:
 - family 1: H(vol)  x  H(grad), placed in the y component of vector output
 - family 2: H(grad) x  H(vol),  placed in the x component of vector output
 */

#ifndef Intrepid2_DerivedBasis_HIV_QUAD_h
#define Intrepid2_DerivedBasis_HIV_QUAD_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family1_QUAD
  : public Basis_TensorBasis<HVOL_LINE, HGRAD_LINE>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis = Basis_TensorBasis<LineHVolBasis,LineGradBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family1_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType)
    :
    TensorBasis(LineHVolBasis(polyOrder_x-1,pointType),
                LineGradBasis(polyOrder_y,pointType))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
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
      // ESEAS implements H(div) as rotated H(curl), which involves weighting family1 with -1.
      // We follow that here to simplify verification tests that involve ESEAS.
      const double weight = -1.0;
      
      Intrepid2::EOperator op1, op2;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        // family 1 goes in the y component; 0 in the x component
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        
        this->TensorBasis::getValues(outputValuesComponent2,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints, weight);
        // place 0 in the y component
        Kokkos::deep_copy(outputValuesComponent1,0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_DIV)
      {
        // family 1 gets a d/dy applied to the second (nonzero) vector component
        // since this is H(VOL)(x) * H(GRAD)(y), this amounts to taking the derivative in the second tensorial component
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_GRAD;
        
        this->TensorBasis::getValues(outputValues,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints, weight);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }

    /** \brief  Fills in coefficients of degrees of freedom for Lagrangian basis on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs have shape (F,D), field dimension matches the cardinality of the basis, and D is the
     basis dimension.

     Degrees of freedom coefficients are such that
     \phi_i(dofCoords_(j)) \cdot dofCoeffs_(j)  = \delta_ij,
     where \phi_i are the basis and \delta_ij the Kronecker delta.
     Note that getDofCoeffs() is supported only for Lagrangian bases.
     */
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),1);
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis::getDofCoeffs(dofCoeffs2);
      //multiply by weight = -1.0
      using ExecutionSpace = typename TensorBasis::ExecutionSpace;
      Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>(0, dofCoeffs2.extent(0)),
          KOKKOS_LAMBDA (const int i){ dofCoeffs2(i) *= -1.0; });
    }

  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family2_QUAD
  : public Basis_TensorBasis<HGRAD_LINE, HVOL_LINE>
  {
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis = Basis_TensorBasis<LineGradBasis,LineHVolBasis>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family2_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType)
    :
    TensorBasis(LineGradBasis(polyOrder_x,pointType),
                LineHVolBasis(polyOrder_y-1,pointType))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
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
        
        // family 2 goes in the x component; 0 in the x component
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
        
        this->TensorBasis::getValues(outputValuesComponent1,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        // place 0 in the y component
        Kokkos::deep_copy(outputValuesComponent2, 0.0);
      }
      else if (operatorType == Intrepid2::OPERATOR_DIV)
      {
        // family 2 gets a d/dx applied to the first (nonzero) vector component
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

    /** \brief  Fills in coefficients of degrees of freedom for Lagrangian basis on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs have shape (F,D), field dimension matches the cardinality of the basis, and D is the
     basis dimension.

     Degrees of freedom coefficients are such that
     \phi_i(dofCoords_(j)) \cdot dofCoeffs_(j)  = \delta_ij,
     where \phi_i are the basis and \delta_ij the Kronecker delta.
     Note that getDofCoeffs() is supported only for Lagrangian bases.
     */
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),1);
      this->TensorBasis::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs2,0.0);
    }

  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_QUAD
  : public Basis_DirectSumBasis<Basis_Derived_HDIV_Family1_QUAD<HGRAD_LINE, HVOL_LINE>,
                                Basis_Derived_HDIV_Family2_QUAD<HGRAD_LINE, HVOL_LINE> >
  {
    using Family1 = Basis_Derived_HDIV_Family1_QUAD<HGRAD_LINE, HVOL_LINE>;
    using Family2 = Basis_Derived_HDIV_Family2_QUAD<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<Family1,Family2>;

    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;

  protected:
    std::string name_;
    ordinal_type order_x_;
    ordinal_type order_y_;
    EPointType pointType_;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Family1(polyOrder_x, polyOrder_y, pointType),
                   Family2(polyOrder_x, polyOrder_y, pointType))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;

      std::ostringstream basisName;
      basisName << "HDIV_QUAD (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_x_ = polyOrder_x;
      order_y_ = polyOrder_y;
      pointType_ = pointType;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_QUAD(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HDIV_QUAD(polyOrder, polyOrder, pointType) {}

    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const {
      return (this->getDofCount(1,0) > 0); //if it has side DOFs, than it needs orientations
    }

    /** \brief  Returns basis name

     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected along normal to the subCell.

        TODO: test this method when different orders are used in different directions
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<ExecutionSpace, OutputValueType, PointValueType>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1) {
        switch(subCellOrd) {
        case 0:
        case 2:
          return Teuchos::rcp( new HVOL_LINE(order_x_-1, pointType_) );
        case 1:
        case 3:
          return Teuchos::rcp( new HVOL_LINE(order_y_-1, pointType_) );
        }
      }

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HIV_QUAD_h */
