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

/** \file   Intrepid2_DerivedBasis_HGRAD_HEX.hpp
    \brief  Implementation of H(grad) basis on the hexahedron that is templated on H(grad) on the line.
    \author Created by N.V. Roberts.
 
 H(grad) on the hexahedron is defined (and implemented) as H(grad) on the quad times H(grad) on the line.
 
 For consistency with the H(div) and H(curl) implementations on the hexahedron, we may in the future define this as
 
 H(grad) x H(grad) x H(grad)
 
 where each component is defined on the line.
 */

#ifndef Intrepid2_DerivedBasis_HGRAD_HEX_h
#define Intrepid2_DerivedBasis_HGRAD_HEX_h

#include "Intrepid2_TensorBasis.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_QUAD.hpp"

namespace Intrepid2
{
  // TODO: make this a subclass of TensorBasis3 instead, following what we've done for H(curl) and H(div)
  template<class HGRAD_LINE>
  class Basis_Derived_HGRAD_HEX
  : public Basis_TensorBasis<typename HGRAD_LINE::ExecutionSpace, typename HGRAD_LINE::OutputValueType, typename HGRAD_LINE::PointValueType>
  {
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineBasis = HGRAD_LINE;
    using QuadBasis = Intrepid2::Basis_Derived_HGRAD_QUAD<HGRAD_LINE>;
    using TensorBasis = Basis_TensorBasis<ExecutionSpace, OutputValueType, PointValueType>;

    std::string name_;
    ordinal_type order_x_;
    ordinal_type order_y_;
    ordinal_type order_z_;
    EPointType pointType_;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new QuadBasis(polyOrder_x,polyOrder_y, pointType)),
                Teuchos::rcp( new LineBasis(polyOrder_z, pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HGRAD;

      std::ostringstream basisName;
      basisName << "HGRAD_HEX (" << this->TensorBasis::getName() << ")";
      name_ = basisName.str();

      order_x_ = polyOrder_x;
      order_y_ = polyOrder_y;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
    }

    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_HEX(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) :
      Basis_Derived_HGRAD_HEX(polyOrder, polyOrder, polyOrder, pointType) {}


    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return (this->getDofCount(1,0) > 1); //if it has more than 1 DOF per edge, than it needs orientations
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      
      if (operatorType == VALUE)
      {
        return OperatorTensorDecomposition(VALUE,VALUE);
      }
      else if (operatorType == GRAD)
      {
        // to evaluate gradient, we need both OP_VALUE and OP_GRAD (thanks to product rule)
        // for quad x line, we will put derivative * value in first component, and value * derivative in second
        
        std::vector< std::vector<EOperator> > ops;
        ops.push_back(std::vector<EOperator>{GRAD,  VALUE});
        ops.push_back(std::vector<EOperator>{VALUE, GRAD});
        
        std::vector<double> weights(ops.size(), 1.0);
        
        return OperatorTensorDecomposition(ops, weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }

    using Basis<ExecutionSpace,OutputValueType,PointValueType>::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x,y dimensions
        \param [in] inputPoints2 - input points in the z dimension
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
        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::make_pair(0,2));
        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        // compute first component -- derivative happens in x and y, and value taken in z
        op1 = Intrepid2::OPERATOR_GRAD;
        op2 = Intrepid2::OPERATOR_VALUE;
        
        this->TensorBasis::getValues(outputValuesComponent1,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        
        // second component -- value in x and y, derivative in z
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

    /** \brief  Returns basis name

     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell should be the restriction to the subCell
        of the bases of the parent cell.
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
        case 4:
        case 6:
          return Teuchos::rcp( new LineBasis(order_x_, pointType_) );
        case 1:
        case 3:
        case 5:
        case 7:
          return Teuchos::rcp( new LineBasis(order_y_, pointType_) );
        case 8:
        case 9:
        case 10:
        case 11:
          return Teuchos::rcp( new LineBasis(order_z_, pointType_) );
        }
      } else if(subCellDim == 2) {
        switch(subCellOrd) {
        case 0:
          return Teuchos::rcp( new QuadBasis(order_x_, order_z_, pointType_) );
        case 1:
          return Teuchos::rcp( new QuadBasis(order_y_,order_z_, pointType_) );
        case 2:
          return Teuchos::rcp( new QuadBasis(order_x_, order_z_, pointType_) );
        case 3:
          return Teuchos::rcp( new QuadBasis(order_z_, order_y_, pointType_) );
        case 4:
          return Teuchos::rcp( new QuadBasis(order_y_, order_x_, pointType_) );
        case 5:
          return Teuchos::rcp( new QuadBasis(order_x_, order_y_, pointType_) );
        }
      }

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_HEX_h */
