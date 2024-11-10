// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HDIV_WEDGE.hpp
    \brief  Implementation of H(div) basis on the wedge that is templated on H(div,tri), H(vol,tri), H(grad,line), and H(vol,line).
    \author Created by N.V. Roberts.
 
 This class constructs the H(curl) space as the direct sum of two families of tensor-product bases on the triangle and line:
 - family 1: H(div,tri) x  H(vol,line),  placed in the x and y components of vector output 
 - family 2: H(vol,tri) x  H(grad,line), placed in the z component of vector output
  
 Our Family 1 corresponds to the following ESEAS entities:
 - quadrilateral faces
 - Interior Family I
 - Interior Family II.
 
 Our Family 2 corresponds to:
 - triangle faces
 - Interior Family III.
 
 See p. 450 in Fuentes et al. (https://doi.org/10.1016/j.camwa.2015.04.027).
 */

#ifndef Intrepid2_DerivedBasis_HDIV_WEDGE_h
#define Intrepid2_DerivedBasis_HDIV_WEDGE_h

#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Sacado.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HDIV_TRI, class HVOL_LINE>
  class Basis_Derived_HDIV_Family1_WEDGE
  : public Basis_TensorBasis<typename HVOL_LINE::BasisBase>
  {

  public:
    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using TriDivBasis   = HDIV_TRI;
    using LineHVolBasis = HVOL_LINE;
    
    using BasisBase   = typename HVOL_LINE::BasisBase;
    using TensorBasis = Basis_TensorBasis<BasisBase>;
    
    using DeviceType      = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family1_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new   TriDivBasis(polyOrder_xy,  pointType) ),
                Teuchos::rcp( new LineHVolBasis(polyOrder_z-1, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
      this->setShardsTopologyAndTags();
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  One-element ops and weights vectors correspond to a single TensorData entry; multiple-element vectors correspond to a VectorData object with axialComponents = false.
    */
    OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      if (operatorType == OPERATOR_DIV)
      {
        // div of (f(x,y)h(z),g(x,y)h(z),0) is (df/dx + dg/dy) * h
        
        std::vector< std::vector<EOperator> > ops(1);
        ops[0] = std::vector<EOperator>{OPERATOR_DIV,OPERATOR_VALUE};
        std::vector<double> weights {1.0};
        OperatorTensorDecomposition opDecomposition(ops, weights);
        return opDecomposition;
      }
      else if (OPERATOR_VALUE == operatorType)
      {
        // family 1 goes in x,y components
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{OPERATOR_VALUE,OPERATOR_VALUE}; // spans (x,y) due to H(div,tri) (first tensorial component)
        ops[1] = std::vector<EOperator>{}; // empty z component
        std::vector<double> weights {1.0,0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
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
                           bool tensorPoints) const override
    {
      EOperator op1, op2;
      if (operatorType == OPERATOR_VALUE)
      {
        op1 = OPERATOR_VALUE;
        op2 = OPERATOR_VALUE;

        // family 1 values goes in (x,y) components, 0 in z
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);

        // place 0 in the z component
        Kokkos::deep_copy(outputValuesComponent3, 0.0);
        this->TensorBasis::getValues(outputValuesComponent12,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);

      }
      else if (operatorType == OPERATOR_DIV)
      {
        // div of (f(x,y)h(z),g(x,y)h(z),0) is (df/dx + dg/dy) * h
        
        op1 = OPERATOR_DIV;
        op2 = OPERATOR_VALUE;

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
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const override {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(), std::make_pair(0,2));
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(), 2);
      this->TensorBasis::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs2,0.0);
    }
  };

  template<class HVOL_TRI, class HGRAD_LINE>
  class Basis_Derived_HDIV_Family2_WEDGE
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using BasisBase = typename HGRAD_LINE::BasisBase;

    using DeviceType      = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    using TriVolBasis   = HVOL_TRI;
    using LineGradBasis = HGRAD_LINE;
    
    using TensorBasis = Basis_TensorBasis<BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family2_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new TriVolBasis(polyOrder_xy-1,pointType)),
                Teuchos::rcp( new LineGradBasis(polyOrder_z,pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
      this->setShardsTopologyAndTags();
    }
    
    using TensorBasis::getValues;
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator & VALUE = OPERATOR_VALUE;
      const EOperator & GRAD  = OPERATOR_GRAD;
      const EOperator & DIV   = OPERATOR_DIV;
      if (operatorType == VALUE)
      {
        // family 2 goes in z component
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{}; // will span x,y because family 1's first component does
        ops[1] = std::vector<EOperator>{VALUE,VALUE};
        std::vector<double> weights {0.,1.0};
        OperatorTensorDecomposition opDecomposition(ops, weights);
        return opDecomposition;
      }
      else if (operatorType == DIV)
      {
        // div of (0, 0, f*g, where f=f(x,y), g=g(z), equals f * g'(z).
        std::vector< std::vector<EOperator> > ops(1);
        ops[0] = std::vector<EOperator>{VALUE,GRAD};
        std::vector<double> weights {1.0};
        OperatorTensorDecomposition opDecomposition(ops, weights);
        return opDecomposition;
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
    }
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const override
    {
      EOperator op1, op2;
      if (operatorType == OPERATOR_VALUE)
      {
        op1 = OPERATOR_VALUE;
        op2 = OPERATOR_VALUE;

        // family 2 values goes in z component
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);

        // place 0 in the x,y components
        Kokkos::deep_copy(outputValuesComponent12, 0.0);
        this->TensorBasis::getValues(outputValuesComponent3,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);

      }
      else if (operatorType == OPERATOR_DIV)
      {
        // div of (0, 0, f*g, where f=f(x,y), g=g(z), equals f * g'(z).
        op1 = OPERATOR_VALUE;
        op2 = OPERATOR_GRAD;

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
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const override {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(), std::make_pair(0,2));
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(), 2);
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis::getDofCoeffs(dofCoeffs2);
    }
  };

  
  template<class HDIV_TRI, class HVOL_TRI, class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_WEDGE
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family1 = Basis_Derived_HDIV_Family1_WEDGE<HDIV_TRI, HVOL_LINE>;
    using Family2 = Basis_Derived_HDIV_Family2_WEDGE<HVOL_TRI, HGRAD_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>;
  public:
    using BasisBase = typename HGRAD_LINE::BasisBase;

  protected:
    std::string name_;
    ordinal_type order_xy_;
    ordinal_type order_z_;
    EPointType pointType_;

  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Teuchos::rcp( new Family1(polyOrder_xy, polyOrder_z, pointType) ),
                   Teuchos::rcp( new Family2(polyOrder_xy, polyOrder_z, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;

      std::ostringstream basisName;
      basisName << "HDIV_WEDGE (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_xy_ = polyOrder_xy;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_WEDGE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HDIV_WEDGE(polyOrder, polyOrder, pointType) {}

    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override
    {
      return true;
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

      The bases of the subCell are the restriction to the subCell
      of the bases of the parent cell.
      TODO: test this method when different orders are used in different directions
      \param [in] subCellDim - dimension of subCell
      \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
      \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
    */
    Teuchos::RCP<BasisBase>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override
    {
      using QuadBasis = Basis_Derived_HVOL_QUAD<HVOL_LINE>;
      using TriBasis  = HVOL_TRI;

      if(subCellDim == 2) {
        switch(subCellOrd) {
        case 0:
          return Teuchos::rcp( new QuadBasis(order_xy_-1, order_z_-1, pointType_) );
        case 1:
          return Teuchos::rcp( new QuadBasis(order_xy_-1, order_z_-1, pointType_) );
        case 2:
          return Teuchos::rcp( new QuadBasis(order_z_-1, order_xy_-1, pointType_) );
        case 3:
          return Teuchos::rcp( new TriBasis(order_xy_-1, pointType_) );
        case 4:
          return Teuchos::rcp( new TriBasis(order_xy_-1, pointType_) );
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"subCellOrd is out of bounds");
        }
      } 
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"subCellDim is out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HDIV_WEDGE<typename HDIV_TRI::HostBasis, typename HVOL_TRI::HostBasis, typename HGRAD_LINE::HostBasis, typename HVOL_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_xy_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HDIV_WEDGE_h */
