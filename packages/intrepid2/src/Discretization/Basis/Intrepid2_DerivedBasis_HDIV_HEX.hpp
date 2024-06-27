// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_Sacado.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family1_HEX
  :
  public Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family1_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis3(Teuchos::rcp( new LineGradBasis(polyOrder_x,pointType)),
                 Teuchos::rcp( new LineHVolBasis(polyOrder_y-1,pointType)),
                 Teuchos::rcp( new LineHVolBasis(polyOrder_z-1,pointType)),
                 true) // true: use shards CellTopology and tags
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, basis2, and basis3.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      const EOperator DIV   = Intrepid2::OPERATOR_DIV;
      
      const double weight = 1.0;
      if (operatorType == VALUE)
      {
        std::vector< std::vector<EOperator> > ops(3);
        ops[0] = std::vector<EOperator>{VALUE,VALUE,VALUE};
        ops[1] = std::vector<EOperator>{};
        ops[2] = std::vector<EOperator>{};
        std::vector<double> weights {weight,0.0,0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == DIV)
      {
        // family 1 is nonzero in the x component, so the div is (GRAD,VALUE,VALUE)
        std::vector< std::vector<EOperator> > ops(1); // scalar value
        ops[0] = std::vector<EOperator>{GRAD,VALUE,VALUE};
        std::vector<double> weights {weight};
        return OperatorTensorDecomposition(ops,weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
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
                           bool tensorPoints) const override
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
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs23 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),std::make_pair(1,3));
      this->TensorBasis3::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs23,0.0);
    }
  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family2_HEX
  :
  public Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family2_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis3(Teuchos::rcp( new LineHVolBasis(polyOrder_x-1,pointType) ),
                 Teuchos::rcp( new LineGradBasis(polyOrder_y,pointType) ),
                 Teuchos::rcp( new LineHVolBasis(polyOrder_z-1,pointType) ),
                 true) // true: use shards CellTopology and tags
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, basis2, and basis3.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      const EOperator DIV   = Intrepid2::OPERATOR_DIV;
      
      const double weight = 1.0;
      if (operatorType == VALUE)
      {
        std::vector< std::vector<EOperator> > ops(3);
        ops[0] = std::vector<EOperator>{};
        ops[1] = std::vector<EOperator>{VALUE,VALUE,VALUE};
        ops[2] = std::vector<EOperator>{};
        std::vector<double> weights {0.0,weight,0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == DIV)
      {
        // family 2 is nonzero in the y component, so the div is (VALUE,GRAD,VALUE)
        std::vector< std::vector<EOperator> > ops(1); // scalar value
        ops[0] = std::vector<EOperator>{VALUE,GRAD,VALUE};
        std::vector<double> weights {weight};
        return OperatorTensorDecomposition(ops,weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
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
                           bool tensorPoints) const override
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
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),1);
      auto dofCoeffs3 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),2);
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis3::getDofCoeffs(dofCoeffs2);
      Kokkos::deep_copy(dofCoeffs3,0.0);
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family3_HEX
  : public Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;
    
    using TensorBasis3 = Basis_TensorBasis3<typename HGRAD_LINE::BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family3_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis3(Teuchos::rcp( new LineHVolBasis(polyOrder_x-1,pointType) ),
                 Teuchos::rcp( new LineHVolBasis(polyOrder_y-1,pointType) ),
                 Teuchos::rcp( new LineGradBasis(polyOrder_z,pointType) ),
                 true) // true: use shards CellTopology and tags
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, basis2, and basis3.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      const EOperator DIV   = Intrepid2::OPERATOR_DIV;
      
      const double weight = 1.0;
      if (operatorType == VALUE)
      {
        std::vector< std::vector<EOperator> > ops(3);
        ops[0] = std::vector<EOperator>{};
        ops[1] = std::vector<EOperator>{};
        ops[2] = std::vector<EOperator>{VALUE,VALUE,VALUE};
        std::vector<double> weights {0.0,0.0,weight};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == DIV)
      {
        // family 3 is nonzero in the z component, so the div is (VALUE,VALUE,GRAD)
        std::vector< std::vector<EOperator> > ops(1); // scalar value
        ops[0] = std::vector<EOperator>{VALUE,VALUE,GRAD};
        std::vector<double> weights {weight};
        return OperatorTensorDecomposition(ops,weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
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
                           bool tensorPoints) const override
    {
      Intrepid2::EOperator op1, op2, op3;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;
        op3 = Intrepid2::OPERATOR_VALUE;
        
        // family 3 goes in the z component; 0 in the x and y components
        auto outputValuesComponent_xy = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::make_pair(0,2));
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
      auto dofCoeffs12 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),std::make_pair(0,2));
      auto dofCoeffs3 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),2);
      Kokkos::deep_copy(dofCoeffs12,0.0);
      this->TensorBasis3::getDofCoeffs(dofCoeffs3);
    }

  };
  
  // ESEAS numbers its H(div) families differently, with the nonzero going in z, x, y for I,II,III.
  // To allow our interior orderings to match that of ESEAS, we put the direct sum in the same order as ESEAS,
  // which is to say that we go 3,1,2.
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_Family3_Family1_HEX
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family3 = Basis_Derived_HDIV_Family3_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family1 = Basis_Derived_HDIV_Family1_HEX<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<typename HGRAD_LINE::BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_Family3_Family1_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType)
    :
    DirectSumBasis(Teuchos::rcp( new Family3(polyOrder_x, polyOrder_y, polyOrder_z, pointType)),
                   Teuchos::rcp( new Family1(polyOrder_x, polyOrder_y, polyOrder_z, pointType))) {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;
    }
  };
  
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HDIV_HEX
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family31 = Basis_Derived_HDIV_Family3_Family1_HEX<HGRAD_LINE, HVOL_LINE>;
    using Family2  = Basis_Derived_HDIV_Family2_HEX        <HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis<typename HGRAD_LINE::BasisBase>;

    std::string name_;
    ordinal_type order_x_;
    ordinal_type order_y_;
    ordinal_type order_z_;
    EPointType pointType_;

  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using BasisBase = typename HGRAD_LINE::BasisBase;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HDIV_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Teuchos::rcp(new Family31(polyOrder_x, polyOrder_y, polyOrder_z, pointType)),
                   Teuchos::rcp(new Family2 (polyOrder_x, polyOrder_y, polyOrder_z, pointType))) {
      this->functionSpace_ = FUNCTION_SPACE_HDIV;

      std::ostringstream basisName;
      basisName << "HDIV_HEX (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_x_ = polyOrder_x;
      order_y_ = polyOrder_y;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
     */
    Basis_Derived_HDIV_HEX(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HDIV_HEX(polyOrder, polyOrder, polyOrder, pointType) {}

    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
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

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected along normal to the subCell.

        TODO: test this method when different orders are used in different directions
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    Teuchos::RCP<BasisBase>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override {

      using QuadBasis = Basis_Derived_HVOL_QUAD<HVOL_LINE>;

      if(subCellDim == 2) {
        switch(subCellOrd) {
        case 0:
          return Teuchos::rcp( new QuadBasis(order_x_-1, order_z_-1, pointType_) );
        case 1:
          return Teuchos::rcp( new QuadBasis(order_y_-1,order_z_-1, pointType_) );
        case 2:
          return Teuchos::rcp( new QuadBasis(order_x_-1, order_z_-1, pointType_) );
        case 3:
          return Teuchos::rcp( new QuadBasis(order_z_-1, order_y_-1, pointType_) );
        case 4:
          return Teuchos::rcp( new QuadBasis(order_y_-1, order_x_-1, pointType_) );
        case 5:
          return Teuchos::rcp( new QuadBasis(order_x_-1, order_y_-1, pointType_) );
        }
      }

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HDIV_HEX<typename HGRAD_LINE::HostBasis, typename HVOL_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_x_, order_y_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HDIV_HEX_h */
