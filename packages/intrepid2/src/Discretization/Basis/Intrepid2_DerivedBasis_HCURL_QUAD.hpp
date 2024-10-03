// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Sacado.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family1_QUAD
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;

    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;

    using BasisBase = typename HGRAD_LINE::BasisBase;

    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;

    using TensorBasis = Basis_TensorBasis<BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_Family1_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new LineHVolBasis(polyOrder_x-1,pointType)),
                Teuchos::rcp( new LineGradBasis(polyOrder_y,pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }

    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      const EOperator CURL  = Intrepid2::OPERATOR_CURL;
      if (operatorType == VALUE)
      {
        // family 1 goes in x component
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{VALUE,VALUE};
        ops[1] = std::vector<EOperator>{};
        std::vector<double> weights {1.0, 0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == CURL)
      {
        // family 1 gets a -d/dy applied to the first (nonzero) vector component
        // since this is H(VOL)(x) * H(GRAD)(y), this amounts to taking the derivative in the second tensorial component
        return OperatorTensorDecomposition(VALUE,GRAD,-1.0);
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
                           const PointViewType inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const override
    {
      Intrepid2::EOperator op1, op2;
      if (operatorType == Intrepid2::OPERATOR_VALUE)
      {
        op1 = Intrepid2::OPERATOR_VALUE;
        op2 = Intrepid2::OPERATOR_VALUE;

        // family 1 goes in the x component; 0 in the y component
        OutputViewType outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
        OutputViewType outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);

        this->TensorBasis::getValues(outputValuesComponent1,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        // place 0 in the y component
        Kokkos::deep_copy(outputValuesComponent2,0);
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
      this->TensorBasis::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs2,0.0);
    }
  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_Family2_QUAD
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {

  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;

    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;

    using LineGradBasis = HGRAD_LINE;
    using LineHVolBasis = HVOL_LINE;

    using BasisBase = typename HGRAD_LINE::BasisBase;

    using TensorBasis = Basis_TensorBasis<BasisBase>;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_Family2_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new LineGradBasis(polyOrder_x,pointType) ),
                Teuchos::rcp( new LineHVolBasis(polyOrder_y-1,pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }

    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      const EOperator GRAD  = Intrepid2::OPERATOR_GRAD;
      const EOperator CURL  = Intrepid2::OPERATOR_CURL;
      if (operatorType == VALUE)
      {
        // family 2 goes in y component
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{};
        ops[1] = std::vector<EOperator>{VALUE,VALUE};
        std::vector<double> weights {0.0, 1.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == CURL)
      {
        // family 2 gets a d/dx applied to the second (nonzero) vector component
        // since this is H(GRAD)(x) * H(VOL)(y), this amounts to taking the derivative in the first tensorial component
        return OperatorTensorDecomposition(GRAD,VALUE,1.0);
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
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis::getDofCoeffs(dofCoeffs2);
    }
  };

  template<class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_QUAD
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family1 = Basis_Derived_HCURL_Family1_QUAD<HGRAD_LINE, HVOL_LINE>;
    using Family2 = Basis_Derived_HCURL_Family2_QUAD<HGRAD_LINE, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>;
  public:
    using BasisBase = typename HGRAD_LINE::BasisBase;

  protected:
    std::string name_;
    ordinal_type order_x_;
    ordinal_type order_y_;
    EPointType pointType_;

  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Teuchos::rcp( new Family1(polyOrder_x, polyOrder_y, pointType) ),
                   Teuchos::rcp( new Family2(polyOrder_x, polyOrder_y, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;

      std::ostringstream basisName;
      basisName << "HCURL_QUAD (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_x_ = polyOrder_x;
      order_y_ = polyOrder_y;
      pointType_ = pointType;
    }

    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_QUAD(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HCURL_QUAD(polyOrder, polyOrder, pointType) {}

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

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected to the subCell line.

        TODO: test this method when different orders are used in different directions
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    Teuchos::RCP<BasisBase>
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

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.

        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HCURL_QUAD<typename HGRAD_LINE::HostBasis, typename HVOL_LINE::HostBasis>;

      auto hostBasis = Teuchos::rcp(new HostBasis(order_x_, order_y_, pointType_));

      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HCURL_QUAD_h */
