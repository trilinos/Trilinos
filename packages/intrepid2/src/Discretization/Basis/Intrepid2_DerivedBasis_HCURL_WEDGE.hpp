// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HCURL_WEDGE.hpp
    \brief  Implementation of H(curl) basis on the wedge that is templated on H(grad,tri), H(curl,tri), H(grad,line), and H(vol,line).
    \author Created by N.V. Roberts.
 
 This class constructs the H(curl) space as the direct sum of two families of tensor-product bases on the triangle and line:
 - family 1: H(curl,tri)  x  H(grad,line), placed in the x and y components of vector output
 - family 2: H(grad,tri) x  H(vol,line),  placed in the z component of vector output
 
 The way that each family decomposes into operators on the component bases cannot be expressed simply with scalar
 weights and EOperators on component bases; instead, a 90-degree rotation is required for the curl evaluations.  This
 motivated the addition of a boolean flag indicating such a rotation in OperatorTensorDecomposition.
 
 Our Family 1 corresponds to the following ESEAS entities:
 - mixed edges
 - triangle faces, Family I and II
 - quadrilateral faces, Family II
 - Interior Family I
 - Interior Family II.
 
 Our Family 2 corresponds to:
 - quadrilateral edges
 - quadrilateral faces, Family I
 - Interior Family III.
 
 See p. 449 in Fuentes et al. (https://doi.org/10.1016/j.camwa.2015.04.027).
 */

#ifndef Intrepid2_DerivedBasis_HCURL_WEDGE_h
#define Intrepid2_DerivedBasis_HCURL_WEDGE_h

#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Sacado.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HCURL_TRI, class HGRAD_LINE>
  class Basis_Derived_HCURL_Family1_WEDGE
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
    
    using TriCurlBasis  = HCURL_TRI;
    using LineGradBasis = HGRAD_LINE;
    
    using TensorBasis = Basis_TensorBasis<BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_Family1_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new TriCurlBasis(polyOrder_xy,pointType)),
                Teuchos::rcp( new LineGradBasis(polyOrder_z,pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }
    
    using TensorBasis::getValues;
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator & VALUE = OPERATOR_VALUE;
      const EOperator & GRAD  = OPERATOR_GRAD;
      const EOperator & CURL  = OPERATOR_CURL;
      if (operatorType == VALUE)
      {
        // family 1 goes in x,y components
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{VALUE,VALUE}; // occupies x,y
        ops[1] = std::vector<EOperator>{}; // zero z
        std::vector<double> weights {1.0, 0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == CURL)
      {
        // curl of (f_x(x,y) g(z), f_y(x,y) g(z), 0), where f is in H(curl,tri), g in H(grad,line)
        // x,y components of curl: rot(f) dg/dz, where rot(f) is a 90-degree rotation of f.
        //   z component  of curl: curl(f) g, where curl(f) is the 2D curl, a scalar.
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{VALUE,GRAD}; // occupies the x,y components
        ops[1] = std::vector<EOperator>{CURL,VALUE}; // z component
        std::vector<double> weights {1.0, 1.0};
        OperatorTensorDecomposition opDecomposition(ops, weights);
        opDecomposition.setRotateXYNinetyDegrees(true);
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

        // family 1 values go in the x,y components; 0 in the z component
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);

        // place 0 in the z component
        Kokkos::deep_copy(outputValuesComponent3, 0.0);
        this->TensorBasis::getValues(outputValuesComponent12,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);

      }
      else if (operatorType == OPERATOR_CURL)
      {
        // curl of (f_x(x,y) g(z), f_y(x,y) g(z), 0), where f is in H(curl,tri), g in H(grad,line)
        // x,y components of curl: rot(f) dg/dz, where rot(f) is a 90-degree rotation of f.
        //   z component  of curl: curl(f) g, where curl(f) is the 2D curl, a scalar.
        
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
        
        op1 = OPERATOR_VALUE;
        op2 = OPERATOR_GRAD;

        this->TensorBasis::getValues(outputValuesComponent12,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{outputValuesComponent12.extent_int(0),outputValuesComponent12.extent_int(1)});
        Kokkos::parallel_for("wedge family 1 curl: rotateXYNinetyDegrees CW", policy,
        KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
          const auto  f_x = outputValuesComponent12(fieldOrdinal,pointOrdinal,0); // copy
          const auto &f_y = outputValuesComponent12(fieldOrdinal,pointOrdinal,1); // reference
          outputValuesComponent12(fieldOrdinal,pointOrdinal,0) = -f_y;
          outputValuesComponent12(fieldOrdinal,pointOrdinal,1) =  f_x;
        });
        
        op1 = OPERATOR_CURL;
        op2 = OPERATOR_VALUE;
        this->TensorBasis::getValues(outputValuesComponent3,
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
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),std::make_pair(0,2));
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),2);
      this->TensorBasis::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs2,0.0);
    }
  };

  template<class HGRAD_TRI, class HVOL_LINE>
  class Basis_Derived_HCURL_Family2_WEDGE
  : public Basis_TensorBasis<typename HVOL_LINE::BasisBase>
  {

  public:
    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using TriGradBasis  = HGRAD_TRI;
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
    Basis_Derived_HCURL_Family2_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new  TriGradBasis(polyOrder_xy,  pointType) ),
                Teuchos::rcp( new LineHVolBasis(polyOrder_z-1, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  One-element ops and weights vectors correspond to a single TensorData entry; multiple-element vectors correspond to a VectorData object with axialComponents = false.
    */
    OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      if (operatorType == OPERATOR_CURL)
      {
        // curl of (0,0,f) is (df/dy, -df/dx, 0)
        
        // this is a rotation of gradient of the triangle H^1 basis, times the H(vol) line value
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{OPERATOR_GRAD,OPERATOR_VALUE}; // occupies the x,y components
        ops[1] = std::vector<EOperator>{};
        std::vector<double> weights {-1.0, 0.0}; // -1 because the rotation goes from (df/dx,df/dy) --> (-df/dy,df/dx), and we want (df/dy,-df/dx)
        OperatorTensorDecomposition opDecomposition(ops, weights);
        opDecomposition.setRotateXYNinetyDegrees(true);
        return opDecomposition;
      }
      else if (OPERATOR_VALUE == operatorType)
      {
        // family 2 goes in z component
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{}; // because family 1 identifies this as spanning (x,y), empty op here will also span (x,y)
        ops[1] = std::vector<EOperator>{OPERATOR_VALUE,OPERATOR_VALUE}; // z component
        std::vector<double> weights {0.0, 1.0};
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

        // family 2 values go in z component, 0 in (x,y)
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);

        // place 0 in the x,y components
        Kokkos::deep_copy(outputValuesComponent12, 0.0);
        this->TensorBasis::getValues(outputValuesComponent3,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);

      }
      else if (operatorType == OPERATOR_CURL)
      {
        // curl of (0,0,f) is (df/dy, -df/dx, 0)
        
        auto outputValuesComponent12 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),std::pair<int,int>{0,2});
        auto outputValuesComponent3  = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),2);
                
        op1 = OPERATOR_GRAD;
        op2 = OPERATOR_VALUE;

        this->TensorBasis::getValues(outputValuesComponent12,
                                     inputPoints1, op1,
                                     inputPoints2, op2, tensorPoints);
        
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{outputValuesComponent12.extent_int(0),outputValuesComponent12.extent_int(1)});
        Kokkos::parallel_for("wedge family 2 curl: rotateXYNinetyDegrees CCW", policy,
        KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
          const auto  f_x = outputValuesComponent12(fieldOrdinal,pointOrdinal,0); // copy
          const auto &f_y = outputValuesComponent12(fieldOrdinal,pointOrdinal,1); // reference
          outputValuesComponent12(fieldOrdinal,pointOrdinal,0) =  f_y;
          outputValuesComponent12(fieldOrdinal,pointOrdinal,1) = -f_x;
        });
        
        Kokkos::deep_copy(outputValuesComponent3, 0.0);
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
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),std::make_pair(0,2));
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),2);
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis::getDofCoeffs(dofCoeffs2);
    }
  };
  
  template<class HGRAD_TRI, class HCURL_TRI, class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_WEDGE
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family1 = Basis_Derived_HCURL_Family1_WEDGE<HCURL_TRI, HGRAD_LINE>;
    using Family2 = Basis_Derived_HCURL_Family2_WEDGE<HGRAD_TRI, HVOL_LINE>;
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
    Basis_Derived_HCURL_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Teuchos::rcp( new Family1(polyOrder_xy, polyOrder_z, pointType) ),
                   Teuchos::rcp( new Family2(polyOrder_xy, polyOrder_z, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;

      std::ostringstream basisName;
      basisName << "HCURL_WEDGE (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_xy_ = polyOrder_xy;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_WEDGE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HCURL_WEDGE(polyOrder, polyOrder, pointType) {}

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
      using LineBasis = HVOL_LINE;
      using TriCurlBasis  = HCURL_TRI;
      using QuadCurlBasis = Basis_Derived_HCURL_QUAD<HGRAD_LINE,HVOL_LINE>;
      if(subCellDim == 1) {
        if(subCellOrd < 6)  //edges associated to basal and upper triagular faces
          return Teuchos::rcp( new LineBasis(order_xy_-1, pointType_) );
        else
          return Teuchos::rcp( new LineBasis(order_z_-1, pointType_) );
      }
      else if(subCellDim == 2) {
        switch(subCellOrd) {
        case 0:
          return Teuchos::rcp( new QuadCurlBasis(order_xy_, order_z_, pointType_) );
        case 1:
          return Teuchos::rcp( new QuadCurlBasis(order_xy_, order_z_, pointType_) );
        case 2:
          return Teuchos::rcp( new QuadCurlBasis(order_z_, order_xy_, pointType_) );
        case 3:
          return Teuchos::rcp( new TriCurlBasis(order_xy_, pointType_) );
        case 4:
          return Teuchos::rcp( new TriCurlBasis(order_xy_, pointType_) );
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
      using HostBasis  = Basis_Derived_HCURL_WEDGE<typename HGRAD_TRI::HostBasis, typename HCURL_TRI::HostBasis, typename HGRAD_LINE::HostBasis, typename HVOL_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_xy_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HCURL_WEDGE_h */
