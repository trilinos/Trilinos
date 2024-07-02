// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HGRAD_WEDGE.hpp
    \brief  Implementation of H(grad) basis on the wedge that is templated on H(grad) on the line, and H(grad) on the triangle.
    \author Created by N.V. Roberts.
 
 H(grad) on the wedge is defined (and implemented) as H(grad, triangle) x H(grad, line).
 */

#ifndef Intrepid2_DerivedBasis_HGRAD_WEDGE_h
#define Intrepid2_DerivedBasis_HGRAD_WEDGE_h

#include "Intrepid2_TensorBasis.hpp"
#include "Intrepid2_DerivedBasis_HGRAD_QUAD.hpp"

namespace Intrepid2
{
  template<class HGRAD_TRI, class HGRAD_LINE>
  class Basis_Derived_HGRAD_WEDGE
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {
  protected:
    std::string name_;
    ordinal_type order_xy_, order_z_;
    EPointType pointType_;
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using BasisBase = typename HGRAD_LINE::BasisBase;
    using TriBasis  = HGRAD_TRI;
    using LineBasis = HGRAD_LINE;
    using TensorBasis = Basis_TensorBasis<BasisBase>;

    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new TriBasis(polyOrder_xy, pointType)),
                Teuchos::rcp( new LineBasis(polyOrder_z, pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HGRAD;

      std::ostringstream basisName;
      basisName << "HGRAD_WEDGE (" << this->TensorBasis::getName() << ")";
      name_ = basisName.str();

      order_xy_= polyOrder_xy;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
      
      this->setShardsTopologyAndTags();
    }

    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in both dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_WEDGE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HGRAD_WEDGE(polyOrder,polyOrder,pointType) {}


    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return (this->getDofCount(1,0) > 1); //if it has more than 1 DOF per edge, than it needs orientations
    }
    
    using BasisBase::getValues;

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
      if(subCellDim == 1) {
        if(subCellOrd < 6)  //edges associated to basal and upper triagular faces
          return Teuchos::rcp( new LineBasis(order_xy_, pointType_) );
        else
          return Teuchos::rcp( new LineBasis(order_z_, pointType_) );
      }
      else if(subCellDim == 2) {
        switch(subCellOrd) {
        case 0:
          return Teuchos::rcp( new Intrepid2::Basis_Derived_HGRAD_QUAD<LineBasis>(order_xy_, order_z_, pointType_) );
        case 1:
          return Teuchos::rcp( new Intrepid2::Basis_Derived_HGRAD_QUAD<LineBasis>(order_xy_, order_z_, pointType_) );
        case 2:
          return Teuchos::rcp( new Intrepid2::Basis_Derived_HGRAD_QUAD<LineBasis>(order_z_, order_xy_, pointType_) );
        case 3:
          return Teuchos::rcp( new TriBasis(order_xy_, pointType_) );
        case 4:
          return Teuchos::rcp( new TriBasis(order_xy_, pointType_) );
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
      using HostBasis  = Basis_Derived_HGRAD_WEDGE<typename HGRAD_TRI::HostBasis, typename HGRAD_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_xy_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_WEDGE_h */
