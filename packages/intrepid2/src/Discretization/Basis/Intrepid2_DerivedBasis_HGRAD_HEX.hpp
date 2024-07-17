// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
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
    using BasisBase = typename HGRAD_LINE::BasisBase;
    using TensorBasis = Basis_TensorBasis<BasisBase>;

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
      
      this->setShardsTopologyAndTags();
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
      return (this->getDofCount(1,0) > 1); //if it has more than 1 DOF per edge, then it needs orientations
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

        The bases of the subCell should be the restriction to the subCell
        of the bases of the parent cell.
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
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HGRAD_HEX<typename HGRAD_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_x_, order_y_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_HEX_h */
