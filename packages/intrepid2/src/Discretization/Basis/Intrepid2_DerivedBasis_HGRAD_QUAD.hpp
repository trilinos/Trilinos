// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HGRAD_QUAD.hpp
    \brief  Implementation of H(grad) basis on the quadrilateral that is templated on H(grad) on the line.
    \author Created by N.V. Roberts.
 
 H(grad) on the quadrilateral is defined (and implemented) as H(grad) x H(grad) on the line.
 */

#ifndef Intrepid2_DerivedBasis_HGRAD_QUAD_h
#define Intrepid2_DerivedBasis_HGRAD_QUAD_h

#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HGRAD_LINE>
  class Basis_Derived_HGRAD_QUAD
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {
  protected:
    std::string name_;
    ordinal_type order_x_, order_y_;
    EPointType pointType_;
  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using BasisBase = typename HGRAD_LINE::BasisBase;
    using LineBasis = HGRAD_LINE;
    using TensorBasis = Basis_TensorBasis<BasisBase>;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_QUAD(int polyOrder_x, int polyOrder_y, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new LineBasis(polyOrder_x, pointType)),
                Teuchos::rcp( new LineBasis(polyOrder_y, pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HGRAD;

      std::ostringstream basisName;
      basisName << "HGRAD_QUAD (" << this->TensorBasis::getName() << ")";
      name_ = basisName.str();

      order_x_= polyOrder_x;
      order_y_ = polyOrder_y;
      pointType_ = pointType;
      
      this->setShardsTopologyAndTags();
    }

    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in both dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HGRAD_QUAD(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HGRAD_QUAD(polyOrder,polyOrder,pointType) {}


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

          The bases of the subCell are the restriction to the subCell
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
          return Teuchos::rcp( new LineBasis(order_x_, pointType_) );
        case 1:
        case 3:
          return Teuchos::rcp( new LineBasis(order_y_, pointType_) );
        }
      }

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HGRAD_QUAD<typename HGRAD_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_x_, order_y_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HGRAD_QUAD_h */
