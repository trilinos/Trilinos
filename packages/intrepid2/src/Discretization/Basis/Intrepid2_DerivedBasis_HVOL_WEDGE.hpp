// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HVOL_WEDGE.hpp
    \brief  Implementation of H(vol) basis on the wedge that is templated on H(grad) on the line, and H(grad) on the triangle.
    \author Created by N.V. Roberts.
 
 H(grad) on the wedge is defined (and implemented) as H(vol, triangle) x H(vol, line).
 */

#ifndef Intrepid2_DerivedBasis_HVOL_WEDGE_h
#define Intrepid2_DerivedBasis_HVOL_WEDGE_h

#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HVOL_TRI, class HVOL_LINE>
  class Basis_Derived_HVOL_WEDGE
  : public Basis_TensorBasis<typename HVOL_LINE::BasisBase>
  {
  protected:
    std::string name_;
    ordinal_type order_xy_, order_z_;
    EPointType pointType_;
  public:
    using ExecutionSpace  = typename HVOL_LINE::ExecutionSpace;
    using OutputValueType = typename HVOL_LINE::OutputValueType;
    using PointValueType  = typename HVOL_LINE::PointValueType;
    
    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using BasisBase = typename HVOL_LINE::BasisBase;
    using TriBasis  = HVOL_TRI;
    using LineBasis = HVOL_LINE;
    using TensorBasis = Basis_TensorBasis<BasisBase>;

    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HVOL_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new TriBasis(polyOrder_xy, pointType)),
                Teuchos::rcp( new LineBasis(polyOrder_z, pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HVOL;

      std::ostringstream basisName;
      basisName << "HVOL_WEDGE (" << this->TensorBasis::getName() << ")";
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
    Basis_Derived_HVOL_WEDGE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HVOL_WEDGE(polyOrder,polyOrder,pointType) {}


    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return false;
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
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HVOL_WEDGE<typename HVOL_TRI::HostBasis, typename HVOL_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_xy_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HVOL_WEDGE_h */
