// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HVOL_HEX.hpp
    \brief  Implementation of H(vol) basis on the hexahedron that is templated on H(vol) on the line.
    \author Created by N.V. Roberts.
 
 H(vol) on the hexahedron is defined (and implemented) as H(vol) on the quad times H(vol) on the line.
 
 For consistency with the H(div) and H(curl) implementations on the hexahedron, we may in the future define this as
 
 H(vol) x H(vol) x H(vol)
 
 where each component is defined on the line.
 */

#ifndef Intrepid2_DerivedBasis_HVOL_HEX_h
#define Intrepid2_DerivedBasis_HVOL_HEX_h

#include "Intrepid2_TensorBasis.hpp"
#include "Intrepid2_DerivedBasis_HVOL_QUAD.hpp"

namespace Intrepid2
{
  /** \class Intrepid2::Basis_Derived_HVOL_HEX
      \brief Implementation of H(vol) basis on the quadrilateral that is templated on H(vol) on the line.
  */
  
  template<class HVOL_LINE>
  class Basis_Derived_HVOL_HEX
  :
  public Basis_TensorBasis<typename HVOL_LINE::BasisBase>
  // TODO: make this a subclass of TensorBasis3 instead, following what we've done for H(curl) and H(div)
  {
    std::string name_;
    ordinal_type polyOrder_x_, polyOrder_y_, polyOrder_z_;
    EPointType pointType_;
  public:
    using ExecutionSpace  = typename HVOL_LINE::ExecutionSpace;
    using OutputValueType = typename HVOL_LINE::OutputValueType;
    using PointValueType  = typename HVOL_LINE::PointValueType;

    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using BasisBase = typename HVOL_LINE::BasisBase;
    
    using LineBasis = HVOL_LINE;
    using QuadBasis = Intrepid2::Basis_Derived_HVOL_QUAD<HVOL_LINE>;
    using TensorBasis = Basis_TensorBasis<BasisBase>;

    /** \brief  Constructor.
        \param [in] polyOrder_x - the polynomial order in the x dimension.
        \param [in] polyOrder_y - the polynomial order in the y dimension.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HVOL_HEX(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp(new QuadBasis(polyOrder_x,polyOrder_y,pointType)),
                Teuchos::rcp(new LineBasis(polyOrder_z,pointType))),
    polyOrder_x_(polyOrder_x),
    polyOrder_y_(polyOrder_y),
    polyOrder_z_(polyOrder_z),
    pointType_(pointType)
    {
      this->functionSpace_ = FUNCTION_SPACE_HVOL;

      std::ostringstream basisName;
      basisName << "HVOL_HEX (" << this->TensorBasis::getName() << ")";
      name_ = basisName.str();
      
      this->setShardsTopologyAndTags();
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HVOL_HEX(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HVOL_HEX(polyOrder, polyOrder, polyOrder,pointType) {}
    
    /** \brief  Returns basis name

     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }

    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return false;
    }

    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator VALUE = Intrepid2::OPERATOR_VALUE;
      
      if (operatorType == VALUE)
      {
        return OperatorTensorDecomposition(VALUE, VALUE);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
      }
    }
    
    using TensorBasis::getValues;
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, typename BasisBase::OutputValueType, typename BasisBase::PointValueType>
    getHostBasis() const override {
      using HostBasisType  = Basis_Derived_HVOL_HEX<typename HVOL_LINE::HostBasis>;
      return Teuchos::rcp( new HostBasisType(polyOrder_x_, polyOrder_y_, polyOrder_z_, pointType_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HVOL_HEX_h */
