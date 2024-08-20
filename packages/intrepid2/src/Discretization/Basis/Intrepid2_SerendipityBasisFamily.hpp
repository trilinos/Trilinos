// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_SerendipityBasisFamily.hpp
    \brief  Stateless classes that act as factories for two families of hierarchical bases.  SerendipityBasisFamily provides bases associated with interface topologies (vertices, edges, and faces), while DGSerendipityBasisFamily associates all members with element interiors.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_SerendipityBasisFamily_h
#define Intrepid2_SerendipityBasisFamily_h

#include "Intrepid2_Basis.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HCURL_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HDIV_QUAD.hpp"
#include "Intrepid2_DerivedBasis_HVOL_QUAD.hpp"

#include "Intrepid2_DerivedBasis_HGRAD_HEX.hpp"
#include "Intrepid2_DerivedBasis_HCURL_HEX.hpp"
#include "Intrepid2_DerivedBasis_HDIV_HEX.hpp"
#include "Intrepid2_DerivedBasis_HVOL_HEX.hpp"

namespace Intrepid2 {
/** \class Intrepid2::SerendipityBasisWrapper
    \brief Helper class that allows SerendipityBasis construction with poly order arguments that are passed to the tensor-basis constructor.  (SerendipityBasis itself requires a BasisPtr at construction.)
*/
  template<class FullBasis, int numPolyOrderArgs>
  class SerendipityBasisWrapper
  :
  public SerendipityBasis<typename FullBasis::BasisBase>
  {
  public:
    using BasisBase = typename FullBasis::BasisBase;
    using BasisPtr  = Teuchos::RCP<BasisBase>;
    using DeviceType = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    using OrdinalTypeArray1D     = typename BasisBase::OrdinalTypeArray1D;
    using OrdinalTypeArray1DHost = typename BasisBase::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename BasisBase::OrdinalTypeArray2DHost;
    using OutputViewType         = typename BasisBase::OutputViewType;
    using PointViewType          = typename BasisBase::PointViewType;
    
    //! single-argument constructor, for isotropic bases.
    SerendipityBasisWrapper(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    SerendipityBasis<typename FullBasis::BasisBase>(Teuchos::rcp(new FullBasis(polyOrder)) )
    {}
    
    //! two-argument constructor; enabled if numPolyOrderArgs is 2.
    template<bool M=(numPolyOrderArgs==2)>
    SerendipityBasisWrapper(int polyOrder_x, int polyOrder_y, const EPointType pointType=POINTTYPE_DEFAULT, typename std::enable_if<M>::type* = 0)
    :
    SerendipityBasis<typename FullBasis::BasisBase>(Teuchos::rcp(new FullBasis(polyOrder_x, polyOrder_y)) )
    {}
    
    //! three-argument constructor; enabled if numPolyOrderArgs is 3.
    template<bool M=(numPolyOrderArgs==3)>
    SerendipityBasisWrapper(int polyOrder_x, int polyOrder_y, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT, typename std::enable_if<M>::type* = 0)
    :
    SerendipityBasis<typename FullBasis::BasisBase>(Teuchos::rcp(new FullBasis(polyOrder_x, polyOrder_y, polyOrder_z)) )
    {}
  };

/** \class Intrepid2::SerendipityBasisFamily
    \brief Serendipity basis family constructed in terms of arbitrary bases on the line, triangle, and tetrahedron.  (These must be hierarchical bases.)
*/
  template<class LineBasisHGRAD, class LineBasisHVOL, class TriangleBasisFamily, class TetrahedronBasisFamily >
  class DerivedSerendipityBasisFamily
  {
  public:
    using ExecutionSpace  = typename LineBasisHGRAD::ExecutionSpace;
    using OutputValueType = typename LineBasisHGRAD::OutputValueType;
    using PointValueType  = typename LineBasisHGRAD::PointValueType;
    
    using Basis      = typename LineBasisHGRAD::BasisBase;
    using BasisPtr   = Teuchos::RCP<Basis>;
    using DeviceType = typename Basis::DeviceType;
    
    // line bases
    using HGRAD_LINE = LineBasisHGRAD;
    using HVOL_LINE  = LineBasisHVOL;
    
    // quadrilateral bases
    using HGRAD_QUAD = SerendipityBasisWrapper<Basis_Derived_HGRAD_QUAD<HGRAD_LINE>,            2>;
    using HCURL_QUAD = SerendipityBasisWrapper<Basis_Derived_HCURL_QUAD<HGRAD_LINE, HVOL_LINE>, 2>;
    using HDIV_QUAD  = SerendipityBasisWrapper<Basis_Derived_HDIV_QUAD <HGRAD_LINE, HVOL_LINE>, 2>;
    using HVOL_QUAD  = SerendipityBasisWrapper<Basis_Derived_HVOL_QUAD <HVOL_LINE>,             2>;
    
    // hexahedron bases
    using HGRAD_HEX = SerendipityBasisWrapper<Basis_Derived_HGRAD_HEX<HGRAD_LINE>,            3>;
    using HCURL_HEX = SerendipityBasisWrapper<Basis_Derived_HCURL_HEX<HGRAD_LINE, HVOL_LINE>, 3>;
    using HDIV_HEX  = SerendipityBasisWrapper<Basis_Derived_HDIV_HEX <HGRAD_LINE, HVOL_LINE>, 3>;
    using HVOL_HEX  = SerendipityBasisWrapper<Basis_Derived_HVOL_HEX <HVOL_LINE>,             3>;
    
    // triangle bases
    using HGRAD_TRI = typename TriangleBasisFamily::HGRAD;
    using HCURL_TRI = typename TriangleBasisFamily::HCURL;
    using HDIV_TRI  = typename TriangleBasisFamily::HDIV;
    using HVOL_TRI  = typename TriangleBasisFamily::HVOL;
    
    // tetrahedron bases
    using HGRAD_TET = typename TetrahedronBasisFamily::HGRAD;
    using HCURL_TET = typename TetrahedronBasisFamily::HCURL;
    using HDIV_TET  = typename TetrahedronBasisFamily::HDIV;
    using HVOL_TET  = typename TetrahedronBasisFamily::HVOL;
  };

/** \class Intrepid2::SerendipityBasisFamily
      \brief Serendipity basis family constructed using the standard (continuous) hierarchical basis family.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using SerendipityBasisFamily = DerivedSerendipityBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar,true>,
                                                                LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>,
                                                                HierarchicalTriangleBasisFamily<DeviceType,OutputScalar,PointScalar,true>,
                                                                HierarchicalTetrahedronBasisFamily<DeviceType,OutputScalar,PointScalar,true>
                                                              >;

  /** \class Intrepid2::DGSerendipityBasisFamily
      \brief Serendipity basis family constructed using the DG hierarchical basis family.
   
   The suitability of this family for DG contexts is primarily due to the fact that the H(grad) basis has a constant member.  Note also that in this family,
   all members are associated with the cell interior; there are no basis functions associated with subcell topologies.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using DGSerendipityBasisFamily = DerivedSerendipityBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar,false>,
                                                                  LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>,
                                                                  HierarchicalTriangleBasisFamily<DeviceType,OutputScalar,PointScalar,false>,
                                                                  HierarchicalTetrahedronBasisFamily<DeviceType,OutputScalar,PointScalar,false>
                                                                >;
}

#endif /* Intrepid2_HierarchicalBasisFamily_h */
