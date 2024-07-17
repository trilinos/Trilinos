// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasisFamily.hpp
    \brief  Stateless classes that act as factories for two families of hierarchical bases.  HierarchicalBasisFamily provides bases associated with interface topologies (vertices, edges, and faces), while DGHierarchicalBasisFamily associates all members with element interiors.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_HierarchicalBasisFamily_h
#define Intrepid2_HierarchicalBasisFamily_h

#include "Intrepid2_DerivedBasisFamily.hpp"

#include "Intrepid2_HierarchicalBasis_HCURL_TRI.hpp"
#include "Intrepid2_HierarchicalBasis_HCURL_TET.hpp"
#include "Intrepid2_HierarchicalBasis_HDIV_TRI.hpp"
#include "Intrepid2_HierarchicalBasis_HDIV_TET.hpp"
#include "Intrepid2_HierarchicalBasis_HDIV_PYR.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_TRI.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_TET.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_PYR.hpp"
#include "Intrepid2_LegendreBasis_HVOL_LINE.hpp"
#include "Intrepid2_LegendreBasis_HVOL_TRI.hpp"
#include "Intrepid2_LegendreBasis_HVOL_TET.hpp"
#include "Intrepid2_LegendreBasis_HVOL_PYR.hpp"

namespace Intrepid2 {
  

//Dummy basis to be temporarily used for Hierarchical bases that have not been implemented yet
  template<typename ExecutionSpace, typename OutputScalar, typename PointScalar>
  class dummyBasis
  : public Basis<ExecutionSpace,OutputScalar,PointScalar> {
  public:
    dummyBasis(int /*order*/, EPointType /*pointType*/= POINTTYPE_DEFAULT) {};
  };

// the following defines a family of hierarchical basis functions that matches the unpermuted ESEAS basis functions
// each basis member is associated with appropriate subcell topologies, making this suitable for continuous Galerkin finite elements.
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true>
  class HierarchicalTriangleBasisFamily
  {
  public:
    // we will fill these in as we implement them
    using HGRAD = IntegratedLegendreBasis_HGRAD_TRI<DeviceType,OutputScalar,PointScalar,defineVertexFunctions>;
    using HCURL = HierarchicalBasis_HCURL_TRI<DeviceType,OutputScalar,PointScalar,defineVertexFunctions>; // last template argument: useCGBasis; corresponds with defineVertexFunctions.
    using HDIV  = HierarchicalBasis_HDIV_TRI<DeviceType,OutputScalar,PointScalar,defineVertexFunctions>; // last template argument: useCGBasis; corresponds with defineVertexFunctions.
    using HVOL  = LegendreBasis_HVOL_TRI<DeviceType,OutputScalar,PointScalar>;
  };
  
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true>
  class HierarchicalTetrahedronBasisFamily
  {
  public:
    // we will fill these in as we implement them
    using HGRAD = IntegratedLegendreBasis_HGRAD_TET<DeviceType,OutputScalar,PointScalar,defineVertexFunctions>;
    using HCURL = HierarchicalBasis_HCURL_TET<DeviceType,OutputScalar,PointScalar>;
    using HDIV  = HierarchicalBasis_HDIV_TET<DeviceType,OutputScalar,PointScalar>;
    using HVOL  = LegendreBasis_HVOL_TET<DeviceType,OutputScalar,PointScalar>;
  };


  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true>
  class HierarchicalPyramidBasisFamily
  {
  public:
    // we will fill these in as we implement them
    using HGRAD = IntegratedLegendreBasis_HGRAD_PYR<DeviceType,OutputScalar,PointScalar,defineVertexFunctions>;
    using HCURL = void;
    using HDIV  = HierarchicalBasis_HDIV_PYR<DeviceType,OutputScalar,PointScalar>;
    using HVOL  = LegendreBasis_HVOL_PYR<DeviceType,OutputScalar,PointScalar>;
  };
  
  /** \class Intrepid2::HierarchicalBasisFamily
      \brief A family of hierarchical basis functions, constructed in a way that follows work by Fuentes et al.
   
   This family is defined using the DerivedBasisFamily, which in turn templates the definitions on the H(grad) and H(vol) Lagrangian bases on the line.
   
   At present, only hypercube topologies (line, quadrilateral, hexahedron) are supported, but other topologies will be supported in the future.
   
   For mathematical details of the construction see:
   Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
   "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
   Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
   https://doi.org/10.1016/j.camwa.2015.04.027.
   
   Our implementation effectively differs from ESEAS only in the fact that we use only the default orientation
   on the cell, because Intrepid2 has other mechanisms for handling basis orientations.  At present, we have implemented
   the basis functions on hypercube topologies.  We plan to add support for simplices, wedges, and pyramids soon.
   
   We have offline tests that verify agreement between our implementation and ESEAS.  We hope to add these to the
   Trilinos continuous integration tests in the future.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using HierarchicalBasisFamily = DerivedBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar,true>,
                                                      LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>,
                                                      HierarchicalTriangleBasisFamily<DeviceType,OutputScalar,PointScalar>,
                                                      HierarchicalTetrahedronBasisFamily<DeviceType,OutputScalar,PointScalar>,
                                                      HierarchicalPyramidBasisFamily<DeviceType,OutputScalar,PointScalar>
                                                      >;
  
  /** \class Intrepid2::HierarchicalBasisFamily
      \brief A family of hierarchical basis functions, constructed in a way that follows work by Fuentes et al., suitable for use in DG contexts.
   
   This family is defined using the DerivedBasisFamily, which in turn templates the definitions on the H(grad) and H(vol) Lagrangian bases on the line.
   
   The suitability of this family for DG contexts is primarily due to the fact that the H(grad) basis has a constant member.  Note also that in this family,
   all members are associated with the cell interior; there are no basis functions associated with subcell topologies.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using DGHierarchicalBasisFamily = DerivedBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar,false>,
                                                        LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>,
                                                        HierarchicalTriangleBasisFamily<DeviceType,OutputScalar,PointScalar,false>,
                                                        HierarchicalTetrahedronBasisFamily<DeviceType,OutputScalar,PointScalar,false>,
                                                        HierarchicalPyramidBasisFamily<DeviceType,OutputScalar,PointScalar,false>
                                                      >;
  
}

#endif /* Intrepid2_HierarchicalBasisFamily_h */
