// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasisFamily.hpp
    \brief  Stateless classes that act as factories for two families of hierarchical bases.  HierarchicalBasisFamily provides bases associated with interface topologies (vertices, edges, and faces), while DGHierarchicalBasisFamily associates all members with element interiors.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_HierarchicalBasisFamily_h
#define Intrepid2_HierarchicalBasisFamily_h

#include "Intrepid2_DerivedBasisFamily.hpp"

#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_LegendreBasis_HVOL_LINE.hpp"

namespace Intrepid2 {
  // the following defines a family of hierarchical basis functions that matches the unpermuted ESEAS basis functions
  // each basis member is associated with appropriate subcell topologies, making this suitable for continuous Galerkin finite elements.
  
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
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using HierarchicalBasisFamily = DerivedBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<ExecutionSpace,OutputScalar,PointScalar,true>,
                                                      LegendreBasis_HVOL_LINE<ExecutionSpace,OutputScalar,PointScalar> >;
  
  /** \class Intrepid2::HierarchicalBasisFamily
      \brief A family of hierarchical basis functions, constructed in a way that follows work by Fuentes et al., suitable for use in DG contexts.
   
   This family is defined using the DerivedBasisFamily, which in turn templates the definitions on the H(grad) and H(vol) Lagrangian bases on the line.
   
   The suitability of this family for DG contexts is primarily due to the fact that the H(grad) basis has a constant member.  Note also that in this family,
   all members are associated with the cell interior; there are no basis functions associated with subcell topologies.
  */
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using DGHierarchicalBasisFamily = DerivedBasisFamily< IntegratedLegendreBasis_HGRAD_LINE<ExecutionSpace,OutputScalar,PointScalar,false>,
                                                        LegendreBasis_HVOL_LINE<ExecutionSpace,OutputScalar,PointScalar> >;
  
}

#endif /* Intrepid2_HierarchicalBasisFamily_h */
