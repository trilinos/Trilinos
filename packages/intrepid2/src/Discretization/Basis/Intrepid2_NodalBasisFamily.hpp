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

/** \file   Intrepid2_NodalBasisFamily.hpp
    \brief  Stateless class that acts as a factory for a family of nodal bases (hypercube topologies only at this point).  NodalBasisFamily is templated on H(vol) and H(grad) bases in a way that is more consistent with the literature and the hierarchical basis family in Intrepid2.  Once we support all standard topologies, we expect to replace the existing high-order nodal basis implementations in Intrepid2 with those from DerivedNodalBasisFamily.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_NodalBasisFamily_h
#define Intrepid2_NodalBasisFamily_h

#include "Intrepid2_DerivedBasisFamily.hpp"

#include <Intrepid2_HGRAD_LINE_Cn_FEM.hpp>
#include <Intrepid2_HVOL_LINE_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_QUAD_Cn_FEM.hpp>
#include <Intrepid2_HCURL_QUAD_In_FEM.hpp>
#include <Intrepid2_HDIV_QUAD_In_FEM.hpp>
#include <Intrepid2_HVOL_QUAD_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_TRI_Cn_FEM.hpp>
#include <Intrepid2_HCURL_TRI_In_FEM.hpp>
#include <Intrepid2_HDIV_TRI_In_FEM.hpp>
#include <Intrepid2_HVOL_TRI_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>
#include <Intrepid2_HCURL_HEX_In_FEM.hpp>
#include <Intrepid2_HDIV_HEX_In_FEM.hpp>
#include <Intrepid2_HVOL_HEX_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>
#include <Intrepid2_HCURL_TET_In_FEM.hpp>
#include <Intrepid2_HDIV_TET_In_FEM.hpp>
#include <Intrepid2_HVOL_TET_Cn_FEM.hpp>

namespace Intrepid2 {

  // the following defines a family of nodal basis functions for the reference Triangle
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true>
  class NodalTriangleBasisFamily
  {
  public:
    using HGRAD = Basis_HGRAD_TRI_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HCURL = Basis_HCURL_TRI_In_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HDIV  = Basis_HDIV_TRI_In_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HVOL  = Basis_HVOL_TRI_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>;
  };

  // the following defines a family of nodal basis functions for the reference Tetrahedron
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
  typename OutputScalar = double,
  typename PointScalar  = double,
  bool defineVertexFunctions = true>
  class NodalTetrahedronBasisFamily
  {
  public:
    using HGRAD = Basis_HGRAD_TET_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HCURL = Basis_HCURL_TET_In_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HDIV  = Basis_HDIV_TET_In_FEM<ExecutionSpace,OutputScalar,PointScalar>;
    using HVOL  = Basis_HVOL_TET_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>;
  };


  // the following defines a family of nodal basis functions, derived from a the standard high-order Intrepid2 bases on the line
  // note that because the standard H(curl) basis uses a lower-order H(grad) basis in place of the H(vol) that is arguably more natural,
  // the following will *not* match the standard nodal basis declared below.  (Similarly for H(div).)
  
  /** \class Intrepid2::DerivedNodalBasisFamily
      \brief A family of nodal basis functions which is related to, but not identical with, the Lagrangian basis family that Intrepid2 has historically supported.
   
   This family is defined using the DerivedBasisFamily, which in turn templates the definitions on the H(grad) and H(vol) Lagrangian bases on the line.
   
   At present, only hypercube topologies (line, quadrilateral, hexahedron) are supported, but other topologies will be supported in the future.
  */
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double>
  using DerivedNodalBasisFamily = DerivedBasisFamily< Basis_HGRAD_LINE_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>,
                                                      Basis_HVOL_LINE_Cn_FEM<ExecutionSpace,OutputScalar,PointScalar>,
                                                      NodalTriangleBasisFamily<ExecutionSpace,OutputScalar,PointScalar>,
                                                      NodalTetrahedronBasisFamily<ExecutionSpace,OutputScalar,PointScalar> >;
  
  /** \class Intrepid2::NodalBasisFamily
      \brief A family of nodal basis functions representing the higher-order Lagrangian basis family that Intrepid2 has historically supported.
   
   This family is defined with reference to the higher-order implementations (the "Cn" and "In" bases).
  */
  template<typename ExecSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double>
  class NodalBasisFamily
  {
  public:
    using ExecutionSpace  = ExecSpace;
    using OutputValueType = OutputScalar;
    using PointValueType  = PointScalar;
    
    using BasisType = Basis<ExecSpace,OutputScalar,PointScalar>;
    using BasisPtr  = Teuchos::RCP<BasisType>;
    
    // line bases
    using HGRAD_LINE = Basis_HGRAD_LINE_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HVOL_LINE  = Basis_HVOL_LINE_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    
    // quadrilateral bases
    using HGRAD_QUAD = Basis_HGRAD_QUAD_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HCURL_QUAD = Basis_HCURL_QUAD_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HDIV_QUAD  = Basis_HDIV_QUAD_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HVOL_QUAD  = Basis_HVOL_QUAD_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    
    // triangle bases
    using HGRAD_TRI = Basis_HGRAD_TRI_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HCURL_TRI = Basis_HCURL_TRI_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HDIV_TRI  = Basis_HDIV_TRI_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HVOL_TRI  = Basis_HVOL_TRI_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    
    // hexahedron bases
    using HGRAD_HEX = Basis_HGRAD_HEX_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HCURL_HEX = Basis_HCURL_HEX_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HDIV_HEX  = Basis_HDIV_HEX_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HVOL_HEX  = Basis_HVOL_HEX_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    
    // tetrahedron bases
    using HGRAD_TET = Basis_HGRAD_TET_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HCURL_TET = Basis_HCURL_TET_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HDIV_TET  = Basis_HDIV_TET_In_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    using HVOL_TET  = Basis_HVOL_TET_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>;
    
  };
}

#endif /* Intrepid2_NodalBasisFamily_h */
