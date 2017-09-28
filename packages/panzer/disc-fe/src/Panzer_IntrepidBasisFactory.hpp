// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_INTREPID_BASIS_FACTORY_H
#define PANZER_INTREPID_BASIS_FACTORY_H

#include <sstream>
#include <string>
#include <map>

#include "Intrepid2_HVOL_C0_FEM.hpp"

#include "Teuchos_RCP.hpp"
#include "Intrepid2_Basis.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"

#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"

#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"

#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"

#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_In_FEM.hpp"

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"

#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"


namespace panzer {


  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName() 
                               after trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename ExecutionSpace, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<ExecutionSpace,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                       const shards::CellTopology & cell_topology)
  {
    // Shards supports extended topologies so the names have a "size"
    // associated with the number of nodes.  We prune the size to
    // avoid combinatorial explosion of checks.
    std::string cell_topology_type = cell_topology.getName();
    std::size_t end_position = 0;
    end_position = cell_topology_type.find("_");
    std::string cell_type = cell_topology_type.substr(0,end_position);
    
    Teuchos::RCP<Intrepid2::Basis<ExecutionSpace,OutputValueType,PointValueType> > basis;
    
    // high order point distribution type; 
    // for now equispaced only; to get a permutation map with different orientation, 
    // orientation coeff matrix should have the same point distribution.
    const Intrepid2::EPointType point_type = Intrepid2::POINTTYPE_EQUISPACED;

    if ( (basis_type == "Const") && (basis_order == 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_C0_FEM<ExecutionSpace,OutputValueType,PointValueType>(cell_topology) );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_C2_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Hexahedron") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    else if ( (basis_type == "HCurl") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_HEX_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Hexahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_HEX_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Hexahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_HEX_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Hexahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_HEX_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Tetrahedron") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TET_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Tetrahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TET_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Tetrahedron") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TET_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> ); 

    else if ( (basis_type == "HDiv") && (cell_type == "Tetrahedron") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TET_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) ); 

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Quadrilateral") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HCurl") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_QUAD_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Quadrilateral") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_QUAD_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HDiv") && (cell_type == "Quadrilateral") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_QUAD_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Quadrilateral") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_QUAD_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );
    
    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order == 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<ExecutionSpace,OutputValueType,PointValueType> );
    
    else if ( (basis_type == "HGrad") && (cell_type == "Triangle") && (basis_order > 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TRI_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    else if ( (basis_type == "HCurl") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TRI_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HCurl") && (cell_type == "Triangle") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HCURL_TRI_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    else if ( (basis_type == "HDiv") && (cell_type == "Triangle") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TRI_I1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HDiv") && (cell_type == "Triangle") && (basis_order > 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HDIV_TRI_In_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HGrad") && (cell_type == "Line") && (basis_order == 1) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_LINE_C1_FEM<ExecutionSpace,OutputValueType,PointValueType> );

    else if ( (basis_type == "HGrad") && (cell_type == "Line") && (basis_order >= 2) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );
    
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
                               "Failed to create the requestedbasis with basis_type=\"" << basis_type << 
                               "\", basis_order=\"" << basis_order << "\", and cell_type=\"" << cell_type << "\"!\n");

    TEUCHOS_TEST_FOR_EXCEPTION(cell_topology!=basis->getBaseCellTopology(),
                               std::runtime_error,
                               "Failed to create basis.  Intrepid2 basis topology does not match mesh cell topology!");

    return basis;
  }

  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName() after
                               trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename ExecutionSpace, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<ExecutionSpace,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                      const Teuchos::RCP<const shards::CellTopology> & cell_topology) 
  {
    return createIntrepid2Basis<ExecutionSpace,OutputValueType,PointValueType>(basis_type,basis_order,*cell_topology);
  }

}


#endif
