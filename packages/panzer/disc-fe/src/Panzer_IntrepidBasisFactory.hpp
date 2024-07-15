// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_INTREPID_BASIS_FACTORY_H
#define PANZER_INTREPID_BASIS_FACTORY_H

#include <sstream>
#include <string>
#include <map>

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"
#include "Intrepid2_HVOL_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_HVOL_TET_Cn_FEM.hpp"

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

    else if ( (basis_type == "HVol") && (basis_order == 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_C0_FEM<ExecutionSpace,OutputValueType,PointValueType>(cell_topology) );

    else if ( (basis_type == "HVol") && (cell_type == "Quadrilateral") && (basis_order > 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_QUAD_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Triangle") && (basis_order > 0) )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_TRI_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Hexahedron") )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_HEX_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

    else if ( (basis_type == "HVol") && (cell_type == "Tetrahedron") )
      basis = Teuchos::rcp( new Intrepid2::Basis_HVOL_TET_Cn_FEM<ExecutionSpace,OutputValueType,PointValueType>(basis_order, point_type) );

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
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_TET_C2_FEM<ExecutionSpace,OutputValueType,PointValueType> );

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

    // we compare that the base topologies are the same
    // we do so using the NAME. This avoids the ugly task of getting the
    // cell topology data and constructing a new cell topology object since you cant
    // just get the baseCellTopology directly from a shards cell topology 
    TEUCHOS_TEST_FOR_EXCEPTION(cell_topology.getBaseName()!=basis->getBaseCellTopology().getName(),
                               std::runtime_error,
                               "Failed to create basis.  Intrepid2 basis base topology does not match mesh cell base topology!");

    return basis;
  }

  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      // TODO BWR Is the cell_topology documentation below correct?
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
