// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#include <percept/PerceptMesh.hpp>
#include <Intrepid2_CellTools.hpp>

#include <percept/element/intrepid/BasisTable.hpp>

  namespace percept {

    BasisTable* BasisTable::instance;
    void BasisTable::setupBasisTable()
    {
      /// these are the cell topologies currently supported for inverse mappings
      /// see Intrepid2::CellTools::mapToReferenceFrame documentation
      /**
         std::vector<shards::CellTopology> supportedTopologies;
         supportedTopologies.push_back(shards::getCellTopologyData<Triangle<3> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Triangle<6> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<4> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<9> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<4> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<10> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<8> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<27> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Wedge<6> >() );
      */

      // FIXME
      //#if !(defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))

      m_basisTable[shards::getCellTopologyData<shards::Line<2> >()-> key]          = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_LINE_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Line<3> >()-> key]          = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_LINE_C2_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Beam<3> >()-> key]          = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_LINE_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::Triangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Triangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<4> >()-> key] = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<8> >()-> key] = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_I2_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<9> >()-> key] = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<8> >()-> key]    = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<20> >()-> key]   = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_I2_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<27> >()-> key]   = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_HEX_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::Tetrahedron<4> >()-> key]   = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TET_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Tetrahedron<10> >()-> key]  = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TET_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::Wedge<6> >()-> key]         = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_WEDGE_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::Wedge<15> >()-> key]        = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_WEDGE_I2_FEM<Kokkos::HostSpace, double, double >() );


      // Shells
      m_basisTable[shards::getCellTopologyData<shards::ShellTriangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C1_FEM<Kokkos::HostSpace, double, double >() );
      m_basisTable[shards::getCellTopologyData<shards::ShellTriangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_TRI_C2_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()-> key] = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<Kokkos::HostSpace, double, double >() );

      m_basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<8> >()-> key] = Teuchos::rcp ( new Intrepid2::Basis_HGRAD_QUAD_I2_FEM<Kokkos::HostSpace, double, double>() );

      //#endif

      // etc....

      // FIXME
    }

    // static
    BasisTable::BasisTypeRCP BasisTable::
    getBasis(shards::CellTopology& topo)
    {
      unsigned key = topo.getKey();
      if (m_basisTable.size() == 0)
        {
          setupBasisTable();
        }
      std::string msg=std::string("No basis available for this topology")+topo.getName();
      if (m_basisTable.find(key) ==  m_basisTable.end()) std::cout << msg << std::endl;
      VERIFY_OP_ON( (m_basisTable.find(key) ==  m_basisTable.end()), == , false, "No basis available for this topology");

      BasisTable::BasisTypeRCP  basis =  m_basisTable[key];

      return basis;
    }
  }
