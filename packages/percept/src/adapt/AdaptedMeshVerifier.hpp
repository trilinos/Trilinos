// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AdaptedMeshVerifier_hpp
#define AdaptedMeshVerifier_hpp

#include <iostream>


#include "Shards_CellTopology.hpp"
//#include "Teuchos_GlobalMPISession.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <percept/PerceptMesh.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

using namespace shards;


  namespace percept
  {

    template<typename T> class Histogram;

    class AdaptedMeshVerifier
    {
      int m_debug;
      double m_initialTotalVolume;
      double m_initialTotalSurfaceArea;
    public:
      AdaptedMeshVerifier(bool debug=false);

      bool isValid(percept::PerceptMesh& eMesh, bool isInitialMesh);

      // for internal use, but could be useful
      double totalVolume(percept::PerceptMesh& eMesh, stk::mesh::EntityRank rank, bool exclude_parents=true);
      bool checkParentChildVol(PerceptMesh& eMesh, bool debug);

      bool hasHangingNodes(percept::PerceptMesh& eMesh);
      void check_mesh_parallel_consistency (const stk::mesh::Ghosting& ghosts, stk::mesh::EntityRank entity_rank, const std::string& msg="");
      bool is_valid_relations_and_entities(percept::PerceptMesh& eMesh, stk::mesh::EntityRank rank, bool exclude_parents);
      static void check_parent_element_field(PerceptMesh& eMesh, const std::string& msg1="", bool debug=false);
      static void checkPolarity(PerceptMesh& eMesh);
    };

  }

#endif
