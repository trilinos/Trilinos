// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_TopologyVerifier_hpp
#define stk_encr_TopologyVerifier_hpp

#include <functional>
#include <set>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include "Edge.hpp"
#include "ShardsInterfaceTable.hpp"


  namespace percept
  {

    class TopologyVerifier // consider: : MeshVerifier, MeshVerifier : Verifier
    {
      typedef std::set<MyEdge<unsigned> > invalid_edge_set_type;
      std::vector<invalid_edge_set_type > m_invalid_edge_set;
      //std::vector<std::set> m_valid_edge_set;
      enum { NELT = interface_table::NUM_ELEM_TYPES };
      void build_invalid_edge_sets();
    public:
      TopologyVerifier();
      bool isTopologyBad( stk::mesh::BulkData& bulk, stk::mesh::Entity elem);
      bool isTopologyBad( stk::mesh::BulkData& bulk);
    };

  }//namespace percept
#endif
