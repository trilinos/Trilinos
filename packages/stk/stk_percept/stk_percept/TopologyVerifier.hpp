#ifndef stk_encr_TopologyVerifier_hpp
#define stk_encr_TopologyVerifier_hpp

#include <functional>
#include <set>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include "Edge.hpp"
#include "ShardsInterfaceTable.hpp"

using namespace stk::mesh;

namespace stk
{
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
      bool isTopologyBad( mesh::Entity  &elem);
      bool isTopologyBad( mesh::BulkData& mesh_bulk_data);
    };

  }//namespace percept
}//namespace stk
#endif
