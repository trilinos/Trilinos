// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <Akri_MasterElementBasis.hpp>
#include <Akri_DiagWriter.hpp>

namespace krino {

const MasterElement&
MasterElementDeterminer::getMasterElement(stk::mesh::Bucket & bucket, FieldRef field)
{
  stk::topology field_topology = get_field_topology(bucket, field);
  return MasterElementDeterminer::getMasterElement(field_topology);
}

stk::topology
MasterElementDeterminer::get_field_topology(const stk::mesh::Bucket & b, const FieldRef field)
{
  // As an optimization, assume there is only 1 field topology active for a given mesh topology and field
  typedef std::map<std::pair<unsigned, stk::topology>, stk::topology> FieldAndMeshTopoToFieldTopoMap;
  static FieldAndMeshTopoToFieldTopoMap field_topo_map;

  const stk::topology mesh_topology = b.topology();
  const unsigned field_ordinal = field.field().mesh_meta_data_ordinal();
  FieldAndMeshTopoToFieldTopoMap::iterator it = field_topo_map.find(std::make_pair(field_ordinal, mesh_topology));
  if (it != field_topo_map.end())
  {
    return it->second;
  }

  stk::mesh::MetaData & meta = field.field().mesh_meta_data();
  stk::topology field_topology = AuxMetaData::get(meta).get_nodal_field_topology(field, b);

  field_topo_map.insert(FieldAndMeshTopoToFieldTopoMap::value_type(std::make_pair(field_ordinal, mesh_topology), field_topology));
  return field_topology;
}

const MasterElement &
MasterElementDeterminer::getMasterElement(stk::topology t)
{
  static std::vector<std::unique_ptr<MasterElement>> all_master_elems(stk::topology::BEGIN_TOPOLOGY + stk::topology::NUM_TOPOLOGIES);
  std::unique_ptr<MasterElement> & master_elem = all_master_elems[t()];
  if (nullptr == master_elem.get())
  {
    std::unique_ptr<Basis> basis;
    switch(t())
    {
    case stk::topology::LINE_2:
        basis = std::make_unique<Basis_LINE_2>();
        break;
    case stk::topology::LINE_3:
        basis = std::make_unique<Basis_LINE_3>();
        break;
    case stk::topology::TRI_3:
    case stk::topology::TRI_3_2D:
        basis = std::make_unique<Basis_TRI_3>();
        break;
    case stk::topology::TRI_6:
    case stk::topology::TRI_6_2D:
        basis = std::make_unique<Basis_TRI_6>();
        break;
    case stk::topology::QUAD_4:
    case stk::topology::QUAD_4_2D:
        basis = std::make_unique<Basis_QUAD_4>();
        break;
    case stk::topology::QUAD_9:
    case stk::topology::QUAD_9_2D:
        basis = std::make_unique<Basis_QUAD_9>();
        break;
    case stk::topology::TET_4:
        basis = std::make_unique<Basis_TET_4>();
        break;
    case stk::topology::TET_10:
        basis = std::make_unique<Basis_TET_10>();
        break;
    case stk::topology::HEX_8:
        basis = std::make_unique<Basis_HEX_8>();
        break;
    case stk::topology::HEX_27:
        basis = std::make_unique<Basis_HEX_27>();
        break;
    case stk::topology::WEDGE_6:
        basis = std::make_unique<Basis_WEDGE_6>();
        break;
    default:
        ThrowRuntimeError("Element topology not found in MasterElementDeterminer::getMasterElement: " << t.name());
        break;
    }
    master_elem = std::make_unique<MasterElement>(t, std::move(basis));
  }
  return *master_elem;
}

} // namespace krino
