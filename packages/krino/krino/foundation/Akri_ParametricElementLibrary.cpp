// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_ParametricElement.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <map>
#include <utility>

namespace krino {

ParametricElementLibrary::ParametricElementLibrary()
  : myParametricElements(stk::topology::BEGIN_TOPOLOGY + stk::topology::NUM_TOPOLOGIES)
{
}

ParametricElementLibrary::~ParametricElementLibrary() = default;

const ParametricElement &
ParametricElementLibrary::get_parametric_element(stk::topology topology) const
{
  std::unique_ptr<ParametricElement> & parametricElem = myParametricElements[topology()];
  if (nullptr == parametricElem.get())
  {
    parametricElem = std::make_unique<ParametricElement>(topology);
  }
  return *parametricElem;
}

const ParametricElement &
ParametricElementLibrary::get_parametric_element(const stk::mesh::Bucket & bucket, const FieldRef field) const
{
  return get_parametric_element(get_field_topology(bucket, field));
}

stk::topology
ParametricElementLibrary::get_field_topology(const stk::mesh::Bucket & b, const FieldRef field) const
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

void
ParametricElementLibrary::clear()
{
  for (auto & parametricElem : myParametricElements)
    parametricElem.reset();
}

}  // namespace krino
