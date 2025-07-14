// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SubElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_io/IossBridge.hpp>

namespace krino {

void
debug_elem_parts_and_relations(const stk::mesh::BulkData & mesh, const Mesh_Element & elem)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Entity elem_obj = elem.entity();

  krinolog << "Elem: id=" << mesh.identifier(elem_obj) << "\n";
  krinolog << "  Mesh parts=";
  for(auto && elem_part : mesh.bucket(elem_obj).supersets())
  {
    krinolog << "\"" << elem_part->name() << "\"" << " ";
  }
  krinolog << "\n";

  krinolog << "  Nodes: ";
  for (auto && elem_node : elem.get_nodes())
  {
    krinolog << elem_node->entityId() << " ";
  }
  krinolog << "\n";

  const stk::topology topology = elem.coord_topology();
  std::vector<stk::mesh::Entity> side_nodes;
  for (unsigned s = 0; s < topology.num_sides(); ++s)
  {
    stk::mesh::Entity side_obj = find_entity_by_ordinal(mesh, elem_obj, meta.side_rank(), s);
    if (!mesh.is_valid(side_obj)) continue;
    krinolog << "  Elem side: side=" << s << ", id=" << mesh.identifier(side_obj);
    krinolog << ", side node ids: ";
    const stk::topology side_topology = topology.side_topology(s);
    const stk::mesh::Entity * elem_nodes_ptr = mesh.begin_nodes(elem_obj);
    side_nodes.resize(side_topology.num_nodes());
    topology.side_nodes(elem_nodes_ptr, s, side_nodes.begin());
    for (unsigned n=0; n<side_topology.num_nodes(); ++n)
    {
      stk::mesh::Entity side_node = side_nodes[n];
      krinolog << mesh.identifier(side_node) << " ";
    }
    krinolog << "\n";
    krinolog << "    Mesh parts=";
    for(auto && side_part : mesh.bucket(side_obj).supersets())
    {
      krinolog << "\"" << side_part->name() << "\"" << " ";
    }
    krinolog << "\n";
  }

  if (elem.have_subelements())
  {
    std::vector<const SubElement *> conformal_subelems;
    elem.get_subelements(conformal_subelems);

    for (auto && subelem : conformal_subelems)
    {
      stk::mesh::Entity subelem_obj = subelem->entity();
      if (!mesh.is_valid(subelem_obj))
      {
        // single conformal subelement
        if (1 == conformal_subelems.size())
        {
          subelem_obj = elem.entity();
        }
        if (!mesh.is_valid(subelem_obj))
        {
          // For debugging intermediate stage before subelement mesh objects are determined.
          krinolog << "  Subelem: id=-1" << "\n";
          continue;
        }
      }

      krinolog << "  Subelem: id=" << mesh.identifier(subelem_obj) << "\n";
      krinolog << "    Mesh parts=";
      for(auto && subelem_part : mesh.bucket(subelem_obj).supersets())
      {
        krinolog << "\"" << subelem_part->name() << "\"" << " ";
      }
      krinolog << "\n";

      krinolog << "    SubElemNodes: ";
      for (auto && sub_node : subelem->get_nodes())
      {
        krinolog << sub_node->entityId() << " ";
      }
      krinolog << "\n";

      for (unsigned s = 0; s < topology.num_sides(); ++s)
      {
        stk::mesh::Entity side_obj = find_entity_by_ordinal(mesh, subelem_obj, meta.side_rank(), s);
        if (!mesh.is_valid(side_obj)) continue;
        krinolog << "    Subelem side: side=" << s << ", id=" << mesh.identifier(side_obj);
        krinolog << ", side node ids: ";
        const stk::topology side_topology = topology.side_topology(s);
        const stk::mesh::Entity * elem_nodes_ptr = mesh.begin_nodes(subelem_obj);
        side_nodes.resize(side_topology.num_nodes());
        topology.side_nodes(elem_nodes_ptr, s, side_nodes.begin());
        for (unsigned n=0; n<side_topology.num_nodes(); ++n)
        {
          stk::mesh::Entity side_node = side_nodes[n];
          krinolog << mesh.identifier(side_node) << " ";
        }
        krinolog << "\n";
        krinolog << "      Mesh parts=";
        const stk::mesh::PartVector & side_parts = mesh.bucket(side_obj).supersets();
        for(stk::mesh::PartVector::const_iterator part_iter = side_parts.begin(); part_iter != side_parts.end(); ++part_iter)
        {
          const stk::mesh::Part * const part = *part_iter;
          krinolog << "\"" << part->name() << "\"" << " ";
        }
        krinolog << "\n";
      }
    }
  }
}

static double filter_negative_zero(const double val)
{
  if (val == 0.) return 0.;
  return val;
}

void
debug_nodal_parts_and_fields(const stk::mesh::BulkData & mesh, const SubElementNode * node)
{
  stk::mesh::Entity node_obj = node->entity();
  if (!mesh.is_valid(node_obj))
  {
    krinolog << "Node: identifier=nullptr\n";
    return;
  }

  std::string type = "UNKNOWN";
  if (nullptr != dynamic_cast<const SubElementSteinerNode *>(node)) type = "internal";
  if (nullptr != dynamic_cast<const SubElementEdgeNode *>(node)) type = "edge";
  if (nullptr != dynamic_cast<const SubElementMeshNode *>(node)) type = "mesh";

  krinolog << "Node: identifier=" << mesh.identifier(node_obj) << ", type=" << type << "\n";
  if (mesh.mesh_meta_data().spatial_dimension() == 2)
  {
    krinolog << "  coords =" << "(" << node->coordinates()[0] << ", " << node->coordinates()[1] << ")" << "\n";
  }
  else if (mesh.mesh_meta_data().spatial_dimension() == 3)
  {
    krinolog << "  coords =" << "(" << node->coordinates()[0] << ", " << node->coordinates()[1] << ", " << node->coordinates()[2] << ")" << "\n";
  }
  krinolog << "  Mesh parts=";
  for(auto && node_part : mesh.bucket(node_obj).supersets())
  {
    krinolog << "\"" << node_part->name() << "\"" << " ";
  }
  krinolog << "\n";

  for ( auto && stk_field : mesh.mesh_meta_data().get_fields() )
  {
    const FieldRef field(stk_field);

    if( field.entity_rank()!=stk::topology::NODE_RANK || !field.type_is<double>() ) continue;

    const unsigned field_length = field.length();

    double * data = field_data<double>(field, node_obj);
    if (nullptr == data)
    {
      krinolog << "  Field: " << field.name() << ", state=" << static_cast<int>(field.state()) << " -> NOT DEFINED." << "\n";
    }
    else
    {
      if (1 == field_length)
      {
        krinolog << "  Field: " << field.name() << ", state=" << static_cast<int>(field.state()) << ", value=" << filter_negative_zero(*data) << "\n";
      }
      else
      {
        krinolog << "  Field: " << field.name() << ", state=" << static_cast<int>(field.state()) << ", values[] = ";
        for (unsigned i = 0; i < field_length; ++i) krinolog << filter_negative_zero(data[i]) << " ";
        krinolog << "\n";
      }
    }
  }
}

void
debug_sides(const stk::mesh::BulkData & mesh, stk::mesh::Part & /*active_part*/)
{
  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_entities( mesh, mesh.mesh_meta_data().side_rank(), sides );

  for ( auto && side : sides )
    krinolog << debug_entity_1line(mesh, side) << "\n";
}

}
