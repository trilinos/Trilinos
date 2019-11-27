// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/SmoothingHelperFunctions.hpp>

namespace percept
{
    static void get_nodes_on_side(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, unsigned element_side_ordinal, std::vector<stk::mesh::Entity>& node_vector)
    {
      stk::mesh::EntityRank side_entity_rank = bulkData.mesh_meta_data().side_rank();

      const CellTopologyData * const element_topo_data = stk::mesh::get_cell_topology(bulkData.bucket(element).topology()).getCellTopologyData();
      shards::CellTopology element_topo(element_topo_data);
      const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );
      const unsigned *  inodes = 0;
      unsigned n_elem_side_nodes = 0;

      if (side_entity_rank == stk::topology::EDGE_RANK)
        {
          inodes = element_topo_data->edge[element_side_ordinal].node;
          n_elem_side_nodes = 2;
        }
      else if (side_entity_rank == stk::topology::FACE_RANK )
        {
          n_elem_side_nodes = element_topo_data->side[element_side_ordinal].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = element_topo_data->side[element_side_ordinal].node;
        }

      for (unsigned knode = 0; knode < n_elem_side_nodes; knode++)
        {
          node_vector.push_back(elem_nodes[inodes[knode]].entity());
        }
    }

  stk::mesh::Part* get_skin_part(stk::mesh::BulkData * bulkData, const std::string& part_name, bool remove_previous_part_nodes)
  {
    stk::mesh::MetaData * metaData = &(bulkData->mesh_meta_data());

    stk::mesh::Part* part = metaData->get_part(part_name);
    VERIFY_OP_ON(part, !=, 0, "Need to call add_inner_skin_part first - no available inner skin part");

    if (remove_previous_part_nodes)
      {
        stk::mesh::Selector on_skin_part(*part);
        std::vector<stk::mesh::Entity> nodes;
        const stk::mesh::BucketVector & buckets = bulkData->buckets(stk::topology::NODE_RANK);
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (bucket.owned()
                && on_skin_part(bucket))
              {
                const unsigned num_nodes_in_bucket = bucket.size();
                for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                  {
                    stk::mesh::Entity node = bucket[iNode];
                    if (bulkData->is_valid(node) && std::find(nodes.begin(),nodes.end(),node)==nodes.end())
                      nodes.push_back(node);
                  }
              }
          }

        std::vector<stk::mesh::PartVector> add_parts(nodes.size()), remove_parts(nodes.size(),stk::mesh::PartVector(1,part));
        bulkData->batch_change_entity_parts(nodes, add_parts, remove_parts);
      }

    stk::mesh::EntitySideVector boundary;

    // select owned
    stk::mesh::Selector owned = metaData->locally_owned_part();

    const stk::mesh::PartVector parts = metaData->get_parts();
    for (unsigned ip=0; ip < parts.size(); ip++)
      {
        bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);
        if (stk_auto) continue;
        unsigned per = parts[ip]->primary_entity_rank();
        if (per == stk::topology::ELEMENT_RANK)
          {
            const stk::topology topology = metaData->get_topology(*parts[ip]);
            if (!topology.is_valid() || topology.dimension() != per)
              {
                std::cout << "Warning: PerceptMesh::get_skin_part: skipping part with dimension < element_rank, part name= " << parts[ip]->name() << std::endl;
                continue;
              }
            //std::cout << "INFO::smoothing: freezing points on boundary: " << parts[ip]->name() << std::endl;
            stk::mesh::EntityVector owned_elements;

            stk::mesh::Selector block(*parts[ip]);
            block = block & owned;
            get_selected_entities( block,
                                   bulkData->buckets(stk::topology::ELEMENT_RANK),
                                   owned_elements);
            //Part * skin_part = 0;
            stk::mesh::EntityVector elements_closure;

            // compute owned closure
            find_closure( *bulkData, owned_elements, elements_closure );

            // compute boundary
            boundary_analysis( *bulkData, elements_closure, stk::topology::ELEMENT_RANK, boundary);
          }
      }

    std::vector<stk::mesh::Entity> nodes;

    std::vector<stk::mesh::Entity> node_vector;
    for (unsigned iesv=0; iesv < boundary.size(); ++iesv)
      {
        stk::mesh::EntitySide& es = boundary[iesv];
        node_vector.resize(0);
        if (bulkData->is_valid(es.inside.entity))
          get_nodes_on_side(*bulkData, es.inside.entity, es.inside.side_ordinal, node_vector);
        if (bulkData->is_valid(es.outside.entity))
          get_nodes_on_side(*bulkData, es.outside.entity, es.outside.side_ordinal, node_vector);
        for (unsigned inv=0; inv < node_vector.size(); inv++)
          {
            if (bulkData->bucket(node_vector[inv]).owned() && std::find(nodes.begin(),nodes.end(),node_vector[inv])==nodes.end())
              nodes.push_back(node_vector[inv]);
          }
      }
    std::vector<stk::mesh::PartVector> add_parts(nodes.size(),stk::mesh::PartVector(1,part)), remove_parts(nodes.size());
    bulkData->batch_change_entity_parts(nodes, add_parts, remove_parts );

    return part;
  }
}
