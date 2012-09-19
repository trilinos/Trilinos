#ifndef STK_SIERRA_IO_MODIFIABLE_MESH_READER_HPP
#define STK_SIERRA_IO_MODIFIABLE_MESH_READER_HPP

#include <Ioss_SubSystem.h>

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/entity_rank.hpp>
#include <sierra/mesh/details/entity_id.hpp>
#include <sierra/mesh/details/part_key.hpp>
#include <sierra/mesh/details/cell_topology.hpp>
#include <sierra/io/ioss_topology.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>

#include <vector>
#include <string>

namespace sierra {
  namespace mesh {
    namespace io {

      typedef details::entity_key entity_key;
      typedef std::vector<entity_key> key_vector;
      typedef std::pair<Ioss::EntityBlock*,details::entity_rank> io_part;
      typedef details::entity_property entity_property;
      typedef details::entity_rank entity_rank;
      typedef details::entity_id entity_id;
      typedef details::relation_position relation_position;
      typedef details::selector selector;

      class ModifiableMeshReader
      {
      public:
        ModifiableMeshReader(Ioss::Region *io_region, modifiable_mesh &mesh);

        void process_mesh();

        template <class Field, typename IossPartIterator>
        void process_field(Field &field, IossPartIterator part_begin, IossPartIterator part_end,
			   const std::string &io_field_name);
      private:
        template <typename T>
           void process_blocks(const std::vector<T*> &blocks, details::entity_rank rank);

        void process_block(Ioss::EntityBlock *block, details::entity_rank rank);

        Ioss::Region *m_ioRegion;

        std::map<entity_rank,key_vector> m_entityKeys;
        modifiable_mesh &m_mesh;

        int m_meshRankProcessed;
      };

      ModifiableMeshReader::ModifiableMeshReader(Ioss::Region *io_region, modifiable_mesh &mesh) :
        m_ioRegion(io_region), m_mesh(mesh), m_meshRankProcessed(-1)
      {}

      void ModifiableMeshReader::process_mesh()
      {
        assert(m_meshRankProcessed == -1);

        m_mesh.set_num_entity_ranks(4); // TODO: Not in generic API yet, remove hardwired value.

        m_entityKeys[m_mesh.node_rank()].resize(m_ioRegion->get_property("node_count").get_int());
        Ioss::NodeBlockContainer node_blocks = m_ioRegion->get_node_blocks();
        process_blocks(node_blocks, m_mesh.node_rank());

        m_entityKeys[m_mesh.element_rank()].resize(m_ioRegion->get_property("element_count").get_int());
        Ioss::ElementBlockContainer element_blocks = m_ioRegion->get_element_blocks();
        process_blocks(element_blocks, m_mesh.element_rank());
      }

      template <typename T>
      void ModifiableMeshReader::process_blocks(const std::vector<T*> &blocks, details::entity_rank rank)
      {
        assert(m_meshRankProcessed < (int)(size_t)rank);
        m_meshRankProcessed = rank;

        for(size_t i=0; i < blocks.size(); i++) {
          T* block = blocks[i];
          process_block(block, rank);
        }
      }

      void ModifiableMeshReader::process_block(Ioss::EntityBlock *block, details::entity_rank rank)
      {
        std::vector<int> global_ids;
        std::vector<int> connectivity;
        block->get_field_data("ids", global_ids);
        block->get_field_data("connectivity", connectivity);

        const Ioss::ElementTopology *topology = block->topology();
        const sierra::mesh::details::cell_topology cell_topo(sierra::mesh::io::map_topology_ioss_to_cell(topology));

        details::part_key parts[3];

        std::ostringstream osstr;
        osstr << "rank_"<<rank<<"_part";

        parts[0] = m_mesh.declare_part(osstr.str(),       rank);
        parts[1] = m_mesh.declare_part(topology->name(),  cell_topo);
        parts[2] = m_mesh.declare_part(block->name(),     io_part(block,rank));

        sierra::mesh::details::part_property& part_prop = m_mesh.get(parts[2]);
        part_prop.add_property(cell_topo);
        part_prop.add_property(rank);

        const int nodes_per_entity = topology->number_nodes();
        const int entity_count = block->get_property("entity_count").get_int();
        const int offset = block->get_offset();

        key_vector &nodes = m_entityKeys[m_mesh.node_rank()];
        key_vector &from  = m_entityKeys[rank];

        for (int i=0; i<entity_count; i++) {
          entity_key entity_from = m_mesh.add_entity(entity_property(rank, entity_id(global_ids[i])));
          assert(offset+i < (int)from.size());
          from[offset+i] = entity_from;
          m_mesh.change_entity_parts(entity_from, parts, parts+3);

          // Setting up connectivity relations from entity to nodes
          if (rank > m_mesh.node_rank()) { // TODO:  See if better way to do this.
            for (int j=0; j<nodes_per_entity; j++) {
              int local_id = connectivity[i*nodes_per_entity + j]-1;
              assert((local_id >= 0) && (local_id < static_cast<int>(nodes.size())));
              entity_key entity_to = nodes[local_id];
              m_mesh.add_relation(entity_from, entity_to, relation_position(m_mesh.node_rank(), j));
              m_mesh.add_relation(entity_to, entity_from, relation_position(rank));

              // move entity_from to blockPart, topologyPart?
              m_mesh.change_entity_parts(entity_to, parts+1, parts+3);
            }
          }
        }
      }

      template <class Field, typename IossPartIterator>
      void ModifiableMeshReader::process_field(Field &field, IossPartIterator part_begin, IossPartIterator part_end, const std::string &io_field_name)
      {
        assert(m_meshRankProcessed == (int)(size_t)m_mesh.element_rank());

        std::vector<details::part_key> io_parts;
        while (part_begin != part_end) {
          io_parts.push_back(m_mesh.find_part((*part_begin)->name()));
          //\TODO should assert here that the returned part-key is valid and is io-part
          part_begin++;
        }

        selector io_part_selector = selectUnion(io_parts.begin(), io_parts.end());
        field.update_from_mesh(io_part_selector, m_mesh);

        BOOST_FOREACH(details::part_key key, io_parts) {
          if (m_mesh[key].has_property<io_part>()) {
            Ioss::EntityBlock *block = m_mesh[key].get_property<io_part>().first;
            details::entity_rank rank = m_mesh[key].get_property<io_part>().second;
            assert(block->field_exists(io_field_name));
            const Ioss::Field &io_field = block->get_fieldref(io_field_name);
            int component_count = io_field.transformed_storage()->component_count();
            assert(component_count == Field::field_length_per_entity); // TODO: Fix

            key_vector &entities  = m_entityKeys[rank];

            std::vector<double> field_data;
            block->get_field_data(io_field_name, field_data);
            const int entity_count = block->get_property("entity_count").get_int();
            const int offset = block->get_offset();

            for (int i=0; i < entity_count; i++) {
              entity_key entity = entities[offset+i];
              typename field_traits<Field>::data_type *data = field[m_mesh.get_bucket_location(entity)];
              int dbegin = (i+offset)*component_count;
              int dend   = dbegin + component_count;
              std::copy(&field_data[dbegin], &field_data[dend],data);
            }
          }
        }
      }

    } // io namespace
  }  // mesh namespace
} // stk namespace

#endif

