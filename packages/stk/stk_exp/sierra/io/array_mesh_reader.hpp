#ifndef STK_SIERRA_IO_ARRAY_MESH_READER_HPP
#define STK_SIERRA_IO_ARRAY_MESH_READER_HPP

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include <sierra/mesh/array_mesh/array_mesh.hpp>
#include <sierra/io/array_mesh_ioss_topology.hpp>

#include <vector>
#include <string>

namespace sierra {
namespace mesh {
namespace io {

      class array_mesh_reader
      {
      public:
    	  /** this constructor doesn't call process_mesh(), it is assumed that the
    	   * calling code will call that...
    	   */
        array_mesh_reader(Ioss::Region *io_region, array_mesh &mesh);

        /** mesh_type can be 'exodus' or 'generated'.
         * if mesh_type is generated, then file_name is mesh-spec e.g. '2x1x1|sideset:xX'.
         *
         * This constructor calls process_mesh() which loads mesh connectivity, sidesets, etc.
         */
        array_mesh_reader(MPI_Comm comm,
        				const std::string& mesh_type, const std::string& file_name,
        				array_mesh& mesh);

        ~array_mesh_reader()
        {
        	if (m_created_io_region) delete m_io_region;
        }

        void process_mesh();

        template <typename IossPartIterator>
        void process_field(std::vector<double>& field, IossPartIterator part_begin, IossPartIterator part_end,
			   const std::string &io_field_name);

        void read_nodal_field(std::vector<double>& field, const std::string& name);
        void read_nodal_field(std::vector<double>& field, const std::string& name, double time);

        /** power-users may want access to the Ioss::Region...
         */
        Ioss::Region* get_ioss_region() { return m_io_region; }

      private:
        template <typename T>
           void process_blocks(const std::vector<T*> &blocks, int rank);

        void process_block(Ioss::EntityBlock *block, int rank);
  
        void process_sidesets();

        void process_sideset(const Ioss::SideSet* sset);

        void process_nodesets();
        void process_nodeset(const Ioss::NodeSet* nset);

        Ioss::Region *m_io_region;
        bool m_created_io_region;

        array_mesh &m_mesh;

        int m_meshRankProcessed;
      };

      inline array_mesh_reader::array_mesh_reader(Ioss::Region *io_region, array_mesh &mesh) :
        m_io_region(io_region), m_created_io_region(false), m_mesh(mesh), m_meshRankProcessed(-1)
      {}

      inline array_mesh_reader::array_mesh_reader(MPI_Comm comm,
    		  	  	  	  	  	  	  	  	  const std::string& mesh_type, const std::string& file_name,
    		  	  	  	  	  	  	  	  	  array_mesh &mesh) :
         m_io_region(NULL), m_created_io_region(true), m_mesh(mesh), m_meshRankProcessed(-1)
       {
    	  Ioss::Init::Initializer init_db;
    	  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, comm);

    	  //Ioss::Region takes ownership of dbi, so we don't need to delete it in this scope.
    	  m_io_region = new Ioss::Region(dbi, "input_model");

    	  process_mesh();
       }

       inline void array_mesh_reader::process_mesh()
      {
        assert(m_meshRankProcessed == -1);
  
        int num_nodes = m_io_region->get_property("node_count").get_int();
        int num_elems = m_io_region->get_property("element_count").get_int();
        m_mesh.reserve(num_nodes, num_elems);

        Ioss::NodeBlockContainer node_blocks = m_io_region->get_node_blocks();
        process_blocks(node_blocks, array_mesh::Node);
  
        Ioss::ElementBlockContainer element_blocks = m_io_region->get_element_blocks();
        process_blocks(element_blocks, array_mesh::Element);

        process_sidesets();
        process_nodesets();
      }

      template <typename T>
      inline
      void array_mesh_reader::process_blocks(const std::vector<T*> &blocks, int rank)
      {
        assert(m_meshRankProcessed < rank);
        m_meshRankProcessed = rank;

        for(size_t i=0; i < blocks.size(); i++) {
          T* block = blocks[i];
          process_block(block, rank);
        }
      }

      inline void array_mesh_reader::process_block(Ioss::EntityBlock *block, int rank)
      {
        std::vector<int> global_ids;
        std::vector<int> connectivity_raw;
        block->get_field_data("ids", global_ids);
        block->get_field_data("connectivity_raw", connectivity_raw);
        int num_elems = global_ids.size();
        int block_id = block->get_property("id").get_int();
  
        if (rank == stk::topology::ELEMENT_RANK) m_mesh.add_element_ids(global_ids.begin(), global_ids.end());
        if (rank == stk::topology::NODE_RANK) m_mesh.add_node_ids(global_ids.begin(), global_ids.end());

        const Ioss::ElementTopology *topology = block->topology();
        const std::string& name = block->name();
        std::string topo_name = topology->name();
        int num_nodes = topology->number_nodes();

//        int array_mesh_topo = map_ioss_topology_to_array_mesh(topo_name, num_nodes);
//
//        stk::topology t;
//        switch(array_mesh_topo) {
//        case sierra::mesh::Tet4::value:
//          t = stk::topology::TET_4;
//          break;
//        case sierra::mesh::Hex8::value:
//          t = stk::topology::HEX_8;
//          break;
//        case sierra::mesh::Node::value:
//          t = stk::topology::NODE;
//          break;
//        default:
//          std::cout << "Unsupported Topology" << std::endl;
//          throw std::runtime_error("array_mesh_reader ERROR unsupported topology");
//        };
//
//        array_mesh::BlockIndex blk = m_mesh.add_block(rank, block_id, num_elems,t, name);

        int offset = 0;

        stk::topology topo = map_ioss_topology_to_array_mesh(topo_name, num_nodes);
        array_mesh::BlockIndex blk = m_mesh.add_block(rank, block_id, num_elems, topo, name);

//Exodus provides 1-based connectivities, but we want 0-based so
//we will now subtract 1 from them.
        if (rank == sierra::mesh::array_mesh::Element) {
          for(size_t i=0; i<connectivity_raw.size(); ++i) connectivity_raw[i] -= 1;
        }

        const int block_offset = block->get_offset();

        for (int i=0; i<num_elems; i++) {
          m_mesh.add_connectivity(blk, block_offset + i, &connectivity_raw[offset], &connectivity_raw[offset]+num_nodes);
          offset += num_nodes;
        } 
      }

      inline void array_mesh_reader::process_sidesets()
      {
        const Ioss::SideSetContainer& side_sets = m_io_region->get_sidesets();

        for(Ioss::SideSetContainer::const_iterator sideset_iterator=side_sets.begin(), sideset_end=side_sets.end(); sideset_iterator != sideset_end; ++sideset_iterator)
        {
          process_sideset(*sideset_iterator);
        }
      }

      inline void array_mesh_reader::process_sideset(const Ioss::SideSet* sset)
      {
        size_t block_count = sset->block_count();
        const std::string& name = sset->name();
        int sideset_id = sset->get_property("id").get_int();
        array_mesh::SidesetIndex sideset_index = m_mesh.add_sideset(sideset_id, name);

        for(size_t i=0; i<block_count; ++i) {
          Ioss::SideBlock* block = sset->get_block(i);
          std::vector<int> elem_side;
          block->get_field_data("element_side_raw", elem_side);

          size_t ii=0;
          while(ii<elem_side.size()) {
        	//exodus provides 1-based numbers, but array_mesh uses 0-based. So subtract 1 here:
            int elem_num = elem_side[ii++] - 1;
            int local_side = elem_side[ii++] - 1;
            m_mesh.add_side(sideset_index, elem_num, local_side);
          }
        }
      }

      inline void array_mesh_reader::process_nodesets()
      {
        const Ioss::NodeSetContainer& node_sets = m_io_region->get_nodesets();

        for(Ioss::NodeSetContainer::const_iterator nodeset_iterator=node_sets.begin(), nodeset_end=node_sets.end(); nodeset_iterator != nodeset_end; ++nodeset_iterator)
        {
          process_nodeset(*nodeset_iterator);
        }
      }

      inline void array_mesh_reader::process_nodeset(const Ioss::NodeSet* nset)
      {
        const std::string& name = nset->name();
        int nodeset_id = nset->get_property("id").get_int();
        array_mesh::NodesetIndex nodeset_index = m_mesh.add_nodeset(nodeset_id, name);

        std::vector<int> node_ids;

        nset->get_field_data("ids_raw", node_ids);
        m_mesh.add_nodes(nodeset_index, node_ids.begin(), node_ids.end());
      }

      template <typename IossPartIterator>
      inline
      void array_mesh_reader::process_field(std::vector<double>& field, IossPartIterator part_iter, IossPartIterator part_end, const std::string &io_field_name)
      {
        assert(m_meshRankProcessed == sierra::mesh::array_mesh::Element);
  
        int offset = 0;
        while (part_iter != part_end) {
          Ioss::EntityBlock *block = *part_iter; ++part_iter;
          assert(block->field_exists(io_field_name));
          const Ioss::Field &io_field = block->get_fieldref(io_field_name);
          int component_count = io_field.transformed_storage()->component_count();

          const size_t entity_count = block->get_property("entity_count").get_int();
          offset += block->get_offset() * component_count;
          size_t data_size = entity_count * component_count;

          if (field.size() < offset+data_size) field.resize(offset+data_size);

          data_size *= sizeof(double);//turn data_size into number-of-bytes

          block->get_field_data(io_field_name, &field[offset], data_size);
        }
      }

      inline
      void array_mesh_reader::read_nodal_field(std::vector<double>& field, const std::string& name)
      {
    	  process_field(field,
    			  m_io_region->get_node_blocks().begin(),
    			  m_io_region->get_node_blocks().end(),
    			  name);
      }

      inline
      void array_mesh_reader::read_nodal_field(std::vector<double>& field, const std::string& name, double time)
      {
    	  // Get state count and all states...
    	  int step_count = m_io_region->get_property("state_count").get_int();

    	  int step = -1;
    	  for (int istep = 1; istep <= step_count; istep++) {
    		  double state_time = m_io_region->get_state_time(istep);
    		  if (std::abs(state_time - time) < 1.e-12) {
    			  step = istep;
    		  }
    	  }

    	  if (step > -1) {
    		  m_io_region->begin_state(step);

    		  process_field(field, m_io_region->get_node_blocks().begin(),
    				  m_io_region->get_node_blocks().end(), name);

    		  m_io_region->end_state(step);
    	  }
    	  else {
    		  std::cout<<"ArrayMeshReader::read_nodal_field WARNING, requested time ("<<time<<") not found in input-mesh-file, skipping this read."<<std::endl;
    	  }
      }
} // io namespace
}  // mesh namespace
} // stk namespace 

#endif

