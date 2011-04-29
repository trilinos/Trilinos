/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <assert.h>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <stk_io/IossBridge.hpp>

#include "PerceptMeshReadWrite.hpp"
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/Util.hpp>

/** \addtogroup stk_io_module
 * \{
 */

/**
 * Example code showing a basic, but complete, mesh to results output
 * coding including subsetting and periodic field input and output.
 * Includes handling of nodeblocks, element blocks, nodesets,
 * and sidesets.  Attribute fields and distribution factor
 * fields are also supported.
 *
 * This example can serve as the basis for adding binary IO support to
 * an application.  The code here uses the Ioss to/from stk::mesh
 * bridge functions in the stk::io namespace defined in IossBridge.hpp
 * include file.
 */
namespace stk {
  namespace percept {
    namespace io_util {

      // same as FEMHelper::declare_element, except the rank is chosen from the part, not defaulted to element_rank()

      stk::mesh::Entity & my_declare_element( stk::mesh::BulkData & mesh ,
                                              stk::mesh::Part & part ,
                                              const stk::mesh::EntityId elem_id ,
                                              const stk::mesh::EntityId node_id[] )
      {
        stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);
        const CellTopologyData * const top = fem_meta.get_cell_topology( part ).getCellTopologyData();

        ThrowErrorMsgIf(top == NULL,
                        "Part " << part.name() << " does not have a local topology");

        stk::mesh::PartVector empty ;
        stk::mesh::PartVector add( 1 ); add[0] = & part ;

        //const stk::mesh::EntityRank entity_rank = fem_meta.element_rank();
        const stk::mesh::EntityRank entity_rank = part.primary_entity_rank();

        stk::mesh::Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

        const stk::mesh::EntityRank node_rank = fem_meta.node_rank();

        for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
          //declare node if it doesn't already exist
          stk::mesh::Entity * node = mesh.get_entity( node_rank , node_id[i]);
          if ( NULL == node) {
            node = & mesh.declare_entity( node_rank , node_id[i], empty );
          }

          mesh.declare_relation( elem , *node , i );
        }
        return elem ;
      }


      //#define USE_LOCAL_PART_PROCESSING !defined(__IBMCPP__)
#define USE_LOCAL_PART_PROCESSING 1
#if USE_LOCAL_PART_PROCESSING
      static bool local_include_entity(Ioss::GroupingEntity *entity)
      {
        //assert(entity);

        std::string name = entity->name();
        if (name.find(PerceptMesh::s_omit_part) != std::string::npos)
          {
            //std::cout << "tmp found omitted part= " << name << std::endl;
            if ( entity->property_exists(std::string("omitted") ) )
              {
                entity->property_erase(std::string("omitted"));
              }
            entity->property_add(Ioss::Property(std::string("omitted"), 1));

            return false;
          }

        // Check whether entity has "omitted" property...
        bool omitted = (entity->property_exists("omitted")) &&
          (entity->get_property("omitted").get_int() == 1);
        //std::cout << "tmp include_entity name = " << name << " omitted= " << omitted << std::endl;

        return !omitted;
      }

      template<class ENTITY>
      void local_internal_part_processing(ENTITY *entity, stk::mesh::fem::FEMMetaData &meta,
                                                 stk::mesh::EntityRank type);

      template<>
      void local_internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::fem::FEMMetaData &meta,
                                                 stk::mesh::EntityRank type)
      {
        if (local_include_entity(entity)) {
          //std::cout << "tmp entity name= " << entity->name() << std::endl;
          stk::mesh::Part & part = meta.declare_part(entity->name(), type);
          stk::io::put_io_part_attribute(part);
        }
      }

      template<>
      void local_internal_part_processing(Ioss::EntityBlock *entity, 
                                          stk::mesh::fem::FEMMetaData &meta,
                                          stk::mesh::EntityRank type)
      {

        if (local_include_entity(entity)) {
          //std::cout << "tmp entity name= " << entity->name() << std::endl;

          //stk::mesh::Part & part = meta.declare_part(entity->name(), type);
          //stk::io::put_io_part_attribute(part);

          const Ioss::ElementTopology *topology = entity->topology();
          const CellTopologyData * const cell_topology = stk::io::map_topology_ioss_to_cell(topology);
          /// \todo IMPLEMENT Determine whether application can work
          /// with this topology type... Perhaps map_topology_ioss_to_cell only
          /// returns a valid topology if the application has registered
          /// that it can handle that specific topology.

          if (cell_topology != NULL) {

#define PERCEPT_ALLOW_PART_RANK_SUBDIM 0

#if PERCEPT_ALLOW_PART_RANK_SUBDIM
            stk::mesh::Part & part = meta.declare_part(entity->name(),  mesh::fem::CellTopology(cell_topology) );
#else
            if (0 && cell_topology->dimension != type)
              {
                std::cout << "tmp type= " << type << " dimension= " << cell_topology->dimension << std::endl;
              }
            stk::mesh::Part & part = meta.declare_part(entity->name(), type);
#endif
            stk::io::put_io_part_attribute(part);

#if PERCEPT_ALLOW_PART_RANK_SUBDIM
#else

            const CellTopologyData* cell_topo = stk::percept::PerceptMesh::get_cell_topology(part);
            std::string ctname="null";
            if (cell_topo)
              {
                ctname=cell_topo->name;
              }
            //std::cout << "tmp part = " << part.name() << " ctname= " << ctname << " cellTopo= " << mesh::fem::CellTopology(cell_topology).getName() << std::endl;
            stk::mesh::fem::set_cell_topology(meta, part,  mesh::fem::CellTopology(cell_topology) );

#endif
          } else {
            /// \todo IMPLEMENT handle cell_topolgy mapping error...
            stk::mesh::Part & part = meta.declare_part(entity->name(), type);
            stk::io::put_io_part_attribute(part);
          }
        }
      }

      template<> void local_internal_part_processing(Ioss::ElementBlock *entity, stk::mesh::fem::FEMMetaData &meta, stk::mesh::EntityRank type) { return local_internal_part_processing((Ioss::EntityBlock*)entity, meta, type);  }
      template<> void local_internal_part_processing(Ioss::SideBlock *entity, stk::mesh::fem::FEMMetaData &meta, stk::mesh::EntityRank type) { return local_internal_part_processing((Ioss::EntityBlock*)entity, meta, type);  }

      template<> void local_internal_part_processing(Ioss::SideSet *entity, stk::mesh::fem::FEMMetaData &meta, stk::mesh::EntityRank type) { return local_internal_part_processing((Ioss::GroupingEntity*)entity, meta, type);  }
      template<> void local_internal_part_processing(Ioss::NodeSet *entity, stk::mesh::fem::FEMMetaData &meta, stk::mesh::EntityRank type) { return local_internal_part_processing((Ioss::GroupingEntity*)entity, meta, type);  }



#endif

      template <typename T>
      static void local_default_part_processing(const std::vector<T*> &entities, stk::mesh::fem::FEMMetaData &meta,
                                                stk::mesh::EntityRank type)
      {
#if !USE_LOCAL_PART_PROCESSING
        stk::io::default_part_processing(entities, meta, type);
#else
        for(size_t i=0; i < entities.size(); i++) {
          T* entity = entities[i];
          local_internal_part_processing(entity, meta, type);
        }
#endif
      }

      void local_default_part_processing(const std::vector<Ioss::SideBlock*> &entities,
					 stk::mesh::fem::FEMMetaData &meta)
      {
#if !USE_LOCAL_PART_PROCESSING
        stk::io::default_part_processing(entities, meta, type);
#else
        for(size_t i=0; i < entities.size(); i++) {
	  Ioss::SideBlock* entity = entities[i];
	  int my_type = entity->topology()->parametric_dimension();
          local_internal_part_processing(entity, meta, (stk::mesh::EntityRank)my_type);
        }
#endif
      }

      void local_default_part_processing(const std::vector<Ioss::SideSet*> &entities,
					 stk::mesh::fem::FEMMetaData &meta)
      {
#if !USE_LOCAL_PART_PROCESSING
        stk::io::default_part_processing(entities, meta, type);
#else
        for(size_t i=0; i < entities.size(); i++) {
	  // A sideset can have entities with multiple parametric dimensions.
	  // Look at all sideblocks in the sideset and get the maximum parametric
	  // dimension...
	  Ioss::SideSet* entity = entities[i];
	  int my_type = entity->max_parametric_dimension();
          local_internal_part_processing(entity, meta, (stk::mesh::EntityRank)my_type);
        }
#endif
      }


      //static int s_spatial_dim = 3;

      // ========================================================================
      void process_read_nodeblocks_meta(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta, int& spatial_dim)
      {
        const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
        assert(node_blocks.size() == 1);

        Ioss::NodeBlock *nb = node_blocks[0];

        assert(nb->field_exists("mesh_model_coordinates"));
        Ioss::Field coordinates = nb->get_field("mesh_model_coordinates");
        spatial_dim = coordinates.transformed_storage()->component_count();
#if 1
        if (!meta.is_FEM_initialized())
          {
            meta.FEM_initialize(spatial_dim);
          }
#endif
        //s_spatial_dim = spatial_dim;
        //std::cout << "PerceptMeshReadWrite::process_read_nodeblocks_meta: spatial_dim= " << spatial_dim << std::endl;

        stk::mesh::Field<double,stk::mesh::Cartesian> & coord_field =
          meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

        stk::mesh::put_field( coord_field, stk::mesh::fem::FEMMetaData::NODE_RANK, meta.universal_part(),
                              spatial_dim);

        /** \todo IMPLEMENT truly handle fields... For this case we are
         * just defining a field for each transient field that is present
         * in the mesh...
         */
        stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk::mesh::fem::FEMMetaData::NODE_RANK);
      }

      // ========================================================================
      void process_read_elementblocks_meta(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta)
      {
        const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
        local_default_part_processing(elem_blocks, meta, meta.element_rank() );

        // Parts were created above, now handle element block specific
        // information (topology, attributes, ...);
        for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
            it != elem_blocks.end(); ++it) {
          Ioss::ElementBlock *entity = *it;

          if (stk::io::include_entity(entity)) {
            stk::mesh::Part* const part = meta.get_part(entity->name());
            assert(part != NULL);

            // Element Block attributes (if any)...
            /** \todo IMPLEMENT truly handle attribute fields... For this
             * case we are just defining a field for each attribute field
             * that is present in the mesh...
             */
            stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE,
                                      *part,
                                      stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );

            /** \todo IMPLEMENT truly handle fields... For this case we
             * are just defining a field for each transient field that is
             * present in the mesh...
             */
            stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                      *part,
                                      stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );

            const CellTopologyData* cell_topo = stk::percept::PerceptMesh::get_cell_topology(*part);
            std::string cell_topo_name = "UNKNOWN";
            if (cell_topo != NULL)
              cell_topo_name = cell_topo->name;

#if 0
            std::cout << entity->type_string() << ": " << entity->name()
                      << " , celltop = " << cell_topo_name
                      << std::endl ;
#endif
          }
        }
      }

      // ========================================================================
      void process_read_nodesets_meta(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta)
      {
        const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
        local_default_part_processing(node_sets, meta, stk::mesh::fem::FEMMetaData::NODE_RANK);

        //!std::cout << "PerceptMeshReadWrite::process_read_nodesets: node_sets size = " << node_sets.size() <<  std::endl;

        /** \todo REFACTOR should "distribution_factor" be a default field
         * that is automatically declared on all objects that it exists
         * on as is done in current framework?
         */
        stk::mesh::Field<double> & distribution_factors_field =
          meta.declare_field<stk::mesh::Field<double> >("distribution_factors");

        /** \todo REFACTOR How to associate distribution_factors field
         * with the nodeset part if a node is a member of multiple
         * nodesets
         */

        for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
            it != node_sets.end(); ++it) {
          Ioss::NodeSet *entity = *it;

          //std::cout << "PerceptMeshReadWrite::process_read_nodesets: node_sets entity= " << entity->name()
          //          << " stk::io::include_entity(entity) = " << stk::io::include_entity(entity) << std::endl;


          if (stk::io::include_entity(entity)) {
            stk::mesh::Part* const part = meta.get_part(entity->name());
            assert(part != NULL);
            assert(entity->field_exists("distribution_factors"));

            stk::mesh::put_field(distribution_factors_field, stk::mesh::fem::FEMMetaData::NODE_RANK, *part);

            /** \todo IMPLEMENT truly handle fields... For this case we
             * are just defining a field for each transient field that is
             * present in the mesh...
             */
            stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                      *part,
                                      stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );
          }
        }
      }

      // ========================================================================
      void process_read_surface_entity_meta(Ioss::SideSet *entity, stk::mesh::fem::FEMMetaData &meta)
      {
        assert(entity->type() == Ioss::SIDESET);
	const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
	local_default_part_processing(blocks, meta);

        stk::mesh::Part* const fs_part = meta.get_part(entity->name());
        assert(fs_part != NULL);

        stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
        bool surface_df_defined = false; // Has the surface df field been defined yet?


        int block_count = entity->block_count();
        for (int i=0; i < block_count; i++) {
          Ioss::SideBlock *fb = entity->get_block(i);
          if (stk::io::include_entity(fb)) {
            //std::cout << "PerceptMeshReadWrite::process_read_surface_entity fb->type_string()= " <<  fb->type_string() << " " << fb->name() << "\n";
            stk::mesh::Part * const fb_part = meta.get_part(fb->name());
            assert(fb_part != NULL);
            meta.declare_part_subset(*fs_part, *fb_part);

            if (fb->field_exists("distribution_factors")) {
              if (!surface_df_defined) {
                std::string field_name = entity->name() + "_df";
                distribution_factors_field =
                  &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(field_name);
                stk::io::set_distribution_factor_field(*fs_part, *distribution_factors_field);
                surface_df_defined = true;
              }
              stk::io::set_distribution_factor_field(*fb_part, *distribution_factors_field);
              int face_node_count = fb->topology()->number_nodes();
              stk::mesh::put_field(*distribution_factors_field,
                                   stk::percept::PerceptMesh::fem_entity_rank( fb_part->primary_entity_rank() ),
                                   *fb_part, face_node_count);
            }

            /** \todo IMPLEMENT truly handle fields... For this case we
             * are just defining a field for each transient field that is
             * present in the mesh...
             */
            stk::io::define_io_fields(fb, Ioss::Field::TRANSIENT,
                                      *fb_part,
                                      stk::percept::PerceptMesh::fem_entity_rank( fb_part->primary_entity_rank() ) );
          }
        }
      }

      // ========================================================================
      void process_read_sidesets_meta(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta)
      {
        const Ioss::SideSetContainer& side_sets = region.get_sidesets();
        local_default_part_processing(side_sets, meta);

        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
          Ioss::SideSet *entity = *it;

          if (stk::io::include_entity(entity)) {
            process_read_surface_entity_meta(entity, meta);
          }
        }
      }

      // ========================================================================
      // Bulk Data
      // ========================================================================

      void add_other_fields(std::vector<stk::mesh::Entity*> nodes, Ioss::NodeBlock *nb, stk::mesh::BulkData &bulk)
      {
        EXCEPTWATCH;
        //std::cout << "PerceptMeshReadWrite::add_other_fields: " << std::endl;
        bool printInfo = true;

        const mesh::fem::FEMMetaData& metaData = mesh::fem::FEMMetaData::get(bulk);

        const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

        unsigned nparts = parts.size();
        if (printInfo) std::cout << "info>    Number of parts = " << nparts << std::endl;

        const mesh::FieldVector & fields =  metaData.get_fields();
        unsigned nfields = fields.size();
        if (printInfo)
          {
            std::cout << "info>    Number of fields = " << fields.size() << std::endl;
            for (unsigned ifld = 0; ifld < nfields; ifld++)
              {
                mesh::FieldBase *field = fields[ifld];
                if (printInfo) std::cout << "info>    Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
                if (printInfo) std::cout << "info>    " << *field << std::endl;

                //if (field->name() == std::string("pressure"))
                //  stk::io::field_data_from_ioss(field, nodes, nb, field->name());
                if (field->name() == std::string("coordinates"))
                  stk::io::field_data_from_ioss(field, nodes, nb, "mesh_model_coordinates");
              }
          }
      }

      void process_read_nodeblocks_bulk(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        // This must be called after the "process_read_element_blocks" call
        // since there may be nodes that exist in the database that are
        // not part of the analysis mesh due to subsetting of the element
        // blocks.

        const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
        assert(node_blocks.size() == 1);

        Ioss::NodeBlock *nb = node_blocks[0];

        std::vector<stk::mesh::Entity*> nodes;
        stk::io::get_entity_list(nb, stk::mesh::fem::FEMMetaData::NODE_RANK, bulk, nodes);

        /** \todo REFACTOR Application would probably store this field
         * (and others) somewhere after the declaration instead of
         * looking it up each time it is needed.
         */
        const stk::mesh::fem::FEMMetaData& meta = mesh::fem::FEMMetaData::get(bulk);
        stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
          meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

        stk::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");

        //add_other_fields(nodes, nb, bulk);

      }

      // ========================================================================
      void process_read_elementblocks_bulk(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();

        for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
            it != elem_blocks.end(); ++it) {
          Ioss::ElementBlock *entity = *it;

          if (stk::io::include_entity(entity)) {
            const std::string &name = entity->name();
            const stk::mesh::fem::FEMMetaData& meta = mesh::fem::FEMMetaData::get(bulk);
            stk::mesh::Part* const part = meta.get_part(name);
            assert(part != NULL);

            const CellTopologyData* cell_topo = stk::percept::PerceptMesh::get_cell_topology(*part);
            if (cell_topo == NULL) {
              std::ostringstream msg ;
              msg << " PerceptMeshReadWrite::process_read_elementblocks_bulk: INTERNAL_ERROR: Part " << part->name() << " returned NULL from stk::percept::PerceptMesh::get_cell_topology()";
              throw std::runtime_error( msg.str() );
            }

            std::vector<int> elem_ids ;
            std::vector<int> connectivity ;

            entity->get_field_data("ids", elem_ids);
            entity->get_field_data("connectivity", connectivity);

            size_t element_count = elem_ids.size();
            int nodes_per_elem = cell_topo->node_count ;
            std::vector<mesh::EntityId> e_connectivity(nodes_per_elem) ;

            std::vector<stk::mesh::Entity*> elements(element_count);
            for(size_t i=0; i<element_count; ++i) {
              int *conn = &connectivity[i*nodes_per_elem];
              for (int j=0; j < nodes_per_elem; j++)
                {
                  e_connectivity[j] = (mesh::EntityId)conn[j];
                }
              mesh::EntityId* e_conn = &(e_connectivity[0]);
              elements[i] = &stk::mesh::fem::declare_element(bulk, *part, elem_ids[i], e_conn);
#if PERCEPT_ALLOW_PART_RANK_SUBDIM
              //elements[i] = &my_declare_element(bulk, *part, elem_ids[i], e_conn);
#endif
            }

            // For this example, we are just taking all attribute fields
            // found on the io database and populating fields on the
            // corresponding mesh part.  In practice, would probably be
            // selective about which attributes to use...
            Ioss::NameList names;
            entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
            for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
              if (*I == "attribute" && names.size() > 1)
                continue;
              stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase>(*I);
              stk::io::field_data_from_ioss(field, elements, entity, *I);

            }
          }
        }
      }

      // ========================================================================
      void process_read_nodesets_bulk(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        // Should only process nodes that have already been defined via the element
        // blocks connectivity lists.
        const Ioss::NodeSetContainer& node_sets = region.get_nodesets();

        for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
            it != node_sets.end(); ++it) {
          Ioss::NodeSet *entity = *it;

          if (stk::io::include_entity(entity)) {
            const std::string & name = entity->name();
            const stk::mesh::fem::FEMMetaData& meta = mesh::fem::FEMMetaData::get(bulk);
            stk::mesh::Part* const part = meta.get_part(name);
            assert(part != NULL);
            stk::mesh::PartVector add_parts( 1 , part );

            std::vector<int> node_ids ;
            int node_count = entity->get_field_data("ids", node_ids);

            std::vector<stk::mesh::Entity*> nodes(node_count);
            for(int i=0; i<node_count; ++i) {
              nodes[i] = bulk.get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i] );
              if (nodes[i] != NULL)
                bulk.declare_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i], add_parts );
            }


            /** \todo REFACTOR Application would probably store this field
             * (and others) somewhere after the declaration instead of
             * looking it up each time it is needed.
             */
            stk::mesh::Field<double> *df_field =
              meta.get_field<stk::mesh::Field<double> >("distribution_factors");

            if (df_field != NULL) {
              stk::io::field_data_from_ioss(df_field, nodes, entity, "distribution_factors");
            }
          }
        }
      }

      // ========================================================================
      void process_read_surface_entity_bulk(const Ioss::SideSet* io ,
                                            stk::mesh::BulkData & bulk)
      {
        assert(io->type() == Ioss::SIDESET);
        const stk::mesh::fem::FEMMetaData& fem_meta = mesh::fem::FEMMetaData::get(bulk);
        //const stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get ( meta );

        int block_count = io->block_count();
        for (int i=0; i < block_count; i++) {
          Ioss::SideBlock *block = io->get_block(i);
          if (stk::io::include_entity(block)) {
            std::vector<int> side_ids ;
            std::vector<int> elem_side ;

            stk::mesh::Part * const fb_part = fem_meta.get_part(block->name());

            block->get_field_data("ids", side_ids);
            block->get_field_data("element_side", elem_side);

            assert(side_ids.size() * 2 == elem_side.size());
            stk::mesh::PartVector add_parts( 1 , fb_part );

            size_t side_count = side_ids.size();
            std::vector<stk::mesh::Entity*> sides(side_count);
            for(size_t is=0; is<side_count; ++is) {


#if PERCEPT_ALLOW_PART_RANK_SUBDIM
              if (0 && fem_meta.element_rank() != fb_part->primary_entity_rank())
                {
                  std::cout << "tmp fem_meta.element_rank()= " << fem_meta.element_rank() <<
                    " fb_part->primary_entity_rank() = " <<  fb_part->primary_entity_rank() << std::endl;
                }
              stk::mesh::Entity* const elem = bulk.get_entity(fb_part->primary_entity_rank(), elem_side[is*2]);
#else
              stk::mesh::Entity* const elem = bulk.get_entity(fem_meta.element_rank(), elem_side[is*2]);
#endif
              // If NULL, then the element was probably assigned to an
              // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
              // Hence the '-1' in the following line.
              int side_ordinal = elem_side[is*2+1] - 1 ;

              // element block that appears in the database, but was
              // subsetted out of the analysis mesh. Only process if
              // non-null.
              if (elem != NULL) {
                stk::mesh::Entity& side =
                  stk::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);

                bulk.change_entity_parts( side, add_parts );
                sides[is] = &side;
              } else {
                sides[is] = NULL;
              }
            }

            const stk::mesh::Field<double, stk::mesh::ElementNode> *df_field =
              stk::io::get_distribution_factor_field(*fb_part);
            if (df_field != NULL) {
              stk::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
            }
          }
        }
      }

      // ========================================================================
      void process_read_sidesets_bulk(Ioss::Region &region, stk::mesh::BulkData &bulk)
      {
        const Ioss::SideSetContainer& side_sets = region.get_sidesets();

        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
          Ioss::SideSet *entity = *it;

          if (stk::io::include_entity(entity)) {
            process_read_surface_entity_bulk(entity, bulk);
          }
        }
      }

      // ========================================================================
      // ========================================================================
      void get_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                          stk::mesh::EntityRank part_type,
                          Ioss::GroupingEntity *io_entity,
                          Ioss::Field::RoleType filter_role)
      {
        //std::cout << "PerceptMeshReadWrite::get_field_data: ... " << std::endl;
        std::vector<stk::mesh::Entity*> entities;
        stk::io::get_entity_list(io_entity, part_type, bulk, entities);

        stk::mesh::fem::FEMMetaData & meta = mesh::fem::FEMMetaData::get(part);
        stk::mesh::Part &universal = meta.universal_part();
        const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

        std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const stk::mesh::FieldBase *f = *I; ++I;
          //std::cout << "PerceptMeshReadWrite::get_field_data: f->name()= " << f->name() << " part.name()= " << part.name() << std::endl;

          if (stk::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
            stk::io::field_data_from_ioss(f, entities, io_entity, f->name());
          }
        }
      }

      void process_read_input_request(Ioss::Region &region,
                                      stk::mesh::BulkData &bulk,
                                      int step)
      {
        //std::cout << "PerceptMeshReadWrite::process_read_input_request begin_state("<<step<<")" << std::endl;
        if (step> 0)
          {
            region.begin_state(step);
          }

        // Special processing for nodeblock (all nodes in model)...
        const stk::mesh::fem::FEMMetaData & meta = mesh::fem::FEMMetaData::get(bulk);

        // ??? Get field data from nodeblock...
        get_field_data(bulk, meta.universal_part(), stk::mesh::fem::FEMMetaData::NODE_RANK,
                       region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

        const stk::mesh::PartVector & all_parts = meta.get_parts();
        for ( stk::mesh::PartVector::const_iterator
                ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

          stk::mesh::Part * const part = *ip;

          // Check whether this part should be output to results database.
          if (stk::io::is_part_io_part(*part)) {
            // Get Ioss::GroupingEntity corresponding to this part...
            Ioss::GroupingEntity *entity = region.get_entity(part->name());
            if (entity != NULL) {
              if (entity->type() == Ioss::SIDESET) {
		Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
                int block_count = sset->block_count();
                for (int i=0; i < block_count; i++) {
                  Ioss::SideBlock *fb = sset->get_block(i);
                  /// \todo REFACTOR Need filtering mechanism.
                  get_field_data(bulk, *part,
                                 stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                                 fb, Ioss::Field::TRANSIENT);
                }
              } else {
                get_field_data(bulk, *part,
                               stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                               entity, Ioss::Field::TRANSIENT);
              }
            } else {
              /// \todo IMPLEMENT handle error... Possibly an assert since
              /// I think the corresponding entity should always exist...
            }
          }
        }

        //std::cout << "PerceptMeshReadWrite::process_read_input_request end_state("<<step<<")" << std::endl;
        if (step > 0)
          {
            region.end_state(step);
          }
      }

      void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                          stk::mesh::EntityRank part_type,
                          Ioss::GroupingEntity *io_entity,
                          Ioss::Field::RoleType filter_role)
      {
        std::vector<stk::mesh::Entity*> entities;
        stk::io::get_entity_list(io_entity, part_type, bulk, entities);

        stk::mesh::fem::FEMMetaData & meta = mesh::fem::FEMMetaData::get(part);
        stk::mesh::Part &universal = meta.universal_part();
        const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

        std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const stk::mesh::FieldBase *f = *I; ++I;
          //std::cout << "PerceptMeshReadWrite::put_field_data: f->name()= " << f->name() <<  " part.name()= " << part.name() << std::endl;
          if (stk::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
            stk::io::field_data_to_ioss(f, entities, io_entity, f->name(), filter_role);
          }
        }
      }

      void process_output_request(Ioss::Region &region,
                                  stk::mesh::BulkData &bulk,
                                  int step)
      {
        //std::cout << "PerceptMeshReadWrite::process_output_request begin_state("<<step<<")" << std::endl;
        if (step > 0)
          {
            region.begin_state(step);
          }
        // Special processing for nodeblock (all nodes in model)...
        const stk::mesh::fem::FEMMetaData & meta = mesh::fem::FEMMetaData::get(bulk);

        put_field_data(bulk, meta.universal_part(), stk::mesh::fem::FEMMetaData::NODE_RANK,
                       region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

        const stk::mesh::PartVector & all_parts = meta.get_parts();
        for ( stk::mesh::PartVector::const_iterator
                ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

          stk::mesh::Part * const part = *ip;

          //           std::cout << "process_output_request: part = " << part->name() << "is_part_io_part = "
          //                     << stk::io::is_part_io_part(*part) << std::endl;
          //           Util::pause(true);


          // Check whether this part should be output to results database.
          if (stk::io::is_part_io_part(*part)) {

            // Get Ioss::GroupingEntity corresponding to this part...
            Ioss::GroupingEntity *entity = region.get_entity(part->name());
            if (entity != NULL) {

              if (entity->type() == Ioss::SIDESET) {
		Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
                int block_count = sset->block_count();

                for (int i=0; i < block_count; i++) {
                  Ioss::SideBlock *fb = sset->get_block(i);
                  /// \todo REFACTOR Need filtering mechanism.
                  put_field_data(bulk, *part,
                                 stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                                 fb, Ioss::Field::TRANSIENT);
                }
              }
              else
                {
                  put_field_data(bulk, *part,
                                 stk::percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                                 entity, Ioss::Field::TRANSIENT);
                }
            }
            else
              {
                /// \todo IMPLEMENT handle error... Possibly an assert since
                /// I think the corresponding entity should always exist...
              }
          }
        }
        //std::cout << "PerceptMeshReadWrite::process_output_request end_state("<<step<<")" << std::endl;
        if (step > 0)
          {
            region.end_state(step);
          }
      }

    }
  }
}
