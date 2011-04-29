/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_IO_UTIL_USECASE_MESH_HPP
#define STK_IO_UTIL_USECASE_MESH_HPP
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/fem/FEMMetaData.hpp>

// *************************************************************************
// NOTE: The functionality provided in this file is deprecated and
// NOTE: will not be supported in the future.  You should instead
// NOTE: be using the functionality provided by ../MeshReadWriteUtils
// *************************************************************************
namespace stk {
  namespace mesh {
    class MetaData;
    class Part;
  }
  namespace io {
    namespace util {
      class Gear;

      class MeshData {
	// Used to maintain state between the meta data and bulk data
	// portions of the mesh generation process for use cases.
      public:
	MeshData() : m_region(NULL),
		     m_generateSkinFaces(false)
	{}

	~MeshData();
	Ioss::Region *m_region;
	std::vector<stk::io::util::Gear*> m_gears;
	bool m_generateSkinFaces;

      private:
	MeshData(const MeshData&); // Do not implement
	MeshData& operator=(const MeshData&); // Do not implement

      };

      /** Output a help message showing the valid options for the mesh
	  read/generation including option for generated mesh and
	  gears mesh.
      */
      void show_mesh_help();


      Ioss::Region *create_output_mesh(const std::string &mesh_filename,
				       const std::string &mesh_extension,
				       const std::string &working_directory,
				       MPI_Comm comm,
				       stk::mesh::BulkData &bulk_data,
				       const Ioss::Region *in_region,
				       stk::mesh::MetaData &meta_data,
				       bool add_transient = true,
				       bool add_all_fields = false);


      void create_input_mesh(const std::string &mesh_type,
			     const std::string &mesh_filename,
			     const std::string &working_directory,
			     MPI_Comm comm,
			     stk::mesh::fem::FEMMetaData &meta_data,
			     stk::io::util::MeshData &mesh_data,
			     bool hex_only = false) ;
      void create_output_mesh(const std::string &mesh_filename,
			      const std::string &mesh_extension,
			      const std::string &working_directory,
			      MPI_Comm comm,
			      stk::mesh::BulkData &bulk_data,
			      stk::mesh::fem::FEMMetaData &meta_data,
                              MeshData &mesh_data,
			      bool add_transient = true,
			      bool add_all_fields = false);
      Ioss::Region *create_output_mesh(const std::string &mesh_filename,
				       const std::string &mesh_extension,
				       const std::string &working_directory,
				       MPI_Comm comm,
				       stk::mesh::BulkData &bulk_data,
				       const Ioss::Region *in_region,
				       stk::mesh::fem::FEMMetaData &meta_data,
				       bool add_transient = true,
				       bool add_all_fields = false) ;

      void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta);
      void process_nodeblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
      void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk);

      void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta);
      void process_elementblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
      void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk);


      void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk);
      void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta);
      void process_sidesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;

      void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk);
      void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta);
      void process_nodesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;

      void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
			  stk::mesh::EntityRank part_type,
			  Ioss::GroupingEntity *io_entity,
			  Ioss::Field::RoleType filter_role,
			  bool add_all = false);


      int process_output_request(MeshData &mesh_data,
                                 stk::mesh::BulkData &bulk,
                                 double time, bool add_all_fields = false);

      void process_output_request(Ioss::Region &region,
				  stk::mesh::BulkData &bulk,
				  int step, bool add_all_fields = false);

      void populate_bulk_data(stk::mesh::BulkData &bulk_data,
			      stk::io::util::MeshData &mesh_data,
			      const std::string &mesh_type,
                              int step=-1);
    }
  }
}
#endif
