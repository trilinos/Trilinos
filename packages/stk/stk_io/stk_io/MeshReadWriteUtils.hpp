/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_IO_MESHREADWRITEUTILS_HPP
#define STK_IO_MESHREADWRITEUTILS_HPP
#include <stk_mesh/base/BulkData.hpp>
#include <Ioss_SubSystem.h>

namespace stk {
  namespace mesh {
    class Part;
    namespace fem {
      class FEMMetaData;
    }
  }
  namespace io {
    class MeshData {
      // Used to maintain state between the meta data and bulk data
      // portions of the mesh generation process for use cases.
    public:
      MeshData() : m_input_region(NULL), m_output_region(NULL)
      {}

      ~MeshData();
      Ioss::Region *m_input_region;
      Ioss::Region *m_output_region;

    private:
      MeshData(const MeshData&); // Do not implement
      MeshData& operator=(const MeshData&); // Do not implement

    };

    /** Output a help message showing the valid options for the mesh
	read/generation including option for generated mesh.
    */
    void show_mesh_help();

    void create_input_mesh(const std::string &mesh_type,
			   const std::string &mesh_filename,
			   MPI_Comm comm,
			   stk::mesh::fem::FEMMetaData &fem_meta,
			   MeshData &mesh_data);

    void populate_bulk_data(stk::mesh::BulkData &bulk_data, stk::io::MeshData &mesh_data);

    void define_input_fields(MeshData &mesh_data, stk::mesh::fem::FEMMetaData &fem_meta);
    
    int process_input_request(MeshData &mesh_data,
			      stk::mesh::BulkData &bulk,
			      double time);

    // ???? STORE Ioss::Region or mesh_data as attribute on FEMMeta?		   stk::io::MeshData &mesh_data);
    // ???? IF do that, then need function for getting it back...

    void create_output_mesh(const std::string &mesh_filename,
			    MPI_Comm comm,
			    stk::mesh::BulkData &bulk_data,
			    MeshData &mesh_data);

    void define_output_fields(MeshData &mesh_data,
			      stk::mesh::fem::FEMMetaData &fem_meta,
			      bool add_all_fields = false);
    
    int process_output_request(MeshData &mesh_data,
			       stk::mesh::BulkData &bulk,
			       double time);

    void process_nodeblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
    void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk);

    void process_elementblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
    void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk);

    void process_edgesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
    void process_edgesets(Ioss::Region &region, stk::mesh::BulkData &bulk);

    void process_facesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
    void process_facesets(Ioss::Region &region, stk::mesh::BulkData &bulk);

    void process_nodesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta) ;
    void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk);

    void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
			stk::mesh::EntityRank part_type,
			Ioss::GroupingEntity *io_entity,
			Ioss::Field::RoleType filter_role);

  }
}
#endif
