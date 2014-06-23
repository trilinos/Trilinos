#ifndef STK_SIERRA_ARRAY_MESH_FIXTURE_HPP
#define STK_SIERRA_ARRAY_MESH_FIXTURE_HPP

#include <sierra/io/array_mesh_reader.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <init/Ionit_Initializer.h>

namespace sierra {
namespace mesh {
namespace io {

class array_mesh_fixture {
 public:
  /** Constructor uses the mesh specified by mesh_type and file_name
     to populate the modifiable_mesh 'mesh'.
  */
  array_mesh_fixture(stk::ParallelMachine comm,
             const std::string& mesh_type, const std::string& file_name)
   : m_mesh(true),
     m_coords()
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, comm);
  
    Ioss::Region *region = new Ioss::Region(dbi, "input_model");
  
    sierra::mesh::io::array_mesh_reader reader(region, m_mesh);
    reader.process_mesh();
  
    // Get Ioss::Region to query fields, ...
    Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
  
    // Read in the coordinates field.
    reader.process_field(m_coords, node_blocks.begin(), node_blocks.end(), "mesh_model_coordinates");

    delete region;
  }

  ~array_mesh_fixture(){}

  array_mesh m_mesh;
  std::vector<double> m_coords;
};//class array_mesh_fixture

}//namespace io
}//namespace mesh
}//namespace sierra

#endif

