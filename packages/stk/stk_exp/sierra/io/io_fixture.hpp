#ifndef STK_SIERRA_IO_FIXTURE_HPP
#define STK_SIERRA_IO_FIXTURE_HPP

#include <sierra/io/modifiable_mesh_reader.hpp>
#include <sierra/mesh/csr/csr_mesh.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <init/Ionit_Initializer.h>

namespace sierra {
namespace mesh {
namespace io {

class io_fixture {
 public:
  typedef sierra::mesh::details::constant_size_field<double,3> CoordinateField;

  /** Constructor uses the mesh specified by mesh_type and file_name
     to populate the modifiable_mesh 'mmesh'.
  */
  io_fixture(stk::ParallelMachine comm,
             const std::string& mesh_type, const std::string& file_name)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, comm);
  
    Ioss::Region region(dbi, "input_model");
  
    sierra::mesh::io::ModifiableMeshReader reader(&region, mmesh);
    reader.process_mesh();
  
    // Get Ioss::Region to query fields, ...
    Ioss::NodeBlockContainer node_blocks = region.get_node_blocks();
  
    // Read in the coordinates field.
    reader.process_field(coords, node_blocks.begin(), node_blocks.end(), "mesh_model_coordinates");
  }

  ~io_fixture(){}

  modifiable_mesh mmesh;
  CoordinateField coords;
};//class io_fixture

}//namespace io
}//namespace mesh
}//namespace sierra

#endif

