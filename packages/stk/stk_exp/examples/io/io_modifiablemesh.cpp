#ifndef __IBMCPP__
#include <Ioss_SubSystem.h>

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/io/modifiable_mesh_reader.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>
#include <sierra/io/io_fixture.hpp>

#include <vector>
#include <string>

#include <init/Ionit_Initializer.h>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

void get_cmd_arg(int argc, char** argv, const std::string& flag, std::string& arg_val)
{
  arg_val.clear();

  for(int i=0; i<argc; ++i) {
    if (argv[i] == NULL) continue;
    if (flag != argv[i]) continue;
    if (i < argc-1) {
      arg_val = argv[i+1];
    }
  }
}

int main(int argc, char **argv)
{
  std::string mesh_type, file_name;
  get_cmd_arg(argc, argv, "-mesh_type", mesh_type);
  get_cmd_arg(argc, argv, "-mesh_file", file_name);

  if (mesh_type.empty() || file_name.empty()) {
    std::cout << "-mesh_type or -mesh_file not specified, using generated mesh." << std::endl;
    mesh_type = "generated";
    file_name = "10x20x30";
  }

  Ioss::Init::Initializer init_db;
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL,
						  use_case_environment.m_comm);

  Ioss::Region *region = new Ioss::Region(dbi, "input_model");

  sierra::mesh::modifiable_mesh mesh;

  sierra::mesh::io::ModifiableMeshReader reader(region, mesh);
  reader.process_mesh();

  // Declare coordinate field
  sierra::mesh::details::constant_size_field<double, 3> coordinates;

  // Get Ioss::Region to query fields, ...
  Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();

  // Read in the coordinates field.
  reader.process_field(coordinates, node_blocks.begin(), node_blocks.end(), "mesh_model_coordinates");

  sierra::mesh::modifiable_mesh::entity_descriptor_range entity_range = mesh.get_entities();
  size_t total_num_entities = boost::distance(entity_range);

  std::cout << "Total number of entities: " << total_num_entities << std::endl;

  delete region;

  {
    sierra::mesh::io::io_fixture fixture(use_case_environment.m_comm, mesh_type, file_name);

    sierra::mesh::modifiable_mesh::entity_descriptor_range entity_range2 = fixture.mmesh.get_entities();
    size_t total_num_entities2 = boost::distance(entity_range2);

    std::cout << "Total number of entities in fixture mesh: " << total_num_entities2 << std::endl;
  }

  return 0;
}
#else
int main(int argc, char **argv)
{
  return 0;
}
#endif

