
#include <Ioss_SubSystem.h>
#include <samba_io/mesh_reader.hpp>
#include <samba_io/sidesets_manager.hpp>
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
  std::string mesh_type, file_name; // Should add dimension to this list.
  get_cmd_arg(argc, argv, "-mesh_type", mesh_type);
  get_cmd_arg(argc, argv, "-mesh_file", file_name);

  if (mesh_type.empty() || file_name.empty()) {
    std::cout << "-mesh_type or -mesh_file not specified, using generated mesh." << std::endl;
    mesh_type = "generated";
    file_name = "1x2x3|sideset:xXyYzZ";
  }

  Ioss::Init::Initializer init_db;
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL,
                                                  use_case_environment.m_comm);
  Ioss::Region *region = new Ioss::Region(dbi, "input_model");

  // Create the mesh to be populated.
  samba::mesh mesh(samba::spatial_dimension::create(3));

  // Read in the entities and connectivity.
  samba::io::mesh_reader reader(mesh);

  // Right now, samba seems buggy in handling fields constructed after a mesh is populated,
  // we constuct these before calling reader.process_mesh(.).
  samba::io::mesh_reader::desc_idx_coordinates_field_type
    node_coordinates_by_eDescriptor(mesh, samba::entity_rank::node());
  samba::io::mesh_reader::key_idx_coordinates_field_type
    node_coordinates_by_eKey(mesh, samba::entity_topology::node());

  // Read in the entities and connectivity.
  reader.process_mesh(region);

  // Get Ioss::Region to query fields, ...
  Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();

  // Read in the coordinates field to descriptor-indexed version.
  reader.process_ent_desc_field(node_coordinates_by_eDescriptor, node_blocks.begin(), node_blocks.end(),
                                "mesh_model_coordinates");

  // Read in the coordinates field to key-indexed version.
  reader.process_ent_key_field(node_coordinates_by_eKey, node_blocks.begin(), node_blocks.end(),
                                "mesh_model_coordinates");


  // Read in sidesets.
  samba::io::sidesets_manager sidesets_mgr(mesh);
  sidesets_mgr.process_region(region, reader);


  delete region;

  return 0;
}
