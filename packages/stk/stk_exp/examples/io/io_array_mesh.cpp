#include <Ioss_SubSystem.h>

#include <sierra/mesh/array_mesh/array_mesh.hpp>
#include <sierra/io/array_mesh_reader.hpp>
#include <sierra/io/array_mesh_writer.hpp>

#include <mpi.h>
#include <vector>
#include <iterator>
#include <string>

#include <init/Ionit_Initializer.h>

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

void print_mesh_info(sierra::mesh::array_mesh& mesh, const std::vector<double>& coordinates);

int main(int argc, char **argv)
{
  std::string mesh_type, file_name;
  get_cmd_arg(argc, argv, "-mesh_type", mesh_type);
  get_cmd_arg(argc, argv, "-mesh_file", file_name);

  if (mesh_type.empty() || file_name.empty()) {
    std::cout << "-mesh_type or -mesh_file not specified, using generated mesh." << std::endl;
    mesh_type = "generated";
    file_name = "2x1x1|sideset:xXyY";
  }

  MPI_Init(&argc, &argv);

  std::cout << "mesh spec: "<<file_name<<std::endl;

  sierra::mesh::array_mesh mesh;
  sierra::mesh::io::array_mesh_reader reader(MPI_COMM_WORLD, mesh_type, file_name, mesh);

  // Get coordinate field
  std::vector<double> coordinates;
  reader.read_nodal_field(coordinates, "mesh_model_coordinates");

  print_mesh_info(mesh, coordinates);

  std::vector<double> temperature;
  reader.read_nodal_field(temperature, "TND", 0.0);

  std::string out_file_name = "out_"+file_name;
  if (mesh_type == "generated") {
	  out_file_name = "out_gmesh.exo";
  }

  sierra::mesh::io::array_mesh_writer writer(MPI_COMM_WORLD, out_file_name, mesh, coordinates);

  writer.write_nodal_field("TND", temperature, 0.0);


  MPI_Finalize();

  return 0;
}

void print_mesh_info(sierra::mesh::array_mesh& mesh, const std::vector<double>& coordinates)
{
 std::cout << "Number of elements: " << mesh.get_num_elements() << std::endl;
 std::cout << "Number of nodes: " << mesh.get_num_nodes() << std::endl;
 sierra::mesh::array_mesh::ConstIntRange node_ids = mesh.get_node_ids();
 const int* node_ids_ptr = &*node_ids.first;
 sierra::mesh::array_mesh::SidesetRange sidesets = mesh.get_sidesets();
 std::cout << "Number of sidesets: " << std::distance(sidesets.first,sidesets.second) << std::endl;
 for(sierra::mesh::array_mesh::SidesetIterator sideset_iterator=sidesets.first, sideset_end=sidesets.second; sideset_iterator!=sideset_end; ++sideset_iterator) {
   std::cout<<"  Number of sides in sideset '"<<mesh.get_name(*sideset_iterator)
   		<<"' id="<<mesh.get_id(*sideset_iterator)<<": "<<mesh.get_sideset_elem_nums(*sideset_iterator).size()<<std::endl;
 }
 std::cout << "Element-node Connectivity:"<<std::endl;
 sierra::mesh::array_mesh::BlockRange blocks = mesh.get_blocks();
 for(sierra::mesh::array_mesh::BlockIterator blk_iter=blocks.first, blk_end=blocks.second; blk_iter != blk_end; ++blk_iter) {
   if (mesh.get_topology(*blk_iter) == stk::topology::NODE) continue;
   size_t num_elems = mesh.get_num_elems(*blk_iter);
   size_t nodes_per_elem = mesh.get_num_nodes_per_elem(*blk_iter);
   std::cout << "Block '"<<mesh.get_name(*blk_iter)<<"', id="<<mesh.get_id(*blk_iter)<<": "<<std::endl;
   const std::vector<int>& conn = mesh.get_block_connectivity(*blk_iter);
   size_t offset = 0;
   for(size_t elem=0; elem<num_elems; ++elem) {
     for(size_t node=0; node<nodes_per_elem; ++node) {
       int node_num = conn[offset++];
        std::cout<<node_ids_ptr[node_num]<<" ";
     }
     offset -= nodes_per_elem;
     std::cout<<std::endl;
     for(size_t node=0; node<nodes_per_elem; ++node) {
       int node_num = conn[offset++];
       int idx = node_num*3;
       std::cout<<"("<<coordinates[idx]<<","<<coordinates[idx+1]<<","<<coordinates[idx+2]<<") ";
     }
     std::cout<<std::endl;
   }
 }
}
