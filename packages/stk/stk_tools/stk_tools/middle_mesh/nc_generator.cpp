#include <iostream>

#include "mpi.h"
#include "stk_middle_mesh/parser.hpp"
#include "stk_middle_mesh_util/stk_interface.hpp"
#include "stk_util/parallel/Parallel.hpp"

int main(int argc, char* argv[])
{
  // handle options
  if (argc != 5)
  {
    std::cerr << argc - 1 << " arguments were provided but 4 are required" << std::endl;
    std::cerr << "Usage: " << argv[0] << "input_fname sideset_name_pairs output_fname output_fname2" << std::endl;
    return 1;
  }

  stk::initialize(&argc, &argv);

  stk::middle_mesh::mesh::impl::MeshInput input;
  input.fnameIn    = argv[1];
  input.fnameOut   = argv[3];
  input.fnameOut2  = argv[4];
  input.interfaces = stk::middle_mesh::utils::impl::parse_sideset_names(argv[2]);

  stk::middle_mesh::stk_interface::impl::StkInterface stkInterface(input);

  stk::middle_mesh::impl::NonconformalOpts opts;
  opts.enableSnapAndQuality = false; // TODO
  opts.enableVolumeSnap     = true;  // TODO
  stkInterface.compute_interface(opts);
  stkInterface.write_output();

  stk::finalize();

  return 0;
}
