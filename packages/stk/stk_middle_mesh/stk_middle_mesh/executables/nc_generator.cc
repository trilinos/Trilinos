#include <iostream>

#include "mpi.h"
#include "parser.h"
#include "stk/stk_interface.h"

int main(int argc, char* argv[])
{
  // handle options
  if (argc != 5)
  {
    std::cerr << argc - 1 << " arguments were provided but 4 are required" << std::endl;
    std::cerr << "Usage: " << argv[0] << "input_fname sideset_name_pairs output_fname output_fname2" << std::endl;
    return 1;
  }

  MPI_Init(&argc, &argv);

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

  MPI_Finalize();

  return 0;
}