// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include <mpi.h>                          // for MPI_Comm_rank, MPI_Comm_size
#include <stddef.h>                       // for size_t
#include <makeparfiles.H>                 // for MakeParFile
#include <string>                         // for string
#include "stk_io/FillMesh.hpp"            // for fill_mesh
#include "stk_mesh/base/BulkData.hpp"     // for BulkData
#include "stk_mesh/base/GetEntities.hpp"  // for count_selected_entities
#include "stk_mesh/base/MetaData.hpp"     // for MetaData
#include "stk_topology/topology.hpp"      // for topology, topology::rank_t:...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace
{

TEST(PMESH, a)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    int num_procs = -1, my_proc_id = -1;
    MPI_Comm_rank(comm, &my_proc_id);
    MPI_Comm_size(comm, &num_procs);

    // ncut_x * ncuts_y * ncuts_z is number of processor
    int ncuts_x = 1, ncuts_y = 1, ncuts_z = 1;
    size_t num_elem_per_edge = 10;
    double lenx = 1.0, leny = 1.0, lenzi = 1.0;
    int numraid = 1;
    const std::string rootdir="";
    const std::string subdir="";

    stk_tools::MakeParFile(my_proc_id, num_procs, ncuts_x, ncuts_y, ncuts_z, num_elem_per_edge, lenx, leny, lenzi, numraid, rootdir.c_str(), subdir.c_str());

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);
    stk::io::fill_mesh("1/cube.par.1.0", bulk);

    size_t numElements = stk::mesh::count_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(num_elem_per_edge*num_elem_per_edge*num_elem_per_edge, numElements);
}

}
