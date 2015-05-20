#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include "test_macrodef.hpp"

#include "test_ichol_by_blocks_graphviz.hpp"
#include "test_tri_solve_by_blocks_graphviz.hpp"
#include "test_ichol_tri_solve_by_blocks_graphviz.hpp"

/// \file test_serial.hpp
/// \brief Test serial execution space
/// \author Kyungjoo Kim (kyukim@sandia.gov)

using namespace std;
using namespace Example;

int g_funct_counter = 0;

int main(int argc, char *argv[]) {
  int r_val = 0;

  string file_input = "mm_crs_input.mtx";
  int nrhs = 1;
  if (argc == 3) {
    file_input = argv[1];
    nrhs = atoi(argv[2]);
  }

  Kokkos::initialize();

  r_val += testICholByBlocksGraphviz
    <double,int,unsigned int,Kokkos::Serial,void>(file_input,
                                                  "ichol_by_blocks.gv");

  r_val += testTriSolveByBlocksGraphviz
    <double,int,unsigned int,Kokkos::Serial,void>(file_input,
                                                  1, nrhs,
                                                  "tri_solve_by_blocks.gv");

  r_val += testICholTriSolveByBlocksGraphviz
    <double,int,unsigned int,Kokkos::Serial,void>(file_input,
                                                  1, nrhs,
                                                  "ichol_tri_solve_by_blocks.gv");

  Kokkos::finalize();

  string eval;
  __EVAL_STRING__(r_val, eval);
  cout << "Testing Graphviz::" << eval << endl;

  return r_val;
}
