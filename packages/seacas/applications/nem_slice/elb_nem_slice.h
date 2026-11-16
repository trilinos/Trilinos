/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "elb.h" // for LB_Description<INT>, get_time, etc

template <typename INT> class NemSlice
{
public:
  NemSlice(int argc, char *argv[]) : time1(0), time2(0), argc(argc), argv(argv) {}

  int run()
  {
    init_mpi();
    int result = setup();
    if (result != 0)
      return result;

    result = find_graph();
    if (result != 0)
      return result;

    result = find_load_balance();
    if (result != 0)
      return result;

    result = write_partition();
    if (result != 0)
      return result;
    finalize_mpi();

    return 0;
  }

  void init_mpi();
  void finalize_mpi();
  int  setup();
  int  find_graph();
  int  find_load_balance();
  int  write_partition();

  Machine_Description    machine;
  LB_Description<INT>    lb;
  Problem_Description    problem;
  Solver_Description     solver;
  Weight_Description     weight;
  Mesh_Description<INT>  mesh;
  Sphere_Info            sphere;
  Graph_Description<INT> graph;

private:
  std::string exoII_inp_file;
  std::string ascii_inp_file;
  std::string nemI_out_file;

  double time1;
  double time2;

  int    argc;
  char **argv;

  void print_input(Machine_Description * /*machine*/, LB_Description<INT> * /*lb*/,
                   Problem_Description * /*prob*/, Solver_Description * /*solver*/,
                   Weight_Description * /*weight*/);
};
