// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include "MeshDifference.hpp"

#include <percept/PerceptMesh.hpp>

#include <percept/RunEnvironment.hpp>

  namespace percept
  {

    void MeshDifference::process_options(RunEnvironment& re)
    {
      re.clp.setDocString("MeshDifference options");

      re.clp.setOption("mesh1",      &mesh_opt_1,        "mesh file #1." );
      re.clp.setOption("mesh2",      &mesh_opt_2,        "mesh file #2." );
      re.clp.setOption("print_all_field_diffs",      &print_all_field_diffs,        "print all field diffs" );

    }

    void MeshDifference::run(int argc, char** argv)
    {
      std::string file_name1, file_name2;

      bool debug_re=true;
      ParallelMachineFinalize pmf(true);
      RunEnvironment run_environment(&argc, &argv, debug_re);
      process_options(run_environment);
      processCommandLine(run_environment.clp, argc, argv);

      {
        if (run_environment.help_opt) {
          printHelp(run_environment.clp);
          std::exit(EXIT_SUCCESS);
        }

        std::string working_directory = "";
        std::string in_filename1 = "";
        std::string in_filename2 = "";

        if (true)
          {
    
            working_directory = run_environment.directory_opt;

            in_filename1 = mesh_opt_1;
            file_name1 = in_filename1;
            if (in_filename1[0] != '/' && !working_directory.empty() > 0) {
              file_name1 = working_directory + in_filename1;
            }

            in_filename2 = mesh_opt_2;
            file_name2 = in_filename2;
            if (in_filename2[0] != '/' && !working_directory.empty() > 0) {
              file_name2 = working_directory + in_filename2;
            }

          } 

        if (file_name1.length() == 0 || file_name2.length() == 0)
          {
            std::cout << "OPTION ERROR: The '--mesh1 <filename> --mesh2 <filename>' option is required.\n";
            std::exit(EXIT_FAILURE);
          }

      }

      PerceptMesh mesh1(3, run_environment.m_comm);
      PerceptMesh mesh2(3, run_environment.m_comm);
      mesh1.open_read_only(file_name1);
      std::cout << "read in file1: " << file_name1 << std::endl;
      mesh2.open_read_only(file_name2);
      std::cout << "read in file2: " << file_name2 << std::endl;

      PerceptMesh::mesh_difference(mesh1, mesh2, "diff", true, print_all_field_diffs);
    }

      
  }//namespace percept


