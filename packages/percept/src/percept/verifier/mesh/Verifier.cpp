// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include "Verifier.hpp"

#include <percept/PerceptMesh.hpp>

#include <percept/RunEnvironment.hpp>

  namespace percept
  {


    void Verifier::process_options(RunEnvironment& re)
    {
      re.clp.setDocString("verifier options");

      re.clp.setOption("mesh",       &mesh_opt,          "mesh file." );
      re.clp.setOption("printTable", &printTable_opt,    "print a table of min/max jacobian and quality measures " );
      re.clp.setOption("fullPrint",  &fullPrint_opt,     "print every element's jacobian" );


    }

    void Verifier::verify(int argc, char** argv)
    {
      //stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);

      bool printTable = false;
      bool fullJacobianPrints = false;
      std::string file_name;
      

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
        std::string in_filename = "";

        if (true)
          {
            in_filename = mesh_opt;
    
            working_directory = run_environment.directory_opt;

            file_name = in_filename;
            if (in_filename[0] != '/' && !working_directory.empty() > 0) {
              file_name = working_directory + in_filename;
            }

          } 
        else 
          {
            std::cout << "OPTION ERROR: The '--mesh <filename>' option is required.\n";
            std::exit(EXIT_FAILURE);
          }

        if (fullPrint_opt)
          {
            fullJacobianPrints = true;
          }
        if (printTable_opt)
          {
            printTable = true;
          }
      }

      PerceptMesh mesh(3, run_environment.m_comm);
      mesh.open_read_only(file_name);
      std::cout << "read in file: " << file_name << std::endl;

      TopologyVerifier topoVerifier;
      bool isBad= topoVerifier.isTopologyBad( *mesh.get_bulk_data() );
      if (isBad)
        {
          std::cout << "found bad topology" << std::endl;
        }

      GeometryVerifier geomVerifier(fullJacobianPrints, 0.0);
      isBad= geomVerifier.isGeometryBad( *mesh.get_bulk_data(), printTable );
      if (isBad)
        {
          std::cout << "found bad jacobian" << std::endl;
        }

    }

    //void MeshGeometryVerifier

      
  }//namespace percept
