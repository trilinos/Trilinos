
#include "Verifier.hpp"

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/RunEnvironment.hpp>

namespace stk
{
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
      RunEnvironment run_environment(&argc, &argv, debug_re);
      process_options(run_environment);
      run_environment.processCommandLine(&argc, &argv);

      {
        if (run_environment.help_opt) {
          run_environment.printHelp();
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

      RunEnvironment::doLoadBalance(run_environment.m_comm, file_name);

      PerceptMesh mesh(3, run_environment.m_comm);
      mesh.openReadOnly(file_name);
      std::cout << "read in file: " << file_name << std::endl;

      TopologyVerifier topoVerifier;
      bool isBad= topoVerifier.isTopologyBad( *mesh.get_bulkData() );
      if (isBad)
        {
          std::cout << "found bad topology" << std::endl;
        }

      GeometryVerifier geomVerifier(fullJacobianPrints);
      isBad= geomVerifier.isGeometryBad( *mesh.get_bulkData(), printTable );
      if (isBad)
        {
          std::cout << "found bad jacobian" << std::endl;
        }

    }

    //void MeshGeometryVerifier

      
  }//namespace percept
}//namespace stk


