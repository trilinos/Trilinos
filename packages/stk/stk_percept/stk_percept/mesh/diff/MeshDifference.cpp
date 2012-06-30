
#include "MeshDifference.hpp"

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/RunEnvironment.hpp>

namespace stk
{
  namespace percept
  {

    void MeshDifference::process_options(RunEnvironment& re)
    {
      re.clp.setDocString("MeshDifference options");

      re.clp.setOption("mesh1",      &mesh_opt_1,        "mesh file #1." );
      re.clp.setOption("mesh2",      &mesh_opt_2,        "mesh file #2." );

    }

    void MeshDifference::run(int argc, char** argv)
    {
      //stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);

      std::string file_name1, file_name2;

      bool debug_re=true;
      ParallelMachineFinalize pmf(true);
      RunEnvironment run_environment(&argc, &argv, debug_re);
      process_options(run_environment);
      run_environment.processCommandLine();

      {
        if (run_environment.help_opt) {
          run_environment.printHelp();
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

      //RunEnvironment::doLoadBalance(run_environment.m_comm, file_name);

      PerceptMesh mesh1(3, run_environment.m_comm);
      PerceptMesh mesh2(3, run_environment.m_comm);
      mesh1.open_read_only(file_name1);
      std::cout << "read in file1: " << file_name1 << std::endl;
      mesh2.open_read_only(file_name2);
      std::cout << "read in file2: " << file_name2 << std::endl;

      PerceptMesh::mesh_difference(mesh1, mesh2, "diff", true, true);
    }

      
  }//namespace percept
}//namespace stk


