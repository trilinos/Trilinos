#include <gtest/gtest.h>

#include "Verifier.hpp"

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/RunEnvironment.hpp>

namespace stk
{
  namespace percept
  {

    TEST(topo, test1)
    {
      EXPECT_TRUE(1);
    }

    static void process_options()
    {
      bopt::options_description desc("verifier options", 132);


      desc.add_options()
        ("help,h",        "produce help message")
        ("mesh,m",        bopt::value<std::string>(), "mesh file." )
        //("directory,d",   bopt::value<std::string>(), "working directory with trailing '/'" )
        ("printTable,p" ,   "print a table of min/max jacobian and quality measures " )
        ("fullPrint,f" ,    "print every element's jacobian" )

        // bogus - need to strip off the sierra.sh options ??
        (",o",             "un-used - for sierra.sh internal use only")
        (",d",             "un-used - for sierra.sh internal use only")
        (",i",             "un-used - for sierra.sh internal use only")
        ;

      stk::get_options_description().add(desc);

    }

    void Verifier::verify(int argc, char** argv)
    {
      //stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);

      bool printTable = false;
      bool fullJacobianPrints = false;
      std::string file_name;
      
      //bopt::variables_map vm = process_options(parallel_machine, argc, argv);
      process_options();

      RunEnvironment run_environment(&argc, &argv);

      bopt::variables_map& vm = stk::get_variables_map();


      {
        //bopt::options_description &desc = stk::get_options_description();
        if (vm.count("exit")) {
          exit(1);
        }

        if (vm.count("help")) {
          run_environment.printHelp();
          std::exit(EXIT_SUCCESS);
        }

        std::string working_directory = "";
        std::string in_filename = "";

        if (vm.count("mesh")) 
          {
            in_filename = boost::any_cast<std::string>(vm["mesh"].value());;
    
            if (vm.count("directory")) {
              working_directory = boost::any_cast<std::string>(vm["directory"].value());
            }

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

        if (vm.count("fullPrint"))
          {
            fullJacobianPrints = true;
          }
        if (vm.count("printTable"))
          {
            printTable = true;
          }
      }

      RunEnvironment::doLoadBalance(run_environment.m_comm, file_name);

      PerceptMesh mesh(run_environment.m_comm);
      mesh.openReadOnly(file_name);
      std::cout << "read in file: " << file_name << std::endl;

      TopologyVerifier topoVerifier;
      bool isBad= topoVerifier.isTopologyBad( *mesh.getBulkData() );
      if (isBad)
        {
          std::cout << "found bad topology" << std::endl;
        }

      GeometryVerifier geomVerifier(fullJacobianPrints);
      isBad= geomVerifier.isGeometryBad( *mesh.getBulkData(), printTable );
      if (isBad)
        {
          std::cout << "found bad jacobian" << std::endl;
        }

    }

    //void MeshGeometryVerifier

      
  }//namespace percept
}//namespace stk


