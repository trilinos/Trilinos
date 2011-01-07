
/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <mpi.h>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_adapt/UniformRefiner.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

#define ALLOW_MEM_TEST 1

extern double s_timers[10]; // = {0,0,0,0,0,0,0,0,0,0};


namespace stk { 

  namespace adapt {

    //extern void test_memory(int, int);
    void test_memory(percept::PerceptMesh& eMesh, int n_elements, int n_nodes)
    {
      vector<Entity *> new_elements;
      vector<Entity *> new_nodes;
      
      eMesh.getBulkData()->modification_begin();

      std::cout << "creating " << n_elements << " elements..." <<std::endl;
      eMesh.createEntities( mesh::Element, n_elements, new_elements);
      std::cout << "... done creating " << n_elements << " elements" << std::endl;

      std::cout << "creating " << n_nodes << " nodes..." <<std::endl;
      eMesh.createEntities( mesh::Node, n_nodes, new_nodes);
      std::cout << "... done creating " << n_nodes << " nodes" << std::endl;

      int num_prints = 100;
      int print_mod = n_elements/num_prints;
      int i_node = 0;
      int n_node_per_element = 4; 
      for (int i_element = 0; i_element < n_elements; i_element++)
        {
          if (!i_element || (i_element % print_mod == 0))
            {
              std::cout << "declare_relation for i_element = " << i_element << " [" << n_elements << "] = " << ((double)i_element)/((double)n_elements)*100 << "%"
                        << std::endl;
            }
          //Entity& element = eMesh.getBulkData()->get_entity( stk::mesh::Element , i_element+1);
          Entity& element = *new_elements[i_element];

          for (int j_node = 0; j_node < n_node_per_element; j_node++)
            {
              Entity& node = *new_nodes[i_node];

              eMesh.getBulkData()->declare_relation(element, node, j_node);
              
              i_node++;
              if (i_node >= n_nodes-1)
                i_node = 0;
            }
        }

      std::cout << " doing modification_end ... " << std::endl;
      eMesh.getBulkData()->modification_end();
      std::cout << " done modification_end ... " << std::endl;
      

    }

    using namespace boost::program_options;

    static void 
    print(boost::program_options::options_description& options)
    {
      const std::vector< boost::shared_ptr<option_description> >& m_options = options.options();

      /* Find the maximum width of the option column */
      unsigned width(23);
      unsigned i; // vc6 has broken for loop scoping
      for (i = 0; i < m_options.size(); ++i)
        {
          const option_description& opt = *m_options[i];
          std::stringstream ss;
          ss << "  " << opt.format_name() << ' ' << opt.format_parameter();
          width = std::max(width, static_cast<unsigned>(ss.str().size()));
        }
        
      /* add an additional space to improve readability */
      ++width;
            
      if(0) std::cout << "width= " << width << std::endl;
      //exit(1);
    }

    static void checkInput(std::string option, std::string value, std::string allowed_values, RunEnvironment& run_environment)
    {
      //if (value.length() == 0) return;
      std::vector<std::string> vals = Util::split(allowed_values, ", ");
      for (unsigned i = 0; i < vals.size(); i++)
        {
          if (vals[i] == value)
            return;
        }
      std::ostringstream oss;
      oss << "AdaptMain::checkInput bad option for " << option << " (= " << value << ") \n allowed values = " << allowed_values;
      //throw std::runtime_error(oss.str());
      {
        std::cout << oss.str() << std::endl;
        run_environment.printHelp();
        exit(1);
      }

    }

    int adapt_main(int argc, char **argv) 
    { 
      EXCEPTWATCH;
      //dw_option_mask.mask("search", use_case::LOG_SEARCH, "log search diagnostics");

      bopt::options_description desc("stk_adapt options", 132);
    
      // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
      std::string input_mesh;
      std::string output_mesh;
      std::string block_name;
      std::string convert;
      std::string refine;
      std::string enrich;
      bool doConvert = true;
      bool doLoadBal = true;
      std::string convert_Hex8_Tet4_24 = "Hex8_Tet4_24";      
      int printInfo=0;
      bool remove_original_elements = true;
      unsigned numRefines = 1;
      bool addProcRankField = false;

      //  Hex8_Tet4_24 (default), Quad4_Quad4_4, Qu
      std::string block_name_desc = 
        "block name(s) to convert: there are 4 options\n"
        "  (1) empty string or option not specified: convert all blocks in the input mesh file\n"
        "  (2) file:my_filename.my_ext (e.g. file:filelist.dat) which will read input block names\n"
        "            from the given file\n"
        "  (3) list:block_name_1,block_name_2,...,block_name_n \n"
        "  (4) a single input block name (e.g. block_3) to be converted";

      std::string convert_options = UniformRefinerPatternBase::s_convert_options;
      std::string refine_options  = UniformRefinerPatternBase::s_refine_options;
      std::string enrich_options  = UniformRefinerPatternBase::s_enrich_options;

      std::string def1= Util::split(convert_options, ", ")[0] ;
      if (0) std::cout << def1 << "tmp split = " << Util::split(convert_options, ", ") << std::endl;
      int test_memory_elements = 0;
      int test_memory_nodes = 0;
      
      desc.add_options()
        //         ("convert",    bopt::value<std::string>(&convert)->default_value(Util::split(convert_options, ", ")[0]), convert_options.c_str())
        //         ("refine",     bopt::value<std::string>(&refine)->default_value(Util::split(refine_options, ", ")[0]), refine_options.c_str())
        //         ("enrich",     bopt::value<std::string>(&enrich)->default_value(Util::split(enrich_options, ", ")[0]), enrich_options.c_str())
        ("convert",    bopt::value<std::string>(&convert), convert_options.c_str())
        ("refine",     bopt::value<std::string>(&refine), refine_options.c_str())
        ("enrich",     bopt::value<std::string>(&enrich), enrich_options.c_str())
        ("input_mesh",    bopt::value<std::string>(&input_mesh), "input mesh name")
        ("output_mesh",    bopt::value<std::string>(&output_mesh), "output mesh name")
        ("number_refines",    bopt::value<unsigned>(&numRefines)->default_value(1),   "number of refinement passes")
        ("block_name",    bopt::value<std::string>(&block_name)->default_value(""), block_name_desc.c_str())
        ("print_info",    bopt::value<int>(&printInfo)->default_value(0),                 ">= 0  (higher values print more info)")
        ("load_balance",    bopt::value<bool>(&doLoadBal)->default_value(true), " load balance (slice/spread) input mesh file")
#if ALLOW_MEM_TEST
        ("test_memory_elements",  bopt::value<int>(&test_memory_elements)->default_value(0), " give a number of elements")
        ("test_memory_nodes",  bopt::value<int>(&test_memory_nodes)->default_value(0), " give a number of nodes")
#endif
        ("proc_rank_field",    bopt::value<bool>(&addProcRankField)->default_value(false), " add an element field to show processor rank")
        ("remove_original_elements",    bopt::value<bool>(&remove_original_elements)->default_value(true), " remove original (converted) elements (default=true)")
        ;

      stk::get_options_description().add(desc);
      print(stk::get_options_description());

      int result = 0;
      unsigned failed_proc_rank = 0u;

      RunEnvironment run_environment(&argc, &argv);
      unsigned p_rank = stk::parallel_machine_rank(run_environment.m_comm);


      try {

        //RunEnvironment run_environment(&argc, &argv);


        bopt::variables_map &vm = stk::get_variables_map();  

        if (vm.count("input_mesh"))
          input_mesh = vm["input_mesh"].as<std::string>();
        if (vm.count("output_mesh"))
          output_mesh = vm["output_mesh"].as<std::string>();

        if (vm.count("convert"))
          {
            convert = vm["convert"].as<std::string>();
            checkInput("convert", convert, convert_options, run_environment);
          }
        if (vm.count("refine"))
          {
            refine = vm["refine"].as<std::string>();
            checkInput("refine", refine, refine_options, run_environment);
          }
        if (vm.count("enrich"))
          {
            enrich = vm["enrich"].as<std::string>();
            checkInput("enrich", enrich, enrich_options, run_environment);
          }

        if (printInfo)
          {
            doConvert = false;
          }

        if (input_mesh.length() == 0 
            || output_mesh.length() == 0
            || (convert.length() == 0 && refine.length()==0 && enrich.length()==0)
            //||  not (convert == "Hex8_Tet4_24" || convert == "Quad4_Quad4_4" || convert == "Quad4_Tri3_6")
            )
          {
            run_environment.printHelp();
            exit(1);
          }

        //try {
        
        EXCEPTWATCH;
        MPI_Barrier( MPI_COMM_WORLD );

        if (doLoadBal)
          {
            RunEnvironment::doLoadBalance(run_environment.m_comm, input_mesh);
          }

        percept::PerceptMesh eMesh;

        //unsigned p_size = eMesh.getParallelSize();
        
        Util::setRank(eMesh.getRank());

        eMesh.open(input_mesh);

        Teuchos::RCP<UniformRefinerPatternBase> pattern;

        if (doConvert)
          {
            BlockNamesType block_names = UniformRefiner::getBlockNames(block_name);

            pattern = UniformRefinerPatternBase::createPattern(refine, enrich, convert, eMesh, block_names);

            if (0)
              {
                run_environment.printHelp();
                exit(1);
              }
          }

        if (0)
          {
            std::cout << "tmp convert = " << convert << std::endl;
            std::cout << "tmp refine = " << refine << std::endl;
            std::cout << "tmp enrich = " << enrich << std::endl;
          }

        int scalarDimension = 0; // a scalar

        FieldBase* proc_rank_field = 0;
        if (addProcRankField)
          eMesh.addField("proc_rank", mesh::Element, scalarDimension);

        eMesh.commit();

        if (test_memory_nodes && test_memory_elements)
          {
            test_memory(eMesh, test_memory_elements, test_memory_nodes);
            return 0;
          }

        // FIXME
        if (0)
          {
            eMesh.saveAs("outtmp.e");
            exit(1);
          }
        
        if (printInfo)
          {
            eMesh.printInfo("convert", printInfo);
          }

        // FIXME
        if (0)
          {
            eMesh.printInfo("before convert", 2);
            exit(1);
          }

        if (doConvert)
          {
            UniformRefiner breaker(eMesh, *pattern, proc_rank_field);
            breaker.setRemoveOldElements(remove_original_elements);

            for (unsigned iBreak = 0; iBreak < numRefines; iBreak++)
              {
                if (!eMesh.getRank())
                  {
                    std::cout << "Refinement pass # " << (iBreak+1) << std::endl;
                  }
                breaker.doBreak();
              }

            // FIXME
            //eMesh.printInfo("after convert", 2);

            eMesh.saveAs(output_mesh);
          }

#if 0
        for (int itime=0; itime < 10; itime++)
          {
            std::cout << "tmp timer[" << itime << "]= " << s_timers[itime] << " " << s_timers[itime]/s_timers[3]*100 << " %" << std::endl;
          }
#endif
        // tmp
        //if (p_rank == 1) throw std::runtime_error("test");
        
      }
      catch ( const std::exception * X ) {
        std::cout << "AdaptMain::  unexpected exception POINTER: " << X->what() << std::endl;
        failed_proc_rank = p_rank+1u;
        //exit(1);
      }
      catch ( const std::exception & X ) {
        std::cout << "AdaptMain:: unexpected exception: " << X.what() << std::endl;
        failed_proc_rank = p_rank+1u;
        //exit(1);
      }
      catch( ... ) {
        std::cout << "AdaptMain::  ... exception" << std::endl;
        failed_proc_rank = p_rank+1u;
        //exit(1);
      }

      stk::all_reduce( run_environment.m_comm, stk::ReduceSum<1>( &failed_proc_rank ) );
      if (failed_proc_rank)
        {
          stk::percept::pout() << "P[" << p_rank << "]  exception found on processor " << (failed_proc_rank-1) << "\n";
          exit(1);
        }
      return result;
    }

  }}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) { 

  return stk::adapt::adapt_main(argc, argv);
}
//#endif
