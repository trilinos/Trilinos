
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
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

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/ProgressMeter.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_adapt/RefinerUtil.hpp>

#include <stk_adapt/UniformRefiner.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <boost/shared_ptr.hpp>

#define ALLOW_MEM_TEST 1

extern double s_timers[10]; // = {0,0,0,0,0,0,0,0,0,0};


namespace stk { 

  namespace adapt {

    //extern void test_memory(int, int);
    void test_memory(percept::PerceptMesh& eMesh, int n_elements, int n_nodes)
    {
      vector<stk::mesh::Entity *> new_elements;
      vector<stk::mesh::Entity *> new_nodes;
      
      eMesh.getBulkData()->modification_begin();

      std::cout << "creating " << n_elements << " elements..." <<std::endl;
      eMesh.createEntities( eMesh.element_rank(), n_elements, new_elements);
      std::cout << "... done creating " << n_elements << " elements" << std::endl;

      std::cout << "creating " << n_nodes << " nodes..." <<std::endl;
      eMesh.createEntities( stk::mesh::fem::FEMMetaData::NODE_RANK, n_nodes, new_nodes);
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
          stk::mesh::Entity& element = *new_elements[i_element];

          for (int j_node = 0; j_node < n_node_per_element; j_node++)
            {
              stk::mesh::Entity& node = *new_nodes[i_node];

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

    static void print_simple_usage(int argc, char **argv)
    {
      std::cout << "AdaptMain::print_simple_usage number of arguments = " << argc << std::endl;
      for (int i = 0; i < argc; i++)
      {
        std::cout << "AdaptMain::print_simple_usage arg[" << i << "]= " << argv[i] << std::endl;
      }
      std::cout << "usage: exe_name [convert|enrich|refine] input_file_name [output_file_name] [number_refines]" << std::endl;
    }

    int adapt_main(int argc, char **argv) ;
    int adapt_main_full_options(int argc, char **argv) ;

    // FIXME
    static int check_for_simple_options(int argc, char **argv) 
    {
      int simple = 0;
      for (int i = 1; i < argc; i++)
        {
          if (std::string(argv[i]) == "refine" || std::string(argv[i]) == "enrich" || std::string(argv[i]) == "convert")  
            return i;
        }
      return simple;
    }

    static bool debug = false;
    int adapt_main_simple_options(int argc_in, char **argv_in) 
    {

      // format: exe_name [convert|enrich|refine] input_file_name
      int simple_options_index = check_for_simple_options(argc_in, argv_in);

      //       if (!simple_options_index)
      //         {
      //           print_simple_usage();
      //           return 1;
      //         }

      int argc = argc_in;
      if (argc != 2 + simple_options_index && argc != 3 + simple_options_index && argc != 4 + simple_options_index )
        {
          print_simple_usage(argc_in, argv_in);
          return 1;
        }

      std::vector<std::string> argv(argc);
      for (int i = 0; i < argc; i++)
        {
          argv[i] = (const char *)argv_in[i];
        }
      if (debug)
        std::cout << "argc = " << argc << " argv= \n" << argv << std::endl;

      std::string exe_name        = argv[0];
      std::string option          = argv[0+simple_options_index];
      if (option != "refine" && option != "enrich" && option != "convert")
        {
          print_simple_usage(argc_in, argv_in);
          return 1;
        }
      std::string input_mesh = argv[1+simple_options_index];
      std::string number_refines = (argc == 3+simple_options_index? argv[2+simple_options_index] : "1");
      int nref=0;

      bool isInt = false;
      try {
         nref = boost::lexical_cast<int>(number_refines);
        isInt = true;
      }
      catch( ... ) 
        {
        }

      std::string output_mesh = input_mesh;
      std::string extension = input_mesh.substr(input_mesh.length()-2,input_mesh.length());
      if (debug) std::cout << " extension= " << extension << std::endl;
      std::string new_ext = "";
      new_ext += "_";
      if (option == "refine")
        new_ext += "refined_"+number_refines+extension;
      else
        new_ext += option+"ed_"+extension;

      if (!isInt && (argc == 3+simple_options_index))
        {
          output_mesh = number_refines;
          number_refines = (argc == 4+simple_options_index? argv[3+simple_options_index] : "1");
          //std::cout << "tmp output_mesh= " << output_mesh << std::endl;
          //std::cout << "tmp number_refines= " << number_refines << std::endl;
        }
      else
        {
          Util::replace(output_mesh, extension, new_ext);
        }
      if (debug) std::cout << " output_mesh= " << output_mesh << std::endl;

      std::vector<std::string> argv_new;
      for (int i = 0; i < simple_options_index; i++)
        argv_new.push_back(argv[i]);
      argv_new.push_back("--input_mesh="+input_mesh);
      argv_new.push_back("--output_mesh="+output_mesh);
      if (option == "refine")
        argv_new.push_back("--refine=DEFAULT");
      else if (option == "enrich")
        argv_new.push_back("--enrich=DEFAULT");
      else
        argv_new.push_back("--convert=Hex8_Tet4_24");
      argv_new.push_back("--load_balance=1");
      argv_new.push_back("--remove_original_elements=1");
      argv_new.push_back("--number_refines="+number_refines);

      if ( debug) std::cout << "new argv = \n" << argv_new << std::endl;
      int argc_new = argv_new.size();
      char **argv_new_cstr = new char*[argc_new];
      for (int i = 0; i < argc_new; i++)
        {
          argv_new_cstr[i] = (char *)argv_new[i].c_str();
        }
      return adapt_main_full_options(argc_new, argv_new_cstr);
    }

    static void dump_args(int argc, char **argv) 
    {
      std::cout << "argc = " << argc << std::endl;
      for (int i = 0; i < argc; i++)
        {
          std::cout << "argv[" << i << "]= " << argv[i] << std::endl;
        }
    }


    int adapt_main(int argc, char **argv) 
    {
      if (debug)
        dump_args(argc, argv);
      // allow positional arguments, etc.
      if (check_for_simple_options(argc, argv))
        return adapt_main_simple_options(argc, argv);
      else
        return adapt_main_full_options(argc, argv);
    }

    int adapt_main_full_options(int argc, char **argv) 
    { 
      EXCEPTWATCH;
      bool debug_re = false;
      RunEnvironment run_environment(&argc, &argv, debug_re);
      unsigned p_rank = stk::parallel_machine_rank(run_environment.m_comm);
      unsigned p_size = stk::parallel_machine_size(run_environment.m_comm);


      std::string options_description_desc = "stk_adapt options";
    
      // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
      std::string input_mesh="";
      std::string input_geometry="";
      std::string output_mesh="";
      std::string block_name_inc = "";
      std::string block_name_exc = "";
      std::string convert="";
      std::string refine="DEFAULT";
      //std::string refine="";
      std::string enrich="";
      bool doRefineMesh = true;
      int load_balance = 1;
      std::string convert_Hex8_Tet4_24 = "Hex8_Tet4_24";      
      int print_info=0;
      int remove_original_elements = 1;
      int number_refines = 1;
      int proc_rank_field = 0;
      int query_only = 0;
      int progress_meter = 0;

      //  Hex8_Tet4_24 (default), Quad4_Quad4_4, Qu
      std::string block_name_desc = 
        "block name(s) to convert: there are 4 options\n"
        "  (1) empty string or option not specified: convert all blocks in the input mesh file\n"
        "  (2) file:my_filename.my_ext (e.g. file:filelist.dat) which will read input block names\n"
        "            from the given file\n"
        "  (3) [+]block_name_1,[+]block_name_2, etc ,block_name_n to include only these blocks, plus sign is optional\n"
        "  (4) a single input block name (e.g. block_3) to be converted \n"
        "  (5) -block_3,-block_5 to exclude blocks from those included (all blocks or include-only blocks), minus sign is mandatory\n"
        "  (6) block_1..block_10 include the range of blocks #1 to #10 \n"
        "  (7) any combination of [+] and - options and range (..) option can be specified \n"
        "Note: wherever you specify block_# this can be replaced with just the #, e.g. \"1,2,4,5\" ";

      std::string convert_options = UniformRefinerPatternBase::s_convert_options;
      std::string refine_options  = UniformRefinerPatternBase::s_refine_options;
      std::string enrich_options  = UniformRefinerPatternBase::s_enrich_options;

      std::string def1= Util::split(convert_options, ", ")[0] ;
      if (0) std::cout << def1 << "tmp split = " << Util::split(convert_options, ", ") << std::endl;
      int test_memory_elements = 0;
      int test_memory_nodes = 0;
      
      //convert_options = "DEFAULT or one of "+convert_options;
      //refine_options = "DEFAULT or one of "+refine_options;
      //enrich_options = "DEFAULT or one of "+enrich_options;

      // : if not specified, use input mesh name appended with _{converted,refined_#refines,enriched}");

      std::string block_name_desc_inc = "which blocks to include, specified as: "+block_name_desc;
      std::string block_name_desc_exc = "which blocks to exclude, specified as: "+block_name_desc;

      run_environment.clp.setOption("convert"                  , &convert                  , convert_options.c_str());
      run_environment.clp.setOption("refine"                   , &refine                   , refine_options.c_str());
      run_environment.clp.setOption("enrich"                   , &enrich                   , enrich_options.c_str());
      run_environment.clp.setOption("input_mesh"               , &input_mesh               , "input mesh name");
      run_environment.clp.setOption("output_mesh"              , &output_mesh              , "output mesh name");

      run_environment.clp.setOption("query_only"               , &query_only               , "query only, no refinement done");
      run_environment.clp.setOption("progress_meter"           , &progress_meter           , "progress meter on or off");

      run_environment.clp.setOption("number_refines"           , &number_refines           , "number of refinement passes");
      run_environment.clp.setOption("block_name"               , &block_name_inc           , block_name_desc_inc.c_str());
      //run_environment.clp.setOption("exclude"                  , &block_name_exc           , block_name_desc_exc.c_str());
      run_environment.clp.setOption("print_info"               , &print_info               , ">= 0  (higher values print more info)");
      run_environment.clp.setOption("load_balance"             , &load_balance             , " load balance (slice/spread) input mesh file");

#if ALLOW_MEM_TEST
      run_environment.clp.setOption("test_memory_elements"     , &test_memory_elements     , " give a number of elements");
      run_environment.clp.setOption("test_memory_nodes"        , &test_memory_nodes        , " give a number of nodes");
#endif
      run_environment.clp.setOption("proc_rank_field"          , &proc_rank_field          , " add an element field to show processor rank");
      run_environment.clp.setOption("remove_original_elements" , &remove_original_elements , " remove original (converted) elements (default=true)");
      run_environment.clp.setOption("input_geometry"           , &input_geometry           , "input geometry name");

      run_environment.processCommandLine(&argc, &argv);

      int result = 0;
      unsigned failed_proc_rank = 0u;

      double t0   = 0.0;
      double t1   = 0.0;
      double cpu0 = 0.0;
      double cpu1 = 0.0;

      try {

        if (convert.length())
          checkInput("convert", convert, convert_options, run_environment);

        if (refine.length())
          checkInput("refine", refine, refine_options, run_environment);

        if (enrich.length())
          checkInput("enrich", enrich, enrich_options, run_environment);

        if (print_info)
          {
            doRefineMesh = false;
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

#if defined( STK_HAS_MPI )
        MPI_Barrier( MPI_COMM_WORLD );
#endif

        if (load_balance)
          {
            RunEnvironment::doLoadBalance(run_environment.m_comm, input_mesh);
          }

        percept::PerceptMesh eMesh(0);  // FIXME
        //percept::PerceptMesh eMesh;  // FIXME

        //unsigned p_size = eMesh.getParallelSize();
        
        eMesh.open(input_mesh);

        Util::setRank(eMesh.getRank());

        Teuchos::RCP<UniformRefinerPatternBase> pattern;

        if (doRefineMesh)
          {
            // FIXME move this next block of code to a method on UniformRefiner
            BlockNamesType block_names(stk::percept::EntityRankEnd+1u);
            if (block_name_inc.length())
              {
                block_names = RefinerUtil::getBlockNames(block_name_inc, eMesh.getRank(), eMesh);
                if (1)
                  {
                    eMesh.commit();
                    block_names = RefinerUtil::correctBlockNamesForPartPartConsistency(eMesh, block_names);
                    eMesh.close();
                    eMesh.open(input_mesh);
                  }
              }
            
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

        stk::mesh::FieldBase* proc_rank_field_ptr = 0;
        if (proc_rank_field)
        {
          proc_rank_field_ptr = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
        }

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
        
        if (print_info)
          {
            eMesh.printInfo("convert", print_info);
          }

        // FIXME
        if (0)
          {
            eMesh.printInfo("before convert", 2);
            exit(1);
          }


        if (doRefineMesh)
          {
            
            t0 =  stk::wall_time(); 
            cpu0 = stk::cpu_time();


            UniformRefiner breaker(eMesh, *pattern, proc_rank_field_ptr);

            ProgressMeter pm(breaker);
            //pm.setActive(true);

            //std::cout << "P[" << p_rank << ", " << p_size << "] input_geometry = " << input_geometry << std::endl; 

            if (input_geometry != "")
                breaker.setGeometryFile(input_geometry);
            breaker.setRemoveOldElements(remove_original_elements);
            breaker.setQueryPassOnly(query_only == 1);
            breaker.setDoProgressMeter(progress_meter == 1);
            //breaker.setIgnoreSideSets(true);

            for (int iBreak = 0; iBreak < number_refines; iBreak++)
              {
                if (!eMesh.getRank())
                  {
                    std::cout << "Refinement pass # " << (iBreak+1) << " start..." << std::endl;
                  }
                //breaker.setPassNumber(iBreak);
                breaker.doBreak();
                if (!eMesh.getRank())
                  {
                    std::cout << std::endl;
                    int ib = iBreak;
                    if (!query_only) ib = 0;
                    RefinementInfoByType::printTable(std::cout, breaker.getRefinementInfoByType(), ib , true);
                    std::cout << std::endl;
                  }
                
              }
            t1 =  stk::wall_time(); 
            cpu1 = stk::cpu_time();

            eMesh.saveAs(output_mesh);
          }

#if 0
        for (int itime=0; itime < 10; itime++)
          {
            std::cout << "tmp timer[" << itime << "]= " << s_timers[itime] << " " << s_timers[itime]/s_timers[3]*100 << " %" << std::endl;
          }
#endif
        
      }
      catch ( const std::exception * X ) {
        std::cout << "AdaptMain::  unexpected exception POINTER: " << X->what() << std::endl;
        failed_proc_rank = p_rank+1u;
      }
      catch ( const std::exception & X ) {
        std::cout << "AdaptMain:: unexpected exception: " << X.what() << std::endl;
        failed_proc_rank = p_rank+1u;
      }
      catch( ... ) {
        std::cout << "AdaptMain::  ... exception" << std::endl;
        failed_proc_rank = p_rank+1u;
      }

      stk::all_reduce( run_environment.m_comm, stk::ReduceSum<1>( &failed_proc_rank ) );
      if (failed_proc_rank)
        {
          stk::percept::pout() << "P[" << p_rank << "]  exception found on processor " << (failed_proc_rank-1) << "\n";
          exit(1);
        }

      stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  wall clock time on processor [" << p_rank << ", " << p_size << "]= " << (t1-t0) << " (sec) "
                           << " cpu time= " << (cpu1 - cpu0) << " (sec)\n";
      std::cout << "P[" << p_rank << ", " << p_size << "]  wall clock time on processor [" << p_rank << ", " << p_size << "]= " << (t1-t0) << " (sec) "
                << " cpu time= " << (cpu1 - cpu0) << " (sec) " << std::endl;

      double cpuMax = (cpu1-cpu0);
      double wallMax = (t1-t0);
      double cpuSum = (cpu1-cpu0);

      stk::all_reduce( run_environment.m_comm, stk::ReduceSum<1>( &cpuSum ) );
      stk::all_reduce( run_environment.m_comm, stk::ReduceMax<1>( &cpuMax ) );
      stk::all_reduce( run_environment.m_comm, stk::ReduceMax<1>( &wallMax ) );

      if (0 == p_rank)
        {
          stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)\n";
          stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)\n";
          stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)\n";
          std::cout << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)" << std::endl;
          std::cout << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)" << std::endl;
          std::cout << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)" << std::endl;
        }

      return result;
    }

  }
}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) { 

  return stk::adapt::adapt_main(argc, argv);
}
//#endif
