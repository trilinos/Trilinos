
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
#include <stdint.h>
#include <map>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/ProgressMeter.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/MallocUsed.h>

#include <stk_mesh/base/MemoryUsage.hpp>

#include <stk_adapt/RefinerUtil.hpp>

#include <stk_adapt/UniformRefiner.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <Ioss_Utils.h>

#include <stk_adapt/SerializeNodeRegistry.hpp>

#define ALLOW_MEM_TEST 1

extern double s_timers[10]; // = {0,0,0,0,0,0,0,0,0,0};


namespace stk { 

  namespace adapt {

    static int s_spatialDim=0;

    typedef uint64_t MemorySizeType;

    BOOST_STATIC_ASSERT(sizeof(uint64_t) == sizeof(size_t));

    static double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }


    struct MemoryMultipliers
    {
      MemorySizeType num_hex8;
      MemorySizeType num_tet4;
      MemorySizeType num_nodes;

      typedef MemorySizeType MemMultType;
      MemMultType mult_hex8;
      MemMultType mult_tet4;
      MemMultType mult_nodes;
      
      MemoryMultipliers(MemMultType mult_hex8=1490, MemMultType mult_tet4=702, MemMultType mult_nodes=0):
        num_hex8(0ul),
        num_tet4(0ul),
        num_nodes(0ul),
        mult_hex8(mult_hex8), 
        mult_tet4(mult_tet4), 
        mult_nodes(mult_nodes)
      {
      }

      void read_simple(std::string file_name)
      {
        std::ifstream file(file_name.c_str());
        if (file.good())
          file >> mult_hex8 >> mult_tet4 >> mult_nodes;
        //std::string line1;
        //file >> line1;
        //std::cout << "mult_hex8= " << mult_hex8 << " mult_tet4= " << mult_tet4 << " mult_nodes=" << mult_nodes << std::endl;
      }

      MemorySizeType estimate_memory()
      {
        return mult_nodes*num_nodes + mult_hex8*num_hex8 + mult_tet4*num_tet4;
      }

      MemorySizeType estimate_memory(std::vector<RefinementInfoByType>& refInfo)
      {
        num_hex8=0ul;
        num_tet4=0ul;
        num_nodes=0ul;

        for (unsigned i = 0; i < refInfo.size(); i++)
          {
            num_nodes= refInfo[0].m_numNewNodes;
            //std::cout << "irank, rank, m_numNewNodes, m_numNewElems= " << i << " " << refInfo[i].m_rank << " " << refInfo[i].m_numNewNodes 
            //<< " " << refInfo[i].m_numNewElemsLast
            //<< std::endl;
            
            //             if (refInfo[i].m_rank == 0)
            //               {
            //                 num_nodes += refInfo[i].m_numNewNodes;
            //               }
            //             else
              {
                switch(refInfo[i].m_topology.getKey())
                  {
                  case shards::Hexahedron<8>::key:
                    num_hex8 += refInfo[i].m_numNewElemsLast;
                    break;
                  case shards::Tetrahedron<4>::key:
                    num_tet4 += refInfo[i].m_numNewElemsLast;
                    break;
                  default:
                    break;
                  }
              }
          }

        return estimate_memory();
      }
      
      static void process_estimate(MemorySizeType tot_mem, PerceptMesh& eMesh, std::vector<RefinementInfoByType>& refInfo, std::string memory_multipliers_file, std::string input_file)
      {
        //const stk::ParallelMachine& comm = eMesh.get_bulk_data()->parallel();

        // this is a data gather pass
        if (tot_mem)
          {
            /*
              mesh::Selector sel_locally_owned(eMesh.get_fem_meta_data()->locally_owned_part());
              mesh::Selector sel_globally_shared(eMesh.get_fem_meta_data()->globally_shared_part());
              mesh::Selector sel_universal(eMesh.get_fem_meta_data()->universal_part());
          
              std::vector<unsigned> count ;
              stk::mesh::count_entities( sel_universal, *eMesh.get_bulk_data(), count );
      
              unsigned nnodes = count[0];

              stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
              stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );
            */
            MemoryMultipliers memMults;
            // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
            if (memory_multipliers_file.size())
              memMults.read_simple(memory_multipliers_file);
            RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
            MemorySizeType estMem = memMults.estimate_memory(refInfo);
            //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
            if (eMesh.get_rank() == 0)
              {
                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory[MB]= " << MegaByte(tot_mem) 
                          << " estMem[MB]= " << MegaByte(estMem) 
                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

                std::cout << "(*MemEstMM: " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem) 
                          << ", " << MegaByte(estMem) << "}" << std::endl;
              }

          }
        else
          {
            // this is an estimate multipliers pass (computes memory using current multipliers)
            MemoryMultipliers memMults;
            // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
            if (memory_multipliers_file.size())
              memMults.read_simple(memory_multipliers_file);
            RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
            MemorySizeType estMem = memMults.estimate_memory(refInfo);
            //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
            if (eMesh.get_rank() == 0)
              {
                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory[MB]= " << MegaByte(tot_mem) 
                          << " estMem[MB]= " << MegaByte(estMem) 
                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

                std::cout << "(*MemEstMM: " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem) 
                          << ", " << MegaByte(estMem) << "}" << std::endl;
              }

          }
      }
    };


    static MemorySizeType memory_dump(int dump_level, const stk::ParallelMachine& comm, stk::mesh::BulkData& bulkData, NodeRegistry* node_reg, std::string msg)
    {
      MemorySizeType returned_total_memory;
      MemorySizeType malloc_used_0 = malloc_used();
      MemorySizeType malloc_footprint_0 = malloc_footprint();
      MemorySizeType malloc_max_footprint_0 = malloc_max_footprint();
      MemorySizeType MB = 1024*1024;

      stk::all_reduce( comm, stk::ReduceSum<1>( &malloc_used_0 ) );
      stk::all_reduce( comm, stk::ReduceSum<1>( &malloc_footprint_0 ) );
      stk::all_reduce( comm, stk::ReduceSum<1>( &malloc_max_footprint_0 ) );

      stk::mesh::MemoryUsage mem_usage;
      stk::mesh::compute_memory_usage(bulkData, mem_usage);

      MemorySizeType node_reg_mem = 0;
      if (node_reg) 
        node_reg_mem = node_reg->get_memory_usage();
      MemorySizeType node_reg_mem_sum = node_reg_mem;
      MemorySizeType node_reg_mem_max = node_reg_mem;
      stk::all_reduce( comm, stk::ReduceSum<1>( &node_reg_mem_sum ) );
      stk::all_reduce( comm, stk::ReduceMax<1>( &node_reg_mem_max ) );

      MemorySizeType mem_total_bytes_sum = mem_usage.total_bytes;
      MemorySizeType mem_total_bytes_max = mem_usage.total_bytes;
      stk::all_reduce( comm, stk::ReduceSum<1>( &mem_total_bytes_sum ) );
      stk::all_reduce( comm, stk::ReduceMax<1>( &mem_total_bytes_max ) );
      if (bulkData.parallel_rank() == 0)
        {
#if !defined(SIERRA_PTMALLOC3_ALLOCATOR) && !defined(SIERRA_PTMALLOC2_ALLOCATOR)
          std::cout << "WARNING: ptmalloc2|3 not compiled in so malloc_used info unavailable.  Recompile with e.g. 'bake allocator=ptmalloc2 (or 3)'.  Printing zeros..." << std::endl;

          //  intelcray

#else

#endif

          if (dump_level > 1)

            {
              std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump stk_mesh counted memory usage at stage [" << msg << "] " 
                " parallel sum, max memory [MB]= " << ((double)mem_total_bytes_sum)/MB << " , " << ((double)mem_total_bytes_max)/MB << std::endl;
              if (dump_level > 2)
                stk::mesh::print_memory_usage(mem_usage, std::cout);

              std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump malloc_used total (sum all proc) at stage [" << msg << "] = " 
                        << " malloc_used malloc_footprint malloc_max_footprint [MB]= " << ((double)malloc_used_0)/MB 
                        << " , " << ((double)malloc_footprint_0)/MB 
                        << " , " << ((double)malloc_max_footprint_0)/MB 
                        << std::endl;
            }
          
          {
            std::cout << "AdaptMain::memory_dump summary for " << msg << " : stk_mesh [sum], NodeRegistry [sum], ptmalloc[sum] [MB] "
                      << ((double)mem_total_bytes_sum)/MB << " , " 
                      << ((double)node_reg_mem_sum)/MB << " , " 
                      << ((double)malloc_used_0)/MB 
                      << std::endl;
          }

        }
      if (malloc_used_0)
        returned_total_memory = malloc_used_0;
      else
        returned_total_memory = mem_total_bytes_sum+node_reg_mem_sum;

      return returned_total_memory;
    }

    //extern void test_memory(int, int);
    void test_memory(percept::PerceptMesh& eMesh, MemorySizeType n_elements, MemorySizeType n_nodes)
    {
      vector<stk::mesh::Entity *> new_elements;
      vector<stk::mesh::Entity *> new_nodes;
      
      eMesh.get_bulk_data()->modification_begin();

      std::cout << "creating " << n_elements << " elements..." <<std::endl;
      eMesh.createEntities( eMesh.element_rank(), n_elements, new_elements);
      std::cout << "... done creating " << n_elements << " elements" << std::endl;

      std::cout << "creating " << n_nodes << " nodes..." <<std::endl;
      eMesh.createEntities( stk::mesh::fem::FEMMetaData::NODE_RANK, n_nodes, new_nodes);
      std::cout << "... done creating " << n_nodes << " nodes" << std::endl;

      MemorySizeType num_prints = std::min(static_cast<MemorySizeType>(100ul), n_elements);
      MemorySizeType print_mod = n_elements/num_prints;
      MemorySizeType i_node = 0;
      int n_node_per_element = 8; 
      for (MemorySizeType i_element = 0; i_element < n_elements; i_element++)
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

              eMesh.get_bulk_data()->declare_relation(element, node, j_node);
              
              i_node++;
              if (i_node >= n_nodes-1)
                i_node = 0;
            }
        }

      std::cout << " doing modification_end ... " << std::endl;
      eMesh.get_bulk_data()->modification_end();
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
        (void)nref;
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
      int ret_val = adapt_main_full_options(argc_new, argv_new_cstr);
      delete[] argv_new_cstr;
      return ret_val;
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
      std::string refine="";
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
      int smooth_geometry = 0;
      int remove_geometry_blocks = 0;
      int sync_io_regions = 1;
      int delete_parents = 1;
      int print_memory_usage = 0;
      // a list of comma-separated names like Entity, Relation, Field, etc.
      std::string memory_multipliers_file="";
      int estimate_memory_usage=0;
      int streaming_size=0;
      int streaming_rank=0;
      int streaming_pass_start= -2;  // FIXME - change to not start from -1 below
      int streaming_pass_end= -2;
      //std::string streaming_instruction="";
      int streaming_W = 0;
      int streaming_iW = 0;

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

      int help = 0;

      run_environment.clp.setOption("help"                     , &help                     , "print this usage message");
      run_environment.clp.setOption("convert"                  , &convert                  , convert_options.c_str());
      run_environment.clp.setOption("refine"                   , &refine                   , refine_options.c_str());
      run_environment.clp.setOption("enrich"                   , &enrich                   , enrich_options.c_str());
      run_environment.clp.setOption("input_mesh"               , &input_mesh               , "input mesh name");
      run_environment.clp.setOption("output_mesh"              , &output_mesh              , "output mesh name");

      run_environment.clp.setOption("query_only"               , &query_only               , "query only, no refinement done");
      run_environment.clp.setOption("progress_meter"           , &progress_meter           , "progress meter on or off");
      run_environment.clp.setOption("smooth_geometry"          , &smooth_geometry          , "smooth geometry - moves nodes after geometry projection to try to avoid bad meshes");
      run_environment.clp.setOption("remove_geometry_blocks"   , &remove_geometry_blocks   , "remove geometry blocks from output Exodus file after refinement/geometry projection");
      run_environment.clp.setOption("sync_io_regions"          , &sync_io_regions          , "synchronize input/output region's Exodus id's");
      run_environment.clp.setOption("delete_parents"           , &delete_parents           , "DEBUG: delete parents from a nested, multi-refine mesh - used for debugging");

      run_environment.clp.setOption("number_refines"           , &number_refines           , "number of refinement passes");
      run_environment.clp.setOption("block_name"               , &block_name_inc           , block_name_desc_inc.c_str());
      //run_environment.clp.setOption("exclude"                  , &block_name_exc           , block_name_desc_exc.c_str());
      run_environment.clp.setOption("print_info"               , &print_info               , ">= 0  (higher values print more info)");
      run_environment.clp.setOption("load_balance"             , &load_balance             , " load balance (slice/spread) input mesh file");

      run_environment.clp.setOption("memory_multipliers_file"  , &memory_multipliers_file  ,
                                    "  filename with 3 space-separated entries, with estimate for bytes-per-hex8 tet4 and nodes, e.g. 300 280 200");
                                    //" filename with a comma-separated string of class names and memory multipliers, one per line, \n  eg Node,80\nField,10\nHexahedron<8>,90");
      run_environment.clp.setOption("estimate_memory_usage"    , &estimate_memory_usage    ,
                                    " if query_only=1, use multipliers from memory_multipliers_file to estimate memory to be used in refinements.\n  If query_only=0, gather memory data and write to the file.");

#if ALLOW_MEM_TEST
      run_environment.clp.setOption("test_memory_elements"     , &test_memory_elements     , " give a number of elements");
      run_environment.clp.setOption("test_memory_nodes"        , &test_memory_nodes        , " give a number of nodes");
#endif
      run_environment.clp.setOption("proc_rank_field"          , &proc_rank_field          , " add an element field to show processor rank");
      run_environment.clp.setOption("remove_original_elements" , &remove_original_elements , " remove original (converted) elements (default=true)");
      run_environment.clp.setOption("input_geometry"           , &input_geometry           , "input geometry name");
      run_environment.clp.setOption("streaming_size"           , &streaming_size      , 
                                    "INTERNAL use only by python script streaming refinement interface:\n"
                                    "  run in streaming mode - this number specifies how many virtual procs the mesh is split into\n"
                                    "    i.e. we expect to see files like file.e.N.iN where N = streaming_size iN=0..N");
      run_environment.clp.setOption("streaming_rank"           , &streaming_rank     , 
                                    "INTERNAL use only by python script streaming refinement interface:\n"
                                    "  run in streaming mode - this number specifies which virtual proc this is.");
      run_environment.clp.setOption("streaming_pass_start"     , &streaming_pass_start           , 
                                    "INTERNAL use only by python script streaming refinement interface:\n");
      run_environment.clp.setOption("streaming_pass_end"       , &streaming_pass_end           , 
                                    "INTERNAL use only by python script streaming refinement interface:\n");
      run_environment.clp.setOption("streaming_W"              , &streaming_W           , 
                                    "INTERNAL use only by python script streaming refinement interface:\n");
      run_environment.clp.setOption("streaming_iW"             , &streaming_iW          , 
                                    "INTERNAL use only by python script streaming refinement interface:\n");
      run_environment.clp.setOption("print_memory_usage"       , &print_memory_usage       , "print memory usage");

      int err_clp = run_environment.processCommandLine();
      if (err_clp) return err_clp;

      int result = 0;
      unsigned failed_proc_rank = 0u;

      double t0   = 0.0;
      double t1   = 0.0;
      double cpu0 = 0.0;
      double cpu1 = 0.0;

      if (convert.length())
        checkInput("convert", convert, convert_options, run_environment);

      if (enrich.length())
        checkInput("enrich", enrich, enrich_options, run_environment);

      if (refine.length())
        checkInput("refine", refine, refine_options, run_environment);

      if (print_info)
        {
          doRefineMesh = false;
        }

      if (help
          || input_mesh.length() == 0 
          || output_mesh.length() == 0
          || (convert.length() == 0 && refine.length()==0 && enrich.length()==0)
          //||  not (convert == "Hex8_Tet4_24" || convert == "Quad4_Quad4_4" || convert == "Quad4_Tri3_6")
          )
        {
          run_environment.printHelp();
          return 1;
        }

#if defined( STK_HAS_MPI )
      MPI_Barrier( MPI_COMM_WORLD );
#endif

      std::string input_mesh_save = input_mesh;
      std::string output_mesh_save = output_mesh;

      // streaming
      int m_M = 1;
      int m_W = 1;
      int m_iW = 0;
      int m_M_0 = 0;
      int m_M_1 = 0;
#if STK_ADAPT_HAVE_YAML_CPP
      if (streaming_size)
        {
          m_M = streaming_size;
          m_W = streaming_W ? streaming_W : 1;
          m_iW = streaming_W ? streaming_iW : 0;
          SerializeNodeRegistry::getStreamingPiece(m_M, m_W, m_iW, m_M_0, m_M_1);
        }
      //std::cout << "tmp srk AdaptMain: " << PERCEPT_OUT(streaming_size) << PERCEPT_OUT(streaming_W) << PERCEPT_OUT(streaming_iW) << PERCEPT_OUT(m_M) << PERCEPT_OUT(m_W) << PERCEPT_OUT(m_iW) << PERCEPT_OUT(m_M_0) << PERCEPT_OUT(m_M_1) << std::endl;
#endif

      // FIXME - starting from -1 pass is bogus
#if STK_ADAPT_HAVE_YAML_CPP
      int remove_original_elements_save = remove_original_elements;
      int delete_parents_save = delete_parents;

      if ((streaming_pass_start == -2 && streaming_pass_end != -2) ||
          (streaming_pass_start != -2 && streaming_pass_end == -2))
        {
          throw std::runtime_error("must specify both streaming_pass_start and streaming_pass_end");
        }

      if (streaming_pass_start != -2 && (streaming_pass_start < -1  || streaming_pass_start > SerializeNodeRegistry::MaxPass))
        {
          throw std::runtime_error("streaming_pass_start bad value");
        }
      if (streaming_pass_end != -2 && (streaming_pass_end < -1  || streaming_pass_end > SerializeNodeRegistry::MaxPass))
        {
          throw std::runtime_error("streaming_pass_end bad value");
        }
      if (streaming_pass_end != -2 && (streaming_pass_end < streaming_pass_start))
        {
          throw std::runtime_error("streaming_pass_start > streaming_pass_end");
        }

      if (streaming_pass_end != -2)
        {
          std::cout << "\n\nWARNING: running passes from command line: streaming_pass_start,end= [" << streaming_pass_start << ", " << streaming_pass_end << "]\n\n";
        }

      // allow for driving this from a script
      if (streaming_pass_start == -2 && streaming_pass_end == -2)
        {
          streaming_pass_start = streaming_size ? -1 : 0;
          streaming_pass_end = streaming_size ? SerializeNodeRegistry::MaxPass : 0;
        }
#else
      streaming_pass_start = 0;
      streaming_pass_end = 0;
#endif

      bool do_normal_pass = true;
      if (m_W > 1 && (m_iW == -1 || m_iW == m_W))
        {
          do_normal_pass = false;
        }

      if (streaming_size && s_spatialDim == 0)
        {
          PerceptMesh eMesh(0);
          std::string mesh_name = Ioss::Utils::decode_filename(input_mesh_save, 0, m_M);
          eMesh.open(mesh_name);
          if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
          s_spatialDim = eMesh.get_spatial_dim();
          VERIFY_OP_ON(s_spatialDim, >=, 2, "AdaptMain bad spatial_dim");
        }
      
      for (int i_pass=streaming_pass_start; i_pass <= streaming_pass_end; i_pass++)
        {

#if STK_ADAPT_HAVE_YAML_CPP
          if (streaming_size) 
            {

              // init pass
              if (m_W > 1 && m_iW == -1)
                {
                  PerceptMesh eMesh(3);
                  eMesh.openEmpty();
                  SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, m_M_0,  m_W, m_iW, m_M_0, m_M_1);
                  snr.pass_init(i_pass);
                  VERIFY_OP_ON(streaming_pass_start, ==, streaming_pass_end, "must have one pass at a time with W > 1");
                  continue;
                }

              remove_original_elements = 0;
              delete_parents = 0;
              if (i_pass == 1)
                {
                  remove_original_elements = remove_original_elements_save;
                  delete_parents = delete_parents_save;
                }
            }
          //std::cout << "tmp srk i_pass= " << i_pass << " delete_parents= " << delete_parents << " remove_original_elements= " << remove_original_elements << std::endl;
#endif

          if (do_normal_pass)
            {
              if (streaming_size)
                {
                  // special cases

                  // read all files, get global parts and node info
                  if (i_pass == -1)
                    {
#if STK_ADAPT_HAVE_YAML_CPP
                      PerceptMesh eMesh(s_spatialDim);
                      eMesh.openEmpty();
                      SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                      snr.pass(i_pass);
#else
                      throw std::runtime_error("must have YAML for streaming refine");
#endif
                      continue;
                    }

                  if (i_pass == 2)
                    {
#if STK_ADAPT_HAVE_YAML_CPP
                      if (m_W == 1)
                        {
                          //  no exodus files i/o
                          PerceptMesh eMesh(s_spatialDim);
                          eMesh.openEmpty();
                          SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                          snr.pass(i_pass);
                        }

#else
                      throw std::runtime_error("must have YAML for streaming refine");
#endif
                      continue;
                    }
                }

              for (int m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
                {
                  if (streaming_size)
                    {
                      input_mesh = Ioss::Utils::decode_filename(input_mesh_save, m_iM, m_M);
                      output_mesh = Ioss::Utils::decode_filename(output_mesh_save, m_iM, m_M);

                      if (i_pass == 1)
                        {
                          output_mesh = output_mesh_save+"-pass1";
                          output_mesh = Ioss::Utils::decode_filename(output_mesh, m_iM, m_M);
                        }
                    }

                  try {

                    if (load_balance)
                      {
                        if (streaming_size) throw std::runtime_error("can't load balance and stream");
                        RunEnvironment::doLoadBalance(run_environment.m_comm, input_mesh);
                      }

                    percept::PerceptMesh eMesh(0);  // FIXME
                    //percept::PerceptMesh eMesh;  // FIXME

                    if (streaming_size)
                      {
                        // special cases
                        if (i_pass == 3)
                          {
#if STK_ADAPT_HAVE_YAML_CPP
                            // fini
                            //if ((m_W > 1 && m_iW == m_W) || (m_W == 1 && m_iM == m_M_1))
                            if (m_iM == m_M_1)
                              {
                                NodeRegistry *some_nr = 0;
                                PerceptMesh eMeshEmpty(0);
                                eMeshEmpty.openEmpty();
                                SerializeNodeRegistry snr(eMeshEmpty, some_nr, input_mesh_save, output_mesh_save, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                                if (remove_geometry_blocks) snr.set_geometry_file(input_geometry);
                                snr.pass3_new();
                              }
#else
                            throw std::runtime_error("must have YAML for streaming refine");
#endif
                            continue;
                          }
                      }

                    // ==============  START  ================================================================== 
        
                    if (do_normal_pass)
                      {
                        eMesh.open(input_mesh);
                        if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
                        if (!sync_io_regions) eMesh.set_sync_io_regions(false);
                        if (!s_spatialDim) s_spatialDim = eMesh.get_spatial_dim();

                        Util::setRank(eMesh.get_rank());

                        Teuchos::RCP<UniformRefinerPatternBase> pattern;

                        if (doRefineMesh)
                          {
                            // FIXME move this next block of code to a method on UniformRefiner
                            BlockNamesType block_names(stk::percept::EntityRankEnd+1u);
                            if (block_name_inc.length())
                              {
                                block_names = RefinerUtil::getBlockNames(block_name_inc, eMesh.get_rank(), eMesh);
                                if (1)
                                  {
                                    eMesh.commit();
                                    block_names = RefinerUtil::correctBlockNamesForPartPartConsistency(eMesh, block_names);

                                    eMesh.close();
                                    eMesh.open(input_mesh);
                                    if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
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
                            proc_rank_field_ptr = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
                          }

#if STK_ADAPT_HAVE_YAML_CPP
                        // FIXME - this is this needed? see above
                        // add global parts not found in this file
                        if (streaming_size)
                          {
                            NodeRegistry *some_nr = 0;
                            SerializeNodeRegistry snr(eMesh, some_nr, input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                            snr.declareGlobalParts();
                          }
#endif
                        eMesh.commit();

                        if (print_memory_usage)
                          memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after file open");

                        if (test_memory_nodes && test_memory_elements)
                          {
                            std::cout << "test_memory_elements and nodes are nonzero, will not refine exodus files." << std::endl;

                            test_memory(eMesh, test_memory_elements, test_memory_nodes);

                            if (print_memory_usage)
                              memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after test memory");

                            if (estimate_memory_usage && !query_only)
                              {
                                MemorySizeType tot_mem = memory_dump(false, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after test memory");

                                //std::cout << "MemEst: num_nodes= " << test_memory_nodes << " num_tet4=0 hum_hex8= " << test_memory_elements << " memory= " << MegaByte(tot_mem) << std::endl;
                                //MemoryMultipliers::process_estimate(tot_mem, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file);
                                MemoryMultipliers memMults;
                                if (memory_multipliers_file.size())
                                  memMults.read_simple(memory_multipliers_file);
                                memMults.num_hex8=test_memory_elements;
                                memMults.num_nodes=test_memory_nodes;
                                MemorySizeType estMem = memMults.estimate_memory();
                                //                 std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory= " << MegaByte(tot_mem) 
                                //                           << " estMem= " << MegaByte(estMem) << std::endl;
                                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory[MB]= " << MegaByte(tot_mem) 
                                          << " estMem[MB]= " << MegaByte(estMem) 
                                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;
                                std::cout << "(*MemEstMM: " << input_mesh << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem) 
                                          << ", " << MegaByte(estMem) << "}" << std::endl;

                              }
                            if (estimate_memory_usage && query_only)
                              {
                                //MemoryMultipliers::process_estimate(0, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                              }

                            return 0;
                          }

                        // FIXME
                        if (0)
                          {
                            eMesh.save_as("outtmp.e");
                            exit(1);
                          }
        
                        if (print_info)
                          {
                            eMesh.print_info("convert", print_info);
                          }

                        // FIXME
                        if (0)
                          {
                            eMesh.print_info("before convert", 2);
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
                              {
                                breaker.setGeometryFile(input_geometry);
                                breaker.setSmoothGeometry(smooth_geometry == 1);
                                //breaker.setRemoveGeometryBlocks(remove_geometry_blocks == 1);
                              }
                            breaker.setRemoveOldElements(remove_original_elements);
                            breaker.setQueryPassOnly(query_only == 1);
                            breaker.setDoProgressMeter(progress_meter == 1 && 0 == p_rank);
                            //breaker.setIgnoreSideSets(true);

                            for (int iBreak = 0; iBreak < number_refines; iBreak++)
                              {
                                if (!eMesh.get_rank())
                                  {
                                    std::cout << "Refinement pass # " << (iBreak+1) << " start..." << std::endl;
                                  }
                                //breaker.setPassNumber(iBreak);
                                breaker.doBreak();
                                //RefinementInfoByType::countCurrentNodes(eMesh, breaker.getRefinementInfoByType());
                                if (!eMesh.get_rank())
                                  {
                                    std::cout << std::endl;
                                    int ib = iBreak;
                                    if (!query_only) ib = 0;
                                    RefinementInfoByType::printTable(std::cout, breaker.getRefinementInfoByType(), ib , true);
                                    std::cout << std::endl;
                                  }
                                if (print_memory_usage)
                                  {
                                    memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(),
                                                std::string("after refine pass: ")+toString(iBreak));
                                  }

                                if (estimate_memory_usage && !query_only)
                                  {
                                    MemorySizeType tot_mem = memory_dump(false, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(),
                                                                         std::string("after refine pass: ")+toString(iBreak));
                                    std::cout << "P[" << p_rank << "] tmp srk tot_mem= " << MegaByte(tot_mem) << std::endl;
                                    MemoryMultipliers::process_estimate(tot_mem, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                                  }
                                if (estimate_memory_usage && query_only)
                                  {
                                    RefinementInfoByType::estimateNew(breaker.getRefinementInfoByType(), iBreak);
                                    MemoryMultipliers::process_estimate(0, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                                  }

                              } // iBreak

                            if (streaming_size)
                              {
#if STK_ADAPT_HAVE_YAML_CPP
                                {
                                  SerializeNodeRegistry snr(eMesh, &breaker.getNodeRegistry(), input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                                  snr.pass(i_pass);
                                }

                                if (i_pass == 0)
                                  {
                                    NodeRegistry *some_nr = 0;
                                    SerializeNodeRegistry snr(eMesh, some_nr, input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                                    snr.getGlobalPartMap();
                                    snr.passM1_mergeGlobalParts();
                                    snr.setGlobalPartMap();
                                  }
#else
                                throw std::runtime_error("must have YAML for streaming refine");
#endif
                              }

                            if (delete_parents)
                              breaker.deleteParentElements();

                            t1 =  stk::wall_time(); 
                            cpu1 = stk::cpu_time();

                            stk::percept::pout() << "P[" << p_rank << "] AdaptMain::  saving mesh... \n";
                            std::cout << "P[" << p_rank << "]  AdaptMain:: saving mesh... " << std::endl;
                            if (streaming_size) eMesh.setStreamingSize(m_M);
                            if (remove_geometry_blocks) eMesh.remove_geometry_blocks_on_output(input_geometry);
                            eMesh.save_as(output_mesh);
                            stk::percept::pout() << "P[" << p_rank << "] AdaptMain:: ... mesh saved\n";
                            std::cout << "P[" << p_rank << "]  AdaptMain:: mesh saved" << std::endl;

                            if (print_memory_usage)
                              memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(), "after final save mesh");

                          } // doRefineMesh
                      } // do_normal_pass

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
                } // m_iM

            } // do_normal_pass

#if STK_ADAPT_HAVE_YAML_CPP
          // fini pass
          if (m_W > 1 && m_iW == m_W)
            {
              VERIFY_OP_ON(s_spatialDim, >=, 2, "AdaptMain bad spatial_dim");
              PerceptMesh eMesh(s_spatialDim);
              eMesh.openEmpty();
              SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0,  m_W, m_iW, m_M_0, m_M_1);
              snr.pass_final(i_pass);
            }
#endif

        } // i_pass

      return result;
    }

  }
}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) { 

  int res=0;
  res = stk::adapt::adapt_main(argc, argv);
  stk::adapt::ParallelMachineFinalize pm(true);
  return res;
}
//#endif
