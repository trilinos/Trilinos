// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "ArgvBuilder.h"
#include "ElementPartition.h"
#include "FileUtils.h"
#include "MeshFixture.h"
#include "NodePartition.h"

#include "scopeguard.h"
#include <copy_string_cpp.h>

#include "add_to_log.h" // for add_to_log
#include "alloc.h"
#include "exodusII.h" // for ex_opts, ex_int64_status, etc
#include "fmt/ostream.h"
#include "nem_spread.h"     // for NemSpread, second, etc
#include "ps_pario_const.h" // for Parallel_IO
#include "rf_allo.h"        // for safe_free
#include "vector_data.h"
#include <cstdint> // for int64_t
#include <cstdio>  // for stderr, etc
#include <cstdlib> // for exit
#include <cstring>
#include <random>
#include <unistd.h> // for getopt, optarg, optind

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include "XGetopt.h"
#else
#include "getopt.h" // for getopt
#endif

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_Enumerate.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

using suplib_cpp::array_alloc;
using suplib_cpp::copy_string;

template <typename T, typename INT>
int nem_spread(NemSpread<T, INT> &spreader, const char *salsa_cmd_file, int subcycles, int cycle);

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

namespace {
  enum PartitionType { LINEAR, ROUND_ROBIN, RANDOM };

  inline std::ostream &operator<<(std::ostream &os, PartitionType type)
  {
    switch (type) {
    case LINEAR: os << "LINEAR"; break;
    case ROUND_ROBIN: os << "ROUND_ROBIN"; break;
    case RANDOM: os << "RANDOM"; break;
    default: os << "UNKNOWN"; break;
    }
    return os;
  }

  template <typename INT>
  void internal_fill_processor_id(const Ioss::ElementBlock                  *eb,
                                  const std::vector<utest_util::EntityProc> &procAssign,
                                  std::vector<double>                       &proc_id)
  {
    proc_id.clear();
    proc_id.reserve(eb->entity_count());

    std::vector<INT> elem_ids;

    eb->get_field_data("ids", elem_ids);

    for (INT id : elem_ids) {
      auto lowerBound =
          std::lower_bound(procAssign.begin(), procAssign.end(), id, utest_util::EntityProcLess());
      bool found = !(lowerBound == procAssign.end() || lowerBound->id != id);

      EXPECT_TRUE(found) << "element " << id << " was not assigned a processor";
      proc_id.push_back(lowerBound->proc);
    }
  }

  void fill_processor_id(const Ioss::ElementBlock                  *eb,
                         const std::vector<utest_util::EntityProc> &procAssign,
                         std::vector<double>                       &proc_id)
  {
    if (eb->get_database()->int_byte_size_api() == 8) {
      internal_fill_processor_id<int64_t>(eb, procAssign, proc_id);
    }
    else {
      internal_fill_processor_id<int>(eb, procAssign, proc_id);
    }
  }

  void add_processor_id_field(Ioss::Region                              *region,
                              const std::vector<utest_util::EntityProc> &procAssign)
  {
    region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    const auto &eblocks = region->get_element_blocks();
    for (const auto &eb : eblocks) {
      eb->field_add(
          Ioss::Field("processor_id", Ioss::Field::REAL, "scalar", Ioss::Field::TRANSIENT));
    }
    region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);

    region->begin_mode(Ioss::STATE_TRANSIENT);

    auto step = region->add_state(0.0);
    region->begin_state(step);

    for (const auto &eb : eblocks) {
      std::vector<double> proc_id;
      fill_processor_id(eb, procAssign, proc_id);
      eb->put_field_data("processor_id", proc_id);
    }

    region->end_state(step);
    region->end_mode(Ioss::STATE_TRANSIENT);
  }

  class NemSpreadTester : public utest_util::MeshFixture
  {
  public:
    NemSpreadTester() : utest_util::MeshFixture(3)
    {
      Ioss::Init::Initializer io;
      setup_empty_mesh();
    }

    ~NemSpreadTester() {}

  protected:
    std::string m_uniqueBaseName;
    std::string m_outputFile;
    std::string m_propertyName{"MATERIAL_PROPERTY"};
    std::string m_propertyValue{"KRYPTONITE"};

    std::string m_nemesisFile;
    std::string m_pexFile;
    std::string m_salsaCmdFile;
    bool        m_force64Bit        = false;
    int         m_startProc         = 0;
    int         m_numProcs          = 0;
    int         m_subcycles         = 0;
    int         m_cycle             = -1;
    int         m_selectedChangeSet = 0;

  protected:
    void parse_options(int argc, char *argv[])
    {
      int c;

      // Reset for new scan
      optind = 1;

      while ((c = getopt(argc, argv, "64Vhp:r:s:n:S:c:C:")) != -1) {
        switch (c) {
        case 'h':
          fmt::print(stderr, " usage:\n");
          fmt::print(stderr,
                     "\tnem_spread  [-s <start_proc>] [-n <num_proc>] [-S <subcycles> -c <cycle>] "
                     "[-C change_set_#] [command_file]\n");
          fmt::print(stderr,
                     "\t\tDecompose for processors <start_proc> to <start_proc>+<num_proc>\n");
          fmt::print(stderr, "\t\tDecompose for cycle <cycle> of <subcycle> groups\n");
          fmt::print(stderr, "\t\tRead from change_set `#` (1-based) if specified.\n");
          fmt::print(stderr, "\tnem_spread  [-V] [-h] (show version or usage info)\n");
          fmt::print(stderr, "\tnem_spread  [command file] [<-p Proc> <-r raid #>]\n");
          exit(1);
          break;
        case 'V':
          fmt::print("{} version {}\n", UTIL_NAME, VER_STR);
          exit(0);
          break;
        case 'p': /* Which proc to use? Also for compatibility */ break;
        case 'r': /* raid number.  Seems to be unused; left around for compatibility */ break;
        case 's': /* Start with processor <x> */ sscanf(optarg, "%d", &m_startProc); break;
        case 'n': /* Number of processors to output files for */
          sscanf(optarg, "%d", &m_numProcs);
          break;
        case 'C': /* change_set index <x> */ sscanf(optarg, "%d", &m_selectedChangeSet); break;
        case '6':
        case '4':
          m_force64Bit = true; /* Force storing output mesh using 64bit integers */
          break;
        case 'S': /* Number of subcycles to use (see below) */
          sscanf(optarg, "%d", &m_subcycles);
          break;
        case 'c': /* Which cycle to spread (see below) */ sscanf(optarg, "%d", &m_cycle); break;
        }
      }

      if (optind >= argc) {
        m_salsaCmdFile = "nem_spread.inp";
      }
      else {
        m_salsaCmdFile = argv[optind];
      }
    }

    void create_nemesis_file(const utest_util::Partition &partition)
    {
      int         exoid;
      std::string method1{}, method2{};
      std::string title;

      int cpu_ws = sizeof(float);
      int io_ws  = sizeof(float);

      fmt::print("Outputting load balance to file {}\n", m_nemesisFile);

      /* Create the load balance file */
      /* Attempt to create a netcdf4-format file; if it fails, then assume
         that the netcdf library does not support that mode and fall back
         to classic netcdf3 format.  If that fails, issue an error and
         return failure.
      */
      int mode3 = EX_CLOBBER;
      int mode4 = mode3 | EX_NETCDF4 | EX_NOCLASSIC;

      if (m_force64Bit)
        mode4 |= EX_ALL_INT64_DB;

      if (partition.api_size() == 8) {
        mode4 |= EX_ALL_INT64_API;
      }

      ex_opts(EX_DEFAULT); // Eliminate misleading error if the first ex_create fails, but the
                           // second succeeds.
      if ((exoid = ex_create(m_nemesisFile.c_str(), mode4, &cpu_ws, &io_ws)) < 0) {
        /* If int64api or int64db non-zero, then netcdf-4 format is required, so
           fail now...
        */
        if (m_force64Bit) {
          ASSERT_TRUE(false) << "fatal: failed to create Nemesis netcdf-4 file";
        }
        if ((exoid = ex_create(m_nemesisFile.c_str(), mode3, &cpu_ws, &io_ws)) < 0) {
          ASSERT_TRUE(false) << "fatal: failed to create Nemesis file";
        }
      }
      ON_BLOCK_EXIT(ex_close, exoid);

      /* Set the error reporting value */
      ex_opts(EX_VERBOSE | EX_DEBUG);

      /* Enable compression (if netcdf-4) */
      ex_set_option(exoid, EX_OPT_COMPRESSION_LEVEL, 1);
      ex_set_option(exoid, EX_OPT_COMPRESSION_SHUFFLE, 1);

      /* Create the title */
      method1 = partition.type();

      title = fmt::format("nem_slice {} load balance file", method1);

      method1 = "method1: ";
      method2 = "method2: ";

      method1 += "Linear decomposition";
      method2 += "no refinement";

      method1 += " via bisection";

      /* Output the info records */
      char *info[3];
      info[0] = const_cast<char *>(title.c_str());
      info[1] = const_cast<char *>(method1.c_str());
      info[2] = const_cast<char *>(method2.c_str());

      ASSERT_FALSE(ex_put_info(exoid, 3, info) < 0) << "warning: output of info records failed";

      /* Generate a QA record for the utility */
      auto        now = std::chrono::system_clock::now();
      std::string time =
          fmt::format("{:%T}", std::chrono::time_point_cast<std::chrono::seconds>(now));
      std::string date = fmt::format("{:%Y/%m/%d}", now);

      char qa_date[32];
      char qa_time[32];
      char qa_name[MAX_STR_LENGTH];
      char qa_vers[32];

      copy_string(qa_time, time);
      copy_string(qa_date, date);
      copy_string(qa_name, "nem_spread_unit_tester");
      copy_string(qa_vers, VER_STR);

      char **lqa_record = reinterpret_cast<char **>(array_alloc(1, 4, sizeof(char *)));
      for (int i2 = 0; i2 < 4; i2++) {
        lqa_record[i2] = reinterpret_cast<char *>(array_alloc(1, MAX_STR_LENGTH + 1, sizeof(char)));
      }

      copy_string(lqa_record[0], qa_name, MAX_STR_LENGTH + 1);
      copy_string(lqa_record[1], qa_vers, MAX_STR_LENGTH + 1);
      copy_string(lqa_record[2], qa_date, MAX_STR_LENGTH + 1);
      copy_string(lqa_record[3], qa_time, MAX_STR_LENGTH + 1);

      fmt::print("QA Record:\n");
      for (int i2 = 0; i2 < 4; i2++) {
        fmt::print("\t{}\n", lqa_record[i2]);
      }

      ASSERT_FALSE(ex_put_qa(exoid, 1, reinterpret_cast<char *(*)[4]>(&lqa_record[0])) < 0)
          << "fatal: unable to output QA records";

      /* free up memory */
      for (int i2 = 0; i2 < 4; i2++) {
        free(lqa_record[i2]);
      }

      free(lqa_record);

      partition.write_nemesis_data(exoid);
    }

    void create_pex_file()
    {
      std::string rootdir = "./";

#if !defined(__IOSS_WINDOWS__)
      char *cwd = getcwd(nullptr, 0);
      if (nullptr != cwd) {
        rootdir = cwd;
      }
      free(cwd);

      if (!rootdir.empty() && rootdir.back() != '/') {
        rootdir += "/";
      }
#endif

      std::ofstream os(m_pexFile, std::ofstream::out);
      os << "Input FEM file                   = " << m_outputFile << std::endl;
      os << "LB file                          = " << m_nemesisFile << std::endl;
      os << "Parallel Results File Base Name  = " << m_uniqueBaseName << std::endl;
      os << "File Extension for Spread Files  = " << ".g" << std::endl;
      os << "Number of Processors             = " << m_numProcs << std::endl;
      os << "------------------------------------------------------------\n";
      os << "                Parallel I/O section\n";
      os << "------------------------------------------------------------\n";
      os << "Parallel Disk Info= number=1, offset=1, zeros, nosubdirectory\n";
      os << "Parallel file location = root=" << rootdir << ", subdir=." << std::endl;
    }

    std::vector<utest_util::EntityProc> get_linear_element_partition()
    {
      std::vector<utest_util::EntityProc> procAssign;
      std::vector<unsigned>               procs;

      size_t numElems = get_mesh().get_num_local_elements();
      fill_linear_proc_distribution(numElems, m_numProcs, procs);

      for (size_t i = 0; i < numElems; i++) {
        utest_util::IossElementData elemData = get_mesh().get_local_element(i);
        utest_util::EntityProc      eProc(elemData.id, procs[i]);

        procAssign.push_back(eProc);
      }

      std::sort(procAssign.begin(), procAssign.end(), utest_util::EntityProcLess());
      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_round_robin_element_partition()
    {
      std::vector<utest_util::EntityProc> procAssign;

      size_t numElems = get_mesh().get_num_local_elements();
      for (size_t i = 0; i < numElems; i++) {
        utest_util::IossElementData elemData = get_mesh().get_local_element(i);
        utest_util::EntityProc      eProc(elemData.id, i % m_numProcs);

        procAssign.push_back(eProc);
      }

      std::sort(procAssign.begin(), procAssign.end(), utest_util::EntityProcLess());
      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_random_element_partition()
    {
      std::vector<utest_util::EntityProc> procAssign = get_linear_element_partition();

      unsigned seed = 1234567;
      std::shuffle(procAssign.begin(), procAssign.end(), std::default_random_engine(seed));

      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_element_partition(const PartitionType type)
    {
      switch (type) {
      case LINEAR: return get_linear_element_partition(); break;
      case ROUND_ROBIN: return get_round_robin_element_partition(); break;
      case RANDOM: return get_random_element_partition(); break;
      default:
        std::ostringstream oss;
        oss << "Partition type: " << type << " not supported";
        throw std::runtime_error(oss.str());
      }

      return std::vector<utest_util::EntityProc>{};
    }

    std::vector<utest_util::EntityProc> get_linear_node_partition()
    {
      std::vector<utest_util::EntityProc> procAssign;
      std::vector<unsigned>               procs;

      size_t numNodes = get_mesh().get_num_local_nodes();
      fill_linear_proc_distribution(numNodes, m_numProcs, procs);

      for (size_t i = 0; i < numNodes; i++) {
        utest_util::IossNodeData nodeData = get_mesh().get_local_node(i);
        utest_util::EntityProc   eProc(nodeData.id, procs[i]);

        procAssign.push_back(eProc);
      }

      std::sort(procAssign.begin(), procAssign.end(), utest_util::EntityProcLess());
      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_round_robin_node_partition()
    {
      std::vector<utest_util::EntityProc> procAssign;

      size_t numNodes = get_mesh().get_num_local_nodes();
      for (size_t i = 0; i < numNodes; i++) {
        utest_util::IossNodeData nodeData = get_mesh().get_local_node(i);
        utest_util::EntityProc   eProc(nodeData.id, i % m_numProcs);

        procAssign.push_back(eProc);
      }

      std::sort(procAssign.begin(), procAssign.end(), utest_util::EntityProcLess());
      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_random_node_partition()
    {
      std::vector<utest_util::EntityProc> procAssign = get_linear_node_partition();

      unsigned seed = 1234567;
      std::shuffle(procAssign.begin(), procAssign.end(), std::default_random_engine(seed));

      return procAssign;
    }

    std::vector<utest_util::EntityProc> get_node_partition(const PartitionType type)
    {
      switch (type) {
      case LINEAR: return get_linear_node_partition(); break;
      case ROUND_ROBIN: return get_round_robin_node_partition(); break;
      case RANDOM: return get_random_node_partition(); break;
      default:
        std::ostringstream oss;
        oss << "Partition type: " << type << " not supported";
        throw std::runtime_error(oss.str());
      }

      return std::vector<utest_util::EntityProc>{};
    }

    void initialize_options(const std::string &baseName, int numProcs, bool is64Bit = false)
    {
      m_uniqueBaseName = utest_util::unique_basename(baseName);
      m_outputFile     = utest_util::unique_filename(baseName, "g");
      m_pexFile        = m_outputFile + ".pex";
      m_nemesisFile    = m_outputFile + ".nem";

      std::string numProcString = std::to_string(numProcs);

      utest_util::ArgvBuilder builder;
      builder.addArgument("nem_spread_unit_test");
      builder.addArgument(m_pexFile.c_str());
      builder.addArgument("-n");
      builder.addArgument(numProcString.c_str());

      if (is64Bit) {
        builder.addArgument("-64");
      }

      parse_options(builder.argc(), builder.argv());
    }

    void create_and_verify_input_mesh_file(const std::string &meshDesc)
    {
      setup_mesh(meshDesc);

      // Add a material property to block_1 from textmesh
      auto inputRegion = get_mesh().get_region();
      auto inputDB     = get_mesh().get_database();

      add_material_property_to_element_block(inputRegion, "block_1", m_propertyName,
                                             m_propertyValue);

      // Write the textmesh to file
      Ioss::PropertyManager properties;
      Ioss::MeshCopyOptions options;
      options.verbose           = false;
      options.output_summary    = true;
      options.debug             = false;
      options.ints_64_bit       = m_force64Bit;
      options.data_storage_type = 1;
      options.add_proc_id       = false;
      write_region_to_file(inputRegion, properties, options, m_outputFile);

      // Load each individual file and test for material property
      test_property_from_file(m_outputFile, m_propertyName, m_propertyValue);
    }

    void create_element_partitioning(PartitionType type)
    {
      // Create partitioning
      std::vector<utest_util::EntityProc> procAssign = get_element_partition(type);

      if (m_force64Bit) {
        utest_util::ElementPartition<int64_t> partition(&get_mesh(), procAssign, m_numProcs);
        create_nemesis_file(partition);
      }
      else {
        utest_util::ElementPartition<int> partition(&get_mesh(), procAssign, m_numProcs);
        create_nemesis_file(partition);
      }
    }

    void create_node_partitioning(PartitionType type)
    {
      // Create partitioning
      std::vector<utest_util::EntityProc> procAssign = get_node_partition(type);

      if (m_force64Bit) {
        utest_util::NodePartition<int64_t> partition(&get_mesh(), procAssign, m_numProcs);
        create_nemesis_file(partition);
      }
      else {
        utest_util::NodePartition<int> partition(&get_mesh(), procAssign, m_numProcs);
        create_nemesis_file(partition);
      }
    }

    void setup_input_files_for_element_decompsition(const std::string &meshDesc,
                                                    PartitionType      type = LINEAR)
    {
      // Create mesh file
      create_and_verify_input_mesh_file(meshDesc);

      // Create pex file
      create_pex_file();

      // Create partitioning
      create_element_partitioning(type);
    }

    void setup_input_files_for_node_decompsition(const std::string &meshDesc,
                                                 PartitionType      type = LINEAR)
    {
      // Create mesh file
      create_and_verify_input_mesh_file(meshDesc);

      // Create pex file
      create_pex_file();

      // Create partitioning
      create_node_partitioning(type);
    }

    void run_nem_spread()
    {
      /* initialize some variables */
      ExoFile.clear();
      Exo_LB_File.clear();
      Exo_Res_File.clear();
      Output_File_Base_Name.clear();

      Debug_Flag = -1;

      Num_Nod_Var  = -1;
      Num_Elem_Var = -1;
      Num_Glob_Var = -1;
      Num_Nset_Var = -1;
      Num_Sset_Var = -1;

      PIO_Info.Dsk_List_Cnt        = -1;
      PIO_Info.Num_Dsk_Ctrlrs      = -1;
      PIO_Info.PDsk_Add_Fact       = -1;
      PIO_Info.Zeros               = -1;
      PIO_Info.NoSubdirectory      = 0;
      PIO_Info.Par_Dsk_Root[0]     = '\0';
      PIO_Info.Par_Dsk_SubDirec[0] = '\0';
      PIO_Info.Staged_Writes       = true;

      ExoFile = m_outputFile;

      static char yo[] = "nem_spread";
      // Open the mesh file and determine word sizes...
      int   io_ws  = 0;
      int   cpu_ws = sizeof(float);
      float version;

      int exoid = ex_open(ExoFile.c_str(), EX_READ, &cpu_ws, &io_ws, &version);
      if (exoid <= 0) {
        fmt::print(stderr, "{} ERROR: Could not open the mesh file '{}'!\n", yo, ExoFile);
        FAIL();
      }

      // Determine whether there are any change sets in file.
      // If there are, check whether user specified a specific
      // change set index and set to that one (if valid), or
      // if not specified, set to the first.
      exoid = check_change_sets(exoid, m_selectedChangeSet);
      if (exoid == 0) {
        FAIL();
      }

      // See if any 64-bit integers stored on database...
      int int64api = 0;
      int int64db  = ex_int64_status(exoid) & EX_ALL_INT64_DB;
      if (int64db != 0) {
        int64api     = EX_ALL_INT64_API;
        m_force64Bit = true;
      }

      int status;
      if (io_ws == 4) {
        if (int64api != 0) {
          NemSpread<float, int64_t> spreader;
          spreader.selected_change_set = m_selectedChangeSet;
          spreader.io_ws               = io_ws;
          spreader.int64db             = int64db;
          spreader.int64api            = int64api;
          spreader.force64db           = m_force64Bit;
          spreader.Proc_Info[4]        = m_startProc;
          spreader.Proc_Info[5]        = m_numProcs;
          status = nem_spread(spreader, m_salsaCmdFile.c_str(), m_subcycles, m_cycle);
        }
        else {
          NemSpread<float, int> spreader;
          spreader.selected_change_set = m_selectedChangeSet;
          spreader.io_ws               = io_ws;
          spreader.int64db             = int64db;
          spreader.int64api            = int64api;
          spreader.force64db           = m_force64Bit;
          spreader.Proc_Info[4]        = m_startProc;
          spreader.Proc_Info[5]        = m_numProcs;
          status = nem_spread(spreader, m_salsaCmdFile.c_str(), m_subcycles, m_cycle);
        }
      }
      else {
        if (int64api != 0) {
          NemSpread<double, int64_t> spreader;
          spreader.selected_change_set = m_selectedChangeSet;
          spreader.io_ws               = io_ws;
          spreader.int64db             = int64db;
          spreader.int64api            = int64api;
          spreader.force64db           = m_force64Bit;
          spreader.Proc_Info[4]        = m_startProc;
          spreader.Proc_Info[5]        = m_numProcs;
          status = nem_spread(spreader, m_salsaCmdFile.c_str(), m_subcycles, m_cycle);
        }
        else {
          NemSpread<double, int> spreader;
          spreader.selected_change_set = m_selectedChangeSet;
          spreader.io_ws               = io_ws;
          spreader.int64db             = int64db;
          spreader.int64api            = int64api;
          spreader.force64db           = m_force64Bit;
          spreader.Proc_Info[4]        = m_startProc;
          spreader.Proc_Info[5]        = m_numProcs;
          status = nem_spread(spreader, m_salsaCmdFile.c_str(), m_subcycles, m_cycle);
        }
      }

      ex_close(exoid);
    }

    void verify_material_property_in_nem_spread_output()
    {
      // Test the material property from the output files
      for (int i = 0; i < m_numProcs; i++) {
        std::string decodedFileName = Ioss::Utils::decode_filename(m_outputFile, i, m_numProcs);
        test_property_from_file(Ioss::ParallelUtils::comm_self(), decodedFileName, m_propertyName,
                                m_propertyValue);
      }
    }

    void cleanup()
    {
      // Delete textmesh file
      unlink(m_outputFile.c_str());

      // Delete pex file
      unlink(m_pexFile.c_str());

      // Delete nemesis file
      unlink(m_nemesisFile.c_str());

      // Delete parallel output files
      for (int i = 0; i < m_numProcs; i++) {
        std::string decodedFileName = Ioss::Utils::decode_filename(m_outputFile, i, m_numProcs);
        unlink(decodedFileName.c_str());
      }
    }
  };

  TEST_F(NemSpreadTester, emptyTest) {}

  TEST_F(NemSpreadTester, twoHexMeshWithMaterialProperties)
  {
    if (get_parallel_size() != 1)
      GTEST_SKIP();

    std::string baseName = "twoHexNemSpread";
    unsigned    numProcs = 2;
    initialize_options(baseName, numProcs);

    unsigned    numElems = 2;
    std::string meshDesc =
        "textmesh:" + get_stacked_hex_element_textmesh_desc_with_coordinates(numElems, 1, true);
    setup_input_files_for_element_decompsition(meshDesc);

    run_nem_spread();

    verify_material_property_in_nem_spread_output();

    cleanup();
  }

  TEST_F(NemSpreadTester, twoBeamMeshWithMaterialProperties)
  {
    if (get_parallel_size() != 1)
      GTEST_SKIP();

    std::string baseName = "twoBeamNemSpread";
    unsigned    numProcs = 2;
    initialize_options(baseName, numProcs);

    unsigned    numElems = 2;
    std::string meshDesc =
        "textmesh:" + get_stacked_beam_element_textmesh_desc_with_coordinates(numElems, 1, true);
    setup_input_files_for_element_decompsition(meshDesc);

    run_nem_spread();

    verify_material_property_in_nem_spread_output();

    cleanup();
  }

  // We can't do nodal decomp testing yet
  //  [1] Nodal decomp of a mesh with elements doesn't work because Ioex_DatabaseIO fails on load:
  //  io_info exhibits this [2] Can't use textmesh to generate a pure node mesh

} // namespace
