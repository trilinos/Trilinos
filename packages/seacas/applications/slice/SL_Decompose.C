// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <SL_Decompose.h>
#include <SL_tokenize.h>
#include <random>

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_DecompositionUtils.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Region.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_Utils.h>

#include <cassert>
#include <fstream>

#include <exodusII.h>
#if !defined __NVCC__
#include <fmt/color.h>
#endif
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <init/Ionit_Initializer.h>

#if USE_METIS
#include <metis.h>
#else
using idx_t = int;
#endif

extern int    debug_level;
extern double seacas_timer();
extern void   progress(const std::string &output);

namespace {
  char **get_name_array(size_t count, int size)
  {
    auto *names = new char *[count];
    for (size_t i = 0; i < count; i++) {
      names[i] = new char[size + 1];
      std::memset(names[i], '\0', size + 1);
    }
    return names;
  }

  void delete_name_array(char **names, int count)
  {
    for (int i = 0; i < count; i++) {
      delete[] names[i];
    }
    delete[] names;
  }

  template <typename INT>
  void create_adjacency_list(const Ioss::Region &region, std::vector<idx_t> &pointer,
                             std::vector<idx_t> &adjacency, INT)
  {
    progress(__func__);
    // Size of pointer list is element count + 1;
    // Size of adjacency list is sum of nodes-per-element for each element.
    size_t      sum   = 0;
    size_t      count = 0;
    const auto &ebs   = region.get_element_blocks();
    for (const auto &eb : ebs) {
      size_t element_count = eb->entity_count();
      size_t element_nodes = eb->topology()->number_nodes();
      sum += element_count * element_nodes;
      count += element_count;
    }

    pointer.reserve(count + 1);
    adjacency.reserve(sum);
    fmt::print(stderr, "\tAdjacency Size = {} for {} elements.\n", fmt::group_digits(sum),
               fmt::group_digits(count));

    // Now, iterate the blocks again, get connectivity and build adjacency structure.
    std::vector<INT> connectivity;
    for (const auto &eb : ebs) {
      eb->get_field_data("connectivity_raw", connectivity);
      size_t element_count = eb->entity_count();
      size_t element_nodes = eb->topology()->number_nodes();

      size_t el = 0;
      for (size_t j = 0; j < element_count; j++) {
        pointer.push_back(adjacency.size());
        for (size_t k = 0; k < element_nodes; k++) {
          INT node = connectivity[el++] - 1;
          adjacency.push_back(node);
        }
      }
    }
    pointer.push_back(adjacency.size());
    assert(pointer.size() == count + 1);
    assert(adjacency.size() == sum);
  }

  void exodus_error(int lineno)
  {
    std::ostringstream errmsg;
    fmt::print(
        errmsg,
        "Exodus error ({}) {} at line {} in file Slice.C. Please report to gdsjaar@sandia.gov "
        "if you need help.",
        exerrval, ex_strerror(exerrval), lineno);

    ex_err(nullptr, nullptr, EX_PRTLASTMSG);
    throw std::runtime_error(errmsg.str());
  }

  bool case_compare(const std::string &s1, const std::string &s2)
  {
    return (s1.size() == s2.size()) &&
           std::equal(s1.begin(), s1.end(), s2.begin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  }

#if USE_METIS
  int get_common_node_count(const Ioss::Region &region)
  {
    progress(__func__);
    // Determine number of nodes that elements must share to be
    // considered connected.  A 8-node hex-only mesh would have 4
    // A 3D shell mesh should have 2.  Basically, use the minimum
    // number of nodes per side for all element blocks...  Omit sphere
    // elements; ignore bars(?)...

    int common_nodes = 999;

    const auto &ebs = region.get_element_blocks();
    for (const auto &eb : ebs) {
      const Ioss::ElementTopology *topology = eb->topology();
      const Ioss::ElementTopology *boundary = topology->boundary_type(0);
      if (boundary != nullptr) {
        common_nodes = std::min(common_nodes, boundary->number_boundaries());
      }
      else {
        // Different topologies on some element faces...
        size_t nb = topology->number_boundaries();
        for (size_t bb = 1; bb <= nb; bb++) {
          boundary = topology->boundary_type(bb);
          if (boundary != nullptr) {
            common_nodes = std::min(common_nodes, boundary->number_boundaries());
          }
        }
      }
    }

    common_nodes = std::max(1, common_nodes);
    fmt::print(stderr, "Setting common_nodes to {}\n", common_nodes);
    return common_nodes;
  }

  template <typename INT>
  void decompose_metis(const Ioss::Region &region, SystemInterface &interFace,
                       std::vector<int> &elem_to_proc, IOSS_MAYBE_UNUSED INT dummy)
  {
    size_t element_count = region.get_property("element_count").get_int();

    std::vector<idx_t> pointer;
    std::vector<idx_t> adjacency;

    double start = seacas_timer();
    create_adjacency_list(region, pointer, adjacency, dummy);
    double end = seacas_timer();
    fmt::print(stderr, "\tCreate Adjacency List = {:.5}\n", end - start);

    // Call Metis to get the partition...
    {
      start                         = seacas_timer();
      idx_t              elem_count = element_count;
      idx_t              common     = get_common_node_count(region);
      idx_t              proc_count = interFace.processor_count();
      idx_t              obj_val    = 0;
      std::vector<idx_t> options((METIS_NOPTIONS));
      METIS_SetDefaultOptions(&options[0]);
      if (interFace.decomposition_method() == "kway") {
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
      }
      else {
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
      }

      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
      if (interFace.contiguous_decomposition()) {
        options[METIS_OPTION_CONTIG] = 1;
      }
      options[METIS_OPTION_DBGLVL]  = 2;
      options[METIS_OPTION_MINCONN] = 1;

      idx_t              node_count = region.get_property("node_count").get_int();
      std::vector<idx_t> node_partition(node_count);
      std::vector<idx_t> elem_partition(element_count);

      fmt::print(stderr, "\tCalling METIS Decomposition routine.\n");

      METIS_PartMeshDual(&elem_count, &node_count, &pointer[0], &adjacency[0], nullptr, nullptr,
                         &common, &proc_count, nullptr, &options[0], &obj_val, &elem_partition[0],
                         &node_partition[0]);

      Ioss::Utils::clear(node_partition);
      elem_to_proc.reserve(element_count);
      std::copy(elem_partition.begin(), elem_partition.end(), std::back_inserter(elem_to_proc));

      end = seacas_timer();
      fmt::print(stderr, "\tMETIS Partition = {:.5}\n", end - start);
      fmt::print(stderr, "Objective value = {}\n", obj_val);

      // TODO Check Error...
    }
  }
#endif

  void scale_decomp(std::vector<int> &elem_to_proc, int iscale, size_t num_proc)
  {
    // Do the scaling (integer division...)
    if (iscale == 0) {
      // Auto scaling was asked for.  Determine max entry in `elem_to_proc` and
      // set the scale factor.
      auto max_proc = *std::max_element(elem_to_proc.begin(), elem_to_proc.end());

      iscale = (max_proc + 1) / num_proc;
      fmt::print(" Element Processor Map automatic scaling factor = {}\n", iscale);

      if (iscale == 0) {
        fmt::print(stderr,
                   "ERROR: Max value in element processor map is {} which is\n"
                   "\tless than the processor count ({}). Scaling values is not possible.",
                   max_proc, num_proc);
        exit(EXIT_FAILURE);
      }
    }
    std::transform(elem_to_proc.begin(), elem_to_proc.end(), elem_to_proc.begin(),
                   [iscale](int p) { return p / iscale; });
  }

  std::pair<int, std::string> extract_iscale_name(const std::string &var_name,
                                                  const std::string &var_type)
  {
    if (var_name.empty()) {
      fmt::print(stderr, "\nERROR: No element decomposition {} specified.\n", var_type);
      exit(EXIT_FAILURE);
    }
    // If the "var_name" string contains a comma, then the value
    // following the comma is either an integer "scale" which is
    // divided into each entry in `elem_to_proc`, or it is the
    // string "auto" which will automatically scale all values by
    // the *integer* "max/processorCount"
    //
    // NOTE: integer division with *no* rounding is used.
    int  iscale = 1;
    auto pos    = var_name.find(',');
    if (pos != std::string::npos) {
      // Extract the string following the comma...
      auto scale = var_name.substr(pos + 1);
      if (scale == "AUTO" || scale == "auto") {
        iscale = 0;
      }
      else {
        iscale = std::stoi(scale);
      }
    }
    return {iscale, var_name.substr(0, pos)};
  }
} // namespace

template std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                             const std::vector<float> &weights,
                                             IOSS_MAYBE_UNUSED int     dummy);
template std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                             const std::vector<float> &weights,
                                             IOSS_MAYBE_UNUSED int64_t dummy);

template <typename INT>
std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                    const std::vector<float> &weights, IOSS_MAYBE_UNUSED INT dummy)
{
  progress(__func__);
  // Populate the 'elem_to_proc' vector with a mapping from element to processor.

  size_t element_count = region.get_property("element_count").get_int();
  size_t elem_per_proc = element_count / interFace.processor_count();
  size_t extra         = element_count % interFace.processor_count();

  std::vector<int> elem_to_proc;
  elem_to_proc.reserve(element_count);

  fmt::print(stderr, "\nDecomposing {} elements across {} processors using method '{}'.\n",
             fmt::group_digits(element_count), fmt::group_digits(interFace.processor_count()),
             interFace.decomposition_method());
  if (interFace.lineDecomp_) {
    fmt::print(stderr, "\tDecomposition will be modified to put element lines/chains/columns on "
                       "same processor rank\n");
  }

  if (interFace.outputDecompMap_) {
    fmt::print(stderr, "\tDecomposition will be output to an element map named '{}'.\n",
               interFace.decomposition_variable());
  }
  if (interFace.outputDecompField_) {
    fmt::print(stderr, "\tDecomposition will be output to an element field named '{}'.\n",
               interFace.decomposition_variable());
  }
  fmt::print(stderr, "\n");

  if (interFace.decomposition_method() == "linear") {
    size_t elem_beg = 0;
    for (size_t proc = 0; proc < interFace.processor_count(); proc++) {
      size_t add      = (proc < extra) ? 1 : 0;
      size_t elem_end = elem_beg + elem_per_proc + add;

      for (size_t elem = elem_beg; elem < elem_end; elem++) {
        elem_to_proc.push_back(proc);
      }
      elem_beg = elem_end;
    }
  }
  else if (interFace.decomposition_method() == "scattered" ||
           interFace.decomposition_method() == "random") {
    // Scattered...
    size_t proc = 0;
    for (size_t elem = 0; elem < element_count; elem++) {
      elem_to_proc.push_back(proc++);
      if (proc >= interFace.processor_count()) {
        proc = 0;
      }
    }
    if (interFace.decomposition_method() == "random") {
      // Random...  Use scattered method and then random_shuffle() the vector.
      // Ensures that each processor has correct number of elements, but
      // they are randomly distributed.
      std::random_device rd;
      std::mt19937       g(rd());
      std::shuffle(elem_to_proc.begin(), elem_to_proc.end(), g);
    }
  }

  else if (interFace.decomposition_method() == "rcb" || interFace.decomposition_method() == "rib" ||
           interFace.decomposition_method() == "hsfc") {
    Ioss::DecompUtils::decompose_zoltan(
        region, interFace.processor_count(), interFace.decomposition_method(), elem_to_proc,
        weights, interFace.ignore_x_, interFace.ignore_y_, interFace.ignore_z_, dummy);
  }

  else if (interFace.decomposition_method() == "rb" || interFace.decomposition_method() == "kway") {
#if USE_METIS
    decompose_metis(region, interFace, elem_to_proc, dummy);
#else
    fmt::print(stderr, "ERROR: Metis library not enabled in this version of slice.\n"
                       "       The 'rb' and 'kway' methods are not available.\n\n");
    std::exit(1);
#endif
  }

  else if (interFace.decomposition_method() == "variable") {
    auto [iscale, elem_variable] =
        extract_iscale_name(interFace.decomposition_variable(), "variable");

    // Get all element blocks and cycle through each reading the
    // values for the processor...
    const auto &blocks   = region.get_element_blocks();
    auto       *c_region = (Ioss::Region *)(&region);
    c_region->begin_state(1);
    for (const auto &block : blocks) {
      if (!block->field_exists(elem_variable)) {
        fmt::print(stderr, "\nERROR: Element variable '{}' does not exist on block {}.\n",
                   elem_variable, block->name());
        exit(EXIT_FAILURE);
      }
      std::vector<double> tmp_vals;
      block->get_field_data(elem_variable, tmp_vals);
      auto block_count = block->entity_count();
      for (int64_t i = 0; i < block_count; i++) {
        elem_to_proc.push_back((int)tmp_vals[i]);
      }
    }
    scale_decomp(elem_to_proc, iscale, interFace.processor_count());
  }
  else if (interFace.decomposition_method() == "map") {
    auto [iscale, map_name] = extract_iscale_name(interFace.decomposition_variable(), "map");

    Ioss::DatabaseIO *db    = region.get_database();
    int               exoid = db->get_file_pointer();

    bool map_read  = false;
    int  map_count = ex_inquire_int(exoid, EX_INQ_ELEM_MAP);
    if (map_count > 0) {
      int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      max_name_length     = max_name_length < 32 ? 32 : max_name_length;
      char **names        = get_name_array(map_count, max_name_length);
      int    error        = ex_get_names(exoid, EX_ELEM_MAP, names);
      if (error < 0) {
        exodus_error(__LINE__);
      }

      for (int i = 0; i < map_count; i++) {
        if (case_compare(names[i], map_name)) {
          elem_to_proc.resize(element_count);
          error = ex_get_num_map(exoid, EX_ELEM_MAP, i + 1, Data(elem_to_proc));
          if (error < 0) {
            exodus_error(__LINE__);
          }
          map_read = true;
          break;
        }
      }
      delete_name_array(names, map_count);
    }

    if (!map_read) {
      fmt::print(stderr, "\nERROR: Element decomposition map '{}' could not be read from file.\n",
                 map_name);
      exit(EXIT_FAILURE);
    }

    scale_decomp(elem_to_proc, iscale, interFace.processor_count());
  }
  else if (interFace.decomposition_method() == "file") {
    // Read the element decomposition mapping from a file.  The
    // syntax of the file is an optional element count followed by
    // the processor for this range.  If the element range is
    // omitted, then the processor applies to the next element in
    // the sequence. All elements must be specified or an error will
    // be raised.
    //
    // Example:
    // 0
    // 100 1
    // 0
    //
    // Will assign:
    // * element 1 to processor 0;
    // * followed by the next 100 elements (2 to 101) to processor 1;
    // * followed by the next element (102) to processor 0.
    //
    // The resulting decomposition will have 2 elements (1, 102) on
    // processor 0 and 100 elements (2..101) on processor 1.

    const std::string &filename = interFace.decomposition_file();
    if (filename.empty()) {
      fmt::print(stderr, "\nERROR: No element decomposition file specified.\n");
      exit(EXIT_FAILURE);
    }

    std::ifstream decomp_file(filename, std::ios::in);
    if (!decomp_file.good()) {
      fmt::print(
          stderr,
          "\nERROR: Element decomposition file '{}' does not exist or could not be opened.\n",
          filename);
      exit(EXIT_FAILURE);
    }

    std::string line;
    size_t      line_num = 0;
    while (std::getline(decomp_file, line)) {
      line_num++;
      // See if 1 or 2 tokens on line...
      std::vector<std::string> tokens;
      tokens       = SLIB::tokenize(line, ", \t");
      size_t proc  = 0;
      size_t count = 1;
      if (tokens.empty()) {
        break;
      }
      if (tokens.size() == 1) {
        // Just a processor specification for the next element...
        proc = std::stoi(tokens[0]);
        elem_to_proc.push_back(proc);
      }
      else {
        // Count and processor specified.
        count = std::stoi(tokens[0]);
        proc  = std::stoi(tokens[1]);
      }
      if (proc > interFace.processor_count()) {
        fmt::print(stderr,
                   "\nERROR: Invalid processor {} specified on line {} of decomposition file.\n"
                   "\tValid range is 0..{}\n",
                   fmt::group_digits(proc), fmt::group_digits(line_num),
                   fmt::group_digits(interFace.processor_count() - 1));
        exit(EXIT_FAILURE);
      }

      if (elem_to_proc.size() + count > element_count) {
        fmt::print(stderr,
                   "\nERROR: The processor specification on line {}"
                   " of the decomposition file results in too many elements being specified.\n"
                   "\tThe total number of elements in the model is {}\n"
                   "\tPrior to this line, {} elements were specified.\n"
                   "\tIncluding this line, {} elements will be specified.\n",
                   fmt::group_digits(line_num), fmt::group_digits(element_count),
                   fmt::group_digits(elem_to_proc.size()),
                   fmt::group_digits(elem_to_proc.size() + count));
        exit(EXIT_FAILURE);
      }

      for (size_t i = 0; i < count; i++) {
        elem_to_proc.push_back(proc);
      }
    }
  }
  else {
    fmt::print(stderr, "ERROR: Unrecognized decomposition method '{}'\n\n",
               interFace.decomposition_method());
    exit(EXIT_FAILURE);
  }

  assert(elem_to_proc.size() == element_count);
  return elem_to_proc;
}
