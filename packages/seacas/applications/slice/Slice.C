// Copyright(C) 2016-2017 National Technology & Engineering Solutions of
// Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
//   copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided
//   with the distribution.
//
// * Neither the name of NTESS nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <SL_SystemInterface.h>
#include <SL_tokenize.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_FileInfo.h>
#include <Ioss_Region.h>
#include <Ioss_SubSystem.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Utils.h>
#include <cassert>
#include <exo_fpp/Iofx_DatabaseIO.h>
#include <init/Ionit_Initializer.h>

#include <exodusII.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#if USE_METIS
#include <metis.h>
#else
using idx_t = int;
#endif

#include <sys/stat.h>
#include <sys/types.h>

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

#define OUTPUT std::cerr

// ========================================================================
// TODO(gdsjaar):
//  * Sideset distribution factors
//  * Variables
//  * All entity types
//  * More efficient border-node-processor communication map.
// ========================================================================

extern double seacas_timer();
bool          minimize_open_files = false;
int           debug_level         = 0;

// size_t partial_count = 100000;
size_t partial_count = 1000000000;

namespace {
  void progress(const std::string &output)
  {
    static auto start = std::chrono::high_resolution_clock::now();

    if ((debug_level & 1) != 0) {
      auto                          now  = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = now - start;
      std::cerr << " [" << std::fixed << std::setprecision(2) << diff.count() << " - "
                << Ioss::Utils::get_memory_info() << "]\t" << output << "\n";
    }
  }

  void proc_progress(int p, int proc_count)
  {
    if (((debug_level & 4) != 0) && ((p + 1) % (proc_count / 20) == 0)) {
      progress("\t\tProcessor " + std::to_string(p + 1));
    }
  }

  void filename_substitution(std::string &filename, const SystemInterface &interface);

  template <typename T> struct remove_pointer
  {
    using type = T;
  };
  template <typename T> struct remove_pointer<T *>
  {
    using type = T;
  };

  template <typename INT>
  void slice(Ioss::Region &region, const std::string &nemfile, SystemInterface &interface,
             INT /*dummy*/);

  int get_free_descriptor_count();

  template <typename INT> bool is_sequential(const std::vector<INT> &map)
  {
    progress(__func__);
    for (size_t i = 0; i < map.size(); i++) {
      if (map[i] != i + 1) {
        return false;
      }
    }
    return true;
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

    auto &ebs = region.get_element_blocks();
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
    OUTPUT << "Setting common_nodes to " << common_nodes << "\n";
    return common_nodes;
  }
#endif

  struct my_numpunct : std::numpunct<char>
  {
  protected:
    char        do_thousands_sep() const override { return ','; }
    std::string do_grouping() const override { return "\3"; }
  };
} // namespace
// ========================================================================

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  std::cout.imbue(std::locale(std::locale(), new my_numpunct));
  std::cerr.imbue(std::locale(std::locale(), new my_numpunct));

  double begin = seacas_timer();

  Ioss::Init::Initializer io;
  SystemInterface::show_version();

  SystemInterface interface;
  bool            ok = interface.parse_options(argc, argv);
  if (!ok) {
    OUTPUT << "\nERROR: Problem parsing command line options.\n\n";
    exit(EXIT_FAILURE);
  }

  std::string nem_file = interface.nemesisFile_;
  std::string path     = interface.output_path();
  if (!path.empty()) {
    filename_substitution(path, interface);

    // See if specified path exists.
    Ioss::FileInfo output_path(path);
    if (!output_path.exists()) {
      // Try to create the directory...
      if (mkdir(path.c_str(), 0777) == -1) {
        OUTPUT << "ERROR: Could not create path '" << path << "' for output of decomposed files.\n";
        exit(EXIT_FAILURE);
      }
    }
    else if (!output_path.is_dir()) {
      OUTPUT << "ERROR: Path '" << path << "' is not a directory.\n";
      exit(EXIT_FAILURE);
    }

    // See if the nem_file already has a path prepended to the
    // filename and if so, extract the basename.
    Ioss::FileInfo nemesis(nem_file);
    std::string    sep = "/";
    if (path[path.length() - 1] == '/') {
      sep = "";
    }
    nem_file = path + sep + nemesis.tailname();
  }

  OUTPUT << "Input:    '" << interface.inputFile_ << '\n';
  OUTPUT << "Output:   '" << nem_file << '\n';
  OUTPUT << '\n';

  // Check whether processor count is larger than maximum number of open files...
  size_t max_files    = get_free_descriptor_count();
  minimize_open_files = (interface.processor_count() + 1 > max_files);

  debug_level   = interface.debug();
  partial_count = interface.partial();

  //========================================================================
  // INPUT ...
  // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
  //========================================================================
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(interface.inputFormat_, interface.inputFile_,
                                                  Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD);
  if (dbi == nullptr || !dbi->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  //  dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
  dbi->set_surface_split_type(Ioss::SPLIT_BY_DONT_SPLIT);
  dbi->set_field_separator(0);

  // NOTE: 'region' owns 'db' pointer at this time...
  Ioss::Region region(dbi, "region_1");

  if (dbi->int_byte_size_api() == 4) {
    progress("4-byte slice");
    slice(region, nem_file, interface, 1);
  }
  else {
    progress("8-byte slice");
    slice(region, nem_file, interface, static_cast<int64_t>(1));
  }

#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif
  OUTPUT << "Total execution time = " << seacas_timer() - begin << "\n";
  OUTPUT << "\nSlice execution successful.\n";
  return EXIT_SUCCESS;
}

namespace {

  template <typename INT>
  void create_adjacency_list(const Ioss::Region &region, std::vector<idx_t> &pointer,
                             std::vector<idx_t> &adjacency, INT /*dummy*/)
  {
    progress(__func__);
    // Size of pointer list is element count + 1;
    // Size of adjacency list is sum of nodes-per-element for each element.
    size_t sum   = 0;
    size_t count = 0;
    auto & ebs   = region.get_element_blocks();
    for (const auto &eb : ebs) {
      size_t element_count = eb->entity_count();
      size_t element_nodes = eb->get_property("topology_node_count").get_int();
      sum += element_count * element_nodes;
      count += element_count;
    }

    pointer.reserve(count + 1);
    adjacency.reserve(sum);
    OUTPUT << "\tAdjacency Size = " << sum << " for " << count << " elements.\n";

    // Now, iterate the blocks again, get connectivity and build adjacency structure.
    std::vector<INT> connectivity;
    for (const auto &eb : ebs) {
      eb->get_field_data("connectivity_raw", connectivity);
      size_t element_count = eb->entity_count();
      size_t element_nodes = eb->get_property("topology_node_count").get_int();

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

  void decompose_elements(const Ioss::Region &region, SystemInterface &interface,
                          std::vector<int> &elem_to_proc)
  {
    progress(__func__);
    // Populate the 'elem_to_proc' vector with a mapping from element to processor.

    size_t element_count = region.get_property("element_count").get_int();
    size_t elem_per_proc = element_count / interface.processor_count();
    size_t extra         = element_count % interface.processor_count();

    elem_to_proc.reserve(element_count);

    OUTPUT << "Decomposing " << element_count << " elements across " << interface.processor_count()
           << " processors using method '" << interface.decomposition_method() << "'.\n\n";

    if (interface.decomposition_method() == "linear") {
      size_t elem_beg = 0;
      for (size_t proc = 0; proc < interface.processor_count(); proc++) {
        size_t add      = (proc < extra) ? 1 : 0;
        size_t elem_end = elem_beg + elem_per_proc + add;

        for (size_t elem = elem_beg; elem < elem_end; elem++) {
          elem_to_proc.push_back(proc);
        }
        elem_beg = elem_end;
      }
    }
    else if (interface.decomposition_method() == "scattered") {
      // Scattered...
      size_t proc = 0;
      for (size_t elem = 0; elem < element_count; elem++) {
        elem_to_proc.push_back(proc++);
        if (proc >= interface.processor_count()) {
          proc = 0;
        }
      }
    }

    else if (interface.decomposition_method() == "rb" ||
             interface.decomposition_method() == "kway") {
#if USE_METIS
      std::vector<idx_t> pointer;
      std::vector<idx_t> adjacency;

      double start = seacas_timer();
      create_adjacency_list(region, interface, pointer, adjacency, dummy);
      double end = seacas_timer();
      OUTPUT << "\tCreate Adjacency List = " << end - start << "\n";

      // Call Metis to get the partition...
      {
        start                         = seacas_timer();
        idx_t              elem_count = element_count;
        idx_t              common     = get_common_node_count(region);
        idx_t              proc_count = interface.processor_count();
        idx_t              obj_val    = 0;
        std::vector<idx_t> options(METIS_NOPTIONS);
        METIS_SetDefaultOptions(&options[0]);
        if (interface.decomposition_method() == "kway") {
          options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        }
        else {
          options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
        }

        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
        if (interface.contiguous_decomposition()) {
          options[METIS_OPTION_CONTIG] = 1;
        }
        options[METIS_OPTION_DBGLVL]  = 2;
        options[METIS_OPTION_MINCONN] = 1;

        idx_t              node_count = region.get_property("node_count").get_int();
        std::vector<idx_t> node_partition(node_count);
        std::vector<idx_t> elem_partition(element_count);

        OUTPUT << "\tCalling METIS Decomposition routine.\n";

        METIS_PartMeshDual(&elem_count, &node_count, &pointer[0], &adjacency[0], nullptr, nullptr,
                           &common, &proc_count, nullptr, &options[0], &obj_val, &elem_partition[0],
                           &node_partition[0]);

        Ioss::Utils::clear(node_partition);
        elem_to_proc.reserve(element_count);
        std::copy(elem_partition.begin(), elem_partition.end(), std::back_inserter(elem_to_proc));

        end = seacas_timer();
        OUTPUT << "\tMETIS Partition = " << end - start << "\n";
        OUTPUT << "Objective value = " << obj_val << "\n";

        // TODO Check Error...
      }
#else
      OUTPUT << "ERROR: Metis library not enabled in this version of slice.\n"
             << "       The 'rb' and 'kway' methods are not available.\n\n";
      std::exit(1);
#endif
    }

    else if (interface.decomposition_method() == "random") {
      // Random...  Use linear method and then random_shuffle() the vector.
      // Ensures that each processor has correct number of elements, but
      // they are randomly distributed.
      size_t elem_beg = 0;
      for (size_t proc = 0; proc < interface.processor_count(); proc++) {
        size_t add      = (proc < extra) ? 1 : 0;
        size_t elem_end = elem_beg + elem_per_proc + add;

        for (size_t elem = elem_beg; elem < elem_end; elem++) {
          elem_to_proc.push_back(proc);
        }
        elem_beg = elem_end;
      }
      std::random_device rd;
      std::mt19937       g(rd());
      std::shuffle(elem_to_proc.begin(), elem_to_proc.end(), g);
    }

    else if (interface.decomposition_method() == "file") {
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
      // Will assign element 1 to processor 0; followed by the next
      // 100 elements to processor 1; followed by the next element
      // (102) to processor 0.  The resulting decomposition will have
      // 2 elements (1, 102) on processor 0 and 100 elements (2..101)
      // on processor 1.

      const std::string &filename = interface.decomposition_file();
      if (filename.empty()) {
        OUTPUT << "\nERROR: No element decomposition file specified.\n";
        exit(EXIT_FAILURE);
      }

      std::ifstream decomp_file(filename, std::ios::in);
      if (!decomp_file.good()) {
        OUTPUT << "\nERROR: Element decomposition file '" << filename
               << "' does not exist or could not be opened.\n";
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
        else if (tokens.size() == 1) {
          // Just a processor specification for the next element...
          proc = std::stoi(tokens[0]);
          elem_to_proc.push_back(proc);
        }
        else {
          // Count and processor specified.
          count = std::stoi(tokens[0]);
          proc  = std::stoi(tokens[1]);
        }
        if (proc > interface.processor_count()) {
          OUTPUT << "\nERROR: Invalid processor " << proc << " specified on line " << line_num
                 << " of decomposition file.\n"
                 << "\tValid range is 0.." << interface.processor_count() - 1 << "\n";
          exit(EXIT_FAILURE);
        }

        if (elem_to_proc.size() + count > element_count) {
          OUTPUT << "\nERROR: The processor specification on line " << line_num
                 << " of the decomposition file results in too many elements being specified.\n"
                 << "\tThe total number of elements in the model is " << element_count << "\n"
                 << "\tPrior to this line, " << elem_to_proc.size() << " elements were specified.\n"
                 << "\tIncluding this line, " << elem_to_proc.size() + count
                 << " elements will be specified.\n";
          exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < count; i++) {
          elem_to_proc.push_back(proc);
        }
      }
    }
    assert(elem_to_proc.size() == element_count);
  }

  template <typename INT>
  void free_connectivity_storage(std::vector<std::vector<std::vector<INT>>> &connectivity)
  {
    progress(__func__);
    size_t processor_count = connectivity.size();
    for (size_t p = 0; p < processor_count; p++) {
      size_t block_count = connectivity[p].size();
      for (size_t b = 0; b < block_count; b++) {
        Ioss::Utils::clear(connectivity[p][b]);
      }
      Ioss::Utils::clear(connectivity[p]);
    }
    Ioss::Utils::clear(connectivity);
  }

  template <typename INT>
  void get_sidesets(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                    const std::vector<int> &elem_to_proc, INT /*dummy*/)
  {
    progress(__func__);
    // This routine reads the sidesets in the global database;
    // and defines corresponding sidesets on each processor...
    size_t proc_count = proc_region.size();

    auto & ss        = region.get_sidesets();
    size_t set_count = ss.size();

    for (size_t s = 0; s < set_count; s++) {
      auto *gss     = ss[s];
      auto &ss_name = gss->name();

      std::vector<Ioss::SideSet *> sset(proc_count);
      for (size_t p = 0; p < proc_count; p++) {
        sset[p] = new Ioss::SideSet(proc_region[p]->get_database(), ss_name);
        proc_region[p]->add(sset[p]);
      }

      auto &side_blocks = gss->get_side_blocks();
      for (auto gsb : side_blocks) {
        std::vector<INT> ss_elems;
        gsb->get_field_data("element_side_raw", ss_elems);

        std::vector<INT> pss(proc_count);
        for (size_t i = 0; i < ss_elems.size(); i += 2 /* elem,side pairs */) {
          int64_t elem = ss_elems[i] - 1;
          int     p    = elem_to_proc[elem];
          pss[p]++;
        }

        auto &name      = gsb->name();
        auto &side_type = gsb->topology()->name();
        auto &elem_type = gsb->parent_element_topology()->name();

        for (size_t p = 0; p < proc_count; p++) {
          auto *side_block = new Ioss::SideBlock(proc_region[p]->get_database(), name, side_type,
                                                 elem_type, pss[p]);
          sset[p]->add(side_block);
        }
      }
    }
  }

  template <typename INT>
  void output_sidesets(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                       const std::vector<int> &elem_to_proc, INT /*dummy*/)
  {
    progress(__func__);
    // This routine reads the sidesets in the global database;
    // and outputs the sidesets on each processor...
    size_t proc_count = proc_region.size();

    auto & ss        = region.get_sidesets();
    size_t set_count = ss.size();

    for (size_t s = 0; s < set_count; s++) {
      if (debug_level & 4) {
        progress("\tSideset " + std::to_string(s + 1));
      }
      Ioss::SideSet *gss     = ss[s];
      auto &         ss_name = gss->name();

      std::vector<Ioss::SideSet *> proc_ss(proc_count);
      for (size_t p = 0; p < proc_count; p++) {
        proc_ss[p] = proc_region[p]->get_sideset(ss_name);
      }

      auto &side_blocks = gss->get_side_blocks();
      for (auto gsb : side_blocks) {
        auto &sb_name = gsb->name();

        std::vector<Ioss::SideBlock *> proc_sb(proc_count);
        std::vector<std::vector<INT>>  psb_elems(proc_count);
        std::vector<std::vector<INT>>  psb_sides(proc_count);
        for (size_t p = 0; p < proc_count; p++) {
          proc_sb[p]        = proc_ss[p]->get_side_block(sb_name);
          size_t elem_count = proc_sb[p]->entity_count();
          psb_elems[p].reserve(elem_count * 2);
        }

        std::vector<INT> ss_elems;
        gsb->get_field_data("element_side_raw", ss_elems);

        std::vector<INT> pss(proc_count);
        for (size_t i = 0; i < ss_elems.size(); i += 2 /* elem,side pairs */) {
          int64_t elem = ss_elems[i] - 1;
          int     p    = elem_to_proc[elem];
          psb_elems[p].push_back(elem + 1);
          psb_elems[p].push_back(ss_elems[i + 1]);
        }

        for (size_t p = 0; p < proc_count; p++) {
          Ioss::SideBlock *psb = proc_sb[p];
          psb->put_field_data("element_side", psb_elems[p]);
          if (minimize_open_files) {
            proc_region[p]->get_database()->closeDatabase();
          }
          proc_progress(p, proc_count);
        }
      }
    }
    if (set_count > 0) {
      OUTPUT << "WARNING: Sideset distribution factors not yet handled correctly.\n";
    }
  }

  template <typename INT>
  void output_communication_map(std::vector<Ioss::Region *> &  proc_region,
                                std::vector<std::vector<INT>> &border_node_proc_map)
  {
    progress(__func__);
    size_t proc_count = proc_region.size();
    for (size_t p = 0; p < proc_count; p++) {
      auto &commset = proc_region[p]->get_commsets()[0];
      commset->put_field_data("entity_processor", border_node_proc_map[p]);
      Ioss::Utils::clear(border_node_proc_map[p]);
      if (minimize_open_files) {
        proc_region[p]->get_database()->closeDatabase();
      }
      proc_progress(p, proc_count);
    }
    Ioss::Utils::clear(border_node_proc_map);
  }

  template <typename INT>
  void define_communication_data(const Ioss::Region &           global_region,
                                 std::vector<Ioss::Region *> &  proc_region,
                                 const std::vector<int> &       node_to_proc,
                                 const std::vector<INT> &       node_to_proc_pointer,
                                 std::vector<std::vector<INT>> &border_node_proc_map)
  {
    progress(__func__);
    // This routine categorizes the nodes on a processor as interior
    // or border.
    // TODO(gdsjaar): Categorize elements also. For now, all treated as
    // interior which works for sierra-based applications

    // The node_to_proc_pointer has information about the number of
    // processors that a node is shared with.  If the count is 1, then
    // the node is interior; otherwise, it is border.

    // Allocates:
    // * interior_nodes         INT size - #interior nodes
    // * border_nodes_proc_map  INT size - (proc-node) pair for each border node

    INT              global_node_count = global_region.get_property("node_count").get_int();
    size_t           proc_count        = proc_region.size();
    std::vector<INT> interior_nodes(proc_count);

    // Iterate all nodes and count the number of processors it is on:
    for (INT i = 0; i < global_node_count; i++) {
      size_t node_proc_count = node_to_proc_pointer[i + 1] - node_to_proc_pointer[i];
      if (node_proc_count == 1) {
        size_t proc = node_to_proc[node_to_proc_pointer[i]];
        interior_nodes[proc]++;
      }
      else {
        // Get the <node,proc> pairs for all border nodes on this processor...
        // Not efficient at this time...
        size_t beg = node_to_proc_pointer[i];
        size_t end = node_to_proc_pointer[i + 1];
        for (size_t j = beg; j < end; j++) {
          size_t node = i + 1;
          size_t proc = node_to_proc[j];
          for (size_t k = beg; k < end; k++) {
            if (j == k) {
              continue;
            }
            size_t p = node_to_proc[k];
            border_node_proc_map[p].push_back(node);
            border_node_proc_map[p].push_back(proc);
          }
        }
      }
    }

    INT global_element_count = global_region.get_property("element_count").get_int();

    // Categorize each element as interior...
    // Categorize the remaining nodes as border...
    for (size_t p = 0; p < proc_count; p++) {
      Ioss::Region *region = proc_region[p];

      INT element_count = region->get_property("element_count").get_int();
      INT node_count    = region->get_property("node_count").get_int();
      INT border_nodes  = node_count - interior_nodes[p];

      region->property_add(Ioss::Property("global_node_count", global_node_count));
      region->property_add(Ioss::Property("global_element_count", global_element_count));
      region->property_add(Ioss::Property("processor_count", static_cast<int>(proc_count)));
      region->property_add(Ioss::Property("my_processor", static_cast<int>(p)));

      region->property_add(Ioss::Property("internal_node_count", interior_nodes[p]));
      region->property_add(Ioss::Property("border_node_count", border_nodes));
      region->property_add(Ioss::Property("internal_element_count", element_count));
      region->property_add(Ioss::Property("border_element_count", 0));

      // Add commset data... The length of the commset is the number
      // of <node,proc> pairs for all border nodes.
      //
      // For each node on this processor that isn't an interior node,
      // create the <node,proc> pair...
      auto *commset = new Ioss::CommSet(region->get_database(), "commset_node", "node",
                                        border_node_proc_map[p].size() / 2);
      commset->property_add(Ioss::Property("id", 1));
      region->add(commset);
      if (debug_level & 2) {
        OUTPUT << "Commset for processor " << p << " has " << border_node_proc_map[p].size() / 2
               << " entries.\n";
      }
    }
  }

  template <typename INT>
  void get_nodesets(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                    const std::vector<int> &node_to_proc,
                    const std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    // This routine reads the nodesets in the global database;
    // and defines corresponding nodesets on each processor...
    size_t proc_count = proc_region.size();

    auto & ns        = region.get_nodesets();
    size_t set_count = ns.size();

    for (size_t s = 0; s < set_count; s++) {
      std::vector<INT> pns(proc_count);
      Ioss::NodeSet *  gns = ns[s];

      std::vector<INT> ns_nodes;
      gns->get_field_data("ids_raw", ns_nodes);

      for (size_t i = 0; i < ns_nodes.size(); i++) {
        int64_t node  = ns_nodes[i] - 1;
        size_t  p_beg = node_to_proc_pointer[node];
        size_t  p_end = node_to_proc_pointer[node + 1];
        for (size_t j = p_beg; j < p_end; j++) {
          size_t p = node_to_proc[j];
          pns[p]++;
        }
      }

      auto &name = ns[s]->name();
      if (debug_level & 2) {
        OUTPUT << "\tNodeset " << name << "-- ";
      }
      for (size_t p = 0; p < proc_count; p++) {
        auto *node_set = new Ioss::NodeSet(proc_region[p]->get_database(), name, pns[p]);
        proc_region[p]->add(node_set);
        if (debug_level & 2) {
          OUTPUT << p << ":" << pns[p] << ", ";
        }
      }
      if (debug_level & 2) {
        OUTPUT << "\n";
      }
    }
  }

  template <typename INT>
  void output_nodesets(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                       const std::vector<int> &node_to_proc,
                       const std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    // This routine reads the nodesets in the global database;
    // and defines corresponding nodesets on each processor...
    size_t proc_count = proc_region.size();

    auto & ns        = region.get_nodesets();
    size_t set_count = ns.size();

    for (size_t s = 0; s < set_count; s++) {
      if (debug_level & 4) {
        progress("\tNodeSet " + std::to_string(s + 1));
      }
      std::vector<INT> pns(proc_count);
      Ioss::NodeSet *  gns = ns[s];

      std::vector<INT> ns_nodes;
      gns->get_field_data("ids_raw", ns_nodes);

      std::vector<double> ns_df;
      gns->get_field_data("distribution_factors", ns_df);

      std::vector<std::vector<INT>>    pns_nodes(proc_count);
      std::vector<std::vector<double>> pns_df(proc_count);
      for (size_t p = 0; p < proc_count; p++) {
        size_t node_count = proc_region[p]->get_nodesets()[s]->entity_count();
        pns_nodes[p].reserve(node_count);
        pns_df[p].reserve(node_count);
      }

      for (size_t i = 0; i < ns_nodes.size(); i++) {
        int64_t node  = ns_nodes[i] - 1;
        size_t  p_beg = node_to_proc_pointer[node];
        size_t  p_end = node_to_proc_pointer[node + 1];
        for (size_t j = p_beg; j < p_end; j++) {
          size_t p = node_to_proc[j];
          pns_nodes[p].push_back(node + 1);
          pns_df[p].push_back(ns_df[i]);
        }
      }

      for (size_t p = 0; p < proc_count; p++) {
        Ioss::NodeSet *proc_ns = proc_region[p]->get_nodesets()[s];
        proc_ns->put_field_data("ids", pns_nodes[p]);
        proc_ns->put_field_data("distribution_factors", pns_df[p]);
        if (minimize_open_files) {
          proc_region[p]->get_database()->closeDatabase();
        }
        proc_progress(p, proc_count);
      }
    }
  }

  template <typename INT>
  void output_node_map(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                       const std::vector<int> &node_to_proc,
                       const std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    // This is the processor-local to global-implicit node map...
    // This maps the 1..#node in the global mesh to each processor...
    size_t node_count = region.get_property("node_count").get_int();
    size_t proc_count = proc_region.size();

    std::vector<std::vector<INT>> proc_map(proc_count);
    for (size_t p = 0; p < proc_count; p++) {
      size_t pnode_count = proc_region[p]->get_property("node_count").get_int();
      proc_map[p].reserve(pnode_count);
    }

    for (size_t i = 0; i < node_count; i++) {
      size_t p_beg = node_to_proc_pointer[i];
      size_t p_end = node_to_proc_pointer[i + 1];
      for (size_t j = p_beg; j < p_end; j++) {
        size_t p = node_to_proc[j];
        proc_map[p].push_back(i + 1);
      }
    }

    for (size_t p = 0; p < proc_count; p++) {
      Ioss::NodeBlock *nb = proc_region[p]->get_node_blocks()[0];
      nb->put_field_data("ids", proc_map[p]);
      if (minimize_open_files) {
        proc_region[p]->get_database()->closeDatabase();
      }
      proc_progress(p, proc_count);
    }
  }

  template <typename INT>
  void output_global_node_map(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                              const std::vector<int> &node_to_proc,
                              const std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    // This is the processor-local to global-implicit node map...
    // This maps the node_number map (if it exists) in the global mesh
    // to each processor...
    std::vector<INT> ids;
    Ioss::NodeBlock *gnb = region.get_node_blocks()[0];
    gnb->get_field_data("ids", ids);

    // Check whether the map is sequential (X maps to X);
    bool sequential = is_sequential(ids);
    if (!sequential) {
      OUTPUT << "Node map is not sequential...\n";
    }

    size_t node_count = region.get_property("node_count").get_int();
    size_t proc_count = proc_region.size();

    std::vector<std::vector<INT>> proc_map(proc_count);
    for (size_t p = 0; p < proc_count; p++) {
      size_t pnode_count = proc_region[p]->get_property("node_count").get_int();
      proc_map[p].reserve(pnode_count);
    }

    for (size_t i = 0; i < node_count; i++) {
      size_t p_beg = node_to_proc_pointer[i];
      size_t p_end = node_to_proc_pointer[i + 1];
      for (size_t j = p_beg; j < p_end; j++) {
        size_t p = node_to_proc[j];
        proc_map[p].push_back(ids[i]);
      }
    }

    for (size_t p = 0; p < proc_count; p++) {
      Ioss::NodeBlock *nb = proc_region[p]->get_node_blocks()[0];
      nb->put_field_data("ids", proc_map[p]);
      if (minimize_open_files) {
        proc_region[p]->get_database()->closeDatabase();
      }
      proc_progress(p, proc_count);
    }
  }

  template <typename INT>
  void output_element_map(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                          const std::vector<int> &elem_to_proc, INT /* dummy */)
  {
    progress(__func__);
    // map[p][b] = map for block b on processor p
    size_t proc_count = proc_region.size();

    auto & ebs         = region.get_element_blocks();
    size_t block_count = ebs.size();

    size_t offset = 0;
    for (size_t b = 0; b < block_count; b++) {
      if (debug_level & 4) {
        progress("\tBlock " + std::to_string(b + 1));
      }
#if 0
      std::vector<INT> ids;
      ebs[b]->get_field_data("ids", ids);
#endif

      std::vector<std::vector<INT>> map(proc_count);
      for (size_t p = 0; p < proc_count; p++) {
        auto & proc_ebs           = proc_region[p]->get_element_blocks();
        size_t proc_element_count = proc_ebs[b]->entity_count();
        map[p].reserve(proc_element_count);
      }

      size_t element_count = ebs[b]->entity_count();

      for (size_t j = 0; j < element_count; j++) {
        size_t p = elem_to_proc[offset + j];
#if 0	
	map[p].push_back(ids[j]);
#else
        map[p].push_back(offset + j + 1);
#endif
      }
      offset += element_count;

      for (size_t p = 0; p < proc_count; p++) {
        auto &proc_ebs = proc_region[p]->get_element_blocks();
        proc_ebs[b]->put_field_data("ids", map[p]);
        if (minimize_open_files) {
          proc_region[p]->get_database()->closeDatabase();
        }
        proc_progress(p, proc_count);
      }
    }
  }

  template <typename INT>
  void output_connectivity(std::vector<Ioss::Region *> &               proc_region,
                           std::vector<std::vector<std::vector<INT>>> &connectivity)
  {
    progress(__func__);
    // connectivity[p][b] = map for block b on processor p
    size_t proc_count = proc_region.size();

    for (size_t p = 0; p < proc_count; p++) {
      auto & proc_ebs    = proc_region[p]->get_element_blocks();
      size_t block_count = proc_ebs.size();
      for (size_t b = 0; b < block_count; b++) {
        Ioss::ElementBlock *eb = proc_ebs[b];
        eb->put_field_data("connectivity", connectivity[p][b]);
        if (minimize_open_files) {
          proc_region[p]->get_database()->closeDatabase();
        }
      }
      proc_progress(p, proc_count);
    }
  }

  template <typename INT>
  void output_coordinates(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                          const std::vector<int> &node_to_proc,
                          const std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    std::vector<double> glob_coord_x;
    std::vector<double> glob_coord_y;
    std::vector<double> glob_coord_z;
    Ioss::NodeBlock *   gnb = region.get_node_blocks()[0];

    // Distribute nodal coordinates to each processor...
    // coordinates[p][i] = x,y,z coordinates on processor p
    size_t                           processor_count = proc_region.size();
    std::vector<std::vector<double>> coordinates_x(processor_count);
    std::vector<std::vector<double>> coordinates_y(processor_count);
    std::vector<std::vector<double>> coordinates_z(processor_count);
    for (size_t p = 0; p < processor_count; p++) {
      size_t pnode_count = proc_region[p]->get_property("node_count").get_int();
      coordinates_x[p].reserve(pnode_count);
      coordinates_y[p].reserve(pnode_count);
      coordinates_z[p].reserve(pnode_count);
    }
    progress("\tReserve processor coordinate vectors");

    Ioss::DatabaseIO *db    = region.get_database();
    Iofx::DatabaseIO *ex_db = dynamic_cast<Iofx::DatabaseIO *>(db);

    size_t node_count = region.get_property("node_count").get_int();

    if (ex_db != nullptr && node_count > partial_count) {
      int exoid = ex_db->get_file_pointer();

      glob_coord_x.resize(partial_count);
      glob_coord_y.resize(partial_count);
      glob_coord_z.resize(partial_count);
      for (size_t beg = 1; beg <= node_count; beg += partial_count) {
        size_t count = partial_count;
        if (beg + count - 1 > node_count) {
          count = node_count - beg + 1;
        }

        ex_get_partial_coord(exoid, beg, count, glob_coord_x.data(), glob_coord_y.data(),
                             glob_coord_z.data());
        progress("\tpartial_coord: " + std::to_string(beg) + " " + std::to_string(count));

        for (size_t i = 0; i < count; i++) {
          size_t ii    = beg + i - 1;
          size_t p_beg = node_to_proc_pointer[ii];
          size_t p_end = node_to_proc_pointer[ii + 1];
          for (size_t j = p_beg; j < p_end; j++) {
            size_t p = node_to_proc[j];
            coordinates_x[p].push_back(glob_coord_x[i]);
            coordinates_y[p].push_back(glob_coord_y[i]);
            coordinates_z[p].push_back(glob_coord_z[i]);
          }
        }
      }
    }
    else {
      gnb->get_field_data("mesh_model_coordinates_x", glob_coord_x);
      gnb->get_field_data("mesh_model_coordinates_y", glob_coord_y);
      gnb->get_field_data("mesh_model_coordinates_z", glob_coord_z);
      progress("\tRead global mesh_model_coordinates");

      for (size_t i = 0; i < node_count; i++) {
        size_t p_beg = node_to_proc_pointer[i];
        size_t p_end = node_to_proc_pointer[i + 1];
        for (size_t j = p_beg; j < p_end; j++) {
          size_t p = node_to_proc[j];
          coordinates_x[p].push_back(glob_coord_x[i]);
          coordinates_y[p].push_back(glob_coord_y[i]);
          coordinates_z[p].push_back(glob_coord_z[i]);
        }
      }
    }
    progress("\tPopulate processor coordinate vectors");
    Ioss::Utils::clear(glob_coord_x);
    Ioss::Utils::clear(glob_coord_y);
    Ioss::Utils::clear(glob_coord_z);

    for (size_t p = 0; p < processor_count; p++) {
      Ioss::NodeBlock *nb = proc_region[p]->get_node_blocks()[0];
      nb->put_field_data("mesh_model_coordinates_x", coordinates_x[p]);
      nb->put_field_data("mesh_model_coordinates_y", coordinates_y[p]);
      nb->put_field_data("mesh_model_coordinates_z", coordinates_z[p]);
      if (minimize_open_files) {
        proc_region[p]->get_database()->closeDatabase();
      }
      proc_progress(p, processor_count);
    }
    progress("\tOutput processor coordinate vectors");
  }

  template <typename INT>
  void get_connectivity(const Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                        const std::vector<int> &                    elem_to_proc,
                        std::vector<std::vector<std::vector<INT>>> &connectivity)
  {
    progress(__func__);
    auto & ebs         = region.get_element_blocks();
    size_t block_count = ebs.size();

    size_t processor_count = proc_region.size();
    connectivity.resize(processor_count);

    for (size_t p = 0; p < processor_count; p++) {
      connectivity[p].resize(block_count);
      const auto &pebs = proc_region[p]->get_element_blocks();
      for (size_t b = 0; b < block_count; b++) {
        size_t element_count = pebs[b]->entity_count();
        size_t element_nodes = pebs[b]->get_property("topology_node_count").get_int();
        connectivity[p][b].reserve(element_count * element_nodes); // Use reserve, not resize
      }
    }
    progress("\tAfter memory allocation");

    Ioss::DatabaseIO *db    = region.get_database();
    Iofx::DatabaseIO *ex_db = dynamic_cast<Iofx::DatabaseIO *>(db);

    std::vector<INT> glob_conn;
    size_t           offset = 0;

    for (size_t b = 0; b < block_count; b++) {
      size_t element_count = ebs[b]->entity_count();
      size_t element_nodes = ebs[b]->get_property("topology_node_count").get_int();
      size_t block_id      = ebs[b]->get_property("id").get_int();

      // Do a 'partial_count' elements at a time...
      if (ex_db != nullptr && element_count >= partial_count) {
        int exoid = ex_db->get_file_pointer();

        glob_conn.resize(partial_count * element_nodes);
        for (size_t beg = 1; beg <= element_count; beg += partial_count) {
          size_t count = partial_count;
          if (beg + count - 1 > element_count) {
            count = element_count - beg + 1;
          }

          ex_get_partial_conn(exoid, EX_ELEM_BLOCK, block_id, beg, count, glob_conn.data(), nullptr,
                              nullptr);
          progress("\tpartial_conn: " + std::to_string(beg) + " " + std::to_string(count));

          size_t el = 0;
          for (size_t j = 0; j < count; j++) {
            size_t p = elem_to_proc[offset + j];
            for (size_t k = 0; k < element_nodes; k++) {
              connectivity[p][b].push_back(glob_conn[el++]);
            }
          }
          offset += count;
        }
      }
      else {
        ebs[b]->get_field_data("connectivity_raw", glob_conn);

        size_t el = 0;
        for (size_t j = 0; j < element_count; j++) {
          size_t p = elem_to_proc[offset + j];
          for (size_t k = 0; k < element_nodes; k++) {
            connectivity[p][b].push_back(glob_conn[el++]);
          }
        }
        offset += element_count;
      }
    }
  }

  template <typename INT>
  void get_proc_elem_block_count(const Ioss::Region &region, std::vector<int> &elem_to_proc,
                                 std::vector<std::vector<INT>> &proc_elem_block_cnt)
  {
    progress(__func__);
    auto & ebs         = region.get_element_blocks();
    size_t block_count = ebs.size();
    size_t begin       = 0;
    for (size_t i = 0; i < block_count; i++) {
      size_t end = begin + ebs[i]->entity_count();
      for (size_t j = begin; j < end; j++) {
        size_t processor = elem_to_proc[j];
        proc_elem_block_cnt[i][processor]++;
      }
      begin = end;
    }

    size_t processor_count = proc_elem_block_cnt[0].size();
    for (size_t i = 0; i < processor_count; i++) {
      size_t sum = 0;
      for (size_t j = 0; j < block_count; j++) {
        sum += proc_elem_block_cnt[j][i];
      }
      proc_elem_block_cnt[block_count][i] = sum;
      if (debug_level & 2) {
        OUTPUT << "\tProcessor " << i << " has " << sum << " elements.\n";
      }
    }
  }

  template <typename INT>
  void get_node_to_proc(Ioss::Region &region, std::vector<Ioss::Region *> &proc_region,
                        std::vector<std::vector<std::vector<INT>>> &connectivity,
                        std::vector<int> &node_to_proc, std::vector<INT> &node_to_proc_pointer)
  {
    progress(__func__);
    // Process each element block connectivity to get the node_to_proc mapping.

    // The 'node_to_proc_pointer' vector maps the processor span in the
    // node_to_proc vector.  The processors that node 'node' (0-based)
    // is on are:
    //  * begin = node_to_proc_pointer[node]
    //  * end   = node_to_proc_pointer[node+1]
    //  * proc_list = node_to_proc[begin] .. node_to_proc[end-1]
    //

    size_t proc_count = connectivity.size();

    size_t                        node_count = region.get_property("node_count").get_int();
    std::vector<std::vector<int>> proc_node(node_count);

    // Assume that the majority of nodes will be on 2 or less
    // processors (hopefully, most are on 1).
    // Preallocate the proc_node[node] vector to 2 to minimize
    // resizes... Use 'reserve' instead of 'resize'
    for (size_t i = 0; i < node_count; i++) {
      proc_node[i].reserve(2);
    }
    progress("\tProc_node reserved");

    size_t sum_on_proc_count = 0;
    for (size_t p = 0; p < proc_count; p++) {
      size_t on_proc_count = 0;
      auto & ebs           = proc_region[p]->get_element_blocks();
      size_t block_count   = ebs.size();

      for (size_t b = 0; b < block_count; b++) {
        size_t element_count = ebs[b]->entity_count();
        size_t element_nodes = ebs[b]->get_property("topology_node_count").get_int();
        for (size_t i = 0; i < element_count * element_nodes; i++) {
          INT node = connectivity[p][b][i] - 1;
          if (proc_node[node].empty() ||
              proc_node[node][proc_node[node].size() - 1] != static_cast<int>(p)) {
            proc_node[node].push_back(p);
            on_proc_count++;
          }
        }
      }
      Ioss::NodeBlock *nb =
          new Ioss::NodeBlock(proc_region[p]->get_database(), "node_block1", on_proc_count, 3);
      proc_region[p]->add(nb);
      if (debug_level & 2) {
        OUTPUT << "\tProcessor " << p << " has " << on_proc_count << " nodes.\n";
      }
      sum_on_proc_count += on_proc_count;
    }
    progress("\tProc_node populated");

    // Have data for each node showing which processors it is on...
    // proc_node[node].size() is number of processors for this node...
    node_to_proc_pointer.reserve(node_count + 1);

    std::vector<size_t> proc_histo(17);

    size_t node_to_proc_pointer_size = 0;
    for (size_t i = 0; i < node_count; i++) {
      size_t num_procs = proc_node[i].size();
      if (num_procs == 0) {
        OUTPUT << "WARNING: Node " << i + 1 << " is not connected to any elements.\n";
      }
      else if (num_procs < proc_histo.size()) {
        proc_histo[num_procs]++;
      }
      else {
        proc_histo[0]++;
      }

      node_to_proc_pointer.push_back(node_to_proc_pointer_size);
      node_to_proc_pointer_size += num_procs;
    }
    // Output histogram..
    OUTPUT << "Processor count per node histogram:\n";
    for (size_t i = 0; i < proc_histo.size(); i++) {
      if (proc_histo[i] > 0) {
        OUTPUT << "\tNodes on " << i << " processors = " << proc_histo[i] << "\t("
               << (proc_histo[i] * 100 + node_count / 2) / node_count << ")%\n";
      }
    }
    if (proc_histo[0] > 0) {
      OUTPUT << "\tNodes on " << proc_histo.size() << " or more processors = " << proc_histo[0]
             << "\t(" << (proc_histo[0] * 100 + node_count / 2) / node_count << ")%\n";
    }
    OUTPUT << "\n";

    node_to_proc_pointer.push_back(node_to_proc_pointer_size);
    node_to_proc.reserve(node_to_proc_pointer_size);
    progress("\tNode_to_proc reserved");
    assert(sum_on_proc_count == node_to_proc_pointer_size);

    for (auto pn : proc_node) {
      size_t num_procs = pn.size();
      for (size_t p = 0; p < num_procs; p++) {
        node_to_proc.push_back(pn[p]);
      }
    }
    assert(node_to_proc.size() == node_to_proc_pointer_size);
    progress("\tNode_to_proc populated");
  }

  template <typename INT>
  void slice(Ioss::Region &region, const std::string &nemfile, SystemInterface &interface,
             INT /*dummy*/)
  {
    progress(__func__);
    std::vector<Ioss::Region *> proc_region(interface.processor_count());
    bool                        ints64 = (sizeof(INT) == 8);

    for (size_t i = 0; i < interface.processor_count(); i++) {
      std::string outfile = Ioss::Utils::decode_filename(nemfile, i, interface.processor_count());
      Ioss::DatabaseIO *dbo =
          Ioss::IOFactory::create("exodus", outfile, Ioss::WRITE_RESTART, (MPI_Comm)MPI_COMM_WORLD);
      if (ints64) {
        dbo->set_int_byte_size_api(Ioss::USE_INT64_API);
      }

      proc_region[i] = new Ioss::Region(dbo);
      proc_region[i]->begin_mode(Ioss::STATE_DEFINE_MODEL);
    }

    double           start = seacas_timer();
    std::vector<int> elem_to_proc;
    decompose_elements(region, interface, elem_to_proc);
    double end = seacas_timer();
    OUTPUT << "Decompose elements = " << end - start << "\n";

    start = seacas_timer();
    // Build the proc_elem_block_cnt[i][j] vector.
    // Gives number of elements in block i on processor j
    size_t block_count = region.get_property("element_block_count").get_int();
    std::vector<std::vector<INT>> proc_elem_block_cnt(block_count + 1);
    for (auto &pebc : proc_elem_block_cnt) {
      pebc.resize(interface.processor_count());
    }
    get_proc_elem_block_count(region, elem_to_proc, proc_elem_block_cnt);
    end = seacas_timer();

    OUTPUT << "Calculate elements per element block on each processor = " << end - start << "\n";

    // Create element blocks for each processor...
    for (size_t p = 0; p < interface.processor_count(); p++) {
      auto & ebs = region.get_element_blocks();
      size_t bc  = ebs.size();
      for (size_t b = 0; b < bc; b++) {
        std::string type = ebs[b]->get_property("topology_type").get_string();
        auto *eb = new Ioss::ElementBlock(proc_region[p]->get_database(), ebs[b]->name(), type,
                                          proc_elem_block_cnt[b][p]);
        proc_region[p]->add(eb);
      }
    }

    start = seacas_timer();
    // Now, build the proc-to-element vector which lists the elements
    // on each processor. Know the size of each vector, so can
    // preallocate.  The 'get_proc_elem_block_count()' function also
    // gives us the number of elements in each block on each
    // processor, so it is useful for more just give a preallocation size...
    std::vector<std::vector<INT>> proc_to_elem(interface.processor_count());
    for (size_t i = 0; i < interface.processor_count(); i++) {
      proc_to_elem[i].reserve(
          proc_elem_block_cnt[block_count][i]); // Reserve, not resize so doesn't initialize
    }

    size_t element_count = region.get_property("element_count").get_int();
    for (size_t i = 0; i < element_count; i++) {
      size_t processor = elem_to_proc[i];
      proc_to_elem[processor].push_back(i);
    }

    for (size_t i = 0; i < interface.processor_count(); i++) {
      assert(proc_to_elem[i].size() == (size_t)proc_elem_block_cnt[block_count][i]);
    }

    // Read connectivity and partition to each processor/block.
    // connectvity[p][b] = connectivity for block b on processor p
    std::vector<std::vector<std::vector<INT>>> connectivity(interface.processor_count());
    get_connectivity(region, proc_region, elem_to_proc, connectivity);
    end = seacas_timer();

    OUTPUT << "Get connectivity lists for each element block on each processor = " << end - start
           << "\n";

    // Now that we have the elements on each processor and the element
    // blocks those elements are in, can generate the node to proc list...

    start = seacas_timer();
    std::vector<int> node_to_proc;
    std::vector<INT> node_to_proc_pointer;
    get_node_to_proc(region, proc_region, connectivity, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "Node Categorization Time = " << end - start << "\n";

    // Communication map data -- interior/border nodes
    start = seacas_timer();
    std::vector<std::vector<INT>> border_node_proc_map(interface.processor_count());
    define_communication_data(region, proc_region, node_to_proc, node_to_proc_pointer,
                              border_node_proc_map);
    end = seacas_timer();
    OUTPUT << "Communication Data Definitions = " << end - start << "\n";

    // Determine nodeset distribution to processor regions.
    start = seacas_timer();
    get_nodesets(region, proc_region, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "Get nodeset data = " << end - start << "\n";

    start = seacas_timer();
    get_sidesets(region, proc_region, elem_to_proc, (INT)0);
    end = seacas_timer();
    OUTPUT << "Get sideset data = " << end - start << "\n";

    start             = seacas_timer();
    double start_comb = start;
    OUTPUT << "Begin writing  output files\n";
    size_t proc_count = interface.processor_count();
    for (size_t p = 0; p < proc_count; p++) {
      proc_region[p]->synchronize_id_and_name(&region);
      proc_region[p]->end_mode(Ioss::STATE_DEFINE_MODEL);
      proc_region[p]->begin_mode(Ioss::STATE_MODEL);
      if (minimize_open_files) {
        proc_region[p]->get_database()->closeDatabase();
      }
      proc_progress(p, proc_count);
    }
    end = seacas_timer();
    OUTPUT << "\tDefine output databases = " << end - start << "\n";

// Generate and output node map...
#if 1
    start = seacas_timer();
    output_node_map(region, proc_region, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "\tNode Map Output = " << end - start << "\n";
#else
    start = seacas_timer();
    output_global_node_map(region, proc_region, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "\tGlobal Node Map Output = " << end - start << "\n";
#endif

    start = seacas_timer();
    output_element_map(region, proc_region, elem_to_proc, (INT)1);
    end = seacas_timer();
    OUTPUT << "\tElement Map Output = " << end - start << "\n";

    start = seacas_timer();
    output_communication_map(proc_region, border_node_proc_map);
    end = seacas_timer();
    OUTPUT << "\tCommunication map Output = " << end - start << "\n";

    start = seacas_timer();
    output_connectivity(proc_region, connectivity);
    end = seacas_timer();
    OUTPUT << "\tConnectivity Output = " << end - start << "\n";

    // Free up connectivity space...
    free_connectivity_storage(connectivity);

    start = seacas_timer();
    output_coordinates(region, proc_region, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "\tCoordinates Output = " << end - start << "\n";

    start = seacas_timer();
    output_nodesets(region, proc_region, node_to_proc, node_to_proc_pointer);
    end = seacas_timer();
    OUTPUT << "\tNodeset Output = " << end - start << "\n";

    start = seacas_timer();
    output_sidesets(region, proc_region, elem_to_proc, (INT)0);
    end = seacas_timer();
    OUTPUT << "\tSideset Output = " << end - start << "\n";

    // Close all files...
    start = seacas_timer();
    for (size_t p = 0; p < interface.processor_count(); p++) {
      proc_region[p]->end_mode(Ioss::STATE_MODEL);
      delete proc_region[p];
    }
    end = seacas_timer();
    OUTPUT << "\tClose and finalize all output databases = " << end - start << "\n";
    OUTPUT << "Total time to write output files = " << end - start_comb << " ("
           << (end - start_comb) / interface.processor_count() << " per file)\n";
  }
#if defined(__PUMAGON__)
#include <stdio.h>
#else
#include <climits>
#include <unistd.h>
#endif

  int get_free_descriptor_count()
  {
// Returns maximum number of files that one process can have open
// at one time. (POSIX)
#if defined(__PUMAGON__)
    int fdmax = FOPEN_MAX;
#else
    int fdmax = sysconf(_SC_OPEN_MAX);
    if (fdmax == -1) {
      // POSIX indication that there is no limit on open files...
      fdmax = INT_MAX;
    }
#endif
    // File descriptors are assigned in order (0,1,2,3,...) on a per-process
    // basis.

    // Assume that we have stdin, stdout, stderr, and output exodus
    // file (4 total).

    return fdmax - 4;

    // Could iterate from 0..fdmax and check for the first
    // EBADF (bad file descriptor) error return from fcntl, but that takes
    // too long and may cause other problems.  There is a isastream(filedes)
    // call on Solaris that can be used for system-dependent code.
    //
    // Another possibility is to do an open and see which descriptor is
    // returned -- take that as 1 more than the current count of open files.
    //
  }

  void filename_substitution(std::string &filename, const SystemInterface &interface)
  {
    // See if filename contains "%P" which is replaced by the number of processors...
    // Assumes that %P only occurs once...
    // filename is changed.
    size_t pos = filename.find("%P");
    if (pos != std::string::npos) {
      // Found the characters...  Replace with the processor count...
      size_t      num_proc = interface.processor_count();
      std::string tmp(filename, 0, pos);
      tmp += std::to_string(num_proc);
      tmp += filename.substr(pos + 2);
      filename = tmp;
    }

    // If contains %M, replace with the decomposition method.
    pos = filename.find("%M");
    if (pos != std::string::npos) {
      // Found the characters...  Replace with the input file basename...
      const std::string &method_name = interface.decomposition_method();
      std::string        tmp(filename, 0, pos);
      tmp += method_name;
      tmp += filename.substr(pos + 2);
      filename = tmp;
    }
  }
} // namespace
