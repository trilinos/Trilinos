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

#if USE_ZOLTAN
#include <zoltan.h>     // for Zoltan_Initialize
#include <zoltan_cpp.h> // for Zoltan
#endif

extern int    debug_level;
extern double seacas_timer();
extern void   progress(const std::string &output);

namespace {
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

#if USE_ZOLTAN
  template <typename INT>
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
  get_element_centroid(const Ioss::Region &region, IOSS_MAYBE_UNUSED INT dummy)
  {
    size_t element_count = region.get_property("element_count").get_int();

    // The zoltan methods supported in slice are all geometry based
    // and use the element centroid.
    std::vector<double> x(element_count);
    std::vector<double> y(element_count);
    std::vector<double> z(element_count);

    const auto         *nb = region.get_node_blocks()[0];
    std::vector<double> coor;
    nb->get_field_data("mesh_model_coordinates", coor);

    const auto &blocks = region.get_element_blocks();
    size_t      el     = 0;
    for (auto &eb : blocks) {
      std::vector<INT> connectivity;
      eb->get_field_data("connectivity_raw", connectivity);
      size_t blk_element_count = eb->entity_count();
      size_t blk_element_nodes = eb->topology()->number_nodes();

      for (size_t j = 0; j < blk_element_count; j++) {
        for (size_t k = 0; k < blk_element_nodes; k++) {
          auto node = connectivity[j * blk_element_nodes + k] - 1;
          x[el] += coor[node * 3 + 0];
          y[el] += coor[node * 3 + 1];
          z[el] += coor[node * 3 + 2];
        }
        x[el] /= blk_element_nodes;
        y[el] /= blk_element_nodes;
        z[el] /= blk_element_nodes;
        el++;
      }
    }
    return {x, y, z};
  }
  /*****************************************************************************/
  /***** Global data structure used by Zoltan callbacks.                   *****/
  /***** Could implement Zoltan callbacks without global data structure,   *****/
  /***** but using the global data structure makes implementation quick.   *****/
  struct
  {
    size_t  ndot; /* Length of x, y, z, and part (== # of elements) */
    int    *vwgt; /* vertex weights */
    double *x;    /* x-coordinates */
    double *y;    /* y-coordinates */
    double *z;    /* z-coordinates */
  } Zoltan_Data;

  /*****************************************************************************/
  /***** ZOLTAN CALLBACK FUNCTIONS *****/
  int zoltan_num_dim(void * /*data*/, int *ierr)
  {
    /* Return dimensionality of coordinate data.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    if (Zoltan_Data.z != nullptr) {
      return 3;
    }
    if (Zoltan_Data.y != nullptr) {
      return 2;
    }
    return 1;
  }

  int zoltan_num_obj(void * /*data*/, int *ierr)
  {
    /* Return number of objects.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    return Zoltan_Data.ndot;
  }

  void zoltan_obj_list(void * /*data*/, int /*ngid_ent*/, int /*nlid_ent*/, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR /*lids*/, int wdim, float *wgts, int *ierr)
  {
    /* Return list of object IDs.
     * Return only global IDs; don't need local IDs since running in serial.
     * gids are array indices for coordinate and vwgts arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    std::iota(gids, gids + Zoltan_Data.ndot, 0);
    if (wdim != 0) {
      for (size_t i = 0; i < Zoltan_Data.ndot; i++) {
        wgts[i] = static_cast<float>(Zoltan_Data.vwgt[i]);
      }
    }

    *ierr = ZOLTAN_OK;
  }

  void zoltan_geom(void * /*data*/, int /*ngid_ent*/, int /*nlid_ent*/, int nobj,
                   const ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR /*lids*/, int ndim, double *geom,
                   int *ierr)
  {
    /* Return coordinates for objects.
     * gids are array indices for coordinate arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */

    for (size_t i = 0; i < static_cast<size_t>(nobj); i++) {
      size_t j       = gids[i];
      geom[i * ndim] = Zoltan_Data.x[j];
      if (ndim > 1) {
        geom[i * ndim + 1] = Zoltan_Data.y[j];
      }
      if (ndim > 2) {
        geom[i * ndim + 2] = Zoltan_Data.z[j];
      }
    }

    *ierr = ZOLTAN_OK;
  }

  template <typename INT>
  void decompose_zoltan(const Ioss::Region &region, int ranks, SystemInterface &interFace,
                        std::vector<int> &elem_to_proc, const std::vector<int> &weights,
                        IOSS_MAYBE_UNUSED INT dummy)
  {
    if (ranks == 1) {
      return;
    }

    size_t element_count = region.get_property("element_count").get_int();
    if (element_count != static_cast<size_t>(static_cast<int>(element_count))) {
      fmt::print(stderr, "ERROR: Cannot have a mesh with more than 2.1 Billion elements in a "
                         "Zoltan decomposition.\n");
      exit(EXIT_FAILURE);
    }

    auto [x, y, z] = get_element_centroid(region, dummy);

    // Copy mesh data and pointers into structure accessible from callback fns.
    Zoltan_Data.ndot = element_count;
    Zoltan_Data.vwgt = const_cast<int *>(Data(weights));

    if (interFace.ignore_x_ && interFace.ignore_y_) {
      Zoltan_Data.x = Data(z);
    }
    else if (interFace.ignore_x_ && interFace.ignore_z_) {
      Zoltan_Data.x = Data(y);
    }
    else if (interFace.ignore_y_ && interFace.ignore_z_) {
      Zoltan_Data.x = Data(x);
    }
    else if (interFace.ignore_x_) {
      Zoltan_Data.x = Data(y);
      Zoltan_Data.y = Data(z);
    }
    else if (interFace.ignore_y_) {
      Zoltan_Data.x = Data(x);
      Zoltan_Data.y = Data(z);
    }
    else if (!interFace.ignore_z_) {
      Zoltan_Data.x = Data(x);
      Zoltan_Data.y = Data(y);
    }
    else {
      Zoltan_Data.x = Data(x);
      Zoltan_Data.y = Data(y);
      Zoltan_Data.z = Data(z);
    }

    // Initialize Zoltan
    int    argc = 0;
    char **argv = nullptr;

    float ver = 0.0;
    Zoltan_Initialize(argc, argv, &ver);
    fmt::print("Using Zoltan version {:.2}, method {}\n", static_cast<double>(ver),
               interFace.decomposition_method());

    Zoltan zz(Ioss::ParallelUtils::comm_world());

    // Register Callback functions
    // Using global Zoltan_Data; could register it here instead as data field.
    zz.Set_Num_Obj_Fn(zoltan_num_obj, nullptr);
    zz.Set_Obj_List_Fn(zoltan_obj_list, nullptr);
    zz.Set_Num_Geom_Fn(zoltan_num_dim, nullptr);
    zz.Set_Geom_Multi_Fn(zoltan_geom, nullptr);

    // Set parameters for Zoltan
    zz.Set_Param("DEBUG_LEVEL", "0");
    std::string str = fmt::format("{}", ranks);
    zz.Set_Param("NUM_GLOBAL_PARTS", str);
    zz.Set_Param("OBJ_WEIGHT_DIM", "1");
    zz.Set_Param("LB_METHOD", interFace.decomposition_method());
    zz.Set_Param("NUM_LID_ENTRIES", "0");
    zz.Set_Param("REMAP", "0");
    zz.Set_Param("RETURN_LISTS", "PARTITION_ASSIGNMENTS");
    zz.Set_Param("RCB_RECTILINEAR_BLOCKS", "1");

    int num_global = sizeof(INT) / sizeof(ZOLTAN_ID_TYPE);
    num_global     = num_global < 1 ? 1 : num_global;

    // Call partitioner
    int           changes           = 0;
    int           num_local         = 0;
    int           num_import        = 1;
    int           num_export        = 1;
    ZOLTAN_ID_PTR import_global_ids = nullptr;
    ZOLTAN_ID_PTR import_local_ids  = nullptr;
    ZOLTAN_ID_PTR export_global_ids = nullptr;
    ZOLTAN_ID_PTR export_local_ids  = nullptr;
    int          *import_procs      = nullptr;
    int          *import_to_part    = nullptr;
    int          *export_procs      = nullptr;
    int          *export_to_part    = nullptr;
    int rc = zz.LB_Partition(changes, num_global, num_local, num_import, import_global_ids,
                             import_local_ids, import_procs, import_to_part, num_export,
                             export_global_ids, export_local_ids, export_procs, export_to_part);

    if (rc != ZOLTAN_OK) {
      fmt::print(stderr, "ERROR: Problem during call to Zoltan LB_Partition.\n");
      goto End;
    }

    // Sanity check
    if (element_count != static_cast<size_t>(num_export)) {
      fmt::print(stderr, "Sanity check failed; ndot {} != num_export {}.\n", element_count,
                 static_cast<size_t>(num_export));
      goto End;
    }

    elem_to_proc.resize(element_count);
    for (size_t i = 0; i < element_count; i++) {
      elem_to_proc[i] = export_to_part[i];
    }

  End:
    /* Clean up */
    Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs, &export_to_part);
    Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs, &export_to_part);
  }
#endif

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

  void output_histogram(const std::vector<size_t> &proc_work, size_t avg_work, size_t median)
  {
    fmt::print("Work-per-processor Histogram\n");
    std::array<size_t, 16> histogram{};

    auto wmin = *std::min_element(proc_work.begin(), proc_work.end());
    auto wmax = *std::max_element(proc_work.begin(), proc_work.end());

    size_t hist_size = std::min(size_t(16), (wmax - wmin));
    hist_size        = std::min(hist_size, proc_work.size());

    if (hist_size <= 1) {
      fmt::print("\tWork is the same on all processors; no histogram needed.\n\n");
      return;
    }

    auto delta = double(wmax + 1 - wmin) / hist_size;
    for (const auto &pw : proc_work) {
      auto bin = size_t(double(pw - wmin) / delta);
      SMART_ASSERT(bin < hist_size)(bin)(hist_size);
      histogram[bin]++;
    }

    size_t proc_width = Ioss::Utils::number_width(proc_work.size(), true);
    size_t work_width = Ioss::Utils::number_width(wmax, true);

    fmt::print("\n\t{:^{}} {:^{}}\n", "Work Range", 2 * work_width + 2, "#", proc_width);
    auto hist_max = *std::max_element(histogram.begin(), histogram.end());
    for (size_t i = 0; i < hist_size; i++) {
      int         max_star = 50;
      int         star_cnt = ((double)histogram[i] / hist_max * max_star);
      std::string stars(star_cnt, '*');
      for (int j = 9; j < star_cnt;) {
        stars[j] = '|';
        j += 10;
      }
      if (histogram[i] > 0 && star_cnt == 0) {
        stars = '.';
      }
      size_t      w1 = wmin + size_t(i * delta);
      size_t      w2 = wmin + size_t((i + 1) * delta);
      std::string postfix;
      if (w1 <= avg_work && avg_work < w2) {
        postfix += "average";
      }
      if (w1 <= median && median < w2) {
        if (!postfix.empty()) {
          postfix += ", ";
        }
        postfix += "median";
      }
      fmt::print("\t{:{}}..{:{}} ({:{}}):\t{:{}}  {}\n", fmt::group_digits(w1), work_width,
                 fmt::group_digits(w2), work_width, fmt::group_digits(histogram[i]), proc_width,
                 stars, max_star, postfix);
    }
    fmt::print("\n");
  }

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
                                             const std::vector<int> &weights,
                                             IOSS_MAYBE_UNUSED int   dummy);
template std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                             const std::vector<int>   &weights,
                                             IOSS_MAYBE_UNUSED int64_t dummy);

template <typename INT>
std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                    const std::vector<int> &weights, IOSS_MAYBE_UNUSED INT dummy)
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
#if USE_ZOLTAN
    decompose_zoltan(region, interFace.processor_count(), interFace, elem_to_proc, weights, dummy);
#else
    fmt::print(stderr, "ERROR: Zoltan library not enabled in this version of slice.\n"
                       "       The 'rcb', 'rib', and 'hsfc' methods are not available.\n\n");
    std::exit(1);
#endif
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
      char **names        = Ioss::Utils::get_name_array(map_count, max_name_length);
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
      Ioss::Utils::delete_name_array(names, map_count);
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

template <typename INT>
std::map<INT, std::vector<INT>> string_chains(const Ioss::chain_t<INT> &element_chains)
{
  std::map<INT, std::vector<INT>> chains;

  for (size_t i = 0; i < element_chains.size(); i++) {
    auto &chain_entry = element_chains[i];
    if (chain_entry.link >= 0) {
      chains[chain_entry.element].push_back(i + 1);
    }
  }
  return chains;
}

template std::vector<int> line_decomp_weights(const Ioss::chain_t<int> &element_chains,
                                              size_t                    element_count);
template std::vector<int> line_decomp_weights(const Ioss::chain_t<int64_t> &element_chains,
                                              size_t                        element_count);

template <typename INT>
std::vector<int> line_decomp_weights(const Ioss::chain_t<INT> &element_chains, size_t element_count)
{
  auto chains = string_chains(element_chains);

  if ((debug_level & 16) != 0) {
    for (const auto &[chain_root, chain_elements] : chains) {
      fmt::print("Chain Root: {} contains: {}\n", chain_root, fmt::join(chain_elements, ", "));
    }
  }

  std::vector<int> weights(element_count, 1);
  // Now, for each chain...
  for (const auto &[chain_root, chain_elements] : chains) {
    // * Set the weights of all elements in the chain...
    // * non-root = 0, root = length of chain.
    for (const auto &element : chain_elements) {
      weights[element - 1] = 0;
    }
    weights[chain_root - 1] = static_cast<int>(chain_elements.size());
  }
  return weights;
}

template void line_decomp_modify(const Ioss::chain_t<int> &element_chains,
                                 std::vector<int> &elem_to_proc, int proc_count);
template void line_decomp_modify(const Ioss::chain_t<int64_t> &element_chains,
                                 std::vector<int> &elem_to_proc, int proc_count);

template <typename INT>
void line_decomp_modify(const Ioss::chain_t<INT> &element_chains, std::vector<int> &elem_to_proc,
                        int proc_count)
{
  // Get a map of all chains and the elements in the chains.  Map key will be root.
  auto chains = string_chains(element_chains);

  // Delta: elements added/removed from each processor...
  std::vector<int> delta(proc_count);

  // Now, for each chain...
  for (const auto &[chain_root, chain_elements] : chains) {
    if ((debug_level & 16) != 0) {
      fmt::print("Chain Root: {} contains: {}\n", chain_root, fmt::join(chain_elements, ", "));
    }

    std::vector<INT> chain_proc_count(proc_count);

    // * get processors used by elements in the chain...
    for (const auto &element : chain_elements) {
      auto proc = elem_to_proc[element - 1];
      chain_proc_count[proc]++;
    }

    // * Now, subtract the `delta` from each count
    for (int i = 0; i < proc_count; i++) {
      chain_proc_count[i] -= delta[i];
    }

    // * Assign all elements in the chain to processor at chain root
    // * Update the deltas for all processors that gain/lose elements...
    auto root_proc = elem_to_proc[chain_root - 1];
    for (const auto &element : chain_elements) {
      if (elem_to_proc[element - 1] != root_proc) {
        auto old_proc             = elem_to_proc[element - 1];
        elem_to_proc[element - 1] = root_proc;
        delta[root_proc]++;
        delta[old_proc]--;
      }
    }
  }

  std::vector<INT> proc_element_count(proc_count);
  for (auto proc : elem_to_proc) {
    proc_element_count[proc]++;
  }
  if ((debug_level & 32) != 0) {
    fmt::print("\nElements/Processor: {}\n", fmt::join(proc_element_count, ", "));
    fmt::print("Delta/Processor:    {}\n", fmt::join(delta, ", "));
  }
}

template void output_decomposition_statistics(const std::vector<int> &elem_to_proc, int proc_count,
                                              size_t number_elements);
template void output_decomposition_statistics(const std::vector<int64_t> &elem_to_proc,
                                              int proc_count, size_t number_elements);
template <typename INT>
void output_decomposition_statistics(const std::vector<INT> &elem_to_proc, int proc_count,
                                     size_t number_elements)
{
  // Output histogram of elements / rank...
  std::vector<size_t> elem_per_rank(proc_count);
  for (INT proc : elem_to_proc) {
    elem_per_rank[proc]++;
  }

  size_t proc_width = Ioss::Utils::number_width(proc_count, false);
  size_t work_width = Ioss::Utils::number_width(number_elements, true);

  auto   min_work = *std::min_element(elem_per_rank.begin(), elem_per_rank.end());
  auto   max_work = *std::max_element(elem_per_rank.begin(), elem_per_rank.end());
  size_t median   = 0;
  {
    auto pw_copy(elem_per_rank);
    std::nth_element(pw_copy.begin(), pw_copy.begin() + pw_copy.size() / 2, pw_copy.end());
    median = pw_copy[pw_copy.size() / 2];
    fmt::print("\nElements per processor:\n\tMinimum = {}, Maximum = {}, Median = {}, Ratio = "
               "{:.3}\n\n",
               fmt::group_digits(min_work), fmt::group_digits(max_work), fmt::group_digits(median),
               (double)(max_work) / min_work);
  }
  if (min_work == max_work) {
    fmt::print("\nWork on all processors is {}\n\n", fmt::group_digits(min_work));
  }
  else {
    int max_star = 40;
    int min_star = max_star * ((double)min_work / (double)(max_work));
    min_star     = std::max(1, min_star);
    int delta    = max_star - min_star;

    double avg_work = (double)number_elements / (double)proc_count;
    for (size_t i = 0; i < elem_per_rank.size(); i++) {
      int star_cnt =
          (double)(elem_per_rank[i] - min_work) / (max_work - min_work) * delta + min_star;
      std::string stars(star_cnt, '*');
      std::string format = "\tProcessor {:{}}, work = {:{}}  ({:.2f})\t{}\n";
      if (elem_per_rank[i] == max_work) {
        fmt::print(
#if !defined __NVCC__
            fg(fmt::color::red),
#endif
            format, i, proc_width, fmt::group_digits(elem_per_rank[i]), work_width,
            (double)elem_per_rank[i] / avg_work, stars);
      }
      else if (elem_per_rank[i] == min_work) {
        fmt::print(
#if !defined __NVCC__
            fg(fmt::color::green),
#endif
            format, i, proc_width, fmt::group_digits(elem_per_rank[i]), work_width,
            elem_per_rank[i] / avg_work, stars);
      }
      else {
        fmt::print(format, i, proc_width, fmt::group_digits(elem_per_rank[i]), work_width,
                   elem_per_rank[i] / avg_work, stars);
      }
    }

    // Output Histogram...
    output_histogram(elem_per_rank, (size_t)avg_work, median);
  }
}
