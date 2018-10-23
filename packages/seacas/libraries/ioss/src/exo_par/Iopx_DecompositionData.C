// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <Ioss_CodeTypes.h>
#include <Ioss_ElementTopology.h> // for ElementTopology
#include <Ioss_Field.h>           // for Field, etc
#include <Ioss_Map.h>             // for Map, MapContainer
#include <Ioss_ParallelUtils.h>   // for ParallelUtils, etc
#include <Ioss_PropertyManager.h> // for PropertyManager
#include <Ioss_Sort.h>
#include <Ioss_Utils.h> // for TOPTR, Utils
#include <exo_par/Iopx_DecompositionData.h>
#include <exodus/Ioex_Utils.h>

#include <algorithm> // for lower_bound, copy, etc
#include <cassert>   // for assert
#include <climits>   // for INT_MAX
#include <cmath>
#include <cstdlib> // for exit, EXIT_FAILURE
#include <cstring>
#include <iostream> // for operator<<, ostringstream, etc
#include <iterator> // for distance
#include <map>      // for map
#include <numeric>  // for accumulate
#include <utility>  // for pair, make_pair

#if !defined(NO_PARMETIS_SUPPORT)
#include <parmetis.h> // for ParMETIS_V3_Mesh2Dual, etc
#endif

#if !defined(NO_ZOLTAN_SUPPORT)
#include <zoltan.h>     // for Zoltan_Initialize
#include <zoltan_cpp.h> // for Zoltan
#endif

namespace {
  // ZOLTAN Callback functions...

#if !defined(NO_ZOLTAN_SUPPORT)
  int zoltan_num_dim(void *data, int *ierr)
  {
    // Return dimensionality of coordinate data.
    Iopx::DecompositionDataBase *zdata = reinterpret_cast<Iopx::DecompositionDataBase *>(data);

    *ierr = ZOLTAN_OK;
    return zdata->spatial_dimension();
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    // Return number of objects (element count) on this processor...
    Iopx::DecompositionDataBase *zdata = reinterpret_cast<Iopx::DecompositionDataBase *>(data);

    *ierr = ZOLTAN_OK;
    return zdata->decomp_elem_count();
  }

  void zoltan_obj_list(void *data, int ngid_ent, int /*nlid_ent*/, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR lids, int wdim, float *wgts, int *ierr)
  {
    // Return list of object IDs, both local and global.
    Iopx::DecompositionDataBase *zdata = reinterpret_cast<Iopx::DecompositionDataBase *>(data);

    // At the time this is called, we don't have much information
    // These routines are the ones that are developing that
    // information...
    size_t element_count  = zdata->decomp_elem_count();
    size_t element_offset = zdata->decomp_elem_offset();

    *ierr = ZOLTAN_OK;

    if (lids != nullptr) {
      std::iota(lids, lids + element_count, 0);
    }

    if (wdim != 0) {
      std::fill(wgts, wgts + element_count, 1.0);
    }

    if (ngid_ent == 1) {
      std::iota(gids, gids + element_count, element_offset);
    }
    else if (ngid_ent == 2) {
      int64_t *global_ids = reinterpret_cast<int64_t *>(gids);
      std::iota(global_ids, global_ids + element_count, element_offset);
    }
    else {
      *ierr = ZOLTAN_FATAL;
    }
    return;
  }

  void zoltan_geom(void *data, int /*ngid_ent*/, int /*nlid_ent*/, int /*nobj*/,
                   ZOLTAN_ID_PTR /*gids*/, ZOLTAN_ID_PTR /*lids*/, int /*ndim*/, double *geom,
                   int *ierr)
  {
    // Return coordinates for objects.
    Iopx::DecompositionDataBase *zdata = reinterpret_cast<Iopx::DecompositionDataBase *>(data);

    std::copy(zdata->centroids().begin(), zdata->centroids().end(), &geom[0]);

    *ierr = ZOLTAN_OK;
  }
#endif
} // namespace

namespace Iopx {
  template DecompositionData<int>::DecompositionData(const Ioss::PropertyManager &props,
                                                     MPI_Comm                     communicator);
  template DecompositionData<int64_t>::DecompositionData(const Ioss::PropertyManager &props,
                                                         MPI_Comm                     communicator);

  template <typename INT>
  DecompositionData<INT>::DecompositionData(const Ioss::PropertyManager &props,
                                            MPI_Comm                     communicator)
      : DecompositionDataBase(communicator), m_decomposition(props, communicator)
  {
    MPI_Comm_rank(comm_, &m_processor);
    MPI_Comm_size(comm_, &m_processorCount);
  }

  template <typename INT> void DecompositionData<INT>::decompose_model(int filePtr)
  {
    m_decomposition.show_progress(__func__);
    // Initial decomposition is linear where processor #p contains
    // elements from (#p * #element/#proc) to (#p+1 * #element/#proc)

    ex_init_params info{};
    ex_get_init_ext(filePtr, &info);

    size_t globalElementCount          = info.num_elem;
    size_t globalNodeCount             = info.num_nodes;
    m_decomposition.m_spatialDimension = info.num_dim;
    el_blocks.resize(info.num_elem_blk);

    // Generate element_dist/node_dist --  size proc_count + 1
    // processor p contains all elements/nodes from X_dist[p] .. X_dist[p+1]
    m_decomposition.generate_entity_distributions(globalNodeCount, globalElementCount);

    generate_adjacency_list(filePtr, m_decomposition);

#if IOSS_DEBUG_OUTPUT
    std::cerr << "Processor " << m_processor << " has " << decomp_elem_count()
              << " elements; offset = " << decomp_elem_offset() << "\n";
    std::cerr << "Processor " << m_processor << " has " << decomp_node_count()
              << " nodes; offset = " << decomp_node_offset() << ".\n";
#endif

    if (m_decomposition.needs_centroids()) {
      // Get my coordinate data using direct exodus calls
      size_t size = decomp_node_count();
      if (size == 0) {
        size = 1; // Workaround for ambiguity in ex_get_partial_coord
      }

      std::vector<double> x(size);
      std::vector<double> y;
      std::vector<double> z;
      if (m_decomposition.m_spatialDimension > 1) {
        y.resize(size);
      }
      if (m_decomposition.m_spatialDimension > 2) {
        z.resize(size);
      }

      m_decomposition.show_progress("\tex_get_partial_coord");
      ex_get_partial_coord(filePtr, decomp_node_offset() + 1, decomp_node_count(), TOPTR(x),
                           TOPTR(y), TOPTR(z));

      m_decomposition.calculate_element_centroids(x, y, z);
    }

#if !defined(NO_ZOLTAN_SUPPORT)
    float version = 0.0;
    Zoltan_Initialize(0, nullptr, &version);

    Zoltan zz(comm_);

    // Register Zoltan Callback functions...
    zz.Set_Num_Obj_Fn(zoltan_num_obj, this);
    zz.Set_Obj_List_Fn(zoltan_obj_list, this);
    zz.Set_Num_Geom_Fn(zoltan_num_dim, this);
    zz.Set_Geom_Multi_Fn(zoltan_geom, this);
#endif

    m_decomposition.decompose_model(
#if !defined(NO_ZOLTAN_SUPPORT)
        zz,
#endif
        el_blocks);

    if (info.num_node_sets > 0) {
      get_nodeset_data(filePtr, info.num_node_sets);
    }

    if (info.num_side_sets > 0) {
      // Create elemGTL map which is used for sidesets (also element sets)
      build_global_to_local_elem_map();
      get_sideset_data(filePtr, info.num_side_sets);
    }

    // Have all the decomposition data needed
    // Can now populate the Ioss metadata...
    m_decomposition.show_progress("\tFinished with Iopx::decompose_model");

    if (m_decomposition.m_showHWM || m_decomposition.m_showProgress) {
      int64_t             min, max, avg;
      Ioss::ParallelUtils pu(m_decomposition.m_comm);
      pu.hwm_memory_stats(min, max, avg);
      int64_t MiB = 1024 * 1024;
      if (m_processor == 0) {
        std::cerr << "\n\tHigh Water Memory at end of Decomposition: " << min / MiB << "M  "
                  << max / MiB << "M  " << avg / MiB << "M\n";
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::generate_adjacency_list(int                       filePtr,
                                                       Ioss::Decomposition<INT> &decomposition)
  {
    m_decomposition.show_progress(__func__);
    // Range of elements currently handled by this processor [)
    size_t p_start = decomp_elem_offset();
    size_t p_end   = p_start + decomp_elem_count();

    size_t block_count = el_blocks.size();

    std::vector<ex_block> ebs(block_count);
    std::vector<INT>      ids(block_count);
    ex_get_ids(filePtr, EX_ELEM_BLOCK, TOPTR(ids));

    size_t sum    = 0; // Size of adjacency vector.
    size_t offset = 0;

    // Get the global element block index list at this time also.
    // The global element at index 'I' (0-based) is on block B
    // if global_block_index[B] <= I && global_block_index[B+1] < I
    decomposition.m_fileBlockIndex.resize(block_count + 1);

    for (size_t b = 0; b < block_count; b++) {
      el_blocks[b].id_ = ids[b];
      ebs[b].id        = ids[b];
      ebs[b].type      = EX_ELEM_BLOCK;
      ex_get_block_param(filePtr, &ebs[b]);

      // Range of elements in element block b [)
      size_t b_start = offset; // offset is index of first element in this block...
      offset += ebs[b].num_entry;
      size_t b_end = b_start + ebs[b].num_entry;

      if (b_start < p_end && p_start < b_end) {
        // Some of this blocks elements are on this processor...
        size_t overlap       = std::min(b_end, p_end) - std::max(b_start, p_start);
        size_t element_nodes = ebs[b].num_nodes_per_entry;

        sum += overlap * element_nodes;
      }
      decomposition.m_fileBlockIndex[b + 1] = decomposition.m_fileBlockIndex[b] + ebs[b].num_entry;
      el_blocks[b].topologyType             = ebs[b].topology;
      if (ebs[b].num_entry == 0 && (std::strcmp(ebs[b].topology, "nullptr") == 0)) {
        el_blocks[b].topologyType = "sphere";
      }
      el_blocks[b].globalCount    = ebs[b].num_entry;
      el_blocks[b].nodesPerEntity = ebs[b].num_nodes_per_entry;
      el_blocks[b].attributeCount = ebs[b].num_attribute;
    }

    // Check that the number of elements matches the m_fileBlockIndex[b+1] entry.
    // Reading a corrupt mesh in which there are elements not in an element block
    // can cause hard to track down problems...
    if (decomposition.m_globalElementCount != decomposition.m_fileBlockIndex[block_count]) {
      if (m_processor == 0) {
        std::ostringstream errmsg;
        errmsg << "ERROR: The sum of the element counts in each element block gives a total of "
               << decomposition.m_fileBlockIndex[block_count]
               << " elements.\n       This does not match the total element count of "
               << decomposition.m_globalElementCount
               << " which indicates a corrupt mesh description.\n"
               << "       Contact gdsjaar@sandia.gov for more details.\n";
        std::cerr << errmsg.str();
      }
      exit(EXIT_FAILURE);
    }

    // Make sure 'sum' can fit in INT...
    INT tmp_sum = (INT)sum;
    if ((size_t)tmp_sum != sum) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The decomposition of this mesh requires 64-bit integers, but is being\n"
             << "       run with 32-bit integer code. Please rerun with the property "
                "INTEGER_SIZE_API\n"
             << "       set to 8. The details of how to do this vary with the code that is being "
                "run.\n"
             << "       Contact gdsjaar@sandia.gov for more details.\n";
      std::cerr << errmsg.str();
      exit(EXIT_FAILURE);
    }

    decomposition.m_pointer.reserve(decomp_elem_count() + 1);
    decomposition.m_adjacency.reserve(sum);

    // Now, populate the vectors...
    offset = 0;
    sum    = 0; // Size of adjacency vector.

    for (auto &block : ebs) {
      // Range of elements in element block b [)
      size_t b_start = offset; // offset is index of first element in this block...
      offset += block.num_entry;
      size_t b_end = b_start + block.num_entry;

      m_decomposition.show_progress("\tex_get_partial_conn loop");
      if (b_start < p_end && p_start < b_end) {
        // Some of this blocks elements are on this processor...
        size_t  overlap       = std::min(b_end, p_end) - std::max(b_start, p_start);
        size_t  element_nodes = block.num_nodes_per_entry;
        int64_t id            = block.id;

        // Get the connectivity (raw) for this portion of elements...
        std::vector<INT> connectivity(overlap * element_nodes);
        size_t           blk_start = std::max(b_start, p_start) - b_start + 1;
#if IOSS_DEBUG_OUTPUT
        std::cerr << "Processor " << m_processor << " has " << overlap
                  << " elements on element block " << id << "\n";
#endif
        ex_get_partial_conn(filePtr, EX_ELEM_BLOCK, id, blk_start, overlap, TOPTR(connectivity),
                            nullptr, nullptr);
        size_t el = 0;
        for (size_t elem = 0; elem < overlap; elem++) {
          decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
          for (size_t k = 0; k < element_nodes; k++) {
            INT node = connectivity[el++] - 1; // 0-based node
            decomposition.m_adjacency.push_back(node);
          }
        }
        sum += overlap * element_nodes;
      }
    }
    decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
  }

  template <typename INT>
  void DecompositionData<INT>::get_nodeset_data(int filePtr, size_t set_count)
  {
    // Issues:
    // 1. Large node count in nodeset(s) that could overwhelm a single
    //    processor.  For example, every node in the model is in a
    //    nodeset.
    //    -- Cannot blindly read all nodeset data on proc 0 and
    //       broadcast.
    //
    // 2. Lots of small nodesets.  Communication involving all
    //    processors would result in lots of overhead if done on a
    //    nodeset by nodeset basis.
    //
    // 3. Application code will be requesting nodeset data on a
    //    nodeset-by-nodeset basis.  If communication needed each
    //    time, this could cause lots of slowdown/overhead...
    //    However, if cache the data, it could eat up memory.
    //
    // 4. Many models only need nodeset node list and then do nothing
    //    else with the nodeset data.
    //
    // 5. Only the processors that own any of the nodes in a nodeset
    //    should read that nodeset.
    //    -- What if nodelist is big and non-sorted so that a
    //       processors nodes cannot be read in a single contiguous
    //       read...
    //
    // 6. Alternatively, only processor 0 reads, but only communicates
    //    to the processors that have nodes in the nodeset.
    //    == REMEMBER: nodes are shared, so a node could be sent to
    //       multiple processors.
    //
    // Assumptions:
    // 1. Since we have already read coordinate data, we probably have
    //    some extra memory we can afford to use without increasing
    //    the high-water mark.
    //    -- Read the nodeset node lists in groups of size
    //       (3*globNodeCount/procCount*sizeof(double)/sizeof(INT)) or
    //       less.

    int root = 0; // Root processor that reads all nodeset bulk data (nodelists)

    node_sets.resize(set_count);

    std::vector<std::vector<INT>> set_nodelists(set_count);
    std::vector<ex_set>           sets(set_count);
    std::vector<INT>              ids(set_count);
    ex_get_ids(filePtr, EX_NODE_SET, TOPTR(ids));

    for (size_t i = 0; i < set_count; i++) {
      node_sets[i].id_                 = ids[i];
      sets[i].id                       = ids[i];
      sets[i].type                     = EX_NODE_SET;
      sets[i].entry_list               = nullptr;
      sets[i].extra_list               = nullptr;
      sets[i].distribution_factor_list = nullptr;
    }

    ex_get_sets(filePtr, sets.size(), TOPTR(sets));

    // Get total length of nset nodelists...
    size_t nodelist_size = 0;
    for (size_t i = 0; i < set_count; i++) {
      nodelist_size += sets[i].num_entry;
      node_sets[i].fileCount = sets[i].num_entry;
    }

    // Calculate the max "buffer" size usable for storing nodeset
    // nodelists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "nodeCount/2*2" is to
    // equalize the nodeCount among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = std::max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(INT));

    bool subsetting = false; // nodelist_size > max_size;

    if (subsetting) {
      assert(1 == 0);
    }
    else {
      // Can handle reading all nodeset node lists on a single
      // processor simultaneously.
      std::vector<INT> nodelist(nodelist_size);

      // Read the nodelists on root processor.
      if (m_processor == root) {
        size_t offset = 0;
        for (size_t i = 0; i < set_count; i++) {
          ex_get_set(filePtr, EX_NODE_SET, sets[i].id, &nodelist[offset], nullptr);
          offset += sets[i].num_entry;
        }
        assert(offset == nodelist_size);
      }

      // Broadcast this data to all other processors...
      MPI_Bcast(TOPTR(nodelist), sizeof(INT) * nodelist.size(), MPI_BYTE, root, comm_);

      // Each processor now has a complete list of all nodes in all
      // nodesets.  Determine which of these are owned by the current
      // processor...
      size_t offset = 0;
      for (size_t i = 0; i < set_count; i++) {
        size_t ns_beg = offset;
        size_t ns_end = ns_beg + sets[i].num_entry;

        for (size_t n = ns_beg; n < ns_end; n++) {
          INT node = nodelist[n];
          // See if node owned by this processor...
          if (i_own_node(node)) {
            // Save node in this processors nodelist for this set.
            // The saved data is this nodes location in the global
            // nodelist for this set.
            node_sets[i].entitylist_map.push_back(n - offset);
          }
        }
        offset = ns_end;
      }

      // Each processor knows how many of the nodeset nodes it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero node count is
      // the "root" for this nodeset.
      {
        std::vector<int> has_nodes_local(set_count);
        for (size_t i = 0; i < set_count; i++) {
          has_nodes_local[i] = node_sets[i].entitylist_map.empty() ? 0 : 1;
        }

        std::vector<int> has_nodes(set_count * m_processorCount);
        MPI_Allgather(TOPTR(has_nodes_local), has_nodes_local.size(), MPI_INT, TOPTR(has_nodes),
                      has_nodes_local.size(), MPI_INT, comm_);

        for (size_t i = 0; i < set_count; i++) {
          node_sets[i].hasEntities.resize(m_processorCount);
          node_sets[i].root_ = m_processorCount;
          int count          = 0;
          for (int p = 0; p < m_processorCount; p++) {
            if (p < node_sets[i].root_ && has_nodes[p * set_count + i] != 0) {
              node_sets[i].root_ = p;
            }
            node_sets[i].hasEntities[p] = has_nodes[p * set_count + i];
            count += has_nodes[p * set_count + i];
          }
	  int color = node_sets[i].hasEntities[m_processor] ? 1 : MPI_UNDEFINED;
	  MPI_Comm_split(comm_, color, m_processor, &node_sets[i].setComm_);
        }
      }

      // Check nodeset distribution factors to determine whether they
      // are all constant or if they contain varying values that must
      // be communicated.  If constant or empty, then they can be
      // "read" with no communication.
      std::vector<double> df_valcon(2 * set_count);
      if (m_processor == root) {
        for (size_t i = 0; i < set_count; i++) {
          df_valcon[2 * i + 0] = 1.0;
          df_valcon[2 * i + 1] = 1;
          if (sets[i].num_distribution_factor > 0) {
            std::vector<double> df(sets[i].num_distribution_factor);
            ex_get_set_dist_fact(filePtr, EX_NODE_SET, sets[i].id, TOPTR(df));
            double val       = df[0];
            df_valcon[2 * i] = val;
            for (int64_t j = 1; j < sets[i].num_distribution_factor; j++) {
              if (val != df[j]) {
                df_valcon[2 * i + 1] = 0;
              }
            }
          }
        }
      }

      // Tell other processors
      MPI_Bcast(TOPTR(df_valcon), df_valcon.size(), MPI_DOUBLE, root, comm_);
      for (size_t i = 0; i < set_count; i++) {
        node_sets[i].distributionFactorCount    = node_sets[i].ioss_count();
        node_sets[i].distributionFactorValue    = df_valcon[2 * i + 0];
        node_sets[i].distributionFactorConstant = (df_valcon[2 * i + 1] == 1.0);
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_sideset_data(int filePtr, size_t set_count)
  {
    m_decomposition.show_progress(__func__);
    // Issues:
    // 0. See 'get_nodeset_data' for most issues.

    int root = 0; // Root processor that reads all sideset bulk data (nodelists)

    side_sets.resize(set_count);

    std::vector<std::vector<INT>> set_elemlists(set_count);
    std::vector<ex_set>           sets(set_count);
    std::vector<INT>              ids(set_count);
    ex_get_ids(filePtr, EX_SIDE_SET, TOPTR(ids));

    for (size_t i = 0; i < set_count; i++) {
      side_sets[i].id_                 = ids[i];
      sets[i].id                       = ids[i];
      sets[i].type                     = EX_SIDE_SET;
      sets[i].entry_list               = nullptr;
      sets[i].extra_list               = nullptr;
      sets[i].distribution_factor_list = nullptr;
    }

    ex_get_sets(filePtr, sets.size(), TOPTR(sets));

    // Get total length of sideset elemlists...
    size_t elemlist_size = 0;
    for (size_t i = 0; i < set_count; i++) {
      elemlist_size += sets[i].num_entry;
      side_sets[i].fileCount = sets[i].num_entry;
    }

    // Calculate the max "buffer" size usable for storing sideset
    // elemlists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "nodeCount/2*2" is to
    // equalize the nodeCount among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = std::max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(INT));

    bool subsetting = false; // elemlist_size > max_size;

    if (subsetting) {
      assert(1 == 0);
    }
    else {
      // Can handle reading all sideset elem lists on a single
      // processor simultaneously.
      std::vector<INT> elemlist(elemlist_size);

      // Read the elemlists on root processor.
      if (m_processor == root) {
        size_t offset = 0;
        for (size_t i = 0; i < set_count; i++) {
          ex_get_set(filePtr, EX_SIDE_SET, sets[i].id, &elemlist[offset], nullptr);
          offset += sets[i].num_entry;
        }
        assert(offset == elemlist_size);
      }

      // Broadcast this data to all other processors...
      MPI_Bcast(TOPTR(elemlist), sizeof(INT) * elemlist.size(), MPI_BYTE, root, comm_);

      // Each processor now has a complete list of all elems in all
      // sidesets.
      // Determine which of these are owned by the current
      // processor...
      {
        size_t offset = 0;
        for (size_t i = 0; i < set_count; i++) {
          size_t ss_beg = offset;
          size_t ss_end = ss_beg + sets[i].num_entry;

          for (size_t n = ss_beg; n < ss_end; n++) {
            INT elem = elemlist[n];
            // See if elem owned by this processor...
            if (i_own_elem(elem)) {
              // Save elem in this processors elemlist for this set.
              // The saved data is this elems location in the global
              // elemlist for this set.
              side_sets[i].entitylist_map.push_back(n - offset);
            }
          }
          offset = ss_end;
        }
      }

      // Each processor knows how many of the sideset elems it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero elem count is
      // the "root" for this sideset.
      {
        std::vector<int> has_elems_local(set_count);
        for (size_t i = 0; i < set_count; i++) {
          has_elems_local[i] = side_sets[i].entitylist_map.empty() ? 0 : 1;
        }

        std::vector<int> has_elems(set_count * m_processorCount);
        MPI_Allgather(TOPTR(has_elems_local), has_elems_local.size(), MPI_INT, TOPTR(has_elems),
                      has_elems_local.size(), MPI_INT, comm_);

        for (size_t i = 0; i < set_count; i++) {
          side_sets[i].hasEntities.resize(m_processorCount);
          side_sets[i].root_ = m_processorCount;
          int count          = 0;
          for (int p = 0; p < m_processorCount; p++) {
            if (p < side_sets[i].root_ && has_elems[p * set_count + i] != 0) {
              side_sets[i].root_ = p;
            }
            side_sets[i].hasEntities[p] = has_elems[p * set_count + i];
            count += has_elems[p * set_count + i];
          }
	  int color = side_sets[i].hasEntities[m_processor] ? 1 : MPI_UNDEFINED;
	  MPI_Comm_split(comm_, color, m_processor, &side_sets[i].setComm_);
        }
      }

      // Check sideset distribution factors to determine whether they
      // are all constant or if they contain varying values that must
      // be communicated.  If constant or empty, then they can be
      // "read" with no communication.
      std::vector<double> df_valcon(3 * set_count);
      // df_valcon[3*i + 0] = if df constant, this is the constant value
      // df_valcon[3*i + 1] = 1 if df constant, 0 if variable
      // df_valcon[3*i + 2] = value = nodecount if all faces have same node count;
      //                    = -1 if variable
      //                      (0 if df values are constant)

      if (m_processor == root) {
        for (size_t i = 0; i < set_count; i++) {
          df_valcon[3 * i + 0] = 1.0;
          df_valcon[3 * i + 1] = 1;
          df_valcon[3 * i + 2] = 0;
          if (sets[i].num_distribution_factor > 0) {
            std::vector<double> df(sets[i].num_distribution_factor);
            ex_get_set_dist_fact(filePtr, EX_SIDE_SET, sets[i].id, TOPTR(df));
            double val       = df[0];
            df_valcon[3 * i] = val;
            for (int64_t j = 1; j < sets[i].num_distribution_factor; j++) {
              if (val != df[j]) {
                df_valcon[3 * i + 1] = 0;
                break;
              }
            }
            Ioss::Utils::clear(df);
            if (df_valcon[3 * i + 1] == 1.0) { // df are constant.
              df_valcon[3 * i + 2] = 0.0;
            }
            else {

              // To determine the size of the df field on the
              // ioss-decomp sidesets, need to know how many nodes per
              // side there are for all sides in the sideset.  Here we
              // check to see if it is a constant number to avoid
              // communicating the entire list for all sidesets.  If not
              // constant, then we will have to communicate.
              std::vector<int> nodes_per_face(side_sets[i].file_count());
              ex_get_side_set_node_count(filePtr, sets[i].id, TOPTR(nodes_per_face));
              int nod_per_face = nodes_per_face[0];
              for (size_t j = 1; j < nodes_per_face.size(); j++) {
                if (nodes_per_face[j] != nod_per_face) {
                  nod_per_face = -1;
                  break;
                }
              }
              df_valcon[3 * i + 2] = static_cast<double>(nod_per_face);
            }
          }
        }
      }

      // Tell other processors
      MPI_Bcast(TOPTR(df_valcon), df_valcon.size(), MPI_DOUBLE, root, comm_);
      for (size_t i = 0; i < set_count; i++) {
        side_sets[i].distributionFactorValue         = df_valcon[3 * i + 0];
        side_sets[i].distributionFactorConstant      = (df_valcon[3 * i + 1] == 1.0);
        side_sets[i].distributionFactorValsPerEntity = static_cast<int>(df_valcon[3 * i + 2]);
      }

      // See if need to communicate the nodes_per_side data on any
      // sidesets...  If not, then the size of those sidesets can be
      // set here...
      size_t count = 0;
      for (size_t i = 0; i < set_count; i++) {
        if (side_sets[i].distributionFactorValsPerEntity < 0) {
          count += side_sets[i].file_count();
        }
        else {
          side_sets[i].distributionFactorCount =
              side_sets[i].ioss_count() * side_sets[i].distributionFactorValsPerEntity;
        }
      }

      if (count > 0) {
        // At least 1 sideset has variable number of nodes per side...
        std::vector<int> nodes_per_face(count); // not INT
        if (m_processor == root) {
          size_t offset = 0;
          for (size_t i = 0; i < set_count; i++) {
            if (side_sets[i].distributionFactorValsPerEntity < 0) {
              ex_get_side_set_node_count(filePtr, sets[i].id, &nodes_per_face[offset]);
              offset += side_sets[i].file_count();
            }
          }
        }

        // Broadcast this data to all other processors...
        MPI_Bcast(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, root, comm_);

        // Each processor now has a list of the number of nodes per
        // face for all sidesets that have a variable number. This can
        // be used to determine the df field size on the ioss_decomp.
        size_t offset = 0;
        for (size_t i = 0; i < set_count; i++) {
          if (side_sets[i].distributionFactorValsPerEntity < 0) {
            int *npf = &nodes_per_face[offset];
            offset += side_sets[i].file_count();
            size_t my_count = 0;
            for (size_t j = 0; j < side_sets[i].ioss_count(); j++) {
              my_count += npf[side_sets[i].entitylist_map[j]];
            }
            side_sets[i].distributionFactorCount = my_count;
          }
        }
      }
    }
  }

  template <typename INT>
  int DecompositionData<INT>::get_node_coordinates(int filePtr, double *ioss_data,
                                                   const Ioss::Field &field) const
  {
    m_decomposition.show_progress(__func__);
    std::vector<double> tmp(decomp_node_count());

    int ierr = 0;
    if (field.get_name() == "mesh_model_coordinates_x") {
      m_decomposition.show_progress("\tex_get_partial_coord X");
      ierr = ex_get_partial_coord_component(filePtr, decomp_node_offset() + 1, decomp_node_count(),
                                            1, TOPTR(tmp));
      if (ierr >= 0) {
        communicate_node_data(TOPTR(tmp), ioss_data, 1);
      }
    }

    else if (field.get_name() == "mesh_model_coordinates_y") {
      m_decomposition.show_progress("\tex_get_partial_coord Y");
      ierr = ex_get_partial_coord_component(filePtr, decomp_node_offset() + 1, decomp_node_count(),
                                            2, TOPTR(tmp));
      if (ierr >= 0) {
        communicate_node_data(TOPTR(tmp), ioss_data, 1);
      }
    }

    else if (field.get_name() == "mesh_model_coordinates_z") {
      m_decomposition.show_progress("\tex_get_partial_coord Z");
      ierr = ex_get_partial_coord_component(filePtr, decomp_node_offset() + 1, decomp_node_count(),
                                            3, TOPTR(tmp));
      if (ierr >= 0) {
        communicate_node_data(TOPTR(tmp), ioss_data, 1);
      }
    }

    else if (field.get_name() == "mesh_model_coordinates") {
      // Data required by upper classes store x0, y0, z0, ... xn,
      // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
      // ..., yn, z0, ..., zn so we have to allocate some scratch
      // memory to read in the data and then map into supplied
      // 'data'

      std::vector<double> ioss_tmp(ioss_node_count());

      // This implementation trades off extra communication for
      // reduced memory overhead.
      // * This method uses 'ioss_node_count' extra memory; 3
      // reads; and 3 communicate_node_data calls.
      //
      // * Other method uses 6*ioss_node_count extra memory; 1 read;
      // and 1 communicate_node_data call.
      //
      // * NOTE: The read difference is not real since the ex_get_partial_coord
      // function does 3 reads internally.

      for (int d = 0; d < m_decomposition.m_spatialDimension; d++) {
        m_decomposition.show_progress("\tex_get_partial_coord XYZ");
        ierr = ex_get_partial_coord_component(filePtr, decomp_node_offset() + 1,
                                              decomp_node_count(), d + 1, tmp.data());
        if (ierr < 0) {
          return ierr;
        }

        communicate_node_data(TOPTR(tmp), TOPTR(ioss_tmp), 1);

        size_t index = d;
        for (size_t i = 0; i < ioss_node_count(); i++) {
          ioss_data[index] = ioss_tmp[i];
          index += m_decomposition.m_spatialDimension;
        }
      }
    }
    return ierr;
  }

  template void DecompositionData<int>::get_block_connectivity(int filePtr, int *data, int64_t id,
                                                               size_t blk_seq, size_t nnpe) const;
  template void DecompositionData<int64_t>::get_block_connectivity(int filePtr, int64_t *data,
                                                                   int64_t id, size_t blk_seq,
                                                                   size_t nnpe) const;

  template <typename INT>
  void DecompositionData<INT>::get_block_connectivity(int filePtr, INT *data, int64_t id,
                                                      size_t blk_seq, size_t nnpe) const
  {
    m_decomposition.show_progress(__func__);
    Ioss::BlockDecompositionData blk = el_blocks[blk_seq];

    // Determine number of file decomp elements are in this block and the offset into the block.
    size_t count  = get_block_element_count(blk_seq);
    size_t offset = get_block_element_offset(blk_seq);

    std::vector<INT> file_conn(count * nnpe);
    m_decomposition.show_progress("\tex_get_partial_conn");
    ex_get_partial_conn(filePtr, EX_ELEM_BLOCK, id, offset + 1, count, TOPTR(file_conn), nullptr,
                        nullptr);
    m_decomposition.communicate_block_data(TOPTR(file_conn), data, blk, nnpe);

    for (size_t i = 0; i < blk.iossCount * nnpe; i++) {
      data[i] = node_global_to_local(data[i]);
    }
  }

  template <typename INT>
  int DecompositionData<INT>::get_var(int filePtr, int step, ex_entity_type type, int var_index,
                                      ex_entity_id id, int64_t num_entity,
                                      std::vector<double> &data) const
  {
    m_decomposition.show_progress(__func__);
    if (type == EX_ELEM_BLOCK) {
      return get_elem_var(filePtr, step, var_index, id, num_entity, data);
    }
    else if (type == EX_NODAL) {
      return get_node_var(filePtr, step, var_index, id, num_entity, data);
    }
    else if (type == EX_NODE_SET || type == EX_SIDE_SET) {
      return get_set_var(filePtr, step, var_index, type, id, num_entity, data);
    }
    else {
      assert(1 == 0);
      return -1;
    }
  }

  template <typename INT>
  int DecompositionData<INT>::get_attr(int filePtr, ex_entity_type obj_type, ex_entity_id id,
                                       size_t attr_count, double *attrib) const
  {
    m_decomposition.show_progress(__func__);
    if (attr_count == 1) {
      return get_one_attr(filePtr, obj_type, id, 1, attrib);
    }

    if (obj_type == EX_ELEM_BLOCK) {
      return get_elem_attr(filePtr, id, attr_count, attrib);
    }
    else if (obj_type == EX_NODAL) {
      return get_node_attr(filePtr, id, attr_count, attrib);
    }
    else if (obj_type == EX_NODE_SET || obj_type == EX_SIDE_SET) {
      return get_set_attr(filePtr, obj_type, id, attr_count, attrib);
    }
    else {
      assert(1 == 0);
      return -1;
    }
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_attr(int filePtr, ex_entity_type obj_type, ex_entity_id id,
                                           int attrib_index, double *attrib) const
  {
    m_decomposition.show_progress(__func__);
    if (obj_type == EX_ELEM_BLOCK) {
      return get_one_elem_attr(filePtr, id, attrib_index, attrib);
    }
    else if (obj_type == EX_NODAL) {
      return get_one_node_attr(filePtr, id, attrib_index, attrib);
    }
    else if (obj_type == EX_NODE_SET || obj_type == EX_SIDE_SET) {
      return get_one_set_attr(filePtr, obj_type, id, attrib_index, attrib);
    }
    else {
      assert(1 == 0);
      return -1;
    }
  }

  template void DecompositionDataBase::communicate_node_data(int *file_data, int *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(int64_t *file_data, int64_t *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(double *file_data, double *ioss_data,
                                                             size_t comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_node_data(T *file_data, T *ioss_data,
                                                    size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_node_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_node_data(file_data, ioss_data, comp_count);
    }
  }

  template void DecompositionDataBase::communicate_element_data(int *file_data, int *ioss_data,
                                                                size_t comp_count) const;
  template void DecompositionDataBase::communicate_element_data(int64_t *file_data,
                                                                int64_t *ioss_data,
                                                                size_t   comp_count) const;
  template void DecompositionDataBase::communicate_element_data(double *file_data,
                                                                double *ioss_data,
                                                                size_t  comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_element_data(T *file_data, T *ioss_data,
                                                       size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_element_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_element_data(file_data, ioss_data, comp_count);
    }
  }

  void DecompositionDataBase::get_node_entity_proc_data(void *                    entity_proc,
                                                        const Ioss::MapContainer &node_map,
                                                        bool                      do_map) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->m_decomposition.get_node_entity_proc_data(reinterpret_cast<int *>(entity_proc),
                                                        node_map, do_map);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->m_decomposition.get_node_entity_proc_data(reinterpret_cast<int64_t *>(entity_proc),
                                                        node_map, do_map);
    }
  }

  int DecompositionDataBase::get_set_mesh_double(int filePtr, ex_entity_type type, ex_entity_id id,
                                                 const Ioss::Field &field, double *ioss_data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      return this32->get_set_mesh_var(filePtr, type, id, field, ioss_data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      return this64->get_set_mesh_var(filePtr, type, id, field, ioss_data);
    }
  }

  int DecompositionDataBase::get_set_mesh_var(int filePtr, ex_entity_type type, ex_entity_id id,
                                              const Ioss::Field &field, void *ioss_data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      return this32->get_set_mesh_var(filePtr, type, id, field, (int *)ioss_data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      return this64->get_set_mesh_var(filePtr, type, id, field, (int64_t *)ioss_data);
    }
  }

  void DecompositionDataBase::get_block_connectivity(int filePtr, void *data, int64_t id,
                                                     size_t blk_seq, size_t nnpe) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_block_connectivity(filePtr, reinterpret_cast<int *>(data), id, blk_seq, nnpe);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_block_connectivity(filePtr, reinterpret_cast<int64_t *>(data), id, blk_seq, nnpe);
    }
  }

  const Ioss::SetDecompositionData &DecompositionDataBase::get_decomp_set(ex_entity_type type,
                                                                          ex_entity_id   id) const
  {
    if (type == EX_NODE_SET) {
      for (const auto &node_set : node_sets) {
        if (node_set.id_ == id) {
          return node_set;
        }
      }
    }
    else if (type == EX_SIDE_SET) {
      for (const auto &side_set : side_sets) {
        if (side_set.id_ == id) {
          return side_set;
        }
      }
    }

    std::ostringstream errmsg;
    if (type != EX_NODE_SET && type != EX_SIDE_SET) {
      errmsg << "ERROR: Invalid set type specified in get_decomp_set. Only node set or side set "
                "supported\n";
    }
    else {
      std::string typestr = type == EX_NODE_SET ? "node set" : "side set";
      errmsg << "ERROR: Count not find " << typestr << " " << id << "\n";
    }
    std::cerr << errmsg.str();
    exit(EXIT_FAILURE);
    return node_sets[0];
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_seq(ex_entity_type type, ex_entity_id id) const
  {
    m_decomposition.show_progress(__func__);
    if (type == EX_ELEM_BLOCK) {
      for (size_t i = 0; i < el_blocks.size(); i++) {
        if (el_blocks[i].id_ == id) {
          return i;
        }
      }
    }
    return el_blocks.size();
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_element_count(size_t blk_seq) const
  {
    m_decomposition.show_progress(__func__);
    // Determine number of file decomp elements are in this block;
    size_t bbeg  = std::max(m_decomposition.m_fileBlockIndex[blk_seq], decomp_elem_offset());
    size_t bend  = std::min(m_decomposition.m_fileBlockIndex[blk_seq + 1],
                           decomp_elem_offset() + decomp_elem_count());
    size_t count = 0;
    if (bend > bbeg) {
      count = bend - bbeg;
    }
    return count;
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_element_offset(size_t blk_seq) const
  {
    m_decomposition.show_progress(__func__);
    size_t offset = 0;
    if (decomp_elem_offset() > m_decomposition.m_fileBlockIndex[blk_seq]) {
      offset = decomp_elem_offset() - m_decomposition.m_fileBlockIndex[blk_seq];
    }
    return offset;
  }

  template <typename INT>
  int DecompositionData<INT>::get_set_var(int filePtr, int step, int var_index, ex_entity_type type,
                                          ex_entity_id         id, int64_t /*num_entity*/,
                                          std::vector<double> &ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find set corresponding to the specified id...
    auto &set = get_decomp_set(type, id);

    std::vector<double> file_data;
    int                 ierr = 0;
    if (m_processor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count());
      m_decomposition.show_progress("\tex_get_var (set)");
      ierr = ex_get_var(filePtr, step, type, var_index, id, set.file_count(), TOPTR(file_data));
    }

    if (ierr >= 0) {
      communicate_set_data(TOPTR(file_data), TOPTR(ioss_data), set, 1);
    }

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_set_attr(int filePtr, ex_entity_type type, ex_entity_id id,
                                           size_t comp_count, double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find set corresponding to the specified id...
    auto &set = get_decomp_set(type, id);

    std::vector<double> file_data;
    int                 ierr = 0;
    if (m_processor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count() * comp_count);
      ierr = ex_get_attr(filePtr, type, id, TOPTR(file_data));
    }

    if (ierr >= 0) {
      communicate_set_data(TOPTR(file_data), ioss_data, set, comp_count);
    }

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_set_attr(int filePtr, ex_entity_type type, ex_entity_id id,
                                               int attr_index, double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find set corresponding to the specified id...
    auto &set = get_decomp_set(type, id);

    std::vector<double> file_data;
    int                 ierr = 0;
    if (m_processor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count());
      ierr = ex_get_one_attr(filePtr, type, id, attr_index, TOPTR(file_data));
    }

    if (ierr >= 0) {
      communicate_set_data(TOPTR(file_data), ioss_data, set, 1);
    }

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_node_var(int filePtr, int step, int var_index, ex_entity_id id,
                                           int64_t /*num_entity*/,
                                           std::vector<double> &ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    std::vector<double> file_data(decomp_node_count());
    m_decomposition.show_progress("\tex_get_partial_var");
    int ierr = ex_get_partial_var(filePtr, step, EX_NODAL, var_index, id, decomp_node_offset() + 1,
                                  decomp_node_count(), TOPTR(file_data));

    if (ierr >= 0) {
      communicate_node_data(TOPTR(file_data), TOPTR(ioss_data), 1);
    }
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_node_attr(int filePtr, ex_entity_id id, size_t comp_count,
                                            double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    std::vector<double> file_data(decomp_node_count() * comp_count);
    int                 ierr = ex_get_partial_attr(filePtr, EX_NODAL, id, decomp_node_offset() + 1,
                                   decomp_node_count(), TOPTR(file_data));

    if (ierr >= 0) {
      communicate_node_data(TOPTR(file_data), ioss_data, comp_count);
    }
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_node_attr(int filePtr, ex_entity_id id, int attr_index,
                                                double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    std::vector<double> file_data(decomp_node_count());
    int ierr = ex_get_partial_one_attr(filePtr, EX_NODAL, id, decomp_node_offset() + 1,
                                       decomp_node_count(), attr_index, TOPTR(file_data));

    if (ierr >= 0) {
      communicate_node_data(TOPTR(file_data), ioss_data, 1);
    }
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_elem_var(int filePtr, int step, int var_index, ex_entity_id id,
                                           int64_t /*num_entity*/,
                                           std::vector<double> &ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count   = get_block_element_count(blk_seq);
    size_t offset  = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count);
    m_decomposition.show_progress("\tex_get_partial_var (elem)");
    int ierr = ex_get_partial_var(filePtr, step, EX_ELEM_BLOCK, var_index, id, offset + 1, count,
                                  TOPTR(file_data));

    if (ierr >= 0) {
      m_decomposition.communicate_block_data(TOPTR(file_data), TOPTR(ioss_data), el_blocks[blk_seq],
                                             1);
    }

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_elem_attr(int filePtr, ex_entity_id id, size_t comp_count,
                                            double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count   = get_block_element_count(blk_seq);
    size_t offset  = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count * comp_count);
    int ierr = ex_get_partial_attr(filePtr, EX_ELEM_BLOCK, id, offset + 1, count, TOPTR(file_data));

    if (ierr >= 0) {
      m_decomposition.communicate_block_data(TOPTR(file_data), ioss_data, el_blocks[blk_seq],
                                             comp_count);
    }

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_elem_attr(int filePtr, ex_entity_id id, int attr_index,
                                                double *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count   = get_block_element_count(blk_seq);
    size_t offset  = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count);
    int ierr = ex_get_partial_one_attr(filePtr, EX_ELEM_BLOCK, id, offset + 1, count, attr_index,
                                       TOPTR(file_data));

    if (ierr >= 0) {
      m_decomposition.communicate_block_data(TOPTR(file_data), ioss_data, el_blocks[blk_seq], 1);
    }

    return ierr;
  }

  template int DecompositionData<int>::get_set_mesh_var(int filePtr, ex_entity_type type,
                                                        ex_entity_id id, const Ioss::Field &field,
                                                        int *ioss_data) const;
  template int DecompositionData<int64_t>::get_set_mesh_var(int filePtr, ex_entity_type type,
                                                            ex_entity_id       id,
                                                            const Ioss::Field &field,
                                                            int64_t *          ioss_data) const;
  template int DecompositionData<int>::get_set_mesh_var(int filePtr, ex_entity_type type,
                                                        ex_entity_id id, const Ioss::Field &field,
                                                        double *ioss_data) const;
  template int DecompositionData<int64_t>::get_set_mesh_var(int filePtr, ex_entity_type type,
                                                            ex_entity_id       id,
                                                            const Ioss::Field &field,
                                                            double *           ioss_data) const;

  template <typename INT>
  template <typename T>
  int DecompositionData<INT>::get_set_mesh_var(int filePtr, ex_entity_type type, ex_entity_id id,
                                               const Ioss::Field &field, T *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    // Sideset Distribution Factor data can be very complicated.
    // For some sanity, handle all requests for those in a separate routine...
    if (type == EX_SIDE_SET && field.get_name() == "distribution_factors") {
      return handle_sset_df(filePtr, id, field, ioss_data);
    }

    auto &set = get_decomp_set(type, id);

    std::vector<T> file_data;
    int            ierr = 0;

    // These fields call back into this routine and are handled on all
    // processors; not just the root processor.
    if (field.get_name() == "element_side") {
      // Sideset only...
      if (type == EX_SIDE_SET) {
        // Interleave the "ids" and "sides" fields...
        std::vector<T> tmp(set.ioss_count());
        Ioss::Field    elem_field("ids", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH,
                               tmp.size());
        get_set_mesh_var(filePtr, type, id, elem_field, TOPTR(tmp));
        for (size_t i = 0; i < tmp.size(); i++) {
          ioss_data[2 * i] = tmp[i];
        }
        Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH,
                               tmp.size());
        get_set_mesh_var(filePtr, type, id, side_field, TOPTR(tmp));
        for (size_t i = 0; i < tmp.size(); i++) {
          ioss_data[2 * i + 1] = tmp[i];
        }
      }
      else {
        return -1;
      }
      return ierr;
    }
    else if (field.get_name() == "element_side_raw") {
      // Sideset only...
      if (type == EX_SIDE_SET) {
        // Interleave the "ids" and "sides" fields...
        std::vector<T> tmp(set.ioss_count());
        Ioss::Field    elem_field("ids_raw", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH,
                               tmp.size());
        get_set_mesh_var(filePtr, type, id, elem_field, TOPTR(tmp));
        for (size_t i = 0; i < tmp.size(); i++) {
          ioss_data[2 * i] = tmp[i];
        }
        Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH,
                               tmp.size());
        get_set_mesh_var(filePtr, type, id, side_field, TOPTR(tmp));
        for (size_t i = 0; i < tmp.size(); i++) {
          ioss_data[2 * i + 1] = tmp[i];
        }
      }
      else {
        return -1;
      }
      return ierr;
    }

    // If the requested field is "distribution_factors" see if
    // they are constant and the read/comm can be skipped...
    if (field.get_name() == "distribution_factors" && set.distributionFactorConstant) {
      // Fill in the ioss decomp with the constant value...
      for (size_t i = 0; i < set.distributionFactorCount; i++) {
        ioss_data[i] = set.distributionFactorValue;
      }
      return 0;
    }

    if (m_processor == set.root_) {
      // Read the nodeset data from the file..
      if (field.get_name() == "ids" || field.get_name() == "ids_raw") {
        file_data.resize(set.file_count());
        ierr = ex_get_set(filePtr, type, id, TOPTR(file_data), nullptr);
      }
      else if (field.get_name() == "sides") {
        // Sideset only...
        if (type == EX_SIDE_SET) {
          file_data.resize(set.file_count());
          ierr = ex_get_set(filePtr, type, id, nullptr, TOPTR(file_data));
        }
        else {
          return -1;
        }
      }
      else if (field.get_name() == "distribution_factors") {
        ex_set set_param[1];
        set_param[0].id                       = id;
        set_param[0].type                     = type;
        set_param[0].entry_list               = nullptr;
        set_param[0].extra_list               = nullptr;
        set_param[0].distribution_factor_list = nullptr;
        ierr                                  = ex_get_sets(filePtr, 1, set_param);

        if (set_param[0].num_distribution_factor == 0) {
          // This should have been caught above.
          assert(1 == 0 && "Internal error in get_set_mesh_var");
        }
        else {
          if (type == EX_NODE_SET) {
            file_data.resize(set_param[0].num_distribution_factor);
            set_param[0].distribution_factor_list = TOPTR(file_data);
            ierr                                  = ex_get_sets(filePtr, 1, set_param);
          }
          else {
            assert(1 == 0 && "Internal error -- should not be here -- sset df");
          }
        }
      }
      else {
        assert(1 == 0 && "Unrecognized field name in get_set_mesh_var");
      }
    }

    if (ierr >= 0) {
      communicate_set_data(TOPTR(file_data), ioss_data, set, 1);
    }

    // Map global 0-based index to local 1-based index.
    if (field.get_name() == "ids" || field.get_name() == "ids_raw") {
      if (type == EX_NODE_SET) {
        for (size_t i = 0; i < set.ioss_count(); i++) {
          ioss_data[i] = node_global_to_local(ioss_data[i]);
        }
      }
      else if (type == EX_SIDE_SET) {
        for (size_t i = 0; i < set.ioss_count(); i++) {
          ioss_data[i] = elem_global_to_local(ioss_data[i]);
        }
      }
      else {
        assert(1 == 0);
      }
    }
    return ierr;
  }

  template <typename INT>
  template <typename T>
  int DecompositionData<INT>::handle_sset_df(int filePtr, ex_entity_id id, const Ioss::Field &field,
                                             T *ioss_data) const
  {
    m_decomposition.show_progress(__func__);
    int ierr = 0;

    // Sideset Distribution Factor data can be very complicated.
    // For some sanity, handle all requests for those here.  Only handles sidesets
    // distribution_factors field.
    assert(field.get_name() == "distribution_factors");

    auto &set = get_decomp_set(EX_SIDE_SET, id);

    // See if df are constant and the read/comm can be skipped...
    if (set.distributionFactorConstant) {
      // Fill in the ioss decomp with the constant value...
      for (size_t i = 0; i < set.distributionFactorCount; i++) {
        ioss_data[i] = set.distributionFactorValue;
      }
      return 0;
    }

    // See if this set only exists on a single processor.
    //    In that case, the file_data is the same as the ioss_data
    //    and we can read the data directly into ioss_data and return...
    size_t proc_active = std::accumulate(set.hasEntities.begin(), set.hasEntities.end(), 0);
    if (proc_active == 1) {
      if (m_processor == set.root_) {
        ex_set set_param[1];
        set_param[0].id                       = id;
        set_param[0].type                     = EX_SIDE_SET;
        set_param[0].entry_list               = nullptr;
        set_param[0].extra_list               = nullptr;
        set_param[0].distribution_factor_list = nullptr;
        ex_get_sets(filePtr, 1, set_param);
        if (set_param[0].num_distribution_factor == 0) {
          // This should have been caught above.
          assert(1 == 0 && "Internal error in handle_sset_df");
        }
        else {
          // Read data directly into ioss_data.
          set_param[0].distribution_factor_list = ioss_data;
          ex_get_sets(filePtr, 1, set_param);
        }
      }
      return 0;
    }

    // At this point, we have nonconstant distribution factors on
    // a sideset split among 2 or more processors...
    // Two alternatives exist:
    // 1. Constant face topology in sideset (e.g., all quad or all tri) [EASY, COMMON]
    // 2. Non-constant face topology in sideset (e.g., mix of quad/tri/...) [HARD, RARE?]

    if (set.distributionFactorValsPerEntity > 0) {
      // Constant face topology in sideset
      // Simply read the values in the file decomposition and
      // communicate with a comp count of set.distributionFactorValsPerEntity.
      std::vector<T> file_data;
      if (m_processor == set.root_) {
        file_data.resize(set.distributionFactorValsPerEntity * set.fileCount);

        ex_set set_param[1];
        set_param[0].id                       = id;
        set_param[0].type                     = EX_SIDE_SET;
        set_param[0].entry_list               = nullptr;
        set_param[0].extra_list               = nullptr;
        set_param[0].distribution_factor_list = TOPTR(file_data);
        ierr                                  = ex_get_sets(filePtr, 1, set_param);
      }
      if (ierr >= 0) {
        communicate_set_data(TOPTR(file_data), ioss_data, set, set.distributionFactorValsPerEntity);
      }

      return ierr;
    }

    // non-constant face topology in sideset, non-constant df on faces.
    // Get total number of df on file for this sset...
    size_t df_count = 0;
    if (m_processor == set.root_) {
      ex_set set_param[1];
      set_param[0].id                       = id;
      set_param[0].type                     = EX_SIDE_SET;
      set_param[0].entry_list               = nullptr;
      set_param[0].extra_list               = nullptr;
      set_param[0].distribution_factor_list = nullptr;
      ex_get_sets(filePtr, 1, set_param);
      df_count = set_param[0].num_distribution_factor;
    }

    // Get the node-count-per-face for all faces in this set...
    std::vector<int> nodes_per_face(set.file_count() + 1);
    if (m_processor == set.root_) {
      ex_get_side_set_node_count(filePtr, set.id_, TOPTR(nodes_per_face));
      nodes_per_face[set.file_count()] = df_count;
    }

    // Send this data to the other processors

    // NOTE That a processor either sends or receives, but never both,
    // so this will not cause a deadlock...
    if (m_processor != set.root_ && set.hasEntities[m_processor]) {
      MPI_Status status{};
      int result = MPI_Recv(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, set.root_, 222,
                            comm_, &status);

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Recv error on processor " << m_processor
               << " receiving nodes_per_face sideset data";
        std::cerr << errmsg.str();
      }
      df_count = nodes_per_face.back();
    }

    if (set.root_ == m_processor) {
      // Sending data to other processors...
      for (int i = m_processor + 1; i < m_processorCount; i++) {
        if (set.hasEntities[i]) {
          // Send same data to all active processors...
          MPI_Send(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, i, 222, comm_);
        }
      }
    }

    // Now, read the df on the root processor and send it to the other active
    // processors for this set...
    std::vector<double> file_data;
    if (m_processor == set.root_) {
      file_data.resize(df_count);

      ex_set set_param[1];
      set_param[0].id                       = id;
      set_param[0].type                     = EX_SIDE_SET;
      set_param[0].entry_list               = nullptr;
      set_param[0].extra_list               = nullptr;
      set_param[0].distribution_factor_list = TOPTR(file_data);
      ex_get_sets(filePtr, 1, set_param);
    }

    // Send this data to the other processors

    if (m_processor != set.root_ && set.hasEntities[m_processor]) {
      file_data.resize(df_count);
      MPI_Status status{};
      int        result =
          MPI_Recv(TOPTR(file_data), file_data.size(), MPI_DOUBLE, set.root_, 333, comm_, &status);

      if (result != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "ERROR: MPI_Recv error on processor " << m_processor
               << " receiving nodes_per_face sideset data";
        std::cerr << errmsg.str();
      }
    }

    if (set.root_ == m_processor) {
      // Sending data to other processors...
      for (int i = m_processor + 1; i < m_processorCount; i++) {
        if (set.hasEntities[i]) {
          // Send same data to all active processors...
          MPI_Send(TOPTR(file_data), file_data.size(), MPI_DOUBLE, i, 333, comm_);
        }
      }
    }

    // Now, each active processor for this set needs to step through the df
    // data in file_data and transfer the data it owns to ioss_data.
    if (set.hasEntities[m_processor]) {
      // Convert nodes_per_face into an offset into the df array...
      Ioss::Utils::generate_index(nodes_per_face);

      size_t k = 0;
      for (size_t i = 0; i < set.ioss_count(); i++) {
        size_t index = set.entitylist_map[i];
        size_t beg   = nodes_per_face[index];
        size_t end   = nodes_per_face[index + 1];
        for (size_t j = beg; j < end; j++) {
          ioss_data[k++] = file_data[j];
        }
      }
    }
    return 0;
  }

  template void DecompositionData<int>::create_implicit_global_map(
      const std::vector<int> &owning_proc, std::vector<int64_t> &global_implicit_map,
      Ioss::Map &node_map, int64_t *locally_owned_count, int64_t *processor_offset);
  template void DecompositionData<int64_t>::create_implicit_global_map(
      const std::vector<int> &owning_proc, std::vector<int64_t> &global_implicit_map,
      Ioss::Map &node_map, int64_t *locally_owned_count, int64_t *processor_offset);

  template <typename INT>
  void DecompositionData<INT>::create_implicit_global_map(const std::vector<int> &owning_proc,
                                                          std::vector<int64_t> &global_implicit_map,
                                                          Ioss::Map &           node_map,
                                                          int64_t *             locally_owned_count,
                                                          int64_t *             processor_offset)
  {
    m_decomposition.show_progress(__func__);
    // Used on composed output database...
    // If the node is locally owned, then its position is basically
    // determined by removing all shared nodes from the list and
    // then compressing the list. This location plus the proc_offset
    // gives its location in the global-implicit file.
    //
    // If it is shared, then need to communicate with the owning
    // processor to determine where that processor is putting it.
    //
    // First, iterate nodeOwningProcessor list and set the implicit
    // map for all locally-owned nodes and also determine how many
    // of my nodes are owned by which other processors.

    global_implicit_map.resize(owning_proc.size());

    std::vector<int64_t> snd_count(m_processorCount);
    std::vector<int64_t> rcv_count(m_processorCount);

    size_t position = 0;
    for (size_t i = 0; i < global_implicit_map.size(); i++) {
      snd_count[owning_proc[i]]++;
      if (owning_proc[i] == m_processor) {
        global_implicit_map[i] = position++;
      }
    }
    snd_count[m_processor] = 0;

    // The number of locally-owned nodes on this processor is 'position'
    *locally_owned_count = position;

    MPI_Allgather(locally_owned_count, 1, MPI_LONG_LONG_INT, &rcv_count[0], 1, MPI_LONG_LONG_INT,
                  comm_);
    m_decomposition.show_progress("\tAllgather finished");

    // Determine the offset of the nodes on this processor. The offset is the
    // total number of locally-owned nodes on all processors prior to this processor.
    *processor_offset = 0;
    for (int i = 0; i < m_processor; i++) {
      *processor_offset += rcv_count[i];
    }

    for (auto &i : global_implicit_map) {
      i += *processor_offset + 1;
    }

    // Now, tell the other processors how many nodes I will be sending
    // them (Nodes they own that I share with them)
    MPI_Alltoall(TOPTR(snd_count), 1, MPI_LONG_LONG_INT, TOPTR(rcv_count), 1, MPI_LONG_LONG_INT,
                 comm_);
    m_decomposition.show_progress("\tCommunication 1 finished");

    std::vector<int64_t> snd_offset(snd_count);
    Ioss::Utils::generate_index(snd_offset);
    std::vector<int64_t> snd_list(*snd_offset.rbegin() + *snd_count.rbegin());

    {
      std::vector<int64_t> tmp_disp(snd_offset);
      // Now create the list of nodes to send...
      for (size_t i = 0; i < global_implicit_map.size(); i++) {
        if (owning_proc[i] != m_processor) {
          int64_t global_id                    = node_map.map()[i + 1];
          snd_list[tmp_disp[owning_proc[i]]++] = global_id;
        }
      }
    }

    std::vector<int64_t> rcv_offset(rcv_count);
    Ioss::Utils::generate_index(rcv_offset);
    std::vector<int64_t> rcv_list(*rcv_offset.rbegin() + *rcv_count.rbegin());

    Ioss::MY_Alltoallv(snd_list, snd_count, snd_offset, rcv_list, rcv_count, rcv_offset, comm_);
    m_decomposition.show_progress("\tCommunication 2 finished");

    // Iterate rcv_list and convert global ids to the global-implicit position...
    for (auto &i : rcv_list) {
      int64_t local_id     = node_map.global_to_local(i) - 1;
      int64_t rcv_position = global_implicit_map[local_id];
      i                    = rcv_position;
    }

    // Send the data back now...
    Ioss::MY_Alltoallv(rcv_list, rcv_count, rcv_offset, snd_list, snd_count, snd_offset, comm_);
    m_decomposition.show_progress("\tCommunication 3 finished");

    // Fill in the remaining portions of the global_implicit_map...
    std::vector<int64_t> tmp_disp(snd_offset);
    for (size_t i = 0; i < global_implicit_map.size(); i++) {
      if (owning_proc[i] != m_processor) {
        int64_t implicit       = snd_list[tmp_disp[owning_proc[i]]++];
        global_implicit_map[i] = implicit;
      }
    }
  }
} // namespace Iopx
