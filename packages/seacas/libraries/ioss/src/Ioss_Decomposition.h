/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef IOSS_DECOMPOSITON_H
#define IOSS_DECOMPOSITON_H

#include <Ioss_CodeTypes.h>
#include <Ioss_Map.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_PropertyManager.h>
#include <algorithm>
#include <assert.h>
#include <map>
#include <string>
#include <vector>

#if !defined(NO_PARMETIS_SUPPORT)
#include <parmetis.h>
#endif

#if !defined(NO_ZOLTAN_SUPPORT)
#undef MPICPP
#include <zoltan_cpp.h>
#endif

namespace Ioss {
  class BlockDecompositionData
  {
  public:
    BlockDecompositionData()
        : zone_(0), section_(0), fileSectionOffset(0), id_(0), fileCount(0), iossCount(0),
          globalCount(0), zoneNodeOffset(0), topologyType("unknown"), nodesPerEntity(0),
          attributeCount(0), localIossOffset(0)
    {
    }

    const std::string &name() const { return name_; }
    int                zone() const { return zone_; }
    int                section() const { return section_; }
    int64_t            id() const { return id_; }
    size_t             file_count() const { return fileCount; }
    size_t             ioss_count() const { return iossCount; }
    size_t             global_count() const { return globalCount; }

    std::string name_;
    int         zone_;
    int         section_;

    size_t  fileSectionOffset; // In partial read, where start
    int64_t id_;
    size_t  fileCount;
    size_t  iossCount;
    size_t  globalCount;

    size_t      zoneNodeOffset;
    std::string topologyType;
    int         nodesPerEntity;
    int         attributeCount;

    // maps from file-block data to ioss-block data
    // The local_map.size() elements starting at localIossOffset are local.
    // ioss[localIossOffset+i] = file[local_map[i]];
    size_t           localIossOffset;
    std::vector<int> localMap;

    // Maps from file-block data to export list.
    // export[i] = file[export_map[i]
    std::vector<int> exportMap;
    std::vector<int> exportCount;
    std::vector<int> exportIndex;

    // Maps from import data to ioss-block data.
    // ioss[import_map[i] = local_map[i];
    std::vector<int> importMap;
    std::vector<int> importCount;
    std::vector<int> importIndex;
  };

  class SetDecompositionData
  {
  public:
    SetDecompositionData()
        : id_(0), zone_(0), section_(0), fileCount(0), root_(0), parentBlockIndex(0),
          distributionFactorValsPerEntity(-1), distributionFactorCount(0),
          distributionFactorValue(0.0), distributionFactorConstant(false)
    {
    }

    const std::string &name() const { return name_; }
    int64_t            id() const { return id_; }
    int                zone() const { return zone_; }
    int                section() const { return section_; }
    size_t             file_count() const { return fileCount; }
    size_t             ioss_count() const { return entitylist_map.size(); }
    size_t             df_count() const { return distributionFactorCount; }

    // contains global entity-list positions for all entities in this set on this processor.
    std::vector<size_t> entitylist_map;
    std::vector<bool>   hasEntities; // T/F if this set exists on processor p

    std::string name_;
    int64_t     id_;
    int         zone_;
    int         section_;
    size_t      fileCount; // Number of nodes in nodelist for file decomposition
    int         root_;     // Lowest number processor that has nodes for this nodest
    std::string topologyType;
    size_t      parentBlockIndex;

    int    distributionFactorValsPerEntity; // number of df / element or node. -1 if nonconstant.
    size_t distributionFactorCount;
    double distributionFactorValue;    // If distributionFactorConstant == true, the constant value
    bool   distributionFactorConstant; // T if all distribution factors the same value.
  };

  template <typename INT> class Decomposition
  {
  public:
    Decomposition(const Ioss::PropertyManager &props, MPI_Comm comm);

    size_t global_node_count() const { return m_globalNodeCount; }
    size_t global_elem_count() const { return m_globalElementCount; }
    size_t ioss_node_count() const { return nodeGTL.size(); }
    size_t ioss_elem_count() const { return localElementMap.size() + importElementMap.size(); }
    size_t file_node_count() const { return m_nodeCount; }
    size_t file_elem_count() const { return m_elementCount; }
    size_t file_node_offset() const { return m_nodeOffset; }
    size_t file_elem_offset() const { return m_elementOffset; }

    bool needs_centroids() const;

    void generate_entity_distributions(size_t globalNodeCount, size_t globalElementCount);

    // T/F if node with global index node owned by this processors ioss-decomp.
    bool i_own_node(size_t global_index) const
    {
      // global_index is 1-based index into global list of nodes [1..global_node_count]
      return std::binary_search(nodeGTL.begin(), nodeGTL.end(), global_index);
    }

    // T/F if element with global index elem owned by this processors ioss-decomp.
    bool i_own_elem(size_t global_index) const
    {
      // global_index is 1-based index into global list of elements [1..global_element_count]
      return elemGTL.count(global_index) != 0;
    }

    size_t node_global_to_local(size_t global_index) const
    {
      // global_index is 1-based index into global list of nodes [1..global_node_count]
      // return value is 1-based index into local list of nodes on this
      // processor (ioss-decomposition)
      // Note that for 'int', equivalence and equality are the same, so
      // lower_bound is OK here (EffectiveSTL, Item 19)
      typename std::vector<INT>::const_iterator I =
          lower_bound(nodeGTL.begin(), nodeGTL.end(), global_index);
      assert(I != nodeGTL.end());
      return std::distance(nodeGTL.begin(), I) + 1; // Convert to 1-based index.
    }

    size_t elem_global_to_local(size_t global_index) const
    {
      // global_index is 1-based index into global list of elements [1..global_node_count]
      // return value is 1-based index into local list of elements on this
      // processor (ioss-decomposition)
      typename std::map<INT, INT>::const_iterator I = elemGTL.find(global_index);
      assert(I != elemGTL.end());
      return I->second;
    }

    void show_progress(const std::string &message) const
    {
      if (m_showProgress) {
        // Use the output below for debugging...
        // std::cerr << "[" << m_processor << "].... " << message << "\n";
        Ioss::ParallelUtils pu(m_comm);
        pu.progress(message);
      }
    }

    void decompose_model(
#if !defined(NO_ZOLTAN_SUPPORT)
        Zoltan &zz,
#endif
        std::vector<BlockDecompositionData> &element_blocks);

    void simple_decompose();

    void simple_node_decompose();

    void calculate_element_centroids(const std::vector<double> &x, const std::vector<double> &y,
                                     const std::vector<double> &z);

#if !defined(NO_ZOLTAN_SUPPORT)
    void zoltan_decompose(Zoltan &zz);

    void get_local_element_list(const ZOLTAN_ID_PTR &export_global_ids, size_t export_count);
#endif

#if !defined(NO_PARMETIS_SUPPORT)
    void metis_decompose(idx_t *pointer, idx_t *adjacency,
                         std::vector<BlockDecompositionData> &el_blocks);

    void internal_metis_decompose(std::vector<BlockDecompositionData> &el_blocks,
                                  idx_t *element_dist, idx_t *pointer, idx_t *adjacency,
                                  idx_t *elem_partition);
#endif

    void get_node_entity_proc_data(INT *entity_proc, const Ioss::MapContainer &node_map,
                                   bool do_map) const;

    void get_element_block_communication(std::vector<BlockDecompositionData> &el_blocks);
    void build_global_to_local_elem_map();
    void get_local_node_list();
    void get_shared_node_list();
    template <typename T>
    void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const;

    template <typename T>
    void communicate_set_data(T *file_data, T *ioss_data, const SetDecompositionData &set,
                              size_t comp_count) const;

    template <typename T, typename U>
    void communicate_block_data(T *file_data, U *ioss_data, const BlockDecompositionData &block,
                                size_t comp_count) const;

    template <typename T>
    void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const;

    MPI_Comm    m_comm;
    int         m_processor{};
    int         m_processorCount{};
    std::string m_method;

    // Values for the file decomposition
    int    m_spatialDimension;
    size_t m_globalElementCount;
    size_t m_elementCount;
    size_t m_elementOffset;
    size_t m_importPreLocalElemIndex;

    size_t m_globalNodeCount;
    size_t m_nodeCount;
    size_t m_nodeOffset;
    size_t m_importPreLocalNodeIndex;

    bool m_retainFreeNodes;
    bool m_showProgress;
    bool m_showHWM;

    std::vector<double> m_centroids;
    std::vector<INT>    m_pointer;   // Index into adjacency, processor list for each element...
    std::vector<INT>    m_adjacency; // Size is sum of element connectivity sizes

    std::vector<INT> m_nodeCommMap; // node/processor pair of the
    // nodes I communicate with.  Stored node#,proc,node#,proc, ...

    // The global element at index 'I' (0-based) is on block B in the
    // file decomposition.
    // if m_fileBlockIndex[B] <= I && m_fileBlockIndex[B+1] < I
    std::vector<size_t> m_fileBlockIndex;

  private:
    // This processor "manages" the elements on the exodus mesh file from
    // element_offset to element_offset+count (0-based). This is
    // 'file' data
    //
    // This processor also appears to the Ioss clients to own other
    // element and node data based on the decomposition.  This is the
    // 'ioss' data.
    //
    // The indices in 'local_element_map' are the elements that are
    // common to both the 'file' data and the 'ioss' data.
    // local_element_map[i] contains the location in 'file' data for
    // the 'ioss' data at location 'i+import_pre_local_elem_index'
    //
    // local_element_map[i]+m_elementOffset is the 0-based global index
    //
    // The indices in 'import_element_map' map the data received via
    // mpi communication from other processors into 'ioss' data.
    // if 'ind=import_element_map[i]', then ioss[ind] = comm_recv[i]
    // Note that this is the reverse direction of the
    // local_element_map mapping.
    //
    // The indices in 'export_element_map' are used to pull from
    // 'file' data into the comm_send vector.  if 'ind =
    // export_element_map[i]', then 'comm_send[i] = file[ind]' for i =
    // 0..#exported_elements
    //
    // local_element_map.size() + import_element_map.size() == #
    // ioss elements on this processor.
    //
    // local_element_map.size() + export_element_map.size() == #
    // file elements on this processor.
    //
    // export_element_map and import_element_map are sorted.
    // The primary key is processor order followed by global id.
    // The processor association is via 'export_proc_disp' and
    // 'import_proc_disp' Both are of size '#processors+1' and
    // the elements for processor p range from [X_proc_disp[p] to
    // X_proc_disp[p+1])

    std::vector<INT> localElementMap;

    std::vector<INT> importElementMap;
    std::vector<INT> importElementCount;
    std::vector<INT> importElementIndex;

    std::vector<INT> exportElementMap;
    std::vector<INT> exportElementCount;
    std::vector<INT> exportElementIndex;

    std::vector<INT> nodeIndex;

    std::vector<INT> exportNodeMap;
    std::vector<INT> exportNodeCount;
    std::vector<INT> exportNodeIndex;

    std::vector<INT>
                     importNodeMap; // Where to put each imported nodes data in the list of all data...
    std::vector<INT> importNodeCount;
    std::vector<INT> importNodeIndex;

    std::vector<INT> localNodeMap;

    std::vector<INT> m_elementDist;
    std::vector<INT> m_nodeDist;

    // Note that nodeGTL is a sorted vector.
    std::vector<INT>   nodeGTL; // Convert from global index to local index (1-based)
    std::map<INT, INT> elemGTL; // Convert from global index to local index (1-based)
  };
} // namespace Ioss
#endif
