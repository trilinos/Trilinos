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
#ifndef IOCGNS_DECOMPOSITONDATA_H
#define IOCGNS_DECOMPOSITONDATA_H

#include <string>
#include <unordered_map>
#include <vector>

#include <cstddef>
#include <cstdint>

#include <Ioss_CodeTypes.h>
#include <Ioss_Decomposition.h>
#include <Ioss_Field.h>
#include <Ioss_MeshType.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_StructuredBlock.h>
#include <cgns/Iocgns_StructuredZoneData.h>

#include <cgnslib.h>

#if 0
#if !defined(NO_PARMETIS_SUPPORT)
#include <parmetis.h>
#endif
#endif

#undef MPICPP
#if !defined(NO_ZOLTAN_SUPPORT)
#include <zoltan_cpp.h>
#endif
namespace Ioss {
  class Field;
  template <typename INT> class Decomposition;
} // namespace Ioss

namespace Iocgns {

  class ZoneData
  {
  public:
    std::string m_name;
    size_t      m_nodeOffset;
    size_t      m_nodeCount;
    size_t      m_elementOffset;
  };

  class DecompositionDataBase
  {
  public:
    DecompositionDataBase(MPI_Comm comm) {}

    virtual ~DecompositionDataBase();
    virtual void   decompose_model(int serFilePtr, int filePtr, Ioss::MeshType mesh_type) = 0;
    virtual size_t ioss_node_count() const                                                = 0;
    virtual size_t ioss_elem_count() const                                                = 0;
    virtual int    int_size() const                                                       = 0;

    virtual int    spatial_dimension() const = 0;
    virtual size_t global_node_count() const = 0;
    virtual size_t global_elem_count() const = 0;

    virtual size_t decomp_node_offset() const = 0;
    virtual size_t decomp_node_count() const  = 0;
    virtual size_t decomp_elem_offset() const = 0;
    virtual size_t decomp_elem_count() const  = 0;

    virtual std::vector<double> &centroids() = 0;

    virtual size_t get_commset_node_size() const = 0;

    virtual void get_node_coordinates(int filePtr, double *ioss_data,
                                      const Ioss::Field &field) const = 0;

    void get_block_connectivity(int filePtr, void *data, int blk_seq) const;

    void get_element_field(int filePtr, int solution_index, int blk_seq, int field_index,
                           double *data) const;

    void get_node_field(int filePtr, int solution_index, int field_index, double *data) const;

    void get_node_entity_proc_data(void *entity_proc, const Ioss::MapContainer &node_map,
                                   bool do_map) const;

    template <typename T>
    void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const;

    template <typename T>
    void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const;

    void get_sideset_element_side(int filePtr, const Ioss::SetDecompositionData &sset,
                                  void *data) const;

    std::vector<ZoneData>                     m_zones;
    std::vector<Ioss::BlockDecompositionData> m_elementBlocks;
    std::vector<Ioss::SetDecompositionData>   m_sideSets;
    std::vector<Iocgns::StructuredZoneData *> m_structuredZones;

    // Maps nodes shared between zones.
    // TODO: Currently each processor has same map; need to figure out how to reduce size
    std::unordered_map<cgsize_t, cgsize_t> m_zoneSharedMap;
  };

  template <typename INT> class DecompositionData : public DecompositionDataBase
  {
  public:
    DecompositionData(const Ioss::PropertyManager &props, MPI_Comm communicator);
    ~DecompositionData() {}

    int int_size() const { return sizeof(INT); }

    void decompose_model(int serFilePtr, int filePtr, Ioss::MeshType mesh_type);

    int spatial_dimension() const { return m_decomposition.m_spatialDimension; }

    size_t global_node_count() const { return m_decomposition.global_node_count(); }
    size_t global_elem_count() const { return m_decomposition.global_elem_count(); }

    size_t ioss_node_count() const { return m_decomposition.ioss_node_count(); }
    size_t ioss_elem_count() const { return m_decomposition.ioss_elem_count(); }

    size_t decomp_node_offset() const { return m_decomposition.file_node_offset(); }
    size_t decomp_node_count() const { return m_decomposition.file_node_count(); }
    size_t decomp_elem_offset() const { return m_decomposition.file_elem_offset(); }
    size_t decomp_elem_count() const { return m_decomposition.file_elem_count(); }

    std::vector<double> &centroids() { return m_decomposition.m_centroids; }

    template <typename T>
    void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const
    {
      m_decomposition.communicate_element_data(file_data, ioss_data, comp_count);
    }

    void communicate_set_data(INT *file_data, INT *ioss_data, const Ioss::SetDecompositionData &set,
                              size_t comp_count) const
    {
      m_decomposition.communicate_set_data(file_data, ioss_data, set, comp_count);
    }

    template <typename T>
    void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const
    {
      m_decomposition.communicate_node_data(file_data, ioss_data, comp_count);
    }

    template <typename U, typename T>
    void communicate_block_data(U *file_data, T *ioss_data,
                                const Ioss::BlockDecompositionData &block, size_t comp_count) const
    {
      m_decomposition.communicate_block_data(file_data, ioss_data, block, comp_count);
    }

    void get_block_connectivity(int filePtr, INT *data, int blk_seq) const;

    void get_element_field(int filePtr, int solution_index, int blk_seq, int field_index,
                           double *data) const;

    void get_node_field(int filePtr, int solution_index, int field_index, double *data) const;

    size_t get_commset_node_size() const { return m_decomposition.m_nodeCommMap.size() / 2; }

    void get_sideset_element_side(int filePtr, const Ioss::SetDecompositionData &sset,
                                  INT *data) const;

  private:
    void decompose_structured(int serFilePtr, int filePtr);
    void decompose_unstructured(int filePtr);

    void get_sideset_data(int filePtr);
    void generate_zone_shared_nodes(int filePtr, INT min_node, INT max_node);

    bool i_own_node(size_t node)
        const // T/F if node with global index node owned by this processors ioss-decomp.
    {
      return m_decomposition.i_own_node(node);
    }

    bool i_own_elem(size_t elem)
        const // T/F if node with global index elem owned by this processors ioss-decomp.
    {
      return m_decomposition.i_own_elem(elem);
    }

    // global_index is 1-based index into global list of nodes [1..global_node_count]
    // return value is 1-based index into local list of nodes on this
    // processor (ioss-decomposition)
    size_t node_global_to_local(size_t global_index) const
    {
      return m_decomposition.node_global_to_local(global_index);
    }

    size_t elem_global_to_local(size_t global_index) const
    {
      return m_decomposition.elem_global_to_local(global_index);
    }

    void build_global_to_local_elem_map()
    {
      return m_decomposition.build_global_to_local_elem_map();
    }

    void get_element_block_communication()
    {
      m_decomposition.get_element_block_communication(m_elementBlocks);
    }

    void generate_adjacency_list(int fileId, Ioss::Decomposition<INT> &decomposition);

    void calculate_element_centroids(int filePtr, std::vector<double> &centroids);

    void get_shared_node_list() { m_decomposition.get_shared_node_list(); }

    void get_local_node_list() { m_decomposition.get_local_node_list(); }

    void get_file_node_coordinates(int filePtr, int direction, double *ioss_data) const;
    void get_node_coordinates(int filePtr, double *ioss_data, const Ioss::Field &field) const;

    double      m_loadBalanceThreshold{1.4};
    std::string m_lineDecomposition{};

  public:
    Ioss::Decomposition<INT> m_decomposition;
  };

} // namespace Iocgns
#endif
