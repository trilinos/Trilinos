// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_BlockStructuredGrid_hpp
#define percept_BlockStructuredGrid_hpp

#include <percept/structured/StructuredBlock.hpp>
#include <percept/structured/StructuredCellIndex.hpp>
#include <percept/MeshType.hpp>

#include <array>
#include <unordered_map>

namespace Ioss {
  class Region;
}

namespace percept {

  class BlockStructuredGrid {

    // create parallel VTK files
    void create_pvd(const std::string& file_prefix, bool paraviewBroken=true);
    void create_pvts(const std::string& file_prefix);
    
  public:
    
    static const int *permute_to_unstructured;
    
    stk::ParallelMachine m_comm;
    Ioss::Region *m_region;
    std::vector<std::shared_ptr<StructuredBlock> > m_sblocks;

    BlockStructuredGrid(stk::ParallelMachine comm, Ioss::Region *region=0);

    void print(std::ostream& out, int level);

    void read_cgns();

    void dump_vtk(const std::string& file_prefix);

    std::map<std::string, std::shared_ptr<MTSGridField> > m_fields;

    void register_field(const std::string& field, unsigned fourth_dim);

    // FIXME - need to find out if Kokkos requires including the last index of Array4D
    //    (i.e. the index over the xyz components)
    //    Also, need to generalize for general fields if we need to prolong fields as
    //    well as coordinates.

    // for one of the blocks, index0 =  indx[L0]-input_sizes.node_min[L0] + sizes[L0]*(indx[L1] - input_sizes.node_min[L1]) + sizes[L0]*sizes[L1]*(indx[L2] - input_sizes.node_min[L2]);
    // So, we now have iblock index, so:
    // index_in = iblock + nblocks*index0;

    inline
    void multi_dim_indices_from_index(const uint64_t& index_in, StructuredCellIndex& indx) const
    {
      std::array<unsigned,3> indx0;
      unsigned iblock=0;
      multi_dim_indices_from_index(index_in, indx0, iblock);
      indx[0] = indx0[0];
      indx[1] = indx0[1];
      indx[2] = indx0[2];
      indx[3] = iblock;
    }

    inline
    void multi_dim_indices_from_index(const uint64_t& index_in, std::array<unsigned,3>& indx, unsigned& iblock) const
    {
      size_t nblocks = m_sblocks.size();
      iblock = index_in % nblocks;
      uint64_t index =  index_in / nblocks;
      std::shared_ptr<StructuredBlock> sblock = m_sblocks[iblock];
      const int L0 = sblock->m_loop_ordering[0], L1 = sblock->m_loop_ordering[1], L2 = sblock->m_loop_ordering[2];
      SGridSizes& sb_sizes = sblock->m_sizes;
      const unsigned sizes[3] = {
        1+ sb_sizes.node_max[L0] - sb_sizes.node_min[L0],
        1+ sb_sizes.node_max[L1] - sb_sizes.node_min[L1],
        1+ sb_sizes.node_max[L2] - sb_sizes.node_min[L2]
      };
      indx[L2] = sb_sizes.node_min[L2] + (index / (sizes[L0]*sizes[L1] ));
      indx[L1] = sb_sizes.node_min[L1] + ((index / sizes[L0]) % sizes[L1] );
      indx[L0] = sb_sizes.node_min[L0] + (index % sizes[L0]);

      if (1)
        {
          uint64_t ii = indx[L0]-sb_sizes.node_min[L0] + sizes[L0]*(indx[L1] - sb_sizes.node_min[L1]) + sizes[L0]*sizes[L1]*(indx[L2] - sb_sizes.node_min[L2]);
          VERIFY_OP_ON(ii, ==, index, "bad index");
        }
    }

    inline
    void index_from_multi_dim_indices(const std::array<unsigned,3>& indx, const unsigned& iblock, uint64_t& index) const
    {
      size_t nblocks = m_sblocks.size();
      std::shared_ptr<StructuredBlock> sblock = m_sblocks[iblock];
      const int L0 = sblock->m_loop_ordering[0], L1 = sblock->m_loop_ordering[1], L2 = sblock->m_loop_ordering[2];
      SGridSizes& sb_sizes = sblock->m_sizes;
      const unsigned sizes[3] = {
        1+ sb_sizes.node_max[L0] - sb_sizes.node_min[L0],
        1+ sb_sizes.node_max[L1] - sb_sizes.node_min[L1],
        1+ sb_sizes.node_max[L2] - sb_sizes.node_min[L2]
      };

      uint64_t indx0 = indx[L0]-sb_sizes.node_min[L0] + sizes[L0]*(indx[L1] - sb_sizes.node_min[L1]) + sizes[L0]*sizes[L1]*(indx[L2] - sb_sizes.node_min[L2]);
      index = iblock + indx0*nblocks;
    }

    inline
    void index_from_multi_dim_indices(const StructuredCellIndex& cindx, uint64_t& index_in) const
    {
      std::array<unsigned,3> indx { {cindx[0], cindx[1], cindx[2]} };
      unsigned iblock = indx[3];
      index_from_multi_dim_indices(indx, iblock, index_in);
    }

    void get_nodes(std::vector<StructuredCellIndex>& nodes, unsigned offset = 0);

    void get_element_nodes(const StructuredCellIndex& element, std::vector<StructuredCellIndex>& nodes);

    void get_elements(std::vector<StructuredCellIndex>& elements);

    void copy_field(typename StructuredGrid::MTField* field_dest, typename StructuredGrid::MTField* field_src);

    void get_nodes_of_sb(std::vector<StructuredCellIndex>& nodes, unsigned iBlock, unsigned offset = 0);

    void get_elements_of_sb(std::vector<StructuredCellIndex>&elements, unsigned iBlock);

    /// axpby calculates: y = alpha*x + beta*y
    void nodal_field_axpby(double alpha, typename StructuredGrid::MTField* field_x, double beta, typename StructuredGrid::MTField* field_y);

    /// axpbypgz calculates: z = alpha*x + beta*y + gamma*z
    void nodal_field_axpbypgz(double alpha, typename StructuredGrid::MTField* field_x,
                              double beta, typename StructuredGrid::MTField* field_y,
                              double gamma, typename StructuredGrid::MTField* field_z);

    long double nodal_field_dot(typename StructuredGrid::MTField* field_x, typename StructuredGrid::MTField* field_y);

    /// set field to constant value
    void nodal_field_set_value(typename StructuredGrid::MTField* field_x, double value = 0.0);

    void comm_fields(std::vector<const typename StructuredGrid::MTField*>& fields);
    void sum_fields(std::vector<const typename StructuredGrid::MTField*>& fields);

    static std::shared_ptr<BlockStructuredGrid>
    fixture_1(stk::ParallelMachine comm, std::array<unsigned,3> sizes=std::array<unsigned,3>{{2u,2u,2u}},std::array<double,3> dim_widths=std::array<double,3>{{1.0,1.0,1.0}},std::array<double,3> dim_offsets=std::array<double,3>{{0.0,0.0,0.0}});

    size_t parallel_count_nodes();
  };

  struct DetermineIfLowestRankedOwner
  {   //for a given node on a given block, determine if the given block has the lowest blockID that contains the given node
      std::shared_ptr<StructuredBlock> m_sblock;
      SGridSizes m_sizes;

      Kokkos::View<unsigned**,DataLayout,MemSpace> localBeg; //no_zones X 3 (ijk)
      Kokkos::View<unsigned**,DataLayout,MemSpace> localEnd; //no_zones X 4 (ijkb)
      unsigned noZones;
      unsigned blockID;

      DetermineIfLowestRankedOwner(std::shared_ptr<StructuredBlock> m_sblock_in) : m_sblock(m_sblock_in)
      {
          Kokkos::resize(localBeg, m_sblock->m_zoneConnectivity.size(),3);
          Kokkos::resize(localEnd, m_sblock->m_zoneConnectivity.size(),4);

          Kokkos::View<unsigned**,DataLayout,MemSpace>::HostMirror localBeg_mir = Kokkos::create_mirror_view(localBeg);
          Kokkos::View<unsigned**,DataLayout,MemSpace>::HostMirror localEnd_mir = Kokkos::create_mirror_view(localEnd);

          noZones = m_sblock->m_zoneConnectivity.size();
          for (unsigned izone=0; izone < noZones; izone++)
          {
              Ioss::ZoneConnectivity zoneConnectivity = m_sblock->m_zoneConnectivity[izone];
              unsigned dblockid = zoneConnectivity.m_donorZone - 1;

              Ioss::IJK_t localBeg;
              localBeg[0] = std::min(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]) - 1;
              localBeg[1] = std::min(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]) - 1;
              localBeg[2] = std::min(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]) - 1;

              Ioss::IJK_t localEnd;
              localEnd[0] = std::max(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]) - 1;
              localEnd[1] = std::max(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]) - 1;
              localEnd[2] = std::max(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]) - 1;

              for(unsigned ijk=0;ijk<3;ijk++)
              {
                  localBeg_mir(izone,ijk) = (unsigned)localBeg[ijk];
                  localEnd_mir(izone,ijk) = (unsigned)localEnd[ijk];
              }
              localEnd_mir(izone,3) = dblockid;
          }
          Kokkos::deep_copy(localBeg,localBeg_mir);
          Kokkos::deep_copy(localEnd,localEnd_mir);

          m_sizes = m_sblock->m_sizes;
          blockID = m_sblock->get_block_id();
      }

      KOKKOS_INLINE_FUNCTION
      bool operator()(Kokkos::Array<unsigned, 3>& index) const
      {   //generalized boundary selector, should work with an arbitrary block structured mesh
          bool candidate = (index[0] == 0 || index[0] == m_sizes.node_size[0]-1 ||
                  index[1] == 0 || index[1] == m_sizes.node_size[1]-1 ||
                  index[2] == 0 || index[2] == m_sizes.node_size[2]-1 ); //first verify coordinate is on one of the "ends" of a Structured Block

          bool isAConnectivity = false; //assume it is a part of no connetivity
          unsigned min_block_id = blockID;
          if(candidate){ //if it's on an end
              for (unsigned izone=0; izone < noZones; izone++)
              {
                  //predicated on assumption i = 0, j = 1, k = 2 locally for all cases, for both zone connectivities and their views
                 bool isCurrentInConnectivity = ( ( localBeg(izone,0)<=index[0] && index[0]<=localEnd(izone,0) ) &&
                          ( localBeg(izone,1)<=index[1] && index[1]<=localEnd(izone,1) ) &&
                          ( localBeg(izone,2)<=index[2] && index[2]<=localEnd(izone,2) ) ); //verify whether or not it is part of an zone connectivity between blocks
                 if(isCurrentInConnectivity)
                 {
                     if(min_block_id>localEnd(izone,3))
                         min_block_id=localEnd(izone,3);
                     isAConnectivity= true; //if it is the current connectivity, then it is in a connectivity
                 }
              }//foreach zone
          }
          return isAConnectivity && min_block_id==blockID;
      }
  };


  struct DeviceSafeSGridBoundaryNotInterfaceSelector : public StructuredGrid::MTSelector {
      //NOTE: selects a node if : its is on a boundary AND NOT on an interface
      std::shared_ptr<StructuredBlock> m_sblock;
      SGridSizes m_sizes;

      Kokkos::View<unsigned**,DataLayout,MemSpace> localBeg; //no_zones X 3 (ijk)
      Kokkos::View<unsigned**,DataLayout,MemSpace> localEnd; //no_zones X 3 (ijk)
      unsigned noZones;

      DeviceSafeSGridBoundaryNotInterfaceSelector(std::shared_ptr<StructuredBlock> m_sblock_in) : m_sblock(m_sblock_in)
      {
          Kokkos::resize(localBeg, m_sblock->m_zoneConnectivity.size(),3);
          Kokkos::resize(localEnd, m_sblock->m_zoneConnectivity.size(),3);

          Kokkos::View<unsigned**,DataLayout,MemSpace>::HostMirror localBeg_mir = Kokkos::create_mirror_view(localBeg);
          Kokkos::View<unsigned**,DataLayout,MemSpace>::HostMirror localEnd_mir = Kokkos::create_mirror_view(localEnd);

          noZones = m_sblock->m_zoneConnectivity.size();
          for (unsigned izone=0; izone < noZones; izone++)
          {
              Ioss::ZoneConnectivity zoneConnectivity = m_sblock->m_zoneConnectivity[izone];

              Ioss::IJK_t localBeg;
              localBeg[0] = std::min(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]) - 1;
              localBeg[1] = std::min(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]) - 1;
              localBeg[2] = std::min(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]) - 1;

              Ioss::IJK_t localEnd;
              localEnd[0] = std::max(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]) - 1;
              localEnd[1] = std::max(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]) - 1;
              localEnd[2] = std::max(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]) - 1;

              for(unsigned ijk=0;ijk<3;ijk++)
              {
                  localBeg_mir(izone,ijk) = (unsigned)localBeg[ijk];
                  localEnd_mir(izone,ijk) = (unsigned)localEnd[ijk];
              }
          }
          Kokkos::deep_copy(localBeg,localBeg_mir);
          Kokkos::deep_copy(localEnd,localEnd_mir);

          m_sizes = m_sblock->m_sizes;
      }

      virtual ~DeviceSafeSGridBoundaryNotInterfaceSelector(){}

      KOKKOS_INLINE_FUNCTION
      bool operator()(Kokkos::Array<unsigned, 3>& index) const
      {   //generalized boundary selector, should work with an arbitrary block structured mesh
          bool candidate = (index[0] == 0 || index[0] == m_sizes.node_size[0]-1 ||
                  index[1] == 0 || index[1] == m_sizes.node_size[1]-1 ||
                  index[2] == 0 || index[2] == m_sizes.node_size[2]-1 ); //first verify coordinate is on one of the "ends" of a Structured Block

          if(candidate){ //if it's on an end
              for (unsigned izone=0; izone < noZones; izone++)
              {
                  //predicated on assumption i = 0, j = 1, k = 2 locally for all cases, for both zone connectivities and their views
                 bool isInConnectivity = ( ( localBeg(izone,0)<=index[0] && index[0]<=localEnd(izone,0) ) &&
                          ( localBeg(izone,1)<=index[1] && index[1]<=localEnd(izone,1) ) &&
                          ( localBeg(izone,2)<=index[2] && index[2]<=localEnd(izone,2) ) ); //verify whether or not it is part of an zone connectivity between blocks
                 candidate =  !isInConnectivity;
                 if(!candidate)
                     break;
              }//foreach zone
          }
          return candidate;
      }
  };


}//PERCEPT


#endif
