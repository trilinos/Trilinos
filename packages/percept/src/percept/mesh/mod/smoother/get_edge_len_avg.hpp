// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef get_edge_len_avg_hpp
#define get_edge_len_avg_hpp

#include <percept/MeshType.hpp>
#include <percept/PerceptUtils.hpp>

namespace percept {

KOKKOS_INLINE_FUNCTION
void sgrid_find_connected_cells(const SGridSizes& sizes, Kokkos::Array<unsigned, 3>& node_ijk, unsigned node_elems[8][3], unsigned& num_adj_elems)
{ //doesn't take into account block interfaces
  const int node_sizes[3] = {static_cast<int>(sizes.node_size[0]),
                        static_cast<int>(sizes.node_size[1]),
                        static_cast<int>(sizes.node_size[2])};
  num_adj_elems = 0;
  unsigned elem_cnt = 0;
  for (int i0 = node_ijk[0]-1; i0 <= static_cast<int>(node_ijk[0]); ++i0)
    {
      if (i0 < 0 || i0 > node_sizes[0]-2) continue;
      for (int i1 = node_ijk[1]-1; i1 <= static_cast<int>(node_ijk[1]); ++i1)
        {
          if (i1 < 0 || i1 > node_sizes[1]-2) continue;
          for (int i2 = node_ijk[2]-1; i2 <= static_cast<int>(node_ijk[2]); ++i2)
            {
              if (i2 < 0 || i2 > node_sizes[2]-2) continue;
              node_elems[elem_cnt][0] =static_cast<unsigned>(i0);
              node_elems[elem_cnt][1] =static_cast<unsigned>(i1);
              node_elems[elem_cnt][2] =static_cast<unsigned>(i2);
              num_adj_elems = elem_cnt +1;
              elem_cnt++;
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION
void get_nodes_sgrid(unsigned element[3],unsigned node_elems[8][3])
{
  unsigned node_cnt=0;
  for (unsigned i2 = element[2]; i2 <= element[2]+1; ++i2)
    {
      for (unsigned i1 = element[1]; i1 <= element[1]+1; ++i1)
        {
          for (unsigned i0 = element[0]; i0 <= element[0]+1; ++i0)
            {
              node_elems[node_cnt][0]=i0;
              node_elems[node_cnt][1]=i1;
              node_elems[node_cnt][2]=i2;
              node_cnt++;
            }
        }
    }

}

KOKKOS_INLINE_FUNCTION
double edge_length_avg(unsigned element[3], StructuredGrid::MTField::Array4D coord_field_block, double* min_edge_length_in, double* max_edge_length_in)
{
  unsigned elem_nodes[8][3];
  get_nodes_sgrid(element, elem_nodes);

  double edge_length_ave=0.0;
  double min_edge_length = -1.0;
  double max_edge_length = -1.0;
  unsigned edge_count = 12;
  const int edges[12][2] = {{0,1},{1,3},{3,2},{2,0}, {4,5},{5,7},{7,6},{6,4}, {0,4},{1,5},{2,6},{3,7} };
  for (unsigned iedgeOrd = 0; iedgeOrd < edge_count; ++iedgeOrd)
    {
      unsigned in0 = edges[iedgeOrd][0];
      unsigned in1 = edges[iedgeOrd][1];

      Double node_coord_data_0[3];
      node_coord_data_0[0] = coord_field_block(elem_nodes[in0][0],elem_nodes[in0][1],elem_nodes[in0][2],0);
      node_coord_data_0[1] = coord_field_block(elem_nodes[in0][0],elem_nodes[in0][1],elem_nodes[in0][2],1);
      node_coord_data_0[2] = coord_field_block(elem_nodes[in0][0],elem_nodes[in0][1],elem_nodes[in0][2],2);

      Double node_coord_data_1[3];
      node_coord_data_1[0] = coord_field_block(elem_nodes[in1][0],elem_nodes[in1][1],elem_nodes[in1][2],0);
      node_coord_data_1[1] = coord_field_block(elem_nodes[in1][0],elem_nodes[in1][1],elem_nodes[in1][2],1);
      node_coord_data_1[2] = coord_field_block(elem_nodes[in1][0],elem_nodes[in1][1],elem_nodes[in1][2],2);

      double edge_length = 0.0;
      for (int iSpaceDimOrd = 0; iSpaceDimOrd < 3; iSpaceDimOrd++)
        {
          edge_length +=
            (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd])*
            (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd]);
        }
      edge_length = std::sqrt(edge_length);
      edge_length_ave += edge_length / ((double)edge_count);
      if(iedgeOrd == 0)
        {
          min_edge_length = edge_length;
          max_edge_length = edge_length;
        }
      else
        {
          if(min_edge_length < edge_length)
              min_edge_length=edge_length;
          if(max_edge_length > edge_length)
              max_edge_length=edge_length ;
        }
    }
  if (min_edge_length_in) *min_edge_length_in = min_edge_length;
  if (max_edge_length_in) *max_edge_length_in = max_edge_length;
  return edge_length_ave;
}



KOKKOS_INLINE_FUNCTION
Double
sgrid_nodal_edge_length_ave(Kokkos::Array<unsigned, 3>& ijk, StructuredGrid::MTField::Array4D block, const SGridSizes& sizes,
        StructuredGrid::MTField::Array4D adj_elem_field, bool adj_elem_field_filled)
{
  Double nm=0.0;

  double min=std::numeric_limits<double>::max();
  unsigned node_elems[8][3];//there are up to 8 adjacent elements to each node
  unsigned num_adj_elems=0; //actual number of elems adjacent to node
  sgrid_find_connected_cells(sizes, ijk, node_elems,num_adj_elems);
  if(!adj_elem_field_filled)
      adj_elem_field(ijk[0],ijk[1],ijk[2],0) = (double)num_adj_elems;

  for (unsigned ii=0; ii < num_adj_elems; ++ii)
    {
      unsigned element[3] = {node_elems[ii][0],  node_elems[ii][1],  node_elems[ii][2]};
      double lmin=0,lmax=0;
      double elem_edge_len = edge_length_avg(element, block, &lmin, &lmax);

      if (lmin < min)
          min=lmin;
      nm += elem_edge_len;
    }
  return nm; //will compute the total and average at this node later on
}

struct sGrid_sum_and_propagate_interfaces
{
    std::shared_ptr<BlockStructuredGrid> m_bsg;

    StructuredGrid::MTField::Array4D don_ela;
    StructuredGrid::MTField::Array4D loc_ela;

    StructuredGrid::MTField::Array4D don_nae;
    StructuredGrid::MTField::Array4D loc_nae;

    mutable Kokkos::Array<int, 3> transform_arr;
    mutable Kokkos::Array<int, 3> local_range_beg;
    mutable Kokkos::Array<int, 3> donor_range_beg;
    Kokkos::Array<unsigned, 3> local_sizes;

    bool doSum;
    bool adj_elem_field_already_fixed;

public:

    sGrid_sum_and_propagate_interfaces(std::shared_ptr<BlockStructuredGrid> m_bsg_in) : m_bsg(m_bsg_in){
        doSum = true;
        adj_elem_field_already_fixed = false;
    }
    void run()
    {
        for (unsigned iblock=0; iblock < m_bsg->m_sblocks.size(); ++iblock)
        {
            //            std::cout << "iblock " << iblock << std::endl;
            std::shared_ptr<StructuredBlock> sblock = m_bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = m_bsg->m_sblocks[dblockid];

                //                std::cout << "      izone " << izone << "dblockid " << dblockid<< std::endl;

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                don_ela = *m_bsg->m_fields["cg_edge_length"].get()->m_block_fields[dblockid];
                loc_ela = *m_bsg->m_fields["cg_edge_length"].get()->m_block_fields[iblock];

                don_nae = *m_bsg->m_fields["num_adj_elems"].get()->m_block_fields[dblockid];
                loc_nae = *m_bsg->m_fields["num_adj_elems"].get()->m_block_fields[iblock];

                for(unsigned ijk=0;ijk<3;ijk++)
                    local_sizes[ijk] = 1 + localEnd[ijk]-localBeg[ijk];

                for(unsigned iii=0;iii<3;iii++)
                {
                    transform_arr[iii] = zoneConnectivity.m_transform[iii];
                    local_range_beg[iii] =  zoneConnectivity.m_ownerRangeBeg[iii];
                    donor_range_beg[iii] = zoneConnectivity.m_donorRangeBeg[iii];
                }
                if(iblock<dblockid && doSum)
                    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,local_sizes[0]*local_sizes[1]*local_sizes[2]),*this);
                else if(iblock>dblockid && !doSum) //prop
                    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,local_sizes[0]*local_sizes[1]*local_sizes[2]),*this);
            }
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& iField) const
    {
        Kokkos::Array<unsigned, 3> indx;

        indx[0] = (local_range_beg[0] -1)
                            + ( iField % local_sizes[0] );
        indx[1] = (local_range_beg[1] -1)
                            + (iField / local_sizes[0]) % local_sizes[1];
        indx[2] = (local_range_beg[2] -1)
                            +((iField / local_sizes[0])/local_sizes[1]) % local_sizes[2];

        Kokkos::Array<int, 3> local_indices;
        for(unsigned iii=0;iii<3;iii++)
            local_indices[iii] = indx[iii] + 1;

        Kokkos::Array<int,3 > index_donor = device_safe_transform_block_indices(local_indices,
                transform_arr,
                local_range_beg,
                donor_range_beg);

        if(doSum){ //summing
            loc_ela(indx[0],indx[1],indx[2],0) += don_ela(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,0);
            if(!adj_elem_field_already_fixed)
                loc_nae(indx[0],indx[1],indx[2],0) += don_nae(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,0);
        }
        else { //propagating
            loc_ela(indx[0],indx[1],indx[2],0) = don_ela(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,0);
            if(!adj_elem_field_already_fixed)
                loc_nae(indx[0],indx[1],indx[2],0) = don_nae(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,0);
        }
    }
};

  struct sGrid_GenericAlgorithm_get_edge_lengths
  {
    PerceptMesh *m_eMesh;
    StructuredGrid::MTField *cg_edge_length_field;
    StructuredGrid::MTField *m_coord_field_original;
    StructuredGrid::MTField * num_adj_elems_nodal_field;


    StructuredGrid::MTField::Array4D cg_edge_length_field_iter;
    StructuredGrid::MTField::Array4D m_coord_field_original_iter;
    StructuredGrid::MTField::Array4D num_adj_elems_nodal_field_iter;

    std::vector<SGridSizes> block_sizes;
    std::vector<std::array<unsigned int, 3>> loop_orderings;
    SGridSizes block_sizes_iterator;
    Kokkos::Array<unsigned int, 3> loop_orderings_iterator;

    unsigned m_noBlocks;
    bool adj_elem_field_filled;
    bool adj_elem_field_fixed;
    bool compute_final_average;

    sGrid_GenericAlgorithm_get_edge_lengths(PerceptMesh *eMesh) : m_eMesh(eMesh)
    {
      if (eMesh->get_block_structured_grid()) {
                  std::shared_ptr<BlockStructuredGrid> bsg =
                          m_eMesh->get_block_structured_grid();

                  m_noBlocks = bsg->m_sblocks.size();

                  cg_edge_length_field = bsg->m_fields["cg_edge_length"].get();
                  m_coord_field_original  = bsg->m_fields["coordinates_NM1"].get();
                  num_adj_elems_nodal_field = bsg->m_fields["num_adj_elems"].get();

                  block_sizes.resize(m_noBlocks);
                  loop_orderings.resize(m_noBlocks);

                  for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) { //populate loop orderings and block sizes for all blocks in mesh
                      block_sizes[iBlock] = bsg->m_sblocks[iBlock]->m_sizes;
                      loop_orderings[iBlock] =
                              bsg->m_sblocks[iBlock]->m_loop_ordering;
                  }
              }//if bsg
      adj_elem_field_filled = false;
      adj_elem_field_fixed = false;
      compute_final_average = false;


    }


    void calc_edge_lengths(bool debug_print = false)
    {
        compute_final_average = false;
        for(unsigned iBlock=0;iBlock<m_noBlocks;iBlock++)
        {
            //orient iterators
            cg_edge_length_field_iter =  *cg_edge_length_field->m_block_fields[iBlock];
            m_coord_field_original_iter = *m_coord_field_original->m_block_fields[iBlock];
            num_adj_elems_nodal_field_iter = *num_adj_elems_nodal_field->m_block_fields[iBlock];

            block_sizes_iterator = block_sizes[iBlock];
            loop_orderings_iterator[0] = loop_orderings[iBlock][0];
            loop_orderings_iterator[1] = loop_orderings[iBlock][1];
            loop_orderings_iterator[2] = loop_orderings[iBlock][2];

            unsigned total_nodes_this_block = (1
                    + block_sizes[iBlock].node_max[0]
                    - block_sizes[iBlock].node_min[0])
                    * (1 + block_sizes[iBlock].node_max[1]
                            - block_sizes[iBlock].node_min[1])
                    * (1 + block_sizes[iBlock].node_max[2]
                            - block_sizes[iBlock].node_min[2]);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,total_nodes_this_block),*this);
        }

        adj_elem_field_filled = true;

        if(debug_print)
            std::cout << "      dot product on ela_field mesh before fixup " <<(m_eMesh->nodal_field_dot("cg_edge_length", "cg_edge_length")) << std::endl;
        //run fixup
        sGrid_sum_and_propagate_interfaces fixup(m_eMesh->get_block_structured_grid());
        fixup.adj_elem_field_already_fixed = adj_elem_field_fixed;//only calculate this field once
        fixup.doSum = true;
        fixup.run(); //sum

        fixup.doSum = false;
        fixup.run(); //prop

        adj_elem_field_fixed  = true;

        compute_final_average = true;
        for(unsigned iBlock=0;iBlock<m_noBlocks;iBlock++)
        {
            //orient iterators
            cg_edge_length_field_iter =  *cg_edge_length_field->m_block_fields[iBlock];
            m_coord_field_original_iter = *m_coord_field_original->m_block_fields[iBlock];
            num_adj_elems_nodal_field_iter = *num_adj_elems_nodal_field->m_block_fields[iBlock];

            block_sizes_iterator = block_sizes[iBlock];
            loop_orderings_iterator[0] = loop_orderings[iBlock][0];
            loop_orderings_iterator[1] = loop_orderings[iBlock][1];
            loop_orderings_iterator[2] = loop_orderings[iBlock][2];

            unsigned total_nodes_this_block = (1
                    + block_sizes[iBlock].node_max[0]
                    - block_sizes[iBlock].node_min[0])
                    * (1 + block_sizes[iBlock].node_max[1]
                            - block_sizes[iBlock].node_min[1])
                    * (1 + block_sizes[iBlock].node_max[2]
                            - block_sizes[iBlock].node_min[2]);

            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,total_nodes_this_block),*this);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
        Kokkos::Array<unsigned,3> ijk;
        sgrid_multi_dim_indices_from_index_node(block_sizes_iterator,loop_orderings_iterator,index,ijk);
        if(!compute_final_average)
        {
            Double edge_length_ave = sgrid_nodal_edge_length_ave(ijk,m_coord_field_original_iter,block_sizes_iterator,num_adj_elems_nodal_field_iter,adj_elem_field_filled);
            cg_edge_length_field_iter(ijk[0],ijk[1],ijk[2],0) = edge_length_ave;
        }
        else //if compute_final_average
        {
            cg_edge_length_field_iter(ijk[0],ijk[1],ijk[2],0) /= num_adj_elems_nodal_field_iter(ijk[0],ijk[1],ijk[2],0);
        }
    }
  };
}//percept
#endif
