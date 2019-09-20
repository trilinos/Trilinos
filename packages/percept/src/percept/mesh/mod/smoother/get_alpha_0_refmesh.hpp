// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/MeshType.hpp>
#include <percept/PerceptUtils.hpp>

#ifndef get_alpha_0_refmesh_hpp
#define get_alpha_0_refmesh_hpp

namespace percept
{



template<typename T>
struct min_scanner {   //finds the min value of a one dimensional view
    Kokkos::View<T*, DataLayout, MemSpace> candidates;

    KOKKOS_INLINE_FUNCTION void init(T&interimmin) const {
        interimmin = candidates(0);
    }

    KOKKOS_INLINE_FUNCTION void join(volatile T& dst,
            const volatile T& src) const {
        if (dst > src) {
            dst = src;
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int64_t index, T&interimmin) const {
        if (interimmin > candidates(index))
            interimmin = candidates(index);
    }

    min_scanner(Kokkos::View<T*, DataLayout, MemSpace> candidates_in) {
        candidates = candidates_in;
    }

    T find_min() {
        T min;
        Kokkos::parallel_reduce(
                Kokkos::RangePolicy<ExecSpace>(0, candidates.size()), *this,
                min);

        return min;
    }
};

  struct SGrid_GenericAlgorithm_get_alpha_0
  {
    PerceptMesh *m_eMesh;
    StructuredGrid::MTSelector *m_boundarySelector;

    StructuredGrid::MTField *cg_s_field;
    StructuredGrid::MTField *cg_edge_length_field;
    StructuredGrid::MTField::Array4D cg_s_field_iter;
    StructuredGrid::MTField::Array4D cg_edge_length_field_iter;
    std::vector<Kokkos::View<Double*, DataLayout,MemSpace>> alpha_candidates;
    Kokkos::View<Double*, DataLayout,MemSpace> alpha_candidates_iter;

    Double alpha;

    std::vector<SGridSizes> block_sizes;
    std::vector<std::array<unsigned int, 3>> loop_orderings;
    SGridSizes block_sizes_iterator;
    Kokkos::Array<unsigned int, 3> loop_orderings_iterator;

    unsigned m_noBlocks;

    SGrid_GenericAlgorithm_get_alpha_0(PerceptMesh *eMesh,StructuredGrid::MTSelector *boundarySelector=0): m_eMesh(eMesh), m_boundarySelector(boundarySelector)
    {
        if (eMesh->get_block_structured_grid()) {

                    std::shared_ptr<BlockStructuredGrid> bsg =
                            m_eMesh->get_block_structured_grid();

                    m_noBlocks = bsg->m_sblocks.size();

                    cg_s_field = bsg->m_fields["cg_s"].get();;
                    cg_edge_length_field = bsg->m_fields["cg_edge_length"].get();

                    block_sizes.resize(m_noBlocks);
                    loop_orderings.resize(m_noBlocks);
                    alpha_candidates.resize(m_noBlocks);

                    for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) { //populate loop orderings and block sizes for all blocks in mesh
                        block_sizes[iBlock] = bsg->m_sblocks[iBlock]->m_sizes;
                        loop_orderings[iBlock] =
                                bsg->m_sblocks[iBlock]->m_loop_ordering;
                        unsigned total_nodes_this_block = (1
                                + block_sizes[iBlock].node_max[0]
                                - block_sizes[iBlock].node_min[0])
                                * (1 + block_sizes[iBlock].node_max[1]
                                        - block_sizes[iBlock].node_min[1])
                                * (1 + block_sizes[iBlock].node_max[2]
                                        - block_sizes[iBlock].node_min[2]);
                        Kokkos::resize(alpha_candidates[iBlock],total_nodes_this_block);
                    }
                    alpha=std::numeric_limits<Double>::max();
                }//if bsg
    }

    void reset_alpha(){alpha=std::numeric_limits<Double>::max();}

    void calc_alpha()
    {   //finds minimum(alpha>0)
        for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) {//compute all candidate minimums
            unsigned total_nodes_this_block = (1
                               + block_sizes[iBlock].node_max[0]
                               - block_sizes[iBlock].node_min[0])
                               * (1 + block_sizes[iBlock].node_max[1]
                                       - block_sizes[iBlock].node_min[1])
                               * (1 + block_sizes[iBlock].node_max[2]
                                       - block_sizes[iBlock].node_min[2]);


            cg_s_field_iter = *cg_s_field->m_block_fields[iBlock];
            cg_edge_length_field_iter =  *cg_edge_length_field->m_block_fields[iBlock];

            alpha_candidates_iter=alpha_candidates[iBlock];

            block_sizes_iterator = block_sizes[iBlock];
            loop_orderings_iterator[0] = loop_orderings[iBlock][0];
            loop_orderings_iterator[1] = loop_orderings[iBlock][1];
            loop_orderings_iterator[2] = loop_orderings[iBlock][2];

            stk::diag::Timer root_timer("Get_alpha "+std::to_string(total_nodes_this_block), rootTimerStructured());
            stk::diag::TimeBlock root_block(root_timer);

            {
                stk::diag::Timer pf(
                        std::string("get_alpha_candidates"),
                        root_timer);
                stk::diag::TimeBlock pf_block(pf);
                Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,total_nodes_this_block),*this);
            }//timer scope

            {
                stk::diag::Timer pr(
                        std::string("get_smallest_alpha"),
                        root_timer);
                stk::diag::TimeBlock pr_block(pr);
                min_scanner<Double> minscan(alpha_candidates_iter);
                Double cand = minscan.find_min();
                if( alpha > cand)
                    alpha=cand;
            }//timer scope

        }

//      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & alpha_set ) );
//      if (!alpha_set)
//        alpha = 1.0;

      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMin<1>( & alpha ) );
      VERIFY_OP_ON(alpha, > , 0.0, "bad alpha");

      if( alpha==std::numeric_limits<Double>::max() )//if all the fields were fixed
          alpha=1.0;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
      Kokkos::Array<unsigned,3> ijk;
      sgrid_multi_dim_indices_from_index_node(block_sizes_iterator,loop_orderings_iterator,index,ijk);

      unsigned i = ijk[0];
      unsigned j = ijk[1];
      unsigned k = ijk[2];

      Double edge_length_ave = cg_edge_length_field_iter(i,j,k,0);

      Double sn = 0.0;
      for (int idim=0; idim < 3; idim++)
        {
          sn += cg_s_field_iter(i,j,k,idim)*cg_s_field_iter(i,j,k,idim);
        }
      sn = std::sqrt(sn);
      if (sn > 0.0)
        {
          Double alpha_candidate = edge_length_ave / sn;
          alpha_candidates_iter(index)=alpha_candidate;
        }
      else
          alpha_candidates_iter(index)=std::numeric_limits<Double>::max(); //mark fields fixed along boundary as max Double. Only fields in cg_s that have magnitude of zero are fixed fields
    }//operator
  };
}//percept

#endif
