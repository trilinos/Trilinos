// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef gradient_funtors_hpp
#define gradient_funtors_hpp

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/MeshType.hpp>
#include <percept/PerceptUtils.hpp>
#include <fstream>

namespace percept
{

struct sGrid_comm_shared_field_data_post_grad
{
    std::shared_ptr<BlockStructuredGrid> m_bsg;

    StructuredGrid::MTField::Array4D don_cg_g;
    StructuredGrid::MTField::Array4D loc_cg_g;

    StructuredGrid::MTField::Array4D don_checker;
    StructuredGrid::MTField::Array4D loc_checker;

    mutable Kokkos::Array<int, 3> transform_arr;
    mutable Kokkos::Array<int, 3> local_range_beg;
    mutable Kokkos::Array<int, 3> donor_range_beg;
    Kokkos::Array<unsigned, 3> local_sizes;

    bool doSum;


    sGrid_comm_shared_field_data_post_grad(std::shared_ptr<BlockStructuredGrid> bsg) : m_bsg(bsg)
    {
        doSum=true;
    }

    sGrid_comm_shared_field_data_post_grad() : m_bsg(0)
    {
        doSum=true;
    }

    void print_field_datas()
    {   //NOT DEVICE SAFE
        for (unsigned iblock=0; iblock < m_bsg->m_sblocks.size(); ++iblock)
        {
            std::cout << "iblock " << iblock << std::endl;
            std::shared_ptr<StructuredBlock> sblock = m_bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = m_bsg->m_sblocks[dblockid];
                std::cout << "      izone " << izone << " dblockid " << dblockid<< std::endl;

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                don_cg_g = *m_bsg->m_fields["cg_g"].get()->m_block_fields[dblockid]; //*cg_g_field->m_block_fields[iBlock];
                loc_cg_g = *m_bsg->m_fields["cg_g"].get()->m_block_fields[iblock];
                don_checker = dblock->m_sgrid_coords;
                loc_checker = sblock->m_sgrid_coords;

                for(unsigned ijk=0;ijk<3;ijk++)
                    local_sizes[ijk] = 1 + localEnd[ijk]-localBeg[ijk];

                for(unsigned iii=0;iii<3;iii++)
                {
                    transform_arr[iii] = zoneConnectivity.m_transform[iii];
                    local_range_beg[iii] =  zoneConnectivity.m_ownerRangeBeg[iii];
                    donor_range_beg[iii] = zoneConnectivity.m_donorRangeBeg[iii];
                }
                for(unsigned iNode=0;iNode<local_sizes[0]*local_sizes[1]*local_sizes[2];iNode++)
                {
                    Kokkos::Array<unsigned, 3> indx;

                    indx[0] = (local_range_beg[0] -1)
                                    + ( iNode % local_sizes[0] );
                    indx[1] = (local_range_beg[1] -1)
                                    + (iNode / local_sizes[0]) % local_sizes[1];
                    indx[2] = (local_range_beg[2] -1)
                                    +((iNode / local_sizes[0])/local_sizes[1]) % local_sizes[2];

                    Ioss::IJK_t indexLocal = {{(int)(indx[0] +1),
                            (int)(indx[1] +1),
                            (int)(indx[2] +1)}};

                    Kokkos::Array<int, 3> local_indices;
                    for(unsigned iii=0;iii<3;iii++)
                        local_indices[iii] = indexLocal[iii];

                    Kokkos::Array<int,3 > index_donor = device_safe_transform_block_indices(local_indices,
                            transform_arr,
                            local_range_beg,
                            donor_range_beg);

                    for(unsigned d=0;d<3;d++){
                        if( device_safe_abs_flt((Double)loc_checker(indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d) - (Double)don_checker(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d) ) > 1e-14 )
                            printf("dissonant coord in coord field at %d %d %d %d\n",indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d); //this never prints so I think the coordinates are parallel
                        break;
                    }

                    for(unsigned d=0;d<3;d++){
                        if( device_safe_abs_flt((Double)loc_cg_g(indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d) - (Double)don_cg_g(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d) ) > 1e-14 )
                            printf("dissonant coord in cg_g field at local ijk %d %d %d and donor ijk %d %d %d\n",indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,
                                    index_donor[0]-1,index_donor[1]-1,index_donor[2]-1);
                        break;
                    }

                    double loc_xyz[3];
                    double don_xyz[3];
                    for(unsigned d=0;d<3;d++)
                    {
                        loc_xyz[d] = loc_cg_g(indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d);
                        don_xyz[d] = don_cg_g(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d);
                    }


                    printf("        iNode %d has LOCAL ijk %d %d %d and DONOR ijk %d %d %d with LOCAL xyz %g %g %g and DONOR %g %g %g\n",iNode,indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,
                            index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,
                            loc_xyz[0],loc_xyz[1],loc_xyz[2],
                            don_xyz[0],don_xyz[1],don_xyz[2]);
                }
            }//izone
        }//iblock
    }


    void comm_field_data()
    {
        for (unsigned iblock=0; iblock < m_bsg->m_sblocks.size(); ++iblock)
        {
            std::shared_ptr<StructuredBlock> sblock = m_bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = m_bsg->m_sblocks[dblockid];


                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_ownerRangeBeg[0],zoneConnectivity.m_ownerRangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_ownerRangeBeg[1],zoneConnectivity.m_ownerRangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_ownerRangeBeg[2],zoneConnectivity.m_ownerRangeEnd[2]);

                don_cg_g = *m_bsg->m_fields["cg_g"].get()->m_block_fields[dblockid]; //*cg_g_field->m_block_fields[iBlock];
                loc_cg_g = *m_bsg->m_fields["cg_g"].get()->m_block_fields[iblock];
                don_checker = dblock->m_sgrid_coords;
                loc_checker = sblock->m_sgrid_coords;

                for(unsigned ijk=0;ijk<3;ijk++)
                    local_sizes[ijk] = 1 + localEnd[ijk]-localBeg[ijk];

                for(unsigned iii=0;iii<3;iii++)
                {
                    transform_arr[iii] = zoneConnectivity.m_transform[iii];
                    local_range_beg[iii] =  zoneConnectivity.m_ownerRangeBeg[iii];
                    donor_range_beg[iii] = zoneConnectivity.m_donorRangeBeg[iii];
                }
                if(iblock<dblockid && doSum)
                {
                    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,local_sizes[0]*local_sizes[1]*local_sizes[2]),*this);
                }
                else if(iblock>dblockid && !doSum)
                {
                    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,local_sizes[0]*local_sizes[1]*local_sizes[2]),*this);
                }
            }//izone
        }//iblock
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& iNode) const
    {
        Kokkos::Array<unsigned, 3> indx;

        indx[0] = (local_range_beg[0] -1)
                + ( iNode % local_sizes[0] );
        indx[1] = (local_range_beg[1] -1)
                + (iNode / local_sizes[0]) % local_sizes[1];
        indx[2] = (local_range_beg[2] -1)
                +((iNode / local_sizes[0])/local_sizes[1]) % local_sizes[2];

        Kokkos::Array<int, 3> local_indices;
        for(unsigned iii=0;iii<3;iii++)
            local_indices[iii] = indx[iii] + 1;

        Kokkos::Array<int,3 > index_donor = device_safe_transform_block_indices(local_indices,
                transform_arr,
                local_range_beg,
                donor_range_beg);


        if(doSum)
            for (int d=0; d<3; d++){
                //from based 1 to based 0 indexing
                loc_cg_g(indx[0],indx[1],indx[2],d) += don_cg_g(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d);
            }
        else //doProp
            for (int d=0; d<3; d++)
            //from based 1 to based 0 indexing
            loc_cg_g(indx[0],indx[1],indx[2],d) = don_cg_g(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d);
    }
};

struct sGrid_GenericAlgorithm_get_gradient_1_fixup{
    //zeros out fixed nodes on a gradient field of a particular block
    StructuredGrid::MTField::Array4D cg_g_field_iter;
    Kokkos::View< unsigned**, DataLayout, MemSpace> fixed_nodes_iterator;

    void zero_out_fixed_nodes(unsigned numberFixed)
    {
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,numberFixed),*this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
        unsigned i = fixed_nodes_iterator(index,0);
        unsigned j = fixed_nodes_iterator(index,1);
        unsigned k = fixed_nodes_iterator(index,2);
        for(unsigned iDim=0;iDim<3;iDim++)
            cg_g_field_iter(i,j,k,iDim) = 0.0;
    }

};

struct sGrid_GenericAlgorithm_get_gradient_1 {

    mutable HexMeshSmootherMetric m_metric;
    PerceptMesh * m_eMesh;

    StructuredGrid::MTSelector *m_boundarySelector;

    StructuredGrid::MTField * m_coord_field_original;
    StructuredGrid::MTField * m_coord_field_current;
    StructuredGrid::MTField * cg_g_field;
//    StructuredGrid::MTField * cg_r_field;
//    StructuredGrid::MTField * cg_edge_length_field;

    StructuredGrid::MTField::Array4D m_coords_current_iter;
    StructuredGrid::MTField::Array4D m_coord_field_orig_iter;
    StructuredGrid::MTField::Array4D cg_g_field_iter;
//    StructuredGrid::MTField::Array4D cg_r_field_iter;
//    StructuredGrid::MTField::Array4D cg_edge_length_field_iter;
    std::shared_ptr<BlockStructuredGrid> m_bsg;

    std::vector<Kokkos::View<unsigned**, DataLayout, MemSpace> > fixed_nodes_per_block;
    Kokkos::View<unsigned**, DataLayout, MemSpace> fixed_nodes_iterator;
    std::vector<unsigned> number_fixed_per_block;

    std::vector<SGridSizes> block_sizes;
    std::vector<std::array<unsigned int, 3>> loop_orderings;
    SGridSizes block_sizes_iterator;
    Kokkos::Array<unsigned int, 3> loop_orderings_iterator;

    unsigned m_noBlocks;

    sGrid_GenericAlgorithm_get_gradient_1_fixup m_fixup;
    sGrid_comm_shared_field_data_post_grad m_shared_data_fixup;

    sGrid_GenericAlgorithm_get_gradient_1(PerceptMesh *eMesh,
            StructuredGrid::MTSelector *boundarySelector=0) :
            m_metric(eMesh), m_eMesh(eMesh), m_boundarySelector(
                    boundarySelector) {
        if (eMesh->get_block_structured_grid()) {

            m_bsg = m_eMesh->get_block_structured_grid();
            m_shared_data_fixup.m_bsg = m_bsg;

            m_noBlocks = m_bsg->m_sblocks.size();

            m_coord_field_original = m_bsg->m_fields["coordinates_NM1"].get();
            m_coord_field_current = m_bsg->m_fields["coordinates"].get();
            cg_g_field = m_bsg->m_fields["cg_g"].get();

            block_sizes.resize(m_noBlocks);
            loop_orderings.resize(m_noBlocks);

            for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) { //populate loop orderings and block sizes for all blocks in mesh
                block_sizes[iBlock] = m_bsg->m_sblocks[iBlock]->m_sizes;
                loop_orderings[iBlock] =
                        m_bsg->m_sblocks[iBlock]->m_loop_ordering;

            }
            //filter out unfixed nodes for fixup
            std::vector<StructuredGrid::MTNode> FixedNodes_from_block;
            std::vector<StructuredGrid::MTNode> nodes_from_block;
            std::pair<bool, int> fixed;

            std::vector<
                    Kokkos::View<unsigned**, DataLayout, MemSpace>::HostMirror> nodes_mirror;
            nodes_mirror.resize(m_noBlocks);
            fixed_nodes_per_block.resize(m_noBlocks);
            number_fixed_per_block.resize(m_noBlocks);

            for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) { //for each block, filter out any unfixed nodes
                unsigned numFixed = 0;
                m_bsg->get_nodes_of_sb(nodes_from_block, iBlock);
                FixedNodes_from_block.resize(nodes_from_block.size());

                for (int64_t iNode = 0;
                        iNode < (int64_t) nodes_from_block.size(); ++iNode) {
                    fixed = get_fixed_flag_sgrid(nodes_from_block[iNode],
                            m_boundarySelector);

                    if (fixed.first) {
                        FixedNodes_from_block[numFixed] =
                                nodes_from_block[iNode];
                        numFixed++;
                    } //cluster all fixed nodes into lower indices of FixedNodes
                }
                Kokkos::resize(fixed_nodes_per_block[iBlock],
                        numFixed, 3);
                Kokkos::resize(nodes_mirror[iBlock], numFixed, 3);
                for (unsigned iElem = 0; iElem < numFixed; iElem++) {
                    nodes_mirror[iBlock](iElem, 0) =
                            FixedNodes_from_block[iElem][0];
                    nodes_mirror[iBlock](iElem, 1) =
                            FixedNodes_from_block[iElem][1];
                    nodes_mirror[iBlock](iElem, 2) =
                            FixedNodes_from_block[iElem][2];
                }
                Kokkos::deep_copy(fixed_nodes_per_block[iBlock],
                        nodes_mirror[iBlock]); //sync fixed nodes to device
                number_fixed_per_block[iBlock] = numFixed;

            } //foreach block
        }
    } //sGrid_GenericAlgorithm_get_gradient cnstr

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const {
        double v_i_current[8][3];
        double v_i_org[8][3];
        unsigned indx[3] = { 0, 0, 0 };
        unsigned II[3] = { 0, 0, 0 };
        Kokkos::Array<unsigned, 3> cell_ijk;

        sgrid_multi_dim_indices_from_index_cell(block_sizes_iterator,loop_orderings_iterator,index, cell_ijk);

        const int A0 = 0, A1 = 1, A2 = 2;

        unsigned nodes_of_elem[8][3]; //node,ijk

        unsigned cnt = 0;//gather nodal coordinates
        for (indx[2] = 0; indx[2] < 2; ++indx[2]) {
            II[2] = indx[2] + cell_ijk[2];
            for (indx[1] = 0; indx[1] < 2; ++indx[1]) {
                II[1] = indx[1] + cell_ijk[1];
                for (indx[0] = 0; indx[0] < 2; ++indx[0]) {
                    II[0] = indx[0] + cell_ijk[0];
                    for (unsigned ic = 0; ic < 3; ++ic) {
                        v_i_current[cnt][ic] = m_coords_current_iter(II[A0],
                                II[A1], II[A2], ic);
                        v_i_org[cnt][ic] = m_coord_field_orig_iter(II[A0],
                                II[A1], II[A2], ic);
                    }
                    nodes_of_elem[cnt][0] = II[A0];
                    nodes_of_elem[cnt][1] = II[A1];
                    nodes_of_elem[cnt][2] = II[A2];
                    ++cnt;
                }
            }
        }

        double analytic_grad[8][4];

        {
            bool gmvalid = true;

            m_metric.grad_metric(v_i_current, v_i_org, gmvalid, analytic_grad);

            if (gmvalid || m_metric.m_untangling) {
                for (unsigned inode = 0; inode < 8; inode++) {

                    unsigned node[3];
                    node[0] = nodes_of_elem[inode][0];
                    node[1] = nodes_of_elem[inode][1];
                    node[2] = nodes_of_elem[inode][2];

                    for (int jdim = 0; jdim < 3; jdim++) {
                        Kokkos::atomic_add(
                                &cg_g_field_iter(node[0], node[1], node[2],
                                        jdim), analytic_grad[inode][jdim]);
                    }
                }
            }
        }
    }//operator


    void run() {
        for (unsigned iBlock = 0; iBlock < m_noBlocks; iBlock++) {

            std::string timer_name;
            if(m_metric.m_untangling)
                timer_name="untangle ";
            else
                timer_name="smooth ";

            unsigned total_elems_this_block = (1
                    + block_sizes[iBlock].cell_max[0]
                    - block_sizes[iBlock].cell_min[0])
                    * (1 + block_sizes[iBlock].cell_max[1]
                            - block_sizes[iBlock].cell_min[1])
                    * (1 + block_sizes[iBlock].cell_max[2]
                            - block_sizes[iBlock].cell_min[2]);

            stk::diag::Timer root_timer("GA_gg_1 "+timer_name+std::to_string(total_elems_this_block), rootTimerStructured());
            stk::diag::TimeBlock root_block(root_timer);

            //orient iterators
            m_coords_current_iter =
                    *m_coord_field_current->m_block_fields[iBlock];
            m_coord_field_orig_iter =
                    *m_coord_field_original->m_block_fields[iBlock];
            cg_g_field_iter = *cg_g_field->m_block_fields[iBlock];

            block_sizes_iterator = block_sizes[iBlock];
            loop_orderings_iterator[0] = loop_orderings[iBlock][0];
            loop_orderings_iterator[1] = loop_orderings[iBlock][1];
            loop_orderings_iterator[2] = loop_orderings[iBlock][2];

            {//calc gradient fields
                stk::diag::Timer pf1(
                        std::string("GA_gg_1 calc grad pf 1 ")+timer_name
                        + std::to_string(
                                (int) std::cbrt((double) total_elems_this_block)),
                                root_timer);
                stk::diag::TimeBlock pf1_block(pf1);

                Kokkos::parallel_for(
                        Kokkos::RangePolicy<ExecSpace>(0, total_elems_this_block),
                        *this);
            }//timer scope

            {//zero fixed fields
                stk::diag::Timer pf2(
                        std::string("GA_gg_1 fixup pf 2 ")
                + std::to_string(
                        (int) std::cbrt((double) total_elems_this_block)),
                        root_timer);
                stk::diag::TimeBlock pf2_block(pf2);

                m_fixup.fixed_nodes_iterator = fixed_nodes_per_block[iBlock];
                m_fixup.cg_g_field_iter = *cg_g_field->m_block_fields[iBlock];
                m_fixup.zero_out_fixed_nodes(number_fixed_per_block[iBlock]);
            }//timer scope
        }//foreach block

        {//communicate fields between blocks
            std::string timer_name;
            if(m_metric.m_untangling)
                timer_name="untangle ";
            else
                timer_name="smooth ";

            stk::diag::Timer pf3(
                    std::string("GA_gg_1 interblock comm 1 ")+timer_name, rootTimerStructured());
            stk::diag::TimeBlock pf3_block(pf3);
            m_shared_data_fixup.doSum = true;
            m_shared_data_fixup.comm_field_data();
            m_shared_data_fixup.doSum = false;
            m_shared_data_fixup.comm_field_data();

        }//timer scope

    }//run()
};
}//percept
#endif
