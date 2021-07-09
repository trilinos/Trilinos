// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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


#ifndef STKMESHINTERFACETONEMESIS_HPP_
#define STKMESHINTERFACETONEMESIS_HPP_

#include <ne_nemesisI.h>
#include <test_utils/NemesisInfo.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>

class StkMeshInterfaceToNemesis
{
public:
    StkMeshInterfaceToNemesis(stk::mesh::BulkData& bulkData, MPI_Comm comm) :
            sideRank(bulkData.mesh_meta_data().side_rank()),
            m_comm(comm),
            exoTranslator(bulkData)
    {
    }

    virtual ~StkMeshInterfaceToNemesis()
    {
    }

    NemesisInfo getNemesisInfo()
    {
        set_global_nemesis_data();
        translate_to_nemesis_set_data();
        translate_to_nemesis_block_data();
        return nemesis_info;
    }

private:
    void set_global_nemesis_data()
    {
        set_global_number_elements();
        set_global_number_nodes();
        set_global_number_element_blocks();
        set_global_number_sets();
    }

    void translate_to_nemesis_set_data()
    {
        translate_nodeset_ids();
        translate_nodeset_data();
        translate_sideset_ids();
        translate_sideset_data();
    }

    void translate_to_nemesis_block_data()
    {
        translate_block_ids();
        translate_block_data();
    }

    void set_global_number_elements()
    {
        nemesis_info.globalNumberElements = exoTranslator.get_number_global_elements();
    }

    void set_global_number_nodes()
    {
        nemesis_info.globalNumberNodes = exoTranslator.get_number_global_nodes();
    }

    void set_global_number_element_blocks()
    {
        nemesis_info.globalNumberElementBlocks = exoTranslator.get_number_element_blocks();
    }

    void set_global_number_sets()
    {
        nemesis_info.globalNumberNodeSets = exoTranslator.get_number_node_sets();
        nemesis_info.globalNumberSideSets = exoTranslator.get_number_side_sets();
    }

    void translate_nodeset_ids()
    {
        std::vector<stk::mesh::ExodusTranslator::IdType> nodesetIds;
        exoTranslator.fill_node_set_ids(nodesetIds);
        translate_to_nemesis_ids(nodesetIds, nemesis_info.nem_nodeset_ids);
    }

    NemesisId get_global_number_nodes_in_nodeset(int nodeset_id)
    {
        return exoTranslator.get_global_num_entities_for_id(nodeset_id, stk::topology::NODE_RANK);
    }

    void size_nemesis_nodeset_data()
    {
        nemesis_info.nem_nodeset_ids.resize(nemesis_info.nem_nodeset_ids.size());
        nemesis_info.num_nodes_per_nodeset.resize(nemesis_info.nem_nodeset_ids.size(),0);
    }

    void translate_nodeset_data()
    {
        size_nemesis_nodeset_data();
        for(size_t i=0;i<nemesis_info.nem_nodeset_ids.size();++i)
            nemesis_info.num_nodes_per_nodeset[i] = get_global_number_nodes_in_nodeset(nemesis_info.nem_nodeset_ids[i]);
    }

    void translate_to_nemesis_ids(const std::vector<int64_t> original_ids, std::vector<NemesisId>& nem_ids)
    {
        nem_ids.resize(original_ids.size());
        for(size_t i=0;i<original_ids.size();++i)
            nem_ids[i] = original_ids[i];
    }

    void translate_sideset_ids()
    {
        std::vector<stk::mesh::ExodusTranslator::IdType> sideset_ids;
        exoTranslator.fill_side_set_ids(sideset_ids);
        translate_to_nemesis_ids(sideset_ids, nemesis_info.nem_sideset_ids);
    }

    void get_num_df_per_sideset()
    {
        nemesis_info.num_df_per_sideset.resize(nemesis_info.nem_sideset_ids.size());
        for(size_t i=0;i<nemesis_info.num_df_per_sideset.size();++i)
            nemesis_info.num_df_per_sideset[i] = exoTranslator.get_global_num_distribution_factors_in_side_set(nemesis_info.nem_sideset_ids[i]);
    }

    void get_num_sides_per_sideset()
    {
        nemesis_info.num_sides_per_sideset.resize(nemesis_info.nem_sideset_ids.size());
        for(size_t i=0;i<nemesis_info.nem_sideset_ids.size();++i)
            nemesis_info.num_sides_per_sideset[i] = exoTranslator.get_global_num_entities_for_id(nemesis_info.nem_sideset_ids[i], sideRank);
    }

    void translate_sideset_data()
    {
        get_num_df_per_sideset();
        get_num_sides_per_sideset();
    }

    void translate_block_ids()
    {
        std::vector<stk::mesh::ExodusTranslator::IdType> block_ids;
        exoTranslator.fill_element_block_ids(block_ids);
        translate_to_nemesis_ids(block_ids, nemesis_info.nem_block_ids);
    }

    void translate_block_data()
    {
        nemesis_info.num_elements_per_block.resize(nemesis_info.nem_block_ids.size());
        std::vector<NemesisId> localElemsPerBlock(nemesis_info.nem_block_ids.size());
        for(size_t i=0;i<nemesis_info.nem_block_ids.size();++i)
        {
            localElemsPerBlock[i] = exoTranslator.get_local_num_entities_for_id(nemesis_info.nem_block_ids[i], stk::topology::ELEM_RANK);
        }
        stk::all_reduce_sum(m_comm, localElemsPerBlock.data(), nemesis_info.num_elements_per_block.data(), nemesis_info.nem_block_ids.size());
    }

private:
    stk::mesh::EntityRank sideRank;
    MPI_Comm m_comm;
    NemesisInfo nemesis_info;
    stk::mesh::ExodusTranslator exoTranslator;
};



#endif /* STKMESHINTERFACETONEMESIS_HPP_ */
