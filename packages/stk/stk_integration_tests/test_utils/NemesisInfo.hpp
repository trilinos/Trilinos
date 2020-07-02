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

#ifndef NEMESIS_INFO_HPP
#define NEMESIS_INFO_HPP

#include <exodusII.h>
#include <ne_nemesisI.h>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

typedef int NemesisId;

struct NemesisInfo
{
    size_t globalNumberElements = 0;
    size_t globalNumberNodes = 0;
    size_t globalNumberElementBlocks = 0;
    size_t globalNumberNodeSets = 0;
    size_t globalNumberSideSets = 0;

    std::vector<NemesisId> nem_nodeset_ids;
    std::vector<NemesisId> num_nodes_per_nodeset;
    std::vector<NemesisId> nem_sideset_ids;
    std::vector<NemesisId> num_sides_per_sideset;
    std::vector<NemesisId> num_df_per_sideset;
    std::vector<NemesisId> nem_block_ids;
    std::vector<NemesisId> num_elements_per_block;
    std::vector<int> procs_for_shared_nodes;

    NemesisId num_node_cmaps = 1;


    NemesisId num_border_nodes_per_subdomain = 0;
    NemesisId num_internal_nodes_per_subdomain = 0;
    std::vector<NemesisId> node_cmap_counts;
    std::vector<NemesisId> nem_internal_node_ids;
    std::vector<NemesisId> nem_border_node_ids;
    std::vector<NemesisId> shared_node_ids;
    std::vector<NemesisId> nem_elem_ids;
    std::vector<NemesisId> node_num_map;
    std::vector<NemesisId> elem_num_map;

    void clear_subdomain_data()
    {
        num_border_nodes_per_subdomain = 0;
        num_internal_nodes_per_subdomain = 0;
        node_cmap_counts.clear();
        nem_internal_node_ids.clear();
        shared_node_ids.clear();
        nem_elem_ids.clear();
        node_num_map.clear();
        elem_num_map.clear();
    }

    size_t get_number_subdomain_nodes() const
    {
        ThrowRequireWithSierraHelpMsg(num_border_nodes_per_subdomain+num_internal_nodes_per_subdomain>0);
        return num_border_nodes_per_subdomain + num_internal_nodes_per_subdomain;
    }

    char* get_nemesis_file_type()
    {
        // "p" is for parallel nemesis file, "s" is for serial nemesis file
        static char ftype[10];
        std::string filetype("p");
        strcpy(ftype, filetype.c_str());
        return ftype;
    }

    void fill_node_num_map(int exoid, size_t num_subdomain_nodes)
    {
        node_num_map.resize(num_subdomain_nodes);
        int map_count = ex_inquire_int(exoid, EX_INQ_NODE_MAP);
        if(map_count == 1)
        {
            ex_get_num_map(exoid, EX_NODE_MAP, 1, node_num_map.data());
        }
        else
        {
            ex_get_id_map(exoid, EX_NODE_MAP, node_num_map.data());
        }
    }

    void fill_elem_num_map(int exoid, size_t num_elements_this_subdomain)
    {
        elem_num_map.resize(num_elements_this_subdomain,0);
        ex_get_id_map(exoid, EX_ELEM_MAP, elem_num_map.data());
    }

    void add_nemesis_info_to_exodus_file(int exoid, stk::mesh::Part* this_subdomain, int which_proc, int num_total_procs)
    {
        ne_put_init_global(exoid, globalNumberNodes, globalNumberElements, globalNumberElementBlocks, globalNumberNodeSets, globalNumberSideSets);
        ne_put_init_info(exoid, num_total_procs, 1, get_nemesis_file_type());
        ne_put_ns_param_global(exoid, nem_nodeset_ids.data(), num_nodes_per_nodeset.data(), num_nodes_per_nodeset.data());
        ne_put_ss_param_global(exoid, nem_sideset_ids.data(), num_sides_per_sideset.data(), num_df_per_sideset.data());
        ne_put_eb_info_global(exoid, nem_block_ids.data(), num_elements_per_block.data());

        int num_external_nodes = 0;
        int num_border_elements = 0;

        ne_put_loadbal_param(exoid, num_internal_nodes_per_subdomain, num_border_nodes_per_subdomain,
                num_external_nodes, nem_elem_ids.size(), num_border_elements, num_node_cmaps, 0, which_proc);

        ne_put_elem_map(exoid, nem_elem_ids.data(), 0, which_proc);

        ne_put_node_map(exoid, nem_internal_node_ids.data(), nem_border_node_ids.data(), 0, which_proc);

        std::vector<int> node_cmap_ids(num_node_cmaps,0);
        node_cmap_ids[0] = 1;

        ne_put_cmap_params(exoid, node_cmap_ids.data(), node_cmap_counts.data(), 0, 0, which_proc);

        ne_put_node_cmap(exoid, node_cmap_ids[0], shared_node_ids.data(), procs_for_shared_nodes.data(), which_proc);
    }

    size_t get_index_of_global_id(const NemesisId global_id, const std::vector<NemesisId>& ids)
    {
        ThrowRequireWithSierraHelpMsg(ids.size()>0);
        std::vector<NemesisId>::const_iterator iter = std::find(ids.begin(), ids.end(), global_id);
        ThrowRequireMsg(iter != ids.end(), "Program error. Could not find " << global_id);
        return std::distance(ids.begin(), iter);
    }

    std::vector<NemesisId> get_mapped_analyst_ids_to_one_base_indices(const stk::mesh::EntityVector& nodes, const std::vector<NemesisId>& id_map,
            const stk::mesh::BulkData& bulkData)
    {
        std::vector<NemesisId> mapped_ids(nodes.size(),0);
        for(size_t j=0;j<nodes.size();++j)
        {
            mapped_ids[j] = bulkData.identifier(nodes[j]);
            size_t index = get_index_of_global_id(mapped_ids[j], id_map);
            mapped_ids[j] = index+1;
        }
        return mapped_ids;
    }

    void set_nemesis_node_ids(const stk::mesh::EntityVector& shared_nodes,
            const stk::mesh::EntityVector& internal_nodes, const stk::mesh::EntityVector& unique_shared_nodes, stk::mesh::BulkData& bulkData)
    {
        nem_internal_node_ids = get_mapped_analyst_ids_to_one_base_indices(internal_nodes, node_num_map, bulkData);
        nem_border_node_ids = get_mapped_analyst_ids_to_one_base_indices(unique_shared_nodes, node_num_map, bulkData);
        shared_node_ids = get_mapped_analyst_ids_to_one_base_indices(shared_nodes, node_num_map, bulkData);
    }

    void set_nemesis_elem_ids(stk::mesh::Part* this_subdomain, stk::mesh::BulkData& bulkData)
    {
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(*this_subdomain, bulkData.buckets(stk::topology::ELEMENT_RANK), elements);
        nem_elem_ids = get_mapped_analyst_ids_to_one_base_indices(elements, elem_num_map, bulkData);
    }

    void set_node_communication_map_counts(size_t num_shared_nodes)
    {
        node_cmap_counts.resize(num_node_cmaps,0);
        node_cmap_counts[0]=num_shared_nodes;
    }

    stk::mesh::EntityVector get_subdomain_nodes(stk::mesh::Part* this_subdomain, stk::mesh::BulkData& bulkData)
    {
        stk::mesh::EntityVector subdomain_nodes;
        stk::mesh::get_selected_entities(*this_subdomain, bulkData.buckets(stk::topology::NODE_RANK), subdomain_nodes);
        return subdomain_nodes;
    }

    stk::mesh::EntityVector get_difference_between_subdomain_and_shared_nodes(stk::mesh::EntityVector &subdomain_nodes, stk::mesh::EntityVector &unique_shared_nodes)
    {
        stk::mesh::EntityVector internal_nodes;
        internal_nodes.reserve(subdomain_nodes.size());
        std::set_difference(subdomain_nodes.begin(), subdomain_nodes.end(),
                unique_shared_nodes.begin(), unique_shared_nodes.end(),
                std::inserter(internal_nodes, internal_nodes.begin()));
        return internal_nodes;
    }

    void fill_internal_nodes(stk::mesh::Part* this_subdomain, stk::mesh::EntityVector& unique_shared_nodes, stk::mesh::EntityVector& internal_nodes,
            stk::mesh::BulkData& bulkData)
    {
        stk::mesh::EntityVector subdomain_nodes = get_subdomain_nodes(this_subdomain, bulkData);
        std::sort(subdomain_nodes.begin(), subdomain_nodes.end());
        internal_nodes = get_difference_between_subdomain_and_shared_nodes(subdomain_nodes, unique_shared_nodes);
    }
};

#endif
