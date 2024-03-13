// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_ChildNodeCreator.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <Akri_EntityIdPool.hpp>

namespace krino {

struct ActiveChildNodeRequest
{
    std::vector<stk::mesh::EntityId> m_parents;
    std::vector<int> m_sharing_procs;
    std::vector<std::pair<int, stk::mesh::EntityId> > m_id_proc_pairs_from_all_procs;
    bool m_id_procs_pairs_have_been_sorted;
    stk::mesh::Entity *m_node_entity;

    ActiveChildNodeRequest(const std::vector<stk::mesh::EntityId> & parents, stk::mesh::Entity *entity_place_holder=NULL)
    : m_parents(parents), m_sharing_procs(), m_id_proc_pairs_from_all_procs(),
      m_id_procs_pairs_have_been_sorted(false), m_node_entity(entity_place_holder)
    {
        std::sort(m_parents.begin(), m_parents.end());
    }

    void add_proc_id_pair(int proc_id, stk::mesh::EntityId id)
    {
        m_id_proc_pairs_from_all_procs.push_back(std::make_pair(proc_id, id));
    }

    void calculate_sharing_procs(stk::mesh::BulkData& mesh)
    {
      STK_ThrowRequire(!m_parents.empty());

      stk::mesh::EntityKey key0(stk::topology::NODE_RANK, m_parents[0]);
      mesh.comm_shared_procs(key0, m_sharing_procs);
      std::sort(m_sharing_procs.begin(), m_sharing_procs.end());

      std::vector<int> sharingProcs;
      for (unsigned i=1; i<m_parents.size(); ++i)
      {
        stk::mesh::EntityKey key(stk::topology::NODE_RANK, m_parents[i]);
        mesh.comm_shared_procs(key, sharingProcs);
        std::sort(sharingProcs.begin(), sharingProcs.end());

        std::vector<int> working_set;
        m_sharing_procs.swap(working_set);
        std::set_intersection(working_set.begin(),working_set.end(),sharingProcs.begin(),sharingProcs.end(),std::back_inserter(m_sharing_procs));
      }
    }

    size_t num_sharing_procs() const
    {
        return m_sharing_procs.size();
    }

    int sharing_proc(int index) const
    {
        return m_sharing_procs[index];
    }

    stk::mesh::EntityId suggested_node_id() const
    {
        STK_ThrowRequireMsg(!m_id_procs_pairs_have_been_sorted, "Invalid use of child node calculation. Contact sierra-help");
        return m_id_proc_pairs_from_all_procs[0].second;
    }

    void sort_id_proc_pairs()
    {
        m_id_procs_pairs_have_been_sorted = true;
        std::sort(m_id_proc_pairs_from_all_procs.begin(), m_id_proc_pairs_from_all_procs.end());
    }

    stk::mesh::EntityId get_id_for_child() const
    {
        STK_ThrowRequireMsg(m_id_procs_pairs_have_been_sorted, "Invalid use of child node calculation. Contact sierra-help");
        return m_id_proc_pairs_from_all_procs[0].second;
    }

    void set_node_entity_for_request(stk::mesh::BulkData& mesh, const stk::mesh::PartVector & node_parts, std::vector<stk::mesh::EntityProc> & nodeSharing)
    {
        this->sort_id_proc_pairs();
        stk::mesh::EntityId id_for_child = get_id_for_child();
        *m_node_entity = mesh.declare_node(id_for_child, node_parts);
        for (size_t i=0;i<m_id_proc_pairs_from_all_procs.size();++i)
        {
            if ( m_id_proc_pairs_from_all_procs[i].first != mesh.parallel_rank() )
            {
                nodeSharing.emplace_back(*m_node_entity, m_id_proc_pairs_from_all_procs[i].first);
            }
        }
    }

    bool operator==(const ActiveChildNodeRequest& other) const { return m_parents == other.m_parents; }
    bool operator<(const ActiveChildNodeRequest & other) const { return m_parents < other.m_parents; }

    std::string dump(stk::mesh::BulkData& mesh) const
    {
      std::ostringstream out;
      out << "Child node " << ((m_node_entity && mesh.is_valid(*m_node_entity)) ? mesh.identifier(*m_node_entity) : 0) << " with parents ";
      for (unsigned i=0; i<m_parents.size(); ++i)
      {
        out << m_parents[i] << " ";
      }
      out << " with sharing ";
      for (size_t i=0;i<m_id_proc_pairs_from_all_procs.size();++i)
      {
        out << m_id_proc_pairs_from_all_procs[i].first << " ";
      }
      out << std::endl;

      return out.str();
    }
};

void batch_add_node_sharing(stk::mesh::BulkData& mesh, const std::vector<stk::mesh::EntityProc> & nodeSharing)
{
  for (auto && nodeAndProc : nodeSharing)
    mesh.add_node_sharing(nodeAndProc.first, nodeAndProc.second);
}

void batch_create_child_nodes(stk::mesh::BulkData & mesh, const std::vector< ChildNodeRequest > & child_node_requests, const stk::mesh::PartVector & node_parts, const GenerateNewEntityIdsFunction & generate_new_ids)
{
    std::vector<bool> communicate_request(child_node_requests.size(), false);

    unsigned num_nodes_requested = child_node_requests.size();
    std::vector<stk::mesh::EntityId> available_node_ids;

    generate_new_ids(stk::topology::NODE_RANK, num_nodes_requested, available_node_ids);

    while ( true )
    {
        int more_work_to_be_done = false;

        std::vector<ActiveChildNodeRequest> active_child_node_requests;

        for (unsigned it_req=0; it_req<child_node_requests.size(); ++it_req)
        {
            if (communicate_request[it_req])
            {
              continue;
            }

            const ChildNodeRequest & request = child_node_requests[it_req];
            const std::vector<stk::mesh::Entity*> & request_parents = request.parents;

            bool request_is_ready = true;
            for (auto && request_parent : request_parents)
            {
              if (!mesh.is_valid(*request_parent))
              {
                request_is_ready = false;
              }
            }

            if (request_is_ready)
            {
                stk::mesh::Entity *request_child = request.child;

                communicate_request[it_req] = true;
                more_work_to_be_done = true;

                std::vector<stk::mesh::EntityId> parent_ids(request_parents.size());
                for (size_t parent_index=0; parent_index<request_parents.size(); ++parent_index )
                {
                  parent_ids[parent_index] = mesh.identifier(*request_parents[parent_index]);
                }

                active_child_node_requests.push_back(ActiveChildNodeRequest(parent_ids, request_child));
                active_child_node_requests.back().add_proc_id_pair(mesh.parallel_rank(), available_node_ids[it_req]);
                active_child_node_requests.back().calculate_sharing_procs(mesh);
            }
        }

        int global_more_work_to_be_done = false;
        stk::all_reduce_max(mesh.parallel(), &more_work_to_be_done, &global_more_work_to_be_done, 1);
        if ( global_more_work_to_be_done == 0 ) break;

        std::sort(active_child_node_requests.begin(), active_child_node_requests.end());

        // By this point, I have list of requests that I need to communicate

        stk::CommSparse comm_spec(mesh.parallel());

        for (int phase=0;phase<2;++phase)
        {
            for (size_t request_index=0;request_index<active_child_node_requests.size();++request_index)
            {
                for (size_t proc_index=0;proc_index<active_child_node_requests[request_index].num_sharing_procs();++proc_index)
                {
                    const int other_proc = active_child_node_requests[request_index].sharing_proc(proc_index);
                    if (other_proc != mesh.parallel_rank())
                    {
                      const std::vector<stk::mesh::EntityId> & request_parents = active_child_node_requests[request_index].m_parents;
                      const stk::mesh::EntityId this_procs_suggested_id = active_child_node_requests[request_index].suggested_node_id();
                      const size_t num_parents = request_parents.size();
                      comm_spec.send_buffer(other_proc).pack(num_parents);
                      for (size_t parent_index=0; parent_index<request_parents.size(); ++parent_index )
                      {
                        comm_spec.send_buffer(other_proc).pack(request_parents[parent_index]);
                      }
                      comm_spec.send_buffer(other_proc).pack(this_procs_suggested_id);
                    }
                }
            }

            if ( phase == 0 )
            {
              comm_spec.allocate_buffers();
            }
            else
            {
                comm_spec.communicate();
            }
        }

        for(int i = 0; i < mesh.parallel_size(); ++i)
        {
            if(i != mesh.parallel_rank())
            {
                while(comm_spec.recv_buffer(i).remaining())
                {
                    size_t num_parents;
                    stk::mesh::EntityId suggested_node_id;
                    comm_spec.recv_buffer(i).unpack(num_parents);
                    std::vector<stk::mesh::EntityId> request_parents(num_parents);
                    for (size_t parent_index=0; parent_index<num_parents; ++parent_index )
                    {
                      comm_spec.recv_buffer(i).unpack(request_parents[parent_index]);
                    }
                    comm_spec.recv_buffer(i).unpack(suggested_node_id);

                    ActiveChildNodeRequest from_other_proc(request_parents);
                    std::vector<ActiveChildNodeRequest>::iterator iter = std::lower_bound(active_child_node_requests.begin(), active_child_node_requests.end(), from_other_proc);

                    if ( iter != active_child_node_requests.end() && *iter == from_other_proc)
                    {
                        iter->add_proc_id_pair(i, suggested_node_id);
                    }
                }
            }
        }

        std::vector<stk::mesh::EntityProc> nodeSharing;
        for (size_t request_index=0;request_index<active_child_node_requests.size();++request_index)
        {
            active_child_node_requests[request_index].set_node_entity_for_request(mesh, node_parts, nodeSharing);
        }
        batch_add_node_sharing(mesh, nodeSharing);
    }

    std::vector<bool>::iterator iter = std::find(communicate_request.begin(), communicate_request.end(), false);
    STK_ThrowRequireMsg(iter == communicate_request.end(), "Invalid child node request. Contact sierra-help.");
}

}



