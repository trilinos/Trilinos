#include "ProcessKilledElements.hpp"
#include "ElemElemGraph.hpp"
#include "ElemElemGraphImpl.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk { namespace mesh {

void ensure_fresh_modifiable_state(stk::mesh::BulkData& bulkData)
{
  if(bulkData.in_modifiable_state()) {
    bulkData.modification_end();
  }
  bulkData.modification_begin();
}

class RemoteDeathBoundary
{
public:
    RemoteDeathBoundary(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph,
        const stk::mesh::EntityVector& killedElements, const stk::mesh::PartVector& parts_for_creating_side, stk::mesh::Part& active, const stk::mesh::PartVector* boundary_mesh_parts) :
            m_bulkData(bulkData), m_elementGraph(elementGraph), m_killedElements(killedElements), m_parts_for_creating_side(parts_for_creating_side), m_active(active),
            m_boundary_mesh_parts(boundary_mesh_parts), m_topology_modified(false)
    {}
    ~RemoteDeathBoundary(){}

    void update_death_boundary_for_remotely_killed_elements(std::vector<stk::mesh::sharing_info> &shared_modified,
                                                            stk::mesh::EntityVector& deletedEntities,
                                                            stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector)
    {
        std::vector<impl::GraphEdgeProc> remote_edges = get_remote_edges();
        process_remote_edges(shared_modified, deletedEntities, remoteActiveSelector, remote_edges);
    }

    void set_topology_is_modified()
    {
        m_topology_modified = true;
    }

    bool get_topology_modification_status() const
    {
        return m_topology_modified;
    }

private:

    void process_remote_edges(std::vector<stk::mesh::sharing_info> &shared_modified,
                              stk::mesh::EntityVector& deletedEntities,
                              stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                              std::vector<impl::GraphEdgeProc>& remote_edges)
    {
      for(impl::GraphEdgeProc& re : remote_edges)
      {
          stk::mesh::EntityId remote_id = re.get_remote_element_global_id();
          remoteActiveSelector[-remote_id] = false;
      }

      for(impl::GraphEdgeProc& re : remote_edges)
      {
          stk::mesh::EntityId local_id = re.get_local_element_global_id();
          int local_side = re.get_local_element_side_index();
          stk::mesh::EntityId remote_id = re.get_remote_element_global_id();
          int remote_side = re.get_remote_element_side_index();

          stk::mesh::Entity element = m_bulkData.get_entity(stk::topology::ELEM_RANK, local_id);

          impl::ParallelInfo &parallel_edge_info = m_elementGraph.get_parallel_edge_info(element, local_side, remote_id, remote_side);

          m_topology_modified = true;
          process_dead_remote_parallel_graph_edge(parallel_edge_info, element, local_side, remote_id,
                                                  shared_modified, deletedEntities, remoteActiveSelector);
      }

      if(m_elementGraph.has_shell_elements()) {
        process_remote_edges_for_shell_sides_with_multiple_adjacencies(shared_modified, deletedEntities, remoteActiveSelector, remote_edges);
      }
    }

    void process_remote_edges_for_shell_sides_with_multiple_adjacencies(std::vector<stk::mesh::sharing_info> &shared_modified,
                                                                        stk::mesh::EntityVector& /*deletedEntities*/,
                                                                        stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                                                                        std::vector<impl::GraphEdgeProc>& remote_edges)
    {
      for(impl::GraphEdgeProc& re : remote_edges)
      {
          stk::mesh::EntityId local_id = re.get_local_element_global_id();
          int local_side = re.get_local_element_side_index();
          stk::mesh::EntityId remote_id = re.get_remote_element_global_id();
          int remote_side = re.get_remote_element_side_index();

          stk::mesh::Entity element = m_bulkData.get_entity(stk::topology::ELEM_RANK, local_id);

          stk::mesh::impl::LocalId elemLocalId = m_elementGraph.get_local_element_id(element);
          stk::mesh::GraphEdgesForElement graphEdges = m_elementGraph.get_edges_for_element(elemLocalId);

          for(size_t i = 0; i < graphEdges.size(); ++i) {
            const GraphEdge& graphEdge = m_elementGraph.get_graph().get_edge_for_element(elemLocalId, i);
            int sideOrdinal = graphEdge.side1();
            if(local_side == sideOrdinal) {
              stk::mesh::impl::LocalId elem2LocalId = graphEdge.elem2();
              bool isParallelGraphEdge = !stk::mesh::impl::is_local_element(elem2LocalId);
              if(isParallelGraphEdge) {
                stk::mesh::EntityId elem2GlobalId = m_elementGraph.convert_negative_local_id_to_global_id(elem2LocalId);

                if((elem2GlobalId != remote_id) || (graphEdge.side2() != remote_side)) {
                  auto iter = remoteActiveSelector.find(elem2LocalId);
                  bool isActiveRemoteGraphElement = (iter != remoteActiveSelector.end() ? iter->second : true);
                  if(isActiveRemoteGraphElement) {
                    impl::ParallelInfo &other_parallel_edge_info = m_elementGraph.get_parallel_edge_info(graphEdge);

                    impl::add_side_into_exposed_boundary(m_bulkData,
                                                         other_parallel_edge_info,
                                                         element,
                                                         local_side,
                                                         m_parts_for_creating_side,
                                                         shared_modified,
                                                         remoteActiveSelector,
                                                         m_active,
                                                         m_boundary_mesh_parts);
                  }
                }
              }
            }
          }
      }
    }

    void process_dead_remote_parallel_graph_edge(const impl::ParallelInfo &parallel_edge_info,
                                                 stk::mesh::Entity element,
                                                 int local_side,
                                                 stk::mesh::EntityId /*remote_id*/,
                                                 std::vector<stk::mesh::sharing_info> &shared_modified,
                                                 stk::mesh::EntityVector& deletedEntities,
                                                 stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector)
    {
        bool create_side = m_bulkData.bucket(element).member(m_active);
        if(create_side==true)
        {
            impl::add_side_into_exposed_boundary(m_bulkData,
                                                 parallel_edge_info,
                                                 element,
                                                 local_side,
                                                 m_parts_for_creating_side,
                                                 shared_modified,
                                                 remoteActiveSelector,
                                                 m_active,
                                                 m_boundary_mesh_parts);
        }
        else
        {
            if(!impl::remove_side_from_death_boundary(m_bulkData, element, m_active, deletedEntities, local_side, remoteActiveSelector)) {
              stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(m_bulkData, element, local_side);

              if(m_bulkData.is_valid(side)) {
                // Update shared info for existing side that can't be removed
                int other_proc = parallel_edge_info.get_proc_rank_of_neighbor();
                int owning_proc = std::min(other_proc, m_bulkData.parallel_rank());
                shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, owning_proc));
              } else {
                // Create the side since it is meant to exist
                impl::add_side_into_exposed_boundary(m_bulkData,
                                                     parallel_edge_info,
                                                     element,
                                                     local_side,
                                                     m_parts_for_creating_side,
                                                     shared_modified,
                                                     remoteActiveSelector,
                                                     m_active,
                                                     m_boundary_mesh_parts);
              }
            }
        }
    }

    std::vector<impl::GraphEdgeProc> get_remote_edges() const
    {
        std::vector<impl::GraphEdgeProc> elements_to_comm = get_elements_to_communicate();
        return impl::communicate_killed_entities(m_bulkData.parallel(), elements_to_comm);
    }

    std::vector<impl::GraphEdgeProc> get_elements_to_communicate() const
    {
        std::vector<impl::GraphEdgeProc> elements_to_comm;

        for(stk::mesh::Entity this_element :m_killedElements)
        {
            for(size_t j=0;j<m_elementGraph.get_num_connected_elems(this_element);++j)
            {
                if(impl::does_element_have_side(m_bulkData, this_element) && !m_elementGraph.is_connected_elem_locally_owned(this_element, j))
                {
                    impl::IdViaSidePair idViaSidePair = m_elementGraph.get_connected_remote_id_and_via_side(this_element,j);
                    stk::mesh::EntityId other_element_id = idViaSidePair.id;
                    int side1 = idViaSidePair.side;
                    int side2 = m_elementGraph.get_connected_elements_side(this_element,j);
                    int other_proc = m_elementGraph.get_owning_proc_id_of_remote_element(this_element, j);
                    elements_to_comm.push_back(impl::GraphEdgeProc(m_bulkData.identifier(this_element), side1, other_element_id, side2, other_proc));
                }
            }
        }

        return elements_to_comm;
    }

    stk::mesh::BulkData& m_bulkData;
    ElemElemGraph& m_elementGraph;
    const stk::mesh::EntityVector& m_killedElements;
    const stk::mesh::PartVector& m_parts_for_creating_side;
    stk::mesh::Part& m_active;
    const stk::mesh::PartVector* m_boundary_mesh_parts;
    bool m_topology_modified;
};

bool process_killed_elements(stk::mesh::BulkData& bulkData,
                             const stk::mesh::EntityVector& killedElements,
                             stk::mesh::Part& active,
                             stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                             const stk::mesh::PartVector& parts_for_creating_side,
                             const stk::mesh::PartVector* boundary_mesh_parts,
                             stk::mesh::impl::MeshModification::modification_optimization modEndOpt)
{
    ensure_fresh_modifiable_state(bulkData);
    bulkData.m_bucket_repository.set_remove_mode_tracking();
    impl::create_sides_created_during_death_part(bulkData.mesh_meta_data());

    std::vector<stk::mesh::sharing_info> shared_modified;
    stk::mesh::EntityVector deletedEntities;

    bulkData.initialize_face_adjacent_element_graph();
    ElemElemGraph& elementGraph = bulkData.get_face_adjacent_element_graph();

    RemoteDeathBoundary remote_death_boundary(bulkData, elementGraph, killedElements, parts_for_creating_side, active, boundary_mesh_parts);
    remote_death_boundary.update_death_boundary_for_remotely_killed_elements(shared_modified, deletedEntities, remoteActiveSelector);

    std::vector<impl::ElementSidePair> element_side_pairs;
    element_side_pairs.reserve(impl::get_element_side_multiplier() * killedElements.size());

    for(size_t k = 0; k < killedElements.size(); ++k)
    {
        stk::mesh::Entity this_element = killedElements[k];

        for(size_t j = 0; j < elementGraph.get_num_connected_elems(this_element); ++j)
        {
            if(impl::does_element_have_side(bulkData, this_element))
            {
                remote_death_boundary.set_topology_is_modified();
                if(elementGraph.is_connected_elem_locally_owned(this_element, j))
                {
                    impl::ElementViaSidePair other_element_via_side = elementGraph.get_connected_element_and_via_side(this_element, j);
                    stk::mesh::Entity other_element = other_element_via_side.element;
                    if(impl::does_element_have_side(bulkData, other_element_via_side.element))
                    {
                        int side_id = other_element_via_side.side;
                        STK_ThrowRequireWithSierraHelpMsg(side_id != -1);

                        bool is_other_element_alive = bulkData.bucket(other_element).member(active);
                        if(is_other_element_alive)
                        {
                            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, this_element, side_id);

                            if(bulkData.is_valid(side))
                            {
                                if(bulkData.bucket(side).owned())
                                {
                                    stk::mesh::ConstPartVector parts = impl::get_stk_parts_for_moving_parts_into_death_boundary(boundary_mesh_parts);
                                    bulkData.change_entity_parts(side, parts);
                                }
                            }
                            else if(impl::can_add_side_into_exposed_boundary(bulkData, this_element, side_id, remoteActiveSelector, active))
                            {
                                stk::mesh::PartVector parts = impl::get_parts_for_creating_side(bulkData, parts_for_creating_side, other_element, side_id);

                                // switch elements
                                stk::mesh::Entity element_with_perm_0 = other_element;
                                stk::mesh::Entity element_with_perm_4 = this_element;

                                int side_id_needed = elementGraph.get_connected_elements_side(this_element, j);

                                STK_ThrowRequireMsg(side_id_needed >= 0, "ERROR: proc " << bulkData.parallel_rank() << " found side_id_needed=" << side_id_needed
                                                << " between elem " << bulkData.identifier(element_with_perm_0)<< " and " << bulkData.identifier(element_with_perm_4)
                                                << " in elem-elem-graph");

                                side = bulkData.declare_element_side(element_with_perm_0, side_id_needed, parts);
                            }
                        }
                        else
                        {
                            if(!impl::remove_side_from_death_boundary(bulkData, this_element, active, deletedEntities, side_id, remoteActiveSelector)) {
                                stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, this_element, side_id);
                                bool sideExists = bulkData.is_valid(side);

                                if(!sideExists && impl::can_add_side_into_exposed_boundary(bulkData, this_element, side_id, remoteActiveSelector, active)) {
                                  // Create the side since it is meant to exist
                                  stk::mesh::PartVector parts = impl::get_parts_for_creating_side(bulkData, parts_for_creating_side, other_element, side_id);

                                  // switch elements
                                  stk::mesh::Entity element_with_perm_0 = other_element;
                                  stk::mesh::Entity element_with_perm_4 = this_element;

                                  int side_id_needed = elementGraph.get_connected_elements_side(this_element, j);

                                  STK_ThrowRequireMsg(side_id_needed >= 0, "ERROR: proc " << bulkData.parallel_rank() << " found side_id_needed=" << side_id_needed
                                                  << " between elem " << bulkData.identifier(element_with_perm_0)<< " and " << bulkData.identifier(element_with_perm_4)
                                                  << " in elem-elem-graph");

                                  side = bulkData.declare_element_side(element_with_perm_0, side_id_needed, parts);
                                }
                            }
                        }
                    }
                }
                else
                {
                    impl::IdViaSidePair remote_id_side_pair = elementGraph.get_connected_remote_id_and_via_side(this_element, j);
                    stk::mesh::EntityId remote_id = remote_id_side_pair.id;
                    int remote_side = elementGraph.get_connected_elements_side(this_element, j);
                    impl::ParallelInfo &parallel_edge_info = elementGraph.get_parallel_edge_info(this_element, remote_id_side_pair.side, remote_id, remote_side);
                    bool other_element_active = remoteActiveSelector[-remote_id];
                    bool create_side = other_element_active;

                    if(create_side)
                    {
                        impl::add_side_into_exposed_boundary(bulkData, parallel_edge_info, this_element, remote_id_side_pair.side, parts_for_creating_side,
                                                             shared_modified, remoteActiveSelector, active, boundary_mesh_parts);
                    }
                    else
                    {
                        int side_id = remote_id_side_pair.side;
                        STK_ThrowRequireWithSierraHelpMsg(side_id != -1);
                        if(!impl::remove_side_from_death_boundary(bulkData, this_element, active, deletedEntities, side_id, remoteActiveSelector)) {
                          stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, this_element, side_id);
                          if(bulkData.is_valid(side)) {
                            // Update shared info for existing side that can't be removed
                            int other_proc = parallel_edge_info.get_proc_rank_of_neighbor();
                            int owning_proc = std::min(other_proc, bulkData.parallel_rank());
                            shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, owning_proc));
                          } else {
                            // Create the side since it is meant to exist
                            impl::add_side_into_exposed_boundary(bulkData, parallel_edge_info, this_element, remote_id_side_pair.side, parts_for_creating_side,
                                                            shared_modified, remoteActiveSelector, active, boundary_mesh_parts);
                          }
                        }
                    }
                }
            }
        }
    }

    stk::mesh::impl::delete_entities_and_upward_relations(bulkData, deletedEntities);
    stk::mesh::update_sharing_info_ownership(shared_modified);
    bulkData.make_mesh_parallel_consistent_after_element_death(shared_modified, deletedEntities, elementGraph, killedElements, active, modEndOpt);
    bulkData.m_bucket_repository.set_remove_mode_fill_and_sort();
    return remote_death_boundary.get_topology_modification_status();
}

}} // end namespaces stk mesh

