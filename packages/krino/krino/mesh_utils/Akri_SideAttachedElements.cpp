#include <Akri_SideAttachedElements.hpp>

#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <vector>

namespace krino {

static bool is_element_attached(const std::vector<bool> & isElemAttached, const stk::mesh::Entity elem)
{
  return isElemAttached[elem.local_offset()];
}

static void mark_element_as_attached_and_if_new_add_neighbors_on_stack(const stk::mesh::BulkData &mesh,
  const stk::mesh::Selector & elementSelector,
  const stk::mesh::Entity elem,
  std::vector<bool> & isElemAttached,
  std::vector<stk::mesh::Entity> & attachedElemStack)
{
  const auto elemOffset = elem.local_offset();
  if (!isElemAttached[elemOffset])
  {
    isElemAttached[elemOffset] = true;
    const std::vector<stk::mesh::Entity> nbrs = get_selected_side_attached_elements(mesh, elementSelector, elem);
    for (auto & nbr : nbrs)
      attachedElemStack.push_back(nbr);
  }
}

static stk::mesh::Entity stack_pop(std::vector<stk::mesh::Entity> &elemStack)
{
  const stk::mesh::Entity elem = elemStack.back();
  elemStack.pop_back();
  return elem;
}

static void mark_elements_as_attached_and_recursively_check_neighbors(const stk::mesh::BulkData &mesh,
    const stk::mesh::Selector & elementSelector,
    std::vector<bool> & isElemAttached,
    std::vector<stk::mesh::Entity> & attachedElemStack)
{
  while (!attachedElemStack.empty())
  {
    const stk::mesh::Entity elem = stack_pop(attachedElemStack);
    mark_element_as_attached_and_if_new_add_neighbors_on_stack(mesh, elementSelector, elem, isElemAttached, attachedElemStack);
  }
}

static void parallel_communicate_and_add_any_new_attached_elements_on_stack(const stk::mesh::BulkData &mesh,
    std::vector<bool> & isElemAttached,
    std::vector<stk::mesh::Entity> & attachedElemStack)
{
  std::vector<stk::mesh::Entity> unownedElemsAttachedOnThisProc;
  for ( auto & bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, !mesh.mesh_meta_data().locally_owned_part()))
    for ( auto elem : *bucket )
      if (is_element_attached(isElemAttached, elem))
        unownedElemsAttachedOnThisProc.push_back(elem);
  std::vector<stk::mesh::Entity> ownedElemsAttachedOnOtherProcs;
  communicate_entities_to_owning_proc(mesh, unownedElemsAttachedOnThisProc, ownedElemsAttachedOnOtherProcs);

  for (auto & elem : ownedElemsAttachedOnOtherProcs)
    if (!is_element_attached(isElemAttached, elem))
      attachedElemStack.push_back(elem);
}

static void parallel_sync_is_element_attached(const stk::mesh::BulkData &mesh,
    const stk::mesh::Selector & elementSelector,
    std::vector<bool> & isElemAttached,
    std::vector<stk::mesh::Entity> & attachedElemStack)
{
  parallel_communicate_and_add_any_new_attached_elements_on_stack(mesh, isElemAttached, attachedElemStack);
  while (stk::is_true_on_any_proc(mesh.parallel(), !attachedElemStack.empty()))
  {
    mark_elements_as_attached_and_recursively_check_neighbors(mesh, elementSelector, isElemAttached, attachedElemStack);
    parallel_communicate_and_add_any_new_attached_elements_on_stack(mesh, isElemAttached, attachedElemStack);
  }
}

std::vector<bool> are_elements_side_attached_to_selected_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector)
{
  std::vector<bool> isElemAttached = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, false);

  std::vector<stk::mesh::Entity> attachedElemStack;
  for ( auto & bucket : mesh.get_buckets(mesh.mesh_meta_data().side_rank(), sideSelector & mesh.mesh_meta_data().locally_owned_part()) )
    for ( auto side : *bucket )
      for (auto && elem : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
        if (elementSelector(mesh.bucket(elem)))
          attachedElemStack.push_back(elem);

  mark_elements_as_attached_and_recursively_check_neighbors(mesh, elementSelector, isElemAttached, attachedElemStack);

  parallel_sync_is_element_attached(mesh, elementSelector, isElemAttached, attachedElemStack);

  return isElemAttached;
}

std::vector<stk::mesh::Entity> get_selected_owned_side_unattached_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector)
{
  std::vector<bool> isElemAttached = are_elements_side_attached_to_selected_sides(mesh, elementSelector, sideSelector);

  std::vector<stk::mesh::Entity> unattachedElems;
  for ( auto & bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector & mesh.mesh_meta_data().locally_owned_part()) )
    for ( auto elem : *bucket )
      if (!is_element_attached(isElemAttached, elem))
        unattachedElems.push_back(elem);
  return unattachedElems;
}

static void assign_group_id_and_if_new_add_neighbors_on_stack(const stk::mesh::BulkData &mesh,
  const stk::mesh::Selector & elementSelector,
  const stk::mesh::Entity elem,
  const size_t elementGroupId,
  std::vector<size_t> & elementGroupIds,
  std::vector<stk::mesh::Entity> & attachedElemStack)
{
  const auto elemOffset = elem.local_offset();
  if (elementGroupIds[elemOffset] == 0)
  {
    elementGroupIds[elemOffset] = elementGroupId;
    const std::vector<stk::mesh::Entity> nbrs = get_selected_side_attached_elements(mesh, elementSelector, elem);
    for (auto & nbr : nbrs)
      attachedElemStack.push_back(nbr);
  }
}

static void assign_element_group_id_and_recursively_check_neighbors(const stk::mesh::BulkData &mesh,
    const stk::mesh::Selector & elementSelector,
    const size_t elementGroupId,
    std::vector<size_t> & elementGroupIds,
    std::vector<stk::mesh::Entity> & attachedElemStack)
{
  while (!attachedElemStack.empty())
  {
    const stk::mesh::Entity elem = stack_pop(attachedElemStack);
    assign_group_id_and_if_new_add_neighbors_on_stack(mesh, elementSelector, elem, elementGroupId, elementGroupIds, attachedElemStack);
  }
}

void assign_local_group_id_for_each_element(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    std::vector<size_t> & elementGroupIds)
{
  std::vector<stk::mesh::Entity> attachedElemStack;
  attachedElemStack.reserve(ownedSelectedElements.size());
  for ( auto elem : ownedSelectedElements )
  {
    const auto elemOffset = elem.local_offset();
    if(elementGroupIds[elemOffset] == 0)
    {
      attachedElemStack.assign(1, elem);
      const size_t elementGroupId = mesh.identifier(elem);
      assign_element_group_id_and_recursively_check_neighbors(mesh, elementSelector, elementGroupId, elementGroupIds, attachedElemStack);
    }
  }
}

static std::vector<stk::mesh::Entity> get_elements_not_in_given_group(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    const size_t elementGroupId,
    const std::vector<size_t> & elementGroupIds)
{
  std::vector<stk::mesh::Entity> ownedElementsNotInGroup;
  for ( auto elem : ownedSelectedElements )
  {
    const auto elemOffset = elem.local_offset();
    if(elementGroupIds[elemOffset] != elementGroupId)
      ownedElementsNotInGroup.push_back(elem);
  }
  return ownedElementsNotInGroup;
}

void pack_element_group_ids_for_ghosting_procs(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> &elements,
    const std::vector<size_t> &localElementGroupIds,
    const std::map<size_t,size_t> &localToGlobalGroupIds,
    stk::CommSparse &commSparse)
{
  std::vector<int> elemCommProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto elem : elements)
    {
      if (mesh.bucket(elem).owned())
      {
        mesh.comm_procs(elem, elemCommProcs);
        const auto elemOffset = elem.local_offset();
        for (int procId : elemCommProcs)
        {
          if (procId != commSparse.parallel_rank())
          {
            commSparse.send_buffer(procId).pack(mesh.identifier(elem));
            commSparse.send_buffer(procId).pack(localToGlobalGroupIds.at(localElementGroupIds[elemOffset]));
          }
        }
      }
    }
  });
}

void receive_element_group_ids_from_owners_and_adjust_global_group_ids(const stk::mesh::BulkData &mesh,
    const std::vector<size_t> & localElementGroupIds,
    std::map<size_t,size_t> & localToGlobalGroupIds,
    bool & didSomethingChange,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::mesh::EntityId elemId = 0;
    size_t newGroupId = 0;
    commSparse.recv_buffer(procId).unpack(elemId);
    commSparse.recv_buffer(procId).unpack(newGroupId);

    const stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
    STK_ThrowRequire(mesh.is_valid(elem));
    const auto elemOffset = elem.local_offset();
    const size_t localGroupId = localElementGroupIds[elemOffset];

    if (localGroupId != 0)
    {
      auto iter = localToGlobalGroupIds.find(localGroupId);
      STK_ThrowRequire(iter != localToGlobalGroupIds.end());
      if(newGroupId < iter->second)
      {
        didSomethingChange = true;
        iter->second = newGroupId;
      }
    }
  });
}

void make_local_to_global_group_id_mapping_parallel_consistent(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    const std::vector<size_t> & localElementGroupIds,
    std::map<size_t,size_t> & localToGlobalGroupIds)
{
  bool didSomethingChange = true;
  while(didSomethingChange)
  {
    didSomethingChange = false;
    stk::CommSparse commSparse(mesh.parallel());
    pack_element_group_ids_for_ghosting_procs(mesh, ownedSelectedElements, localElementGroupIds, localToGlobalGroupIds, commSparse);
    receive_element_group_ids_from_owners_and_adjust_global_group_ids(mesh, localElementGroupIds, localToGlobalGroupIds, didSomethingChange, commSparse);
    didSomethingChange = stk::is_true_on_any_proc(mesh.parallel(), didSomethingChange);
  }
}

std::map<size_t,size_t> generate_parallel_consistent_local_to_global_group_id_mapping(const stk::mesh::BulkData &mesh,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    const std::vector<size_t> & localElementGroupIds)
{
  std::map<size_t,size_t> localToGlobalGroupIds;
  for (size_t localElementGroupId : localElementGroupIds)
    if (localElementGroupId != 0)
      localToGlobalGroupIds[localElementGroupId] = localElementGroupId;
  make_local_to_global_group_id_mapping_parallel_consistent(mesh, ownedSelectedElements, localElementGroupIds, localToGlobalGroupIds);
  return localToGlobalGroupIds;
}

static std::map<size_t,size_t> get_local_group_ids_sizes(const std::vector<size_t> & elementGroupIds)
{
  std::map<size_t,size_t> localGroupIdSizes;
  for (auto elementGroupId : elementGroupIds)
    if (elementGroupId != 0)
      ++(localGroupIdSizes[elementGroupId]);
  return localGroupIdSizes;
}

static void parallel_sum_group_id_sizes(std::map<size_t,size_t> & groupIdSizes, stk::ParallelMachine comm)
{
  if (stk::parallel_machine_size(comm) == 1)
    return;
  stk::CommSparse commSparse(comm);
  std::vector<size_t> localGroupSizes;
  localGroupSizes.reserve(2*groupIdSizes.size());
  for (auto & entry : groupIdSizes)
  {
    localGroupSizes.push_back(entry.first);
    localGroupSizes.push_back(entry.second);
  }

  std::vector<size_t> globalGroupSizes;
  stk::parallel_vector_concat(comm, localGroupSizes, globalGroupSizes);

  groupIdSizes.clear();
  for (size_t i=0; i<globalGroupSizes.size(); i+=2)
    groupIdSizes[globalGroupSizes[i]] += globalGroupSizes[i+1];
}

static size_t find_id_of_largest_group(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    const std::vector<size_t> & elementGroupIds)
{
  std::map<size_t,size_t> groupIdSizes = get_local_group_ids_sizes(elementGroupIds);
  parallel_sum_group_id_sizes(groupIdSizes, mesh.parallel());
  const auto iter = std::max_element(groupIdSizes.begin(), groupIdSizes.end(), [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b){ return a.second < b.second; });
  STK_ThrowRequire(iter != groupIdSizes.end());
  return iter->first;
}

void make_group_ids_parallel_consistent(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements,
    std::vector<size_t> & elementGroupIds)
{
  const std::map<size_t,size_t> localToGlobalGroupIds = generate_parallel_consistent_local_to_global_group_id_mapping(mesh, ownedSelectedElements, elementGroupIds);
  for (auto & elementGroupId : elementGroupIds)
    if (elementGroupId != 0)
      elementGroupId = localToGlobalGroupIds.at(elementGroupId);
}

static std::vector<size_t> determine_element_group_ids(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Entity> & ownedSelectedElements)
{
  const size_t initialValue = 0;
  std::vector<size_t> elementGroupIds = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, initialValue);
  assign_local_group_id_for_each_element(mesh, elementSelector, ownedSelectedElements, elementGroupIds);
  make_group_ids_parallel_consistent(mesh, ownedSelectedElements, elementGroupIds);
  return elementGroupIds;
}

std::vector<stk::mesh::Entity> find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector)
{
  std::vector<stk::mesh::Entity> ownedSelectedElements;
  stk::mesh::get_selected_entities( elementSelector & mesh.mesh_meta_data().locally_owned_part(), mesh.buckets( stk::topology::ELEMENT_RANK ), ownedSelectedElements, false );

  const std::vector<size_t> elementGroupIds = determine_element_group_ids(mesh, elementSelector, ownedSelectedElements);

  const size_t groupIdOfLargestGroup = find_id_of_largest_group(mesh, ownedSelectedElements, elementGroupIds);
  return get_elements_not_in_given_group(mesh, ownedSelectedElements, groupIdOfLargestGroup, elementGroupIds);
}

}
