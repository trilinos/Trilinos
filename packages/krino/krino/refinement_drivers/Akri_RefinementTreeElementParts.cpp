// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RefinementTreeElementParts.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_RefinementManager.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <algorithm>
#include <functional>

namespace krino {


static void fill_sorted_owned_leaf_children_and_their_element_part_and_their_unique_parents(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & nextIterationOfParentsToProcess)
{
  const stk::mesh::Selector adaptivityLeafChild = mesh.mesh_meta_data().locally_owned_part() & refinement.child_part() & !refinement.parent_part();
  const auto & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, adaptivityLeafChild);

  for( auto&& bucket : buckets )
  {
    for (auto && elem : *bucket)
    {
      sortedChildIdsAndElementPart.emplace_back(mesh.identifier(elem), find_element_part(mesh, elem).mesh_meta_data_ordinal());
      nextIterationOfParentsToProcess.push_back(refinement.get_parent(elem));
    }
  }
  std::sort(sortedChildIdsAndElementPart.begin(), sortedChildIdsAndElementPart.end());
  stk::util::sort_and_unique(nextIterationOfParentsToProcess);
}

static const unsigned * find_child_element_part_ordinal(const std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart, stk::mesh::EntityId childId)
{
  auto compare = [](const std::pair<stk::mesh::EntityId,unsigned> & idAndElementPart, const stk::mesh::EntityId id) { return idAndElementPart.first < id; };
  auto it = std::lower_bound(sortedChildIdsAndElementPart.begin(), sortedChildIdsAndElementPart.end(), childId, compare);
  if (it != sortedChildIdsAndElementPart.end() && it->first == childId)
    return &(it->second);
  return nullptr;
}

static std::set<unsigned> get_child_element_part_ordinals(const std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    const std::vector<stk::mesh::EntityId> & childElemIds)
{
  std::set<unsigned> childElemPartOrdinals;
  for (auto&& childElemId : childElemIds)
  {
    const unsigned * childElemOrd = find_child_element_part_ordinal(sortedChildIdsAndElementPart, childElemId);
    if (nullptr == childElemOrd)
    {
      // not ready to process, some children not yet available
      childElemPartOrdinals.clear();
      break;
    }
    childElemPartOrdinals.insert(*childElemOrd);
  }
  return childElemPartOrdinals;
}

struct RelocatedChildElemWithElemPart {
  int parentElemOwningProc;
  stk::mesh::EntityId childElemId;
  stk::mesh::EntityId parentElemId;
  unsigned childElemPartOrdinal;
};

static
void pack_child_element_data_for_parent_element_processor(const stk::mesh::BulkData & /*mesh*/,
    const std::vector<RelocatedChildElemWithElemPart> & childElemDataToCommunicate,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto && data : childElemDataToCommunicate)
    {
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.childElemId);
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.parentElemId);
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.childElemPartOrdinal);
    }
  });
}

static
void unpack_child_element_data(const stk::mesh::BulkData & mesh,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & iterationChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId childElemId;
      commSparse.recv_buffer(procId).unpack(childElemId);
      stk::mesh::EntityId parentElemId;
      commSparse.recv_buffer(procId).unpack(parentElemId);
      unsigned childElemPartOrdinal;
      commSparse.recv_buffer(procId).unpack(childElemPartOrdinal);
      iterationChildIdsAndElementPart.emplace_back(childElemId, childElemPartOrdinal);
      const stk::mesh::Entity parentElem = mesh.get_entity(stk::topology::ELEMENT_RANK, parentElemId);
      parentsToProcess.push_back(parentElem);
    }
  });
}

static void communicate_elements_to_parent_element_processor(const stk::mesh::BulkData& mesh,
    const std::vector<RelocatedChildElemWithElemPart> & childElemDataToCommunicate,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & iterationChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_child_element_data_for_parent_element_processor(mesh, childElemDataToCommunicate, commSparse);
  unpack_child_element_data(mesh, iterationChildIdsAndElementPart, parentsToProcess, commSparse);
}

static void update_children_and_their_element_part_and_their_unique_parents(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const std::vector<std::pair<stk::mesh::Entity,unsigned>> & iterationChildrenAndElementPart,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess)
{
  parentsToProcess.clear();

  std::vector<std::pair<stk::mesh::EntityId,unsigned>> iterationChildIdsAndElementPart;
  iterationChildIdsAndElementPart.reserve(iterationChildrenAndElementPart.size());

  std::vector<RelocatedChildElemWithElemPart> childElemDataToCommunicate;

  for (auto && childAndElementPart : iterationChildrenAndElementPart)
  {
    const stk::mesh::Entity child = childAndElementPart.first;
    const stk::mesh::EntityId childId = mesh.identifier(child);
    iterationChildIdsAndElementPart.emplace_back(childId, childAndElementPart.second);
    const auto parentIdAndOwnerRank = refinement.get_parent_id_and_parallel_owner_rank(child);
    if (parentIdAndOwnerRank.second == mesh.parallel_rank())
      parentsToProcess.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK, parentIdAndOwnerRank.first));
    else
      childElemDataToCommunicate.push_back(RelocatedChildElemWithElemPart{parentIdAndOwnerRank.second, childId, parentIdAndOwnerRank.first, childAndElementPart.second});
  }

  communicate_elements_to_parent_element_processor(mesh, childElemDataToCommunicate, iterationChildIdsAndElementPart, parentsToProcess);

  std::sort(iterationChildIdsAndElementPart.begin(), iterationChildIdsAndElementPart.end());
  stk::util::insert_keep_sorted(iterationChildIdsAndElementPart, sortedChildIdsAndElementPart, std::less<std::pair<stk::mesh::EntityId,unsigned>>());
  stk::util::sort_and_unique(parentsToProcess);
}

std::vector<std::pair<stk::mesh::Entity,unsigned>> get_owned_refinement_parents_and_their_element_part(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const std::function<unsigned(const std::set<unsigned> &)> & determineParentElementPartOrdinal)
{
  // This used to be pretty simple.  We just recursed down to the leaves of the tree and checked the block
  // parts of the conforming elements.  But now the tree may be spread across processors, so we need
  // to communicate this information up the tree.

  std::vector<std::pair<stk::mesh::Entity,unsigned>> parentsAndElementPart;

  std::vector<std::pair<stk::mesh::EntityId,unsigned>> sortedChildIdsAndElementPart;
  std::vector<stk::mesh::Entity> parentsToProcess;

  fill_sorted_owned_leaf_children_and_their_element_part_and_their_unique_parents(mesh, refinement, sortedChildIdsAndElementPart, parentsToProcess);

  std::vector<stk::mesh::EntityId> childElemIds;
  std::vector<std::pair<stk::mesh::Entity,unsigned>> iterationChildrenAndElementPart;

  while (stk::is_true_on_any_proc(mesh.parallel(), !parentsToProcess.empty()))
  {
    iterationChildrenAndElementPart.clear();

    for( auto&& parent : parentsToProcess )
    {
      refinement.fill_child_element_ids(parent, childElemIds);
      const std::set<unsigned> childElemPartOrdinals = get_child_element_part_ordinals(sortedChildIdsAndElementPart, childElemIds);
      if (!childElemPartOrdinals.empty())
      {
        const unsigned elementPartOrdinal = determineParentElementPartOrdinal(childElemPartOrdinals);
        parentsAndElementPart.emplace_back(parent, elementPartOrdinal);

        if (refinement.is_child(parent))
          iterationChildrenAndElementPart.emplace_back(parent, elementPartOrdinal);
      }
    }

    update_children_and_their_element_part_and_their_unique_parents(mesh, refinement, iterationChildrenAndElementPart, sortedChildIdsAndElementPart, parentsToProcess);
  }
  return parentsAndElementPart;
}

} // namespace krino
