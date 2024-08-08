/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_PATCH_HPP
#define STK_TRANSFER_UTIL_PATCH_HPP

#include "stk_mesh/base/Types.hpp"         // for EntityRank, etc
#include "stk_topology/topology.hpp"       // for topology, etc
#include "stk_util/util/ReportHandler.hpp" // for STK_ThrowAssert, etc
#include <algorithm>                       // for lower_bound
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <vector> // for vector, vector<>::iterator, etc

namespace stk {
namespace transfer {

/// A Filter object that always returns true.
struct FilterNone {
  bool pass(stk::mesh::Entity, const stk::mesh::BulkData&) const { return true; }
  ~FilterNone() {}
};

namespace {

inline void
insert_unique_entity(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& list, stk::mesh::Entity entity)
{
  std::vector<stk::mesh::Entity>::iterator i =
      lower_bound(list.begin(), list.end(), entity, stk::mesh::EntityLess(bulkData));
  if(i == list.end() || *i != entity) {
    list.insert(i, entity);
  }
}

template <typename FILTER>
inline void patch_if(const stk::mesh::BulkData& bulkData, stk::mesh::Entity patchSeed, const FILTER& filter,
                     const stk::mesh::Selector& selector, std::vector<stk::mesh::Entity>& elements,
                     std::vector<stk::mesh::Entity>& nodes, const stk::mesh::EntityRank patchType)
{
  stk::mesh::EntityRank patchSeedRank = bulkData.entity_rank(patchSeed);
  bool invalid = (patchSeedRank != stk::topology::NODE_RANK && patchType != stk::topology::CONSTRAINT_RANK);

  STK_ThrowRequireMsg(!invalid, "\n\nERROR in " << __FUNCTION__ << "\n"
                                            << "Calling object must be a node or patch type must be CONSTRAINT; "
                                            << "instead, this function was called with a CONSTRAINT patch and \n"
                                            << "an entity of type '" << patchSeedRank << "'\n");

  elements.clear();
  nodes.clear();

  // loop on active or scratch ghosted elements connected to node
  stk::mesh::Entity const* targets = bulkData.begin(patchSeed, patchType);
  const int numConn = bulkData.num_connectivity(patchSeed, patchType);

  for(int eItr = 0; eItr < numConn; ++eItr) {
    stk::mesh::Entity elem = targets[eItr];
    if(selector(bulkData.bucket(elem))) {
      if(filter.pass(elem, bulkData)) insert_unique_entity(bulkData, elements, elem);
    }
  }

  // gather nodes touched by the element patch
  for(auto&& elem : elements) {
    stk::mesh::Entity const* entityRelation = bulkData.begin(elem, patchSeedRank);
    for(int i = 0, ie = bulkData.num_connectivity(elem, patchSeedRank); i < ie; ++i) {
      insert_unique_entity(bulkData, nodes, entityRelation[i]);
    }
  }
}

}

template <typename FILTER>
inline void nodes_patch_if(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::Entity>& nodes,
                           const FILTER& predicate, const stk::mesh::Selector& selector,
                           std::vector<stk::mesh::Entity>& elementsInPatch,
                           std::vector<stk::mesh::Entity>& nodesInPatch,
                           const stk::mesh::EntityRank patchType = stk::topology::ELEMENT_RANK)
{
  std::vector<stk::mesh::Entity> nodalPatch;
  std::vector<stk::mesh::Entity> nodesInNodalPatch;
  for(stk::mesh::Entity node : nodes) {
    if(bulkData.entity_rank(node) != stk::topology::NODE_RANK) continue;

    patch_if(bulkData, node, predicate, selector, nodalPatch, nodesInNodalPatch, patchType);

    // Insert unique entries into object-patch lists
    for(auto&& entity : nodalPatch) {
      insert_unique_entity(bulkData, elementsInPatch, entity);
    }

    for(auto& nodeInNodalPatch : nodesInNodalPatch) {
      insert_unique_entity(bulkData, nodesInPatch, nodeInNodalPatch);
    }
  }
}

/**
 * As entity_patch_if(), but only includes an entity in the patch if predicate(entity) is true.
 *
 * @param patch (returned) set of entities in the patch; for these entities,
 *              predicate(entity) is true
 *
 * @param nodes_in_patch (returned) list of all nodes used by entities in 'patch'
 */

template <typename FILTER>
inline void entity_patch_if(const stk::mesh::BulkData& bulkData, stk::mesh::Entity patchSeed, const FILTER& predicate,
                            const stk::mesh::Selector& selector, std::vector<stk::mesh::Entity>& elementsInPatch,
                            std::vector<stk::mesh::Entity>& nodesInPatch,
                            const stk::mesh::EntityRank patchType = stk::topology::ELEMENT_RANK)
{
  if(bulkData.entity_rank(patchSeed) == stk::topology::NODE_RANK || patchType == stk::topology::CONSTRAINT_RANK) {
    patch_if(bulkData, patchSeed, predicate, selector, elementsInPatch, nodesInPatch, patchType);
    return;
  }

  std::vector<stk::mesh::Entity> nodes(bulkData.begin_nodes(patchSeed), bulkData.end_nodes(patchSeed));
  nodes_patch_if(bulkData, nodes, predicate, selector, elementsInPatch, nodesInPatch, patchType);
}

template <typename FILTER>
class Patch {
 public:
  Patch(const stk::mesh::BulkData& bulk, stk::mesh::Entity seed, const FILTER& predicate, const stk::mesh::Selector& selector)
    : m_bulk(bulk)
    , m_filter(predicate)
    , m_selector(selector)
  {

  }

  virtual ~Patch() {}

  const std::vector<stk::mesh::Entity>& get_patch_nodes() const { return m_patchNodes; }
  const std::vector<stk::mesh::Entity>& get_patch_elements() const { return m_patchElements; }
  const std::vector<stk::mesh::Entity>& get_patch_entities() const
  {
    return m_bulk.entity_rank(m_patchSeed) == stk::topology::NODE_RANK ? m_patchNodes : m_patchElements;
  }
  size_t size() const { return get_patch_entities().size(); }

  stk::mesh::Entity get_patch_seed() const { return m_patchSeed; }
  const stk::mesh::BulkData& get_bulk_data() const { return m_bulk; }

 protected:
  const stk::mesh::BulkData& m_bulk;
  stk::mesh::Entity m_patchSeed;
  const FILTER& m_filter;
  const stk::mesh::Selector m_selector;

  std::vector<stk::mesh::Entity> m_patchElements;
  std::vector<stk::mesh::Entity> m_patchNodes;
};

template <typename FILTER>
class LinearPatch : public Patch<FILTER> {
 public:
  LinearPatch(const stk::mesh::BulkData& bulk, stk::mesh::Entity seed, const FILTER& predicate, const stk::mesh::Selector& selector)
    : Patch<FILTER>(bulk, seed, predicate, selector)
  {
    construct_patch(seed);
  }

  void construct_patch(stk::mesh::Entity seed)
  {
    stk::mesh::EntityRank rank = Patch<FILTER>::m_bulk.entity_rank(seed);
    STK_ThrowRequireMsg(rank == stk::topology::NODE_RANK || rank == stk::topology::ELEM_RANK,
                        "Input patch seed must be an element or node");

    Patch<FILTER>::m_patchSeed = seed;

    entity_patch_if(Patch<FILTER>::m_bulk, seed, Patch<FILTER>::m_filter, Patch<FILTER>::m_selector,
                    Patch<FILTER>::m_patchElements, Patch<FILTER>::m_patchNodes, stk::topology::ELEMENT_RANK);

    STK_ThrowRequireMsg(!Patch<FILTER>::m_patchElements.empty(),
                        "Entity " << Patch<FILTER>::m_bulk.entity_key(seed) << " on processor '"
                        << Patch<FILTER>::m_bulk.parallel_rank() << "' has an empty patch\n");
  }
};

template <typename FILTER>
using QuadraticPatch = LinearPatch<FILTER>;

template <typename FILTER>
class CubicPatch : public Patch<FILTER> {
 public:
  CubicPatch(const stk::mesh::BulkData& bulk, stk::mesh::Entity seed, const FILTER& predicate, const stk::mesh::Selector& selector)
    : Patch<FILTER>(bulk, seed, predicate, selector)
  {
    construct_patch(seed);
  }

  void construct_patch(stk::mesh::Entity seed)
  {
    stk::mesh::EntityRank rank = Patch<FILTER>::m_bulk.entity_rank(seed);
    STK_ThrowRequireMsg(rank == stk::topology::NODE_RANK || rank == stk::topology::ELEM_RANK,
                        "Input patch seed must be an element or node");

    Patch<FILTER>::m_patchSeed = seed;

    entity_patch_if(Patch<FILTER>::m_bulk, seed, Patch<FILTER>::m_filter, Patch<FILTER>::m_selector,
                    Patch<FILTER>::m_patchElements, Patch<FILTER>::m_patchNodes, stk::topology::ELEMENT_RANK);

    std::vector<stk::mesh::Entity> neighborPatchNodes = Patch<FILTER>::m_patchNodes;

    nodes_patch_if(Patch<FILTER>::m_bulk, neighborPatchNodes, Patch<FILTER>::m_filter, Patch<FILTER>::m_selector,
                    Patch<FILTER>::m_patchElements, Patch<FILTER>::m_patchNodes, stk::topology::ELEMENT_RANK);

    STK_ThrowRequireMsg(!Patch<FILTER>::m_patchElements.empty(),
                        "Entity " << Patch<FILTER>::m_bulk.entity_key(seed) << " on processor '"
                        << Patch<FILTER>::m_bulk.parallel_rank() << "' has an empty patch\n");
  }
};

} // namespace transfer
} // namespace stk

#endif
