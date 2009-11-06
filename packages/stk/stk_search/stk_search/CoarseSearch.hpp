#ifndef stk_search_CoarseSearch_hpp
#define stk_search_CoarseSearch_hpp

#include <vector>
#include <utility>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/BihTree.hpp>
#include <stk_search/BihTreeParallelOps.hpp>
#include <stk_search/OctTreeOps.hpp>

namespace stk {
namespace search {

// Structure to hold factory specification
struct FactoryOrder {
  // Enumerate search algorithms
  enum Algorithm {OCTREE, BIHTREE};

  // Define virtual base class to derive holder for input processor cuts.
  class Cuts {
    Cuts(){}
    // assignment operator not needed
    // copy constructor not needed
    virtual ~Cuts();
  };

  Algorithm             m_algorithm;
  stk::ParallelMachine  m_communicator;
  Cuts *                m_cuts;
  unsigned *            m_search_tree_stats;
  unsigned              m_debug_level;
  FactoryOrder() :
    m_algorithm        (OCTREE),
    m_communicator     (MPI_COMM_SELF),
    m_cuts             (NULL),
    m_search_tree_stats(NULL),
    m_debug_level      (0)
  {}
};

template <class RangeBoundingVolume,class DomainBoundingVolume>
bool parallel_bihtree_search(
    std::vector<std::pair<typename DomainBoundingVolume::Key,typename RangeBoundingVolume::Key> > & domain_to_range_keys,
    std::vector<RangeBoundingVolume> &range,
    std::vector<DomainBoundingVolume> &domain,
    const FactoryOrder & order)
{

  //find the global box
  std::vector<float> global_box(6);
  stk::search::box_global_bounds(order.m_communicator,
      domain.size() ,
      &domain[0] ,
      range.size(),
      &range[0],
      &global_box[0]);

  //
  stk::search::oct_tree_bih_tree_proximity_search(order.m_communicator,
      &global_box[0],
      domain.size() ,
      &domain[0] ,
      range.size(),
      &range[0],
      NULL,
      domain_to_range_keys ,
      order.m_search_tree_stats);

  return true;
}


template <class RangeBoundingVolume,class DomainBoundingVolume>
bool coarse_search_bihtree(
    std::vector<std::pair<typename DomainBoundingVolume::Key,typename RangeBoundingVolume::Key> > & domain_to_range_keys,
    std::vector<RangeBoundingVolume> &range,
    std::vector<DomainBoundingVolume> &domain,
    const FactoryOrder & order)
{
  const unsigned p_size = parallel_machine_size( order.m_communicator );
  //const unsigned p_rank = parallel_machine_rank( order.m_communicator );

  if (p_size == 1) {
    bih::BihTree<RangeBoundingVolume> tree(range.begin(),range.end());
    unsigned size = domain.size();
    for (unsigned i = 0; i < size ; ++i) {
      tree.intersect(domain[i],domain_to_range_keys);
    }
  }
  else {
    parallel_bihtree_search(domain_to_range_keys,range,domain,order);
  }
  return true;
}

template <class RangeBoundingVolume,class DomainBoundingVolume>
bool coarse_search_octree(
    std::vector<std::pair<typename DomainBoundingVolume::Key,typename RangeBoundingVolume::Key> > & domain_to_range_keys,
    std::vector<RangeBoundingVolume> &range,
    std::vector<DomainBoundingVolume> &domain,
    const FactoryOrder & order)
{

  std::vector<float> global_box(6);
  stk::search::box_global_bounds(order.m_communicator,
      domain.size() ,
      &domain[0] ,
      range.size(),
      &range[0],
      &global_box[0]);

  stk::search::oct_tree_proximity_search(order.m_communicator,
      &global_box[0],
      domain.size() ,
      &domain[0] ,
      range.size(),
      &range[0],
      NULL,
      domain_to_range_keys ,
      order.m_search_tree_stats);
  return true;
}

template <class RangeBoundingVolume,class DomainBoundingVolume>
bool coarse_search(
    std::vector<std::pair<typename DomainBoundingVolume::Key,typename RangeBoundingVolume::Key> > & domain_to_range_keys,
    std::vector<RangeBoundingVolume> &range,
    std::vector<DomainBoundingVolume> &domain,
    const FactoryOrder & order)
{
  switch (order.m_algorithm) {
    case FactoryOrder::BIHTREE:
      return coarse_search_bihtree(domain_to_range_keys,range,domain,order);
    case FactoryOrder::OCTREE:
      return coarse_search_octree(domain_to_range_keys,range,domain,order);
   break;
  }
  return false;
}


} // namespace search
} // namespace stk

#endif // stk_search_CoarseSearch_hpp
