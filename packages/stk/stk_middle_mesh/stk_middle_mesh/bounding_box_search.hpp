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
// 

#ifndef  STK_MIDDLE_MESH_BOUNDING_BOX_SEARCH_HPP
#define  STK_MIDDLE_MESH_BOUNDING_BOX_SEARCH_HPP

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/StaticAssert.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>
#include "search_mesh_element_bounding_box_base.hpp"
#include "search_mesh_element_bounding_box.hpp"
#include "search_mesh_element_bounding_box_normal.hpp"
#include "search_mesh_vertex.hpp"
#include "bounding_box_search_opts.hpp"

namespace stk {
namespace middle_mesh {
namespace search {
namespace impl {

template <typename T, typename U>
inline void inflate_bounding_box(stk::search::Sphere<T>& s, U const& mult_fact, U const& add_fact)
{
  s.set_radius(s.radius() * mult_fact + add_fact);
}

template <typename T, typename U>
inline void inflate_bounding_box(stk::search::Box<T>& b, U const& mult_fact, U const& add_fact)
{
  const double zero = 0.0;

  STK_ThrowRequire(mult_fact >= zero);
  STK_ThrowRequire(add_fact >= zero);

  stk::search::Point<T>& min_corner = b.min_corner();
  stk::search::Point<T>& max_corner = b.max_corner();

  T diag = 0.0;
  for(int i = 0; i < 3; ++i) {
    const T dd = max_corner[i] - min_corner[i];
    diag += dd * dd;
  }
  diag = std::sqrt(diag);

  //TODO: this isn't quite right: moving the max corner *and* the min corner
  //      by d results in a total change of 2d.
  const T d = mult_fact * diag + add_fact;

  for(int i = 0; i < 3; ++i) {
    min_corner[i] -= d;
    max_corner[i] += d;
  }
}

template <typename ForwardIterator, typename Compare>
bool local_is_sorted(ForwardIterator first, ForwardIterator last, Compare compare)
{
  if (first == last) return true;
  ForwardIterator next = first;
  while (++next!=last) {
      if ( compare(*next, *first) ) return false;
      ++first;
  }
  return true;
}

} // namespace impl

enum class SplitCommColor {
      RECV = 0,
      SEND,
      INVALID
};

template <class FROM, class TO> class BoundingBoxSearchType
{
public:
  using MeshA                     = FROM;
  using MeshB                     = TO;
  using EntityKeyA                = typename MeshA::EntityKey;
  using EntityKeyB                = typename MeshB::EntityKey;
  using EntityProcA               = typename MeshA::EntityProc;
  using EntityProcB               = typename MeshB::EntityProc;
  using EntityProcRelation        = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec     = std::vector<EntityProcRelation>;
};

template <typename SEARCH> class BoundingBoxSearch
{
public :

  using  SearchClass               = SEARCH;
  using  SendMesh                  = typename SearchClass::MeshA;
  using  RecvMesh                  = typename SearchClass::MeshB;

  using  MeshA                     = typename SearchClass::MeshA;
  using  MeshB                     = typename SearchClass::MeshB;
  using  EntityKeyA                = typename SearchClass::EntityKeyA;
  using  EntityKeyB                = typename SearchClass::EntityKeyB;

  using  EntityProcA               = typename SearchClass::EntityProcA;
  using  EntityProcB               = typename SearchClass::EntityProcB;

  using  EntityProcRelation        = typename SearchClass::EntityProcRelation;
  using  EntityProcRelationVec     = typename SearchClass::EntityProcRelationVec;

  using  EntityKeySetA             = typename std::set     <EntityKeyA>;
  using  EntityKeySetB             = typename std::set     <EntityKeyB>;
  using  BoundingBoxA              = typename MeshA::BoundingBox;
  using  BoundingBoxB              = typename MeshB::BoundingBox;

  using  Point                     = stk::search::Point<double>;

  enum {Dimension = 3};

  BoundingBoxSearch(std::shared_ptr<MeshA> sendMesh,
                    std::shared_ptr<MeshB> recvMesh,
                    const std::string &name,
                    stk::ParallelMachine pm,
                    const bool doParallelSearch=true,
                    const BoundingBoxSearchOpts& searchOpts = BoundingBoxSearchOpts(),

                    const stk::search::SearchMethod search_method = stk::search::KDTREE);
  virtual ~BoundingBoxSearch(){};
  void coarse_search() ;

  const std::shared_ptr<SendMesh> send_mesh() const {return m_sendMesh;}
  const std::shared_ptr<RecvMesh> recv_mesh() const {return m_recvMesh;}

  EntityProcRelationVec& get_range_to_domain() { return m_global_range_to_domain; }
  const EntityProcRelationVec& get_range_to_domain() const { return m_global_range_to_domain; }
  const std::vector<BoundingBoxB> get_unpaired_recv_entities() const { return m_unpaired_recv_entities; }

  stk::ParallelMachine get_comm() const {return m_shared_comm; }

private :
  void delete_range_points_found(std::vector<BoundingBoxB>& range_vector, const EntityProcRelationVec& del);

  std::shared_ptr<MeshA>               m_sendMesh;
  std::shared_ptr<MeshB>               m_recvMesh;

  const std::string     m_name;
  const bool            m_doParallelSearch;
  const BoundingBoxSearchOpts m_searchOpts;
  const stk::search::SearchMethod m_search_method;
  stk::ParallelMachine m_shared_comm;

  EntityProcRelationVec m_global_range_to_domain;
  SearchClass           m_search;

  std::vector<BoundingBoxB> m_unpaired_recv_entities;

  struct compare {
    bool operator()(const BoundingBoxB &a, const EntityProcB  &b) const
    {
      return a.second < b;
    }

    bool operator()(const EntityProcB  &a, const BoundingBoxB &b) const
    {
      return a < b.second;
    }
  };

  template <class BoundingBoxType>
  struct BoundingBoxCompare{
    bool operator()(const BoundingBoxType & a, const BoundingBoxType & b) const
    {
      return a.second.id() < b.second.id();
    }
  };
};

template <typename SEARCH> BoundingBoxSearch<SEARCH>::BoundingBoxSearch
(std::shared_ptr<MeshA> sendMesh,
 std::shared_ptr<MeshB> recvMesh,
 const std::string        &name,
 stk::ParallelMachine pm,
 const bool doParallelSearch,
 const BoundingBoxSearchOpts& searchOpts,

 const stk::search::SearchMethod search_method) :
      m_sendMesh(sendMesh), m_recvMesh(recvMesh), m_name(name), m_doParallelSearch(doParallelSearch), 
      m_searchOpts(searchOpts),
      m_search_method(search_method),
      m_shared_comm(pm)
  {
  }

template <class SEARCH> void BoundingBoxSearch<SEARCH>::delete_range_points_found
(std::vector<BoundingBoxB>& range_vector, const EntityProcRelationVec& del)
{
  std::vector<EntityProcB> range_entities_found;
  range_entities_found.reserve(del.size());
  for (auto && i : del) {
    range_entities_found.push_back(i.first);
  }
  {
    std::sort(range_entities_found.begin(), range_entities_found.end());
    const typename std::vector<EntityProcB>::iterator it = std::unique(range_entities_found.begin(), range_entities_found.end());
    range_entities_found.resize(it-range_entities_found.begin());
  }

  std::vector<BoundingBoxB> difference(range_vector.size());
  {
    const auto it =
      std::set_difference(
        range_vector.        begin(), range_vector.        end(),
        range_entities_found.begin(), range_entities_found.end(),
        difference.begin(), compare());
    difference.resize(it-difference.begin());
  }
  swap(difference, range_vector);
}

template <class SEARCH> void BoundingBoxSearch<SEARCH>::coarse_search()
{
  m_global_range_to_domain.clear();
  m_unpaired_recv_entities.clear();

  std::vector<BoundingBoxA> domain_vector;
  std::vector<BoundingBoxB> range_vector;

  if(m_sendMesh) m_sendMesh->fill_bounding_boxes(domain_vector);
  if(m_recvMesh) m_recvMesh->fill_bounding_boxes(range_vector);

  if(!impl::local_is_sorted(domain_vector.begin(), domain_vector.end(), BoundingBoxCompare<BoundingBoxA>()))
    std::sort(domain_vector.begin(), domain_vector.end(), BoundingBoxCompare<BoundingBoxA>());
  if(!impl::local_is_sorted(range_vector.begin(), range_vector.end(), BoundingBoxCompare<BoundingBoxB>()))
    std::sort(range_vector.begin(), range_vector.end(), BoundingBoxCompare<BoundingBoxB>());

  int not_empty_count = 0;

  unsigned range_vector_not_empty = !range_vector.empty();
  stk::all_reduce(m_shared_comm, stk::ReduceSum<1>(&range_vector_not_empty));

  for(BoundingBoxA& i : domain_vector) {
    impl::inflate_bounding_box(i.first, m_searchOpts.initialExpansionFactor, m_searchOpts.initialExpansionSum);
  }

  while(range_vector_not_empty) { // Keep going until all range points are processed.
    // Slightly confusing: coarse_search documentation has domain->range
    // relations sorted by domain key.  We want range->domain type relations
    // sorted on range key. It might appear we have the arguments revered
    // in coarse_search call, but really, this is what we want.
    EntityProcRelationVec rng_to_dom;
    EntityProcRelationVec& rng_to_dom_vec = not_empty_count == 0 ? m_global_range_to_domain : rng_to_dom;
    stk::search::coarse_search(range_vector, domain_vector, m_search_method, m_shared_comm, rng_to_dom_vec, m_doParallelSearch);

    if(not_empty_count > 0) {
      m_global_range_to_domain.insert(m_global_range_to_domain.end(), rng_to_dom_vec.begin(), rng_to_dom_vec.end());
    }

    not_empty_count++;

    delete_range_points_found(range_vector, rng_to_dom_vec);

    range_vector_not_empty = !range_vector.empty();
    stk::all_reduce(m_shared_comm, stk::ReduceSum<1>(&range_vector_not_empty));

    if(m_searchOpts.repeatExpansionFactor <= 1.0) {
      if(range_vector_not_empty) {
        size_t range_vector_size = range_vector.size();
        size_t g_range_vector_size = stk::get_global_sum(m_shared_comm, range_vector_size);
        sierra::Env::outputP0() << m_name
                                << ": BoundingBoxSearch<SEARCH>::coarse_search(): Number of points not found: "
                                << g_range_vector_size << " in initial coarse search" << std::endl;
      }

      break;
    }

    if(range_vector_not_empty) {
      for(BoundingBoxB& i : range_vector) {
        impl::inflate_bounding_box(i.first, m_searchOpts.repeatExpansionFactor, m_searchOpts.repeatExpansionSum);
      }
      for(BoundingBoxA& i : domain_vector) {
        impl::inflate_bounding_box(i.first, m_searchOpts.repeatExpansionFactor, m_searchOpts.repeatExpansionSum);
      }

      size_t range_vector_size = range_vector.size();
      size_t g_range_vector_size = 0;
      stk::all_reduce_max(m_shared_comm, &range_vector_size, &g_range_vector_size, 1);
      sierra::Env::outputP0() << m_name
                              << ": BoundingBoxSearch<SEARCH>::coarse_search(): Number of points not found: "
                              << g_range_vector_size << " after expanding bounding boxes: " << not_empty_count
                              << " time(s)" << std::endl;
      sierra::Env::outputP0()
          << "...will now expand the set of candidate bounding boxes and re-attempt the coarse search" << std::endl;
    }
  }

  m_unpaired_recv_entities.swap(range_vector);
  std::sort(m_global_range_to_domain.begin(), m_global_range_to_domain.end());
}


using ElementToElementBoundingBoxSearch =
    BoundingBoxSearch<BoundingBoxSearchType<mesh::impl::SearchMeshElementBoundingBoxBase,
                                            mesh::impl::SearchMeshElementBoundingBoxBase>>;

using ElementToVertBoundingBoxSearch = 
    BoundingBoxSearch<BoundingBoxSearchType<mesh::impl::SearchMeshElementBoundingBoxBase,
                                            mesh::impl::SearchMeshVertex>>;
}
}
}

#endif
