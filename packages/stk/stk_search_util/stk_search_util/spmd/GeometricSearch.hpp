/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCH_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <ostream>   // for std::ostream
#include <limits>    // for std::numeric_limits
#include <algorithm> // for std::sort

#include "stk_search/CoarseSearch.hpp"
#include "stk_search/FilterCoarseSearch.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include "stk_util/util/ReportHandler.hpp"
#include "stk_search_util/Inspector.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {
namespace impl {

inline void throw_error_on_coarse_expansion_limit(stk::ParallelMachine comm, const std::string local_err_msg)
{
  if (stk::util::get_common_coupling_version() >= 15) {
    std::ostringstream global_err_msg;
    stk::all_write_string(comm, global_err_msg, local_err_msg);
    int err = (int) global_err_msg.str().size();
    int global_err = 0;
    stk::all_reduce_sum(comm, &err, &global_err, 1);
    STK_ThrowRequireMsg(!global_err, global_err_msg.str());
  }
}

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

  const T d = mult_fact * diag + add_fact;

  for(int i = 0; i < 3; ++i) {
    min_corner[i] -= d;
    max_corner[i] += d;
  }
}

template <typename T>
inline stk::search::Point<T> compute_centroid(const stk::search::Sphere<T>& s)
{
  return s.center();
}

template <typename T>
inline stk::search::Point<T> compute_centroid(const stk::search::Box<T>& b)
{
  const stk::search::Point<T>& min_corner = b.min_corner();
  const stk::search::Point<T>& max_corner = b.max_corner();

  stk::search::Point<T> centroid;

  for(int i = 0; i < 3; ++i) {
    centroid[i] = (max_corner[i] + min_corner[i]) / 2;
  }

  return centroid;
}

template <typename RECVMESH>
struct compare {
  using EntityProcB = typename RECVMESH::EntityProc;
  using BoundingBoxB = typename RECVMESH::BoundingBox;

  bool operator()(const BoundingBoxB &a, const EntityProcB  &b) const
  {
    return a.second < b;
  }

  bool operator()(const EntityProcB  &a, const BoundingBoxB &b) const
  {
    return a < b.second;
  }
};

template <typename BoundingBox>
struct BoundingBoxCompare{
  bool operator()(const BoundingBox & a, const BoundingBox & b) const
  {
    return a.second.id() < b.second.id();
  }
};

template <typename ForwardIterator, typename Compare>
bool local_is_sorted(ForwardIterator first, ForwardIterator last, Compare comparison)
{
  if (first == last) return true;
  ForwardIterator next = first;
  while (++next!=last) {
    if ( comparison(*next, *first) ) return false;
    ++first;
  }
  return true;
}
} // namespace impl

class SearchBase {
public :
  SearchBase(){};
  virtual ~SearchBase(){};

  virtual void initialize()    = 0;
  virtual void coarse_search() = 0;
  virtual void local_search()  = 0;
};

template <typename SENDMESH, typename RECVMESH>
class GeometricSearch : public SearchBase {
 public:
  using EntityA = typename SENDMESH::Entity;
  using EntityB = typename RECVMESH::Entity;
  using EntityKeyA = typename SENDMESH::EntityKey;
  using EntityKeyB = typename RECVMESH::EntityKey;
  using EntityProcA = typename SENDMESH::EntityProc;
  using EntityProcB = typename RECVMESH::EntityProc;
  using BoundingBoxA = typename SENDMESH::BoundingBox;
  using BoundingBoxB = typename RECVMESH::BoundingBox;

  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;
  using EntityProcVecA = std::vector<EntityProcA>;
  using EntityProcVecB = std::vector<EntityProcB>;

  using EntityKeyVecA = std::vector<EntityKeyA>;
  using EntityKeyVecB = std::vector<EntityKeyB>;

  GeometricSearch(std::shared_ptr<SENDMESH>& mesha, std::shared_ptr<RECVMESH>& meshb, const std::string& name,
                  stk::ParallelMachine pm, const double expansionFactor = 0.5, const double expansionSum = 0.0,
                  const stk::search::SearchMethod searchMethod = stk::search::KDTREE)
    : m_expansionFactor(expansionFactor)
    , m_expansionSum(expansionSum)
    , m_useNearestNodeForClosestBoundingBox(false)
    , m_useCentroidForGeometricProximity(false)
    , m_printSearchWarnings(true)
    , m_doInitialSearchExpansion(true)
    , m_comm(pm)
    , m_initialized(false)
    , m_mesha(mesha)
    , m_meshb(meshb)
    , m_searchMethod(searchMethod)
    , m_isListenerActive(true)
    , m_syncCount(0)
    , m_name(name)
  {
  }

  virtual void coarse_search() override;
  virtual void local_search() override;

  virtual void initialize() override;
  void reinitialize();

  EntityProcRelationVec& get_range_to_domain() { return m_globalRangeToDomain; }
  const EntityProcRelationVec& get_range_to_domain() const { return m_globalRangeToDomain; }
  const std::vector<BoundingBoxB> get_unpaired_recv_entities() const { return m_unpairedRecvEntities; }

  stk::search::FilterCoarseSearchResult<RECVMESH>& get_search_filter_result() { return m_filteredSearchResult; }
  const stk::search::FilterCoarseSearchResult<RECVMESH>& get_search_filter_result() const { return m_filteredSearchResult; }

  void coarse_search_for_objects_outside_domain(bool includeGhosts=false);
  void coarse_search_for_objects_inside_domain();
  void determine_entities_to_copy(EntityProcVecA& entities_to_copy) const;
  void ghost_from_elements();
  void update_ghosted_keys();
  void post_mesh_modification_event(bool activeFlag);

  void set_closest_bounding_box_using_nearest_node(bool flag) { m_useNearestNodeForClosestBoundingBox = flag; }
  bool get_closest_bounding_box_using_nearest_node() const { return m_useNearestNodeForClosestBoundingBox; }

  void set_use_centroid_for_geometric_proximity(bool flag) { m_useCentroidForGeometricProximity = flag; }
  bool get_use_centroid_for_geometric_proximity() const { return m_useCentroidForGeometricProximity; }

  void set_do_initial_search_expansion(bool flag) { m_doInitialSearchExpansion = flag; }
  bool get_do_initial_search_expansion() const { return m_doInitialSearchExpansion; }

  void set_print_search_warnings(bool flag) { m_printSearchWarnings = flag; }
  bool get_print_search_warnings() const { return m_printSearchWarnings; }

  void set_output_stream(std::ostream& out) { m_outputStream = &out; }
  std::ostream& get_output_stream() { return *m_outputStream; }

  const std::shared_ptr<SENDMESH> send_mesh() const { return m_mesha; }
  const std::shared_ptr<RECVMESH> recv_mesh() const { return m_meshb; }

  void register_inspector(const std::string& fileName, const EntityKeyVecB& rangeEntities);
  void inspect_user_defined_entities();
  void register_output_formatter(std::shared_ptr< stk::search::InspectorOutputFormatterBase<SENDMESH, RECVMESH> > outputFormatter);

  void set_fractional_limit_for_objects_outside_domain(double f) {
    STK_ThrowRequireMsg(f > 0, "Fractional limit: " << f << " must be greater than 0.0");
    m_fractionalLimitForObjectsOutsideDomain = f;
  }
  double get_fractional_limit_for_objects_outside_domain() const { return m_fractionalLimitForObjectsOutsideDomain; }

  void set_expansion_limit(double expansionLimit) {
    STK_ThrowRequireMsg(expansionLimit > 0, "Expansion limit: " << expansionLimit << " must be greater than 0.0");
    m_expansionLimit = expansionLimit;
  }
  double get_expansion_limit() const { return m_expansionLimit; }

  void set_expansion_factor(const double expansionFactor) { m_expansionFactor = expansionFactor; }
  double get_expansion_factor() const { return m_expansionFactor; }

  void set_expansion_padding(const double expansionPadding) { m_expansionSum = expansionPadding; }
  double get_expansion_padding() const { return m_expansionSum; }

  void set_search_method(const stk::search::SearchMethod searchMethod) { m_searchMethod = searchMethod; }
  stk::search::SearchMethod get_search_method() const { return m_searchMethod; }

 protected:
  double m_expansionFactor{0.5};
  double m_expansionSum{0.0};
  double m_expansionLimit{std::numeric_limits<double>::max()};
  bool m_useNearestNodeForClosestBoundingBox{false};
  bool m_useCentroidForGeometricProximity{false};
  bool m_printSearchWarnings{true};
  bool m_doInitialSearchExpansion{true};
  stk::ParallelMachine m_comm;
  EntityProcRelationVec m_globalRangeToDomain;

  std::vector<BoundingBoxB> m_unpairedRecvEntities;

  bool m_initialized{false};

  std::vector< stk::search::Inspector<SENDMESH, RECVMESH> > m_inspectors;
  std::shared_ptr< stk::search::InspectorOutputFormatterBase<SENDMESH, RECVMESH> > m_outputFormatter;

  std::shared_ptr<SENDMESH> m_mesha;
  std::shared_ptr<RECVMESH> m_meshb;
  stk::search::SearchMethod m_searchMethod{stk::search::KDTREE};
  bool m_isListenerActive{true};
  size_t m_syncCount{0};
  const std::string m_name;

  std::ostream* m_outputStream{&stk::outputP0()};

  stk::search::FilterCoarseSearchResultVector<RECVMESH> m_filteredSearchResult;

  double m_fractionalLimitForObjectsOutsideDomain{0.0};

  void coarse_search_for_objects_outside_domain_by_default_filter(bool includeGhosts);
  void coarse_search_for_objects_outside_domain_by_nearest_node(bool includeGhosts);

  void delete_range_points_found(std::vector<BoundingBoxB>& range_vector, const EntityProcRelationVec& del);
  void coarse_search_by_range();

  void print_expansion_warnings(int number_coarse_search_passes, size_t domain_size);
};

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search()
{
   m_globalRangeToDomain.clear();

   coarse_search_by_range();
   ghost_from_elements();
   update_ghosted_keys();
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::local_search()
{
  unsigned numUnpairedRecvEntities = stk::get_global_sum(m_comm, m_unpairedRecvEntities.size());
  if (numUnpairedRecvEntities > 0) {
    const bool includeGhosts = true;
    coarse_search_for_objects_outside_domain(includeGhosts);
    numUnpairedRecvEntities = stk::get_global_sum(m_comm, m_unpairedRecvEntities.size());
  }

  stk::search::FilterCoarseSearchOptions filterOptions(*m_outputStream, m_mesha->get_extrapolate_option());
  filterOptions.m_useNearestNodeForClosestBoundingBox = m_useNearestNodeForClosestBoundingBox;
  filterOptions.m_useCentroidForGeometricProximity = m_useCentroidForGeometricProximity;
  filterOptions.m_verbose = m_printSearchWarnings;

  m_filteredSearchResult.clear();
  stk::search::filter_coarse_search(m_name, m_globalRangeToDomain, *m_mesha, *m_meshb, filterOptions, m_filteredSearchResult);
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::initialize()
{
   m_mesha->initialize();
   m_meshb->initialize();
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::reinitialize()
{

}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search_for_objects_outside_domain_by_default_filter(bool includeGhosts)
{
  std::vector<BoundingBoxA> domain_vector;

  m_mesha->bounding_boxes(domain_vector, includeGhosts);

  std::vector<double> parametricCoords;
  std::vector<double> coords;

  for(const BoundingBoxB& boxB : m_unpairedRecvEntities) {
    const EntityProcB& entityProcB = boxB.second;
    EntityKeyB const& entityKeyB = entityProcB.id();
    m_meshb->coordinates(entityKeyB, coords);

    double bestX = std::numeric_limits<double>::max();
    EntityProcA bestEntityProcA;

    for(const BoundingBoxA& boxA : domain_vector) {
      const EntityProcA& entityProcA = boxA.second;
      EntityKeyA const& entityKeyA = entityProcA.id();

      double geometricDistanceSquared = std::numeric_limits<double>::max();
      double parametricDistance = std::numeric_limits<double>::max();
      bool isWithinParametricTolerance = false;
      bool isWithinGeometricTolerance = false;

      m_mesha->find_parametric_coords(entityKeyA, coords, parametricCoords, parametricDistance, isWithinParametricTolerance);

      if(parametricDistance == std::numeric_limits<double>::max()) {
        continue;
      }

      bool modified = m_mesha->modify_search_outside_parametric_tolerance(entityKeyA, coords, parametricCoords,
                                                                          geometricDistanceSquared, isWithinGeometricTolerance);

      double distance = modified ? geometricDistanceSquared : parametricDistance;
      if(stk::search::less_than_with_centroid_fallback<SENDMESH,RECVMESH>(distance, bestX, bestEntityProcA.id(),
                                                                          entityKeyA, entityKeyB, m_mesha.get(), m_meshb.get()))
      {
        bestX = distance;
        bestEntityProcA = entityProcA;
      }
    }

    if(std::numeric_limits<double>::max() > bestX) {
      m_globalRangeToDomain.emplace_back(entityProcB, bestEntityProcA);
    }
  }

  delete_range_points_found(m_unpairedRecvEntities, m_globalRangeToDomain);

  unsigned numUnpairedRecvEntities = stk::get_global_sum(m_comm, m_unpairedRecvEntities.size());

  if(numUnpairedRecvEntities > 0 && m_printSearchWarnings) {
    std::ostringstream os;
    os << m_name << ": Coarse Search: There are " << numUnpairedRecvEntities
       << " points (mesh object centroids) in the receive region '" << m_meshb->name()
       << "' which do not lie inside a send mesh object and are outside by greater "
       << "than the parametric tolerance of " << m_meshb->get_parametric_tolerance() << "."
       << "  These points will not be handled.";

    (*m_outputStream) << os.str() << std::endl;
  }
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search_for_objects_outside_domain_by_nearest_node(bool includeGhosts)
{
  if (stk::is_true_on_all_procs(m_comm, m_unpairedRecvEntities.empty())) {
    return;
  }

  std::vector<BoundingBoxA> domain_vector;
  std::vector<double> coords;

  m_mesha->bounding_boxes(domain_vector, includeGhosts);

  std::vector< std::pair<size_t, double> > sortedDistances;
  sortedDistances.resize(domain_vector.size());

  for(const BoundingBoxB& boxB : m_unpairedRecvEntities) {
    const EntityProcB& entityProcB = boxB.second;
    EntityKeyB const& entityKeyB = entityProcB.id();
    m_meshb->coordinates(entityKeyB, coords);

    for(size_t i = 0; i < domain_vector.size(); ++i) {
      const EntityProcA& entityProcA = domain_vector[i].second;
      EntityKeyA const& entityKeyA = entityProcA.id();

      sortedDistances[i].first = i;
      sortedDistances[i].second = m_mesha->get_distance_from_nearest_node(entityKeyA, coords);
    }

    std::sort(sortedDistances.begin(), sortedDistances.end(),
              [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) { return a.second < b.second; });

    int64_t cutOffIndex = sortedDistances.size() > 0 ? 0 : -1;
    double cutOffDistance = sortedDistances.size() > 0 ? sortedDistances[0].second * (1.0 + m_fractionalLimitForObjectsOutsideDomain)
                                                       : std::numeric_limits<double>::max();

    for(size_t i = 0; i < sortedDistances.size(); ++i) {
      size_t index = sortedDistances[i].first;
      const EntityProcA& entityProcA = domain_vector[index].second;
      EntityKeyA const& entityKeyA = entityProcA.id();

      const EntityProcA& bestEntityProcA = domain_vector[sortedDistances[0].first].second;

      if(stk::search::less_than_with_centroid_fallback<SENDMESH,RECVMESH>(cutOffDistance, sortedDistances[i].second, bestEntityProcA.id(),
                                                                          entityKeyA, entityKeyB, m_mesha.get(), m_meshb.get())) {
        break;
      }
      cutOffIndex = i;
    }

    for(int64_t i = 0; i <= cutOffIndex; ++i) {
      size_t index = sortedDistances[i].first;
      const EntityProcA& entityProcA = domain_vector[index].second;

      m_globalRangeToDomain.emplace_back(entityProcB, entityProcA);
    }
  }

  delete_range_points_found(m_unpairedRecvEntities, m_globalRangeToDomain);

  unsigned numUnpairedRecvEntities = stk::get_global_sum(m_comm, m_unpairedRecvEntities.size());

  if(numUnpairedRecvEntities > 0 && m_printSearchWarnings) {
    std::ostringstream os;
    os << m_name << ": Coarse Search: There are " << numUnpairedRecvEntities
       << " points (mesh object centroids) in the receive region '" << m_meshb->name()
       << "' which do not lie inside a send mesh object and are outside by greater "
       << "than the parametric tolerance of " << m_meshb->get_parametric_tolerance() << "."
       << "  These points will not be handled.";

    (*m_outputStream) << os.str() << std::endl;
  }
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search_for_objects_outside_domain(bool includeGhosts)
{
  stk::search::ObjectOutsideDomainPolicy extrapolatePolicy = m_mesha->get_extrapolate_option();

  if(extrapolatePolicy != stk::search::ObjectOutsideDomainPolicy::IGNORE) {
    if(m_fractionalLimitForObjectsOutsideDomain > 0.0) {
      coarse_search_for_objects_outside_domain_by_nearest_node(includeGhosts);
    } else {
      coarse_search_for_objects_outside_domain_by_default_filter(includeGhosts);
    }
  }
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::delete_range_points_found(std::vector<BoundingBoxB>& range_vector, const EntityProcRelationVec& del)
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
        difference.begin(), impl::compare<RECVMESH>());
    difference.resize(it-difference.begin());
  }
  swap(difference, range_vector);
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::print_expansion_warnings(int number_coarse_search_passes, size_t range_size)
{
  if(!m_printSearchWarnings) return;

  size_t g_range_vector_size = stk::get_global_max(m_comm, range_size);
  (*m_outputStream) << m_name
                    << ": GeometricSearch::coarse_search(): Number of points not found: "
                    << g_range_vector_size << " after expanding bounding boxes: " << number_coarse_search_passes
                    << " time(s)" << std::endl;
  (*m_outputStream)
      << "...will now expand the set of candidate bounding boxes and re-attempt the coarse search" << std::endl;
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search_for_objects_inside_domain()
{
  std::vector<BoundingBoxA> domain_vector;
  std::vector<BoundingBoxB> range_vector;

  if(m_mesha) m_mesha->bounding_boxes(domain_vector);

  if (stk::is_true_on_all_procs(m_comm, domain_vector.empty())) {
    return;
  }

  if(m_meshb) m_meshb->bounding_boxes(range_vector);

  if(!local_is_sorted(domain_vector.begin(), domain_vector.end(), impl::BoundingBoxCompare<BoundingBoxA>()))
    std::sort(domain_vector.begin(), domain_vector.end(), impl::BoundingBoxCompare<BoundingBoxA>());
  if(!local_is_sorted(range_vector.begin(), range_vector.end(), impl::BoundingBoxCompare<BoundingBoxB>()))
    std::sort(range_vector.begin(), range_vector.end(), impl::BoundingBoxCompare<BoundingBoxB>());

  // keep track of how many times the coarse search needs to exercise the expansion_factor
  int num_coarse_search_passes = 0;

  bool range_vector_empty = stk::is_true_on_all_procs(m_comm, range_vector.empty());

  if(m_doInitialSearchExpansion) {
    // THIS IS DANGEROUS. The expansion factor is (by default) used for an initial box expansion/padding,
    // but will also trigger multiple coarse searches if it is ever changed to be > 1. If you never
    // want it to do multiple coarse searches here, make sure it is <= 1.0. This initial expansion is
    // a legacy of Framework Srch to ensure consistency for Goodyear
    for(BoundingBoxA& i : domain_vector) {
      impl::inflate_bounding_box(i.first, m_expansionFactor, m_expansionSum);
    }
  }

  while(!range_vector_empty) { // Keep going until all range points are processed.
    // Slightly confusing: coarse_search documentation has domain->range
    // relations sorted by domain key.  We want range->domain type relations
    // sorted on range key. It might appear we have the arguments revered
    // in coarse_search call, but really, this is what we want.
    EntityProcRelationVec rng_to_dom;
    EntityProcRelationVec& rng_to_dom_vec = num_coarse_search_passes == 0 ? m_globalRangeToDomain : rng_to_dom;
    stk::search::coarse_search(range_vector, domain_vector, m_searchMethod, m_comm, rng_to_dom_vec);

    if(num_coarse_search_passes > 0) {
      m_globalRangeToDomain.insert(m_globalRangeToDomain.end(), rng_to_dom_vec.begin(), rng_to_dom_vec.end());
    }

    // increment how many times we are within the while loop
    num_coarse_search_passes++;

    delete_range_points_found(range_vector, rng_to_dom_vec);

    range_vector_empty = stk::is_true_on_all_procs(m_comm, range_vector.empty());

    const bool terminate_on_first_pass = (range_vector_empty || m_expansionFactor <= 1.0) && num_coarse_search_passes == 1;

    if(terminate_on_first_pass) {
      if(!range_vector_empty && m_printSearchWarnings) {
        size_t range_vector_size = range_vector.size();
        size_t g_range_vector_size = stk::get_global_max(m_comm, range_vector_size);
        (*m_outputStream) << m_name
                          << ": GeometricSearch::coarse_search(): Number of points not found: "
                          << g_range_vector_size << " in initial coarse search" << std::endl;
      }

      break;
    }

    if(!range_vector_empty) {
      std::string local_err_msg;

      for(BoundingBoxB& box : range_vector) {
        const auto box_scale = stk::search::length_scale(box.first);

        if(box_scale > m_expansionLimit)
        {
          std::ostringstream err;
          err << "Bounding box size: " << box_scale << " exceeds the expansion tolerance: "
              << m_expansionLimit << ".  Box info: " << box.first << std::endl;
          local_err_msg +=  err.str();
        }

        // If points were missed, increase search radius.
        impl::inflate_bounding_box(box.first, m_expansionFactor, m_expansionSum);
      }
      for(BoundingBoxA& box : domain_vector) {
        impl::inflate_bounding_box(box.first, m_expansionFactor, m_expansionSum);
      }

      impl::throw_error_on_coarse_expansion_limit(m_comm, local_err_msg);

      // sum and provide message to user
      print_expansion_warnings(num_coarse_search_passes, range_vector.size());
    }
  }

  m_unpairedRecvEntities.swap(range_vector);
  std::sort(m_globalRangeToDomain.begin(), m_globalRangeToDomain.end());
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::coarse_search_by_range()
{
  m_unpairedRecvEntities.clear();
  coarse_search_for_objects_inside_domain();
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::update_ghosted_keys()
{
  const unsigned my_rank = stk::parallel_machine_rank(m_comm);

  for (auto& i : m_globalRangeToDomain) {
    EntityProcA& sendEntityProc = i.second;
    EntityProcB& recvEntityProc = i.first;

    const unsigned send_owning_rank = sendEntityProc.proc();
    const unsigned recv_owning_rank = recvEntityProc.proc();

    if (send_owning_rank != my_rank && recv_owning_rank == my_rank) {
      EntityKeyA& sendKey = sendEntityProc.id();
      m_mesha->update_ghosted_key(sendKey);
    }
  }

}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::determine_entities_to_copy(EntityProcVecA& entities_to_copy) const {
  entities_to_copy.clear();

  const unsigned my_rank = stk::parallel_machine_rank(m_comm);

  for (auto& i : m_globalRangeToDomain) {
    const unsigned            domain_owning_rank = i.second.proc();
    const unsigned             range_owning_rank = i.first.proc();
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const EntityKeyA entity = i.second.id();
      const EntityProcA ep(entity, range_owning_rank);
      entities_to_copy.push_back(ep);
    }
  }
  std::sort(entities_to_copy.begin(), entities_to_copy.end());
  typename EntityProcVecA::iterator del = std::unique(entities_to_copy.begin(), entities_to_copy.end());
  entities_to_copy.resize(std::distance(entities_to_copy.begin(), del));
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::ghost_from_elements()
{
  EntityProcVecA entity_keys;
  determine_entities_to_copy(entity_keys);

  m_mesha->update_ghosting(entity_keys);

  post_mesh_modification_event(false);
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::post_mesh_modification_event(bool activeFlag)
{
  m_isListenerActive = activeFlag;
  m_mesha->post_mesh_modification_event();
  m_isListenerActive = true;
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::register_inspector(const std::string& fileName, const EntityKeyVecB& rangeEntities)
{
  for(stk::search::Inspector<SENDMESH,RECVMESH>& gadget : m_inspectors) {
    STK_ThrowRequireMsg(fileName != gadget.get_file_name(), "Attempt to register inspector with previously used file name: " << fileName);
  }

  m_inspectors.emplace_back(m_mesha, m_meshb, fileName, rangeEntities);

  if(m_outputFormatter != nullptr) {
    m_inspectors.back().register_output_formatter(m_outputFormatter);
  }
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::inspect_user_defined_entities()
{
  for(stk::search::Inspector<SENDMESH,RECVMESH>& gadget : m_inspectors) {
    gadget.execute(m_globalRangeToDomain);
  }
}

template <typename SENDMESH, typename RECVMESH>
void GeometricSearch<SENDMESH,RECVMESH>::register_output_formatter(std::shared_ptr< stk::search::InspectorOutputFormatterBase<SENDMESH, RECVMESH> > outputFormatter)
{
  m_outputFormatter = outputFormatter;

  for(stk::search::Inspector<SENDMESH,RECVMESH>& gadget : m_inspectors) {
    gadget.register_output_formatter(outputFormatter);
  }
}

} // namespace spmd
} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCH_HPP_ */
