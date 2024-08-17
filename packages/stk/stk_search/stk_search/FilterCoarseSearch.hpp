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

#ifndef STK_SEARCH_FILTER_TO_NEAREST_HPP
#define STK_SEARCH_FILTER_TO_NEAREST_HPP

#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <sstream>
#include <iostream>
#include <ostream>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>

#include "stk_search/DistanceComparison.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk {
namespace search {

template <class SENDMESH, class RECVMESH>
using FilterCoarseSearchProcRelation = std::pair<typename RECVMESH::EntityProc, typename SENDMESH::EntityProc>;

template <class SENDMESH, class RECVMESH>
using FilterCoarseSearchProcRelationVec = std::vector<FilterCoarseSearchProcRelation<SENDMESH, RECVMESH>>;

//BEGINObjectOutsideDomainPolicy 
enum class ObjectOutsideDomainPolicy { IGNORE, EXTRAPOLATE, TRUNCATE, PROJECT, ABORT, UNDEFINED_OBJFLAG = 0xff };
//ENDObjectOutsideDomainPolicy 

inline ObjectOutsideDomainPolicy get_object_outside_domain_policy(const std::string& id)
{
  if(id == "IGNORE") return ObjectOutsideDomainPolicy::IGNORE;
  if(id == "EXTRAPOLATE") return ObjectOutsideDomainPolicy::EXTRAPOLATE;
  if(id == "TRUNCATE") return ObjectOutsideDomainPolicy::TRUNCATE;
  if(id == "PROJECT") return ObjectOutsideDomainPolicy::PROJECT;
  if(id == "ABORT") return ObjectOutsideDomainPolicy::ABORT;

  return ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG;
}

inline std::string get_object_outside_domain_policy(const ObjectOutsideDomainPolicy id)
{
  switch(id) {
  case ObjectOutsideDomainPolicy::IGNORE:
    return "IGNORE";
  case ObjectOutsideDomainPolicy::EXTRAPOLATE:
    return "EXTRAPOLATE";
  case ObjectOutsideDomainPolicy::TRUNCATE:
    return "TRUNCATE";
  case ObjectOutsideDomainPolicy::PROJECT:
    return "PROJECT";
  case ObjectOutsideDomainPolicy::ABORT:
    return "ABORT";
  case ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG:
    return "UNDEFINED";
  }

  return std::string("");
}

//BEGINFilterCoarseSearchOptions
struct FilterCoarseSearchOptions
{
  std::ostream& m_outputStream{std::cout};
  ObjectOutsideDomainPolicy m_extrapolatePolicy{ObjectOutsideDomainPolicy::EXTRAPOLATE};
  bool m_useNearestNodeForClosestBoundingBox{false};
  bool m_useCentroidForGeometricProximity{false};
  bool m_verbose{true};
//ENDFilterCoarseSearchOptions

  FilterCoarseSearchOptions(std::ostream& out)
  : m_outputStream(out) {}
  FilterCoarseSearchOptions(std::ostream& out, ObjectOutsideDomainPolicy extrapolatePolicy)
  : m_outputStream(out), m_extrapolatePolicy(extrapolatePolicy) {}
  FilterCoarseSearchOptions(std::ostream& out,
                         ObjectOutsideDomainPolicy extrapolatePolicy,
                         bool useNearestNodeForClosestBoundingBox,
                         bool useCentroidForGeometricProximity,
                         bool verbose)
  : m_outputStream(out)
  , m_extrapolatePolicy(extrapolatePolicy)
  , m_useNearestNodeForClosestBoundingBox(useNearestNodeForClosestBoundingBox)
  , m_useCentroidForGeometricProximity(useCentroidForGeometricProximity)
  , m_verbose(verbose) {}
};

//BEGINFilterCoarseSearchResult
template <class RECVMESH>
class FilterCoarseSearchResult
{
public:
  using EntityKey = typename RECVMESH::EntityKey;

  virtual void add_search_filter_info(const EntityKey key,
                                      const std::vector<double>&paramCoords,
                                      const double parametricDistance,
                                      const bool isWithinParametricTolerance,
                                      const double geometricDistanceSquared,
                                      const bool isWithinGeometricTolerance) = 0;

  virtual void get_parametric_coordinates(const EntityKey key, std::vector<double>& paramCoords) const = 0;
  virtual void clear() = 0;

  virtual ~FilterCoarseSearchResult() {}
};
//ENDFilterCoarseSearchResult

//BEGINFilterCoarseSearchResultMap
template <class RECVMESH>
class FilterCoarseSearchResultMap : public FilterCoarseSearchResult<RECVMESH>
//ENDFilterCoarseSearchResultMap
{
public:
  using EntityKey = typename RECVMESH::EntityKey;

  void add_search_filter_info(const EntityKey key,
                              const std::vector<double>& paramCoords,
                              const double parametricDistance,
                              const bool isWithinParametricTolerance,
                              const double geometricDistanceSquared,
                              const bool isWithinGeometricTolerance) override
  {
    m_searchFilterInfo[key] = paramCoords;
  }

  void get_parametric_coordinates(const EntityKey key, std::vector<double>& paramCoords) const override
  {
    const std::vector<double>& returnVal = m_searchFilterInfo.at(key);
    paramCoords = returnVal;
  }

  void clear() override { m_searchFilterInfo.clear(); }
  virtual ~FilterCoarseSearchResultMap() {}

private:
  std::map<EntityKey, std::vector<double>> m_searchFilterInfo;
};

//BEGINFilterCoarseSearchResultVector
template <class RECVMESH>
class FilterCoarseSearchResultVector : public FilterCoarseSearchResult<RECVMESH>
//ENDFilterCoarseSearchResultVector
{
public:
  using EntityKey = typename RECVMESH::EntityKey;
  using VectorType = std::pair<EntityKey, std::vector<double> >;

  struct VectorTypeLess
  {
    inline bool operator()(const VectorType& a, const VectorType& b) const { return a.first < b.first; }
    inline bool operator()(const VectorType& a, const EntityKey&  b) const { return a.first < b;       }
    inline bool operator()(const EntityKey&  a, const VectorType& b) const { return a       < b.first; }
    inline bool operator()(const EntityKey&  a, const EntityKey&  b) const { return a       < b;       }
  };

  FilterCoarseSearchResultVector() : m_isSorted(false) {}

  void add_search_filter_info(const EntityKey key,
                              const std::vector<double>& paramCoords,
                              const double parametricDistance,
                              const bool isWithinParametricTolerance,
                              const double geometricDistanceSquared,
                              const bool isWithinGeometricTolerance) override
  {
    m_searchFilterInfo.push_back(std::make_pair(key, paramCoords));
    m_isSorted = false;
  }

  void get_parametric_coordinates(const EntityKey key, std::vector<double>& paramCoords) const override
  {
    if(!m_isSorted) {
      stk::util::sort_and_unique(m_searchFilterInfo, VectorTypeLess());
      m_isSorted = true;
    }

    auto result = std::lower_bound(m_searchFilterInfo.begin(), m_searchFilterInfo.end(), key, VectorTypeLess());
    STK_ThrowRequireMsg((result != m_searchFilterInfo.end()) && ((*result).first == key),
                        "Could not find key in search filter info");

    paramCoords = (*result).second;
  }

  void clear() override { m_searchFilterInfo.clear(); }
  virtual ~FilterCoarseSearchResultVector() {}

private:
  mutable bool m_isSorted{false};
  mutable std::vector<VectorType> m_searchFilterInfo;
};

namespace impl {

struct FilterCoarseSearchStats
{
  unsigned numEntitiesWithinTolerance{0};
  unsigned numEntitiesOutsideTolerance{0};

  double longestParametricExtrapolation{0.0};
  double longestGeometricExtrapolation{0.0};
};

template <class MESH>
inline
double get_distance_squared_from_centroid(MESH& mesh, const typename MESH::EntityKey k, const double* toCoords)
{
  std::vector<double> centroid = mesh.centroid(k);
  return stk::search::distance_sq(centroid.size(), centroid.data(), toCoords);
}

template <class MESH>
inline
double get_distance_from_centroid(MESH& mesh, const typename MESH::EntityKey k, const double* toCoords)
{
  double distanceSquared = get_distance_squared_from_centroid<MESH>(mesh, k, toCoords);
  return std::sqrt(distanceSquared);
}

template <typename EntityProcRelationVec>
void remove_non_local_range_entities(EntityProcRelationVec& rangeToDomain, int localProc)
{
  using ValueType = typename EntityProcRelationVec::value_type;
  rangeToDomain.erase(
    std::remove_if(rangeToDomain.begin(), rangeToDomain.end(),
                   [&](const ValueType& x) { return (static_cast<int>(x.first.proc()) != localProc); }),
    rangeToDomain.end());
}

template <typename EntityProcRelationVec>
typename EntityProcRelationVec::const_iterator find_next_key(typename EntityProcRelationVec::const_iterator current,
                                                             const EntityProcRelationVec& rangeToDomain)
{
  typename EntityProcRelationVec::const_iterator next = current;
  while(next != rangeToDomain.end() && next->first.id() == current->first.id()) {
    ++next;
  }
  return next;
}

template <class SENDMESH, class RECVMESH>
struct FilterResult {
  FilterResult()
    : parametricDistance(std::numeric_limits<double>::max())
    , geometricDistanceSquared(std::numeric_limits<double>::max())
    , isWithinGeometricTolerance(false)
    , isWithinParametricTolerance(false)
    , isAccepted(false)
  {
  }

  double parametricDistance;
  double geometricDistanceSquared;
  std::vector<double> parametricCoords;
  bool isWithinGeometricTolerance;
  bool isWithinParametricTolerance;
  bool isAccepted;
  typename FilterCoarseSearchProcRelationVec<SENDMESH, RECVMESH>::const_iterator nearest;
};

template <class SENDMESH, class RECVMESH>
void accept_candidate(std::vector<double>& parametricCoords,
                      const double parametricDistance,       const bool isWithinParametricTolerance,
                      const double geometricDistanceSquared, const bool isWithinGeometricTolerance,
                      const typename FilterCoarseSearchProcRelationVec<SENDMESH, RECVMESH>::const_iterator& ii,
                      FilterResult<SENDMESH, RECVMESH>& bestCandidate)
{
  bestCandidate.isWithinParametricTolerance = isWithinParametricTolerance;
  bestCandidate.isWithinGeometricTolerance = isWithinGeometricTolerance;
  bestCandidate.isAccepted = true;
  bestCandidate.parametricDistance = parametricDistance;
  bestCandidate.nearest = ii;
  bestCandidate.parametricCoords.swap(parametricCoords);
  bestCandidate.geometricDistanceSquared = geometricDistanceSquared;
}

template <class SENDMESH>
void set_geometric_info(SENDMESH& sendMesh,
                        const typename SENDMESH::EntityKey sendEntity, const double* tocoords,
                        const double searchToleranceSquared, const bool useCentroidForGeometricProximity,
                        double& geometricDistanceSquared, bool& isWithinGeometricTolerance)
{
  if(!isWithinGeometricTolerance) {

    if(useCentroidForGeometricProximity) {
      geometricDistanceSquared = sendMesh.get_distance_squared_from_centroid(sendEntity, tocoords);
    }
    else {
      geometricDistanceSquared = sendMesh.get_closest_geometric_distance_squared(sendEntity, tocoords);
    }

    isWithinGeometricTolerance = geometricDistanceSquared <= searchToleranceSquared;
  }
}

template <class SENDMESH, class RECVMESH>
FilterCoarseSearchStats filter_coarse_search_by_range(FilterCoarseSearchProcRelationVec<SENDMESH, RECVMESH>& rangeToDomain,
                                                      SENDMESH& sendMesh, RECVMESH& recvMesh,
                                                      const bool useNearestNodeForClosestBoundingBox,
                                                      const bool useCentroidForGeometricProximity,
                                                      const ObjectOutsideDomainPolicy extrapolatePolicy,
                                                      FilterCoarseSearchResult<RECVMESH>& filterResult)
{
  using const_iterator = typename FilterCoarseSearchProcRelationVec<SENDMESH,RECVMESH>::const_iterator;
  using iterator = typename FilterCoarseSearchProcRelationVec<SENDMESH,RECVMESH>::iterator;

  FilterCoarseSearchStats stats;

  int thisProc = stk::parallel_machine_rank(recvMesh.comm());

  remove_non_local_range_entities(rangeToDomain, thisProc);

  const_iterator current_key = rangeToDomain.begin();
  const_iterator next_key = find_next_key(current_key, rangeToDomain);
  iterator keep = rangeToDomain.begin();
  double searchTolerance = recvMesh.get_search_tolerance();
  double searchToleranceSquared = searchTolerance * searchTolerance;

  std::vector<double> parametricCoords;

  while(current_key != rangeToDomain.end()) {
    FilterResult<SENDMESH, RECVMESH> bestCandidate;

    double geometricDistanceSquared = std::numeric_limits<double>::max();
    double parametricDistance= std::numeric_limits<double>::max();
    bool isWithinParametricTolerance = false;
    bool isWithinGeometricTolerance = false;

    const typename RECVMESH::EntityKey recvEntity = current_key->first.id();

    const double* tocoords = recvMesh.coord(recvEntity);

    std::pair<const_iterator, const_iterator> keys = std::make_pair(current_key, next_key);
    bestCandidate.nearest = keys.second;

    for(const_iterator ii = keys.first; ii != keys.second; ++ii) {
      const typename SENDMESH::EntityKey sendEntity = ii->second.id();
      sendMesh.find_parametric_coords(sendEntity, tocoords, parametricCoords, parametricDistance, isWithinParametricTolerance);

      if(parametricDistance == std::numeric_limits<double>::max()) {
        continue;
      }

      sendMesh.modify_search_outside_parametric_tolerance(sendEntity, tocoords, parametricCoords,
                                                          geometricDistanceSquared, isWithinGeometricTolerance);

      if(useNearestNodeForClosestBoundingBox) {
        const double* fromcoords = sendMesh.coord(sendEntity);
        double distance = recvMesh.get_distance_from_nearest_node(recvEntity, fromcoords);
        geometricDistanceSquared = std::pow(distance, 2);
        isWithinParametricTolerance = false;
        isWithinGeometricTolerance = distance <= searchTolerance;
      }

      bool accept = false;

      if(isWithinParametricTolerance || bestCandidate.isWithinParametricTolerance) {
        accept = parametricDistance < bestCandidate.parametricDistance;
      }
      else {
        set_geometric_info<SENDMESH>(sendMesh, sendEntity, tocoords, searchToleranceSquared,
                                     useCentroidForGeometricProximity,
                                     geometricDistanceSquared,
                                     isWithinGeometricTolerance);

        accept = geometricDistanceSquared < bestCandidate.geometricDistanceSquared;
      }

      if(accept) {
        accept_candidate<SENDMESH, RECVMESH>(parametricCoords,
                                             parametricDistance,       isWithinParametricTolerance,
                                             geometricDistanceSquared, isWithinGeometricTolerance,
                                             ii, bestCandidate);
      }
    }

    bool ignoredNearest = false;

    if(bestCandidate.isAccepted) {
      if(bestCandidate.isWithinParametricTolerance || bestCandidate.isWithinGeometricTolerance) {
        stats.numEntitiesWithinTolerance++;
        filterResult.add_search_filter_info(recvEntity,
                                            bestCandidate.parametricCoords,
                                            bestCandidate.parametricDistance-1.0,
                                            bestCandidate.isWithinParametricTolerance,
                                            bestCandidate.geometricDistanceSquared,
                                            bestCandidate.isWithinGeometricTolerance);
      }
      else {
        stats.longestGeometricExtrapolation =
            std::max(stats.longestGeometricExtrapolation, bestCandidate.geometricDistanceSquared);
        stats.longestParametricExtrapolation =
            std::max(stats.longestParametricExtrapolation, (bestCandidate.parametricDistance - 1.0));

        if(extrapolatePolicy != ObjectOutsideDomainPolicy::IGNORE) {
          filterResult.add_search_filter_info(recvEntity,
                                              bestCandidate.parametricCoords,
                                              bestCandidate.parametricDistance-1.0,
                                              bestCandidate.isWithinParametricTolerance,
                                              bestCandidate.geometricDistanceSquared,
                                              bestCandidate.isWithinGeometricTolerance);
        }
        else {
          ignoredNearest = true;
        }
        stats.numEntitiesOutsideTolerance++;
      }
    } else {
      stats.numEntitiesOutsideTolerance++;
    }

    if(bestCandidate.isAccepted && !ignoredNearest) {
      *keep = *bestCandidate.nearest;
      ++keep;
    }
    current_key = next_key;
    if (current_key != rangeToDomain.end()) {
      next_key = find_next_key(current_key, rangeToDomain);
    }
  }
  rangeToDomain.resize(std::distance(rangeToDomain.begin(),keep));

  stats.longestGeometricExtrapolation = std::sqrt(stats.longestGeometricExtrapolation);

  return stats;
}

struct OutputLogistic {
  const std::string name;
  unsigned numEntitiesWithinTolerance;
  unsigned numEntitiesOutsideTolerance;
  double longestParametricExtrapolation;
  double longestGeometricExtrapolation;
  const double parametricTol;
  const double geometricTol;
};

template <typename SENDMESH, typename RECVMESH>
void output_summary_outside_tolerance(SENDMESH& sendMesh, RECVMESH& recvMesh, OutputLogistic stat,
                                      const ObjectOutsideDomainPolicy& extrapolatePolicy,
                                      std::ostream& outputStream)
{
  double g_longestExtrapolation[2];
  double l_longestExtrapolation[2] = { stat.longestParametricExtrapolation, stat.longestGeometricExtrapolation };
  stk::all_reduce_max(sendMesh.comm(), l_longestExtrapolation, g_longestExtrapolation, 2);

  std::ostringstream os;
  os << stat.name << ": Fine Search: There are " << stat.numEntitiesWithinTolerance
     << " points (mesh object points) in the receive mesh '" << recvMesh.name() << "' which "
     << " are within the parametric space tolerance of " << stat.parametricTol
     << " and will be treated as inside the receiving mesh."
     << "  Additionally there are " << stat.numEntitiesOutsideTolerance
     << " points (mesh object points) in the receive mesh '" << recvMesh.name()
     << "' which do not lie inside a send mesh object and are outside by greater "
     << "than the parametric tolerance of " << stat.parametricTol << " and/or the geometric tolerance of "
     << stat.geometricTol << "."
     << "  These points will be handled via: " << get_object_outside_domain_policy(extrapolatePolicy)
     << ". The longest parametric extrapolation distance for these points is " << g_longestExtrapolation[0]
     << " (parametric units)"
     << " and the longest geometric extrapolation distance for these points is " << g_longestExtrapolation[1]
     << " (mesh distance units)";

  outputStream << os.str() << std::endl;
}
}

//BEGINfilter_coarse_search_impl
template <class SENDMESH, class RECVMESH>
void filter_coarse_search(const std::string& name,
                          FilterCoarseSearchProcRelationVec<SENDMESH, RECVMESH>& rangeToDomain,
                          SENDMESH& sendMesh, RECVMESH& recvMesh,
                          FilterCoarseSearchOptions& filterOptions,
                          FilterCoarseSearchResult<RECVMESH>& filterResult)
//ENDfilter_coarse_search_impl
{
  const double parametricTol = recvMesh.get_parametric_tolerance();
  const double geometricTol = recvMesh.get_search_tolerance();
  impl::FilterCoarseSearchStats stats = impl::filter_coarse_search_by_range(rangeToDomain, sendMesh, recvMesh,
                                                                            filterOptions.m_useNearestNodeForClosestBoundingBox,
                                                                            filterOptions.m_useCentroidForGeometricProximity,
                                                                            filterOptions.m_extrapolatePolicy,
                                                                            filterResult);

  stats.numEntitiesWithinTolerance = stk::get_global_sum(sendMesh.comm(), stats.numEntitiesWithinTolerance);
  stats.numEntitiesOutsideTolerance = stk::get_global_sum(sendMesh.comm(), stats.numEntitiesOutsideTolerance);

  if(stats.numEntitiesOutsideTolerance) {

    if(filterOptions.m_verbose) {
      impl::OutputLogistic outputSummary{ name,
                                          stats.numEntitiesWithinTolerance,
                                          stats.numEntitiesOutsideTolerance,
                                          stats.longestParametricExtrapolation,
                                          stats.longestGeometricExtrapolation,
                                          parametricTol,
                                          geometricTol };
      impl::output_summary_outside_tolerance(sendMesh, recvMesh, outputSummary,
                                             filterOptions.m_extrapolatePolicy, filterOptions.m_outputStream);
    }

    if(filterOptions.m_extrapolatePolicy == ObjectOutsideDomainPolicy::ABORT) {
      throw std::runtime_error("Aborting search due to user specified option");
    }
  }
}

}}
#endif
