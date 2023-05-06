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

#ifndef STK_STK_SEARCH_STK_SEARCH_FINESEARCH_HPP_
#define STK_STK_SEARCH_STK_SEARCH_FINESEARCH_HPP_

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
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace search {

enum class ObjectOutsideDomainPolicy { IGNORE, EXTRAPOLATE, TRUNCATE, PROJECT, ABORT, UNDEFINED_OBJFLAG = 0xff };

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

struct FilterToNearestOptions
{
  std::ostream& m_outputStream{std::cout};
  ObjectOutsideDomainPolicy m_extrapolatePolicy{ObjectOutsideDomainPolicy::EXTRAPOLATE};
  bool m_useNearestNodeForClosestBoundingBox{false};
  bool m_useCentroidForGeometricProximity{false};
  bool m_verbose{true};

  FilterToNearestOptions(std::ostream& out)
  : m_outputStream(out) {}
  FilterToNearestOptions(std::ostream& out, ObjectOutsideDomainPolicy extrapolatePolicy)
  : m_outputStream(out), m_extrapolatePolicy(extrapolatePolicy) {}
  FilterToNearestOptions(std::ostream& out,
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

template <class RECVMESH>
class FilterToNearestResult
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

  virtual ~FilterToNearestResult() {}
};

template <class RECVMESH>
class FilterToNearestResultMap : public FilterToNearestResult<RECVMESH>
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

  virtual ~FilterToNearestResultMap() {}

private:
  std::map<EntityKey, std::vector<double>> m_searchFilterInfo;
};

template <class RECVMESH>
class FilterToNearestResultVector : public FilterToNearestResult<RECVMESH>
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

  FilterToNearestResultVector() : m_isSorted(false) {}

  void add_search_filter_info(const EntityKey key,
                              const std::vector<double>& paramCoords,
                              const double parametricDistance,
                              const bool isWithinParametricTolerance,
                              const double geometricDistanceSquared,
                              const bool isWithinGeometricTolerance) override
  {
    m_searchFilterInfo.push_back(std::make_pair(key, paramCoords));
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

  virtual ~FilterToNearestResultVector() {}

private:
  mutable bool m_isSorted{false};
  mutable std::vector<VectorType> m_searchFilterInfo;
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

}}

#endif /* STK_STK_SEARCH_STK_SEARCH_FINESEARCH_HPP_ */
