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

#ifndef STK_STK_SEARCH_STK_SEARCH_SEARCHINTERFACE_HPP_
#define STK_STK_SEARCH_STK_SEARCH_SEARCHINTERFACE_HPP_

#include <memory>
#include <vector>
#include <utility>
#include <sstream>
#include <iostream>
#include <ostream>
#include <string>
#include <limits>
#include <cmath>
#include <string>
#include <algorithm>

#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...

namespace stk {
namespace search {

template <typename TOPOLOGY, typename ENTITYKEY>
class EvaluatePointsInterface {
 public:
  virtual size_t num_points(ENTITYKEY, TOPOLOGY) = 0;
  virtual const double* coords(ENTITYKEY, size_t) = 0;

  virtual ~EvaluatePointsInterface() {}
};

template <typename ENTITYKEY>
class ExternalPointHandlerInterface {
 public:
  virtual bool
  handle_point(const ENTITYKEY k, const double* toCoords,
               std::vector<double>& parametricCoords,
               double& geometricDistanceSquared,
               bool& isWithinGeometricTolerance) const = 0;

  virtual ~ExternalPointHandlerInterface() {}
};

template <typename ENTITYKEY>
class FindParametricCoordsInterface {
 public:
  virtual void
  find_parametric_coords(const ENTITYKEY k, const double* toCoords,
                         std::vector<double>& paramCoords, double& paramDistance,
                         bool& isWithinParametricTolerance) const = 0;
  virtual void
  evaluate_parametric_coords(const ENTITYKEY k, const std::vector<double>& paramCoords,
                             std::vector<double>& evalPoint) const = 0;
  virtual ~FindParametricCoordsInterface() {}
};

template <typename TOPOLOGY>
class ProvideMasterElementInterface {
 public:
  virtual void evaluate_field(const TOPOLOGY& topo,
                              const std::vector<double>& paramCoords, // (numParamCoords)
                              const unsigned numFieldComponents,
                              const std::vector<double>& fieldData, // (numFieldComponents x numNodes)
                              std::vector<double>& result) = 0; // (numFieldComponents)

  virtual void find_parametric_coordinates(const TOPOLOGY& topo,
                                           const unsigned numCoordComponents,
                                           const std::vector<double>& elemNodeCoords,  // (numCoordComponents x numNodes)
                                           const std::vector<double>& inputCoords, // (numCoordComponents)
                                           std::vector<double>& paramCoords,
                                           double& paramDistance) = 0;

  virtual void coordinate_center(const TOPOLOGY& topo, std::vector<double>&) = 0;
  virtual unsigned num_parametric_coordinates(const TOPOLOGY& topo) = 0;
  virtual unsigned num_integration_points(const TOPOLOGY& topo) = 0;
  virtual const double* integration_points(const TOPOLOGY& topo) = 0;
  virtual ~ProvideMasterElementInterface() {}
};

template <typename MESH>
struct MeshTraits;

template <typename SENDMESH>
class SourceMeshInterface
{
public:
  using BOUNDINGBOX = typename MeshTraits<SENDMESH>::BoundingBox;
  using ENTITYKEY = typename MeshTraits<SENDMESH>::EntityKey;

  SourceMeshInterface() = default;
  virtual ~SourceMeshInterface() = default;

  virtual stk::ParallelMachine comm() const = 0;

  virtual std::string name() const = 0;

  virtual void bounding_boxes(std::vector<BOUNDINGBOX>& boxes) const = 0;

  virtual void find_parametric_coords(const ENTITYKEY k, const double* toCoords,
                                      std::vector<double>& parametricCoords,
                                      double& parametricDistance,
                                      bool& isWithinParametricTolerance) const = 0;

  virtual bool modify_search_outside_parametric_tolerance(const ENTITYKEY k,
                                                          const double* toCoords,
                                                          std::vector<double>& parametricCoords,
                                                          double& geometricDistanceSquared,
                                                          bool& isWithinGeometricTolerance) const = 0;

  virtual double get_distance_from_nearest_node(const ENTITYKEY k, const double* point) const = 0;

  virtual double get_closest_geometric_distance_squared(const ENTITYKEY k, const double* toCoords) const = 0;

  virtual double get_distance_from_centroid(const ENTITYKEY k, const double* toCoords) const = 0;

  virtual double get_distance_squared_from_centroid(const ENTITYKEY k, const double* toCoords) const = 0;

  virtual void centroid(const ENTITYKEY k, std::vector<double>& centroidVec) const = 0;

  virtual const double* coord(const ENTITYKEY k) const = 0;
};

template <typename RECVMESH>
class DestinationMeshInterface
{
public:
  using BOUNDINGBOX = typename MeshTraits<RECVMESH>::BoundingBox;
  using ENTITYKEY = typename MeshTraits<RECVMESH>::EntityKey;

  DestinationMeshInterface() = default;
  virtual ~DestinationMeshInterface() = default;

  virtual stk::ParallelMachine comm() const = 0;

  virtual std::string name() const = 0;

  virtual void bounding_boxes(std::vector<BOUNDINGBOX>& v) const = 0;

  virtual const double* coord(const ENTITYKEY k) const = 0;
  virtual double get_search_tolerance() const = 0;
  virtual double get_parametric_tolerance() const = 0;

  virtual void centroid(const ENTITYKEY k, std::vector<double>& centroidVec) const = 0;
  virtual double get_distance_from_nearest_node(const ENTITYKEY k, const double* toCoords) const = 0;
};

}}

#endif /* STK_STK_SEARCH_STK_SEARCH_SEARCHINTERFACE_HPP_ */
