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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "stk_search/ObjectOutsideDomainPolicy.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machin...

namespace stk
{
namespace search
{
template <typename TOPOLOGY, typename ENTITYKEY>
class EvaluatePointsInterface
{
 public:
  using Topology  = TOPOLOGY;
  using EntityKey = ENTITYKEY;

  virtual size_t num_points(const ENTITYKEY&, const TOPOLOGY&) = 0;
  virtual void coordinates(const ENTITYKEY&, size_t, std::vector<double>& coords) = 0;

  virtual ~EvaluatePointsInterface() {}
};

template <typename ENTITYKEY>
class ExternalPointHandlerInterface
{
 public:
  using EntityKey = ENTITYKEY;

  virtual bool handle_point(const ENTITYKEY& k,
      const std::vector<double>& point,
      std::vector<double>& parametricCoords,
      double& geometricDistanceSquared,
      bool& isWithinGeometricTolerance) const = 0;

  virtual ~ExternalPointHandlerInterface() {}
};

template <typename ENTITYKEY>
class FindParametricCoordinatesInterface
{
 public:
  using EntityKey = ENTITYKEY;

  virtual void find_parametric_coords(const ENTITYKEY& k,
      const std::vector<double>& point,
      std::vector<double>& paramCoords,
      double& paramDistance,
      bool& isWithinParametricTolerance) const = 0;
  virtual void evaluate_parametric_coords(
      const ENTITYKEY& k, const std::vector<double>& paramCoords, std::vector<double>& evalPoint) const = 0;
  virtual ~FindParametricCoordinatesInterface() {}
};

template <typename TOPOLOGY, typename ENTITYKEY, typename FIELD>
class ProvideMasterElementInterface
{
 public:
  using Topology  = TOPOLOGY;
  using EntityKey = ENTITYKEY;

  virtual void evaluate_field(const TOPOLOGY& topo,
      const std::vector<double>& paramCoords,     // (numParamCoords)
      const unsigned numFieldComponents,
      const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
      std::vector<double>& result) const = 0;     // (numFieldComponents)

  virtual void evaluate_field(const TOPOLOGY& topo,
      const unsigned numEvalPoints,
      const std::vector<double>& paramCoords,     // (numParamCoords x numEvalPoints)
      const unsigned numFieldComponents,
      const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
      std::vector<double>& result) const = 0;     // (numFieldComponents x numEvalPoints)

  virtual void nodal_field_data(const ENTITYKEY& key,
      const FIELD& field,
      unsigned& numFieldComponents,
      unsigned& numNodes,
      std::vector<double>& fieldData) const = 0;  // (numFieldComponents x numTopologyNodes)

  virtual void nodal_field_data(const std::vector<ENTITYKEY>& nodeKeys,
      const FIELD& field,
      unsigned& numFieldComponents,
      std::vector<double>& fieldData) const = 0;  // (numFieldComponents x numTopologyNodes)

  virtual void find_parametric_coordinates(const TOPOLOGY& topo,
      const unsigned numCoordComponents,
      const std::vector<double>& elemNodeCoords,  // (numCoordComponents x numTopologyNodes)
      const std::vector<double>& inputCoords,     // (numCoordComponents)
      std::vector<double>& paramCoords,
      double& paramDistance) const = 0;

  virtual void coordinate_center(const TOPOLOGY& topo, std::vector<double>& coords) const = 0;
  virtual unsigned num_parametric_coordinates(const TOPOLOGY& topo) const = 0;
  virtual unsigned num_integration_points(const TOPOLOGY& topo) const = 0;
  virtual void integration_points(const TOPOLOGY& topo, std::vector<double>& gaussPoints) const = 0;
  virtual ~ProvideMasterElementInterface() {}
};

//BEGINSearch_Interface
template <typename MESH>
struct MeshTraits;

template <typename SENDMESH>
class SourceMeshInterface
{
 public:
  // static polymorphism with CRTP (Curiously Recurring Template Pattern)
  using Entity = typename MeshTraits<SENDMESH>::Entity;
  using EntityVec = typename MeshTraits<SENDMESH>::EntityVec;
  using EntityKey = typename MeshTraits<SENDMESH>::EntityKey;
  using EntityProc = typename MeshTraits<SENDMESH>::EntityProc;
  using EntityProcVec = typename MeshTraits<SENDMESH>::EntityProcVec;
  using Point = typename MeshTraits<SENDMESH>::Point;
  using Box = typename MeshTraits<SENDMESH>::Box;
  using Sphere = typename MeshTraits<SENDMESH>::Sphere;
  using BoundingBox = typename MeshTraits<SENDMESH>::BoundingBox;

  SourceMeshInterface() = default;
  virtual ~SourceMeshInterface() = default;

  virtual stk::ParallelMachine comm() const = 0;

  virtual std::string name() const = 0;

  virtual void set_name(const std::string& meshName) = 0;

  virtual void initialize() = 0;

  virtual ObjectOutsideDomainPolicy get_extrapolate_option() const = 0;

  virtual void update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix = "") = 0;

  virtual void update_ghosted_key(EntityKey& k) = 0;

  virtual void destroy_ghosting() = 0;

  virtual void post_mesh_modification_event() = 0;

  virtual std::vector<std::string> get_part_membership(const EntityKey& k) const = 0;

  virtual void bounding_boxes(std::vector<BoundingBox>& boxes, bool includeGhosts=false) const = 0;

  virtual void find_parametric_coords(
    const EntityKey& k,
    const std::vector<double>& point,
    std::vector<double>& parametricCoords,
    double& parametricDistance,
    bool& isWithinParametricTolerance) const = 0;

  virtual bool modify_search_outside_parametric_tolerance(
    const EntityKey& k,
    const std::vector<double>& point,
    std::vector<double>& parametricCoords,
    double& geometricDistanceSquared,
    bool& isWithinGeometricTolerance) const = 0;

  virtual double get_distance_from_nearest_node(
    const EntityKey& k, const std::vector<double>& point) const = 0;

  virtual double get_closest_geometric_distance_squared(
    const EntityKey& k, const std::vector<double>& point) const = 0;

  virtual double get_distance_from_centroid(
    const EntityKey& k, const std::vector<double>& point) const = 0;

  virtual double get_distance_squared_from_centroid(
    const EntityKey& k, const std::vector<double>& point) const = 0;

  virtual void centroid(const EntityKey& k, std::vector<double>& centroidVec) const = 0;

  virtual void coordinates(const EntityKey& k, std::vector<double>& coords) const = 0;
};

template <typename RECVMESH>
class DestinationMeshInterface
{
 public:
  // static polymorphism with CRTP (Curiously Recurring Template Pattern)
  using Entity = typename MeshTraits<RECVMESH>::Entity;
  using EntityVec = typename MeshTraits<RECVMESH>::EntityVec;
  using EntityKey = typename MeshTraits<RECVMESH>::EntityKey;
  using EntityProc = typename MeshTraits<RECVMESH>::EntityProc;
  using EntityProcVec = typename MeshTraits<RECVMESH>::EntityProcVec;
  using Point = typename MeshTraits<RECVMESH>::Point;
  using Sphere = typename MeshTraits<RECVMESH>::Sphere;
  using Box = typename MeshTraits<RECVMESH>::Box;
  using BoundingBox = typename MeshTraits<RECVMESH>::BoundingBox;

  DestinationMeshInterface() = default;
  virtual ~DestinationMeshInterface() = default;

  virtual stk::ParallelMachine comm() const = 0;

  virtual std::string name() const = 0;

  virtual void set_name(const std::string& meshName) = 0;

  virtual void initialize() = 0;

  virtual void bounding_boxes(std::vector<BoundingBox>& v) const = 0;

  virtual void coordinates(const EntityKey& k, std::vector<double>& coords) const = 0;

  virtual double get_search_tolerance() const = 0;

  virtual double get_parametric_tolerance() const = 0;

  virtual void centroid(const EntityKey& k, std::vector<double>& centroidVec) const = 0;

  virtual double get_distance_from_nearest_node(
    const EntityKey& k, const std::vector<double>& point) const = 0;
};
//ENDSearch_Interface

}  // namespace search
}  // namespace stk

#endif /* STK_STK_SEARCH_STK_SEARCH_SEARCHINTERFACE_HPP_ */
