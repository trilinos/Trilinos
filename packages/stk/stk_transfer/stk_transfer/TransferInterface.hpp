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

#ifndef STK_STK_TRANSFER_STK_TRANSFER_TRANSFERINTERFACE_HPP_
#define STK_STK_TRANSFER_STK_TRANSFER_TRANSFERINTERFACE_HPP_

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

#include "stk_util/parallel/Parallel.hpp"  // for parallel_machin...
#include "stk_transfer/TransferTypes.hpp"

namespace stk
{
namespace transfer
{

template <typename ENTITYKEY>
class FieldInterpolatorInterface {
 public:
  virtual void interpolate_fields(const ENTITYKEY& k, const std::vector<double>& evalPoint,
                                  const std::vector<double>& parametricCoords, InterpolationData& data) const = 0;
  virtual ~FieldInterpolatorInterface() {}
};


//BEGINTransfer_Interface
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

  virtual void update_values() = 0;

  virtual void initialize() = 0;

  virtual void interpolate_fields(const EntityKey& k, std::vector<double>& recvRealCoords,
                                  std::vector<double>& sendParametricCoords, InterpolationData& data) const = 0;

  virtual void destroy_ghosting() = 0;

  virtual void update_ghosting(const EntityProcVec& entity_keys, const std::string& suffix = "") = 0;

  virtual std::vector<std::string> get_part_membership(const EntityKey& k) const = 0;
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

  virtual std::string name() const = 0;

  virtual void set_name(const std::string& meshName) = 0;

  virtual void initialize() = 0;

  virtual double* value(const EntityKey& k, const unsigned fieldIndex) const = 0;

  virtual unsigned value_size(const EntityKey& k, const unsigned fieldIndex) const = 0;

  virtual unsigned num_values(const EntityKey& e) const = 0;

  virtual unsigned max_num_values() const = 0;

  virtual unsigned value_key(const EntityKey& k, const unsigned fieldIndex) const = 0;

  virtual void update_values() = 0;

  virtual unsigned get_index(const unsigned i) const = 0;
};
//ENDTransfer_Interface

}  // namespace transfer
}  // namespace stk



#endif /* STK_STK_TRANSFER_STK_TRANSFER_TRANSFERINTERFACE_HPP_ */
