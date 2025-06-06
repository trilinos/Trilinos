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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERUTILS_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERUTILS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_transfer/TransferBase.hpp>
#include <stk_transfer_util/spmd/ElementRecvMesh.hpp>
#include <stk_transfer_util/spmd/ElementSendMesh.hpp>
#include <stk_transfer_util/spmd/NodeRecvMesh.hpp>
#include <stk_transfer_util/spmd/NodeSendMesh.hpp>
#include <stk_transfer_util/spmd/GeometricInterp.hpp>
#include <stk_transfer_util/spmd/GeometricTransfer.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <memory>
#include <functional>
#include <string>
#include <vector>

namespace stk::mesh { class BulkData; }
namespace stk::mesh { class FieldBase; }
namespace stk::mesh { class MetaData; }
namespace stk::mesh { class Part; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

using NodeLoopFunction = std::function<void(const stk::mesh::BulkData& sendBulk, const stk::mesh::EntityKey& sendEntityKey,
                                            const stk::mesh::BulkData& recvBulk, const stk::mesh::EntityKey& recvEntityKey)>;

using ElemLoopFunction = std::function<void(const stk::mesh::BulkData& sendBulk, const stk::mesh::EntityKey& sendEntityKey,
                                            const stk::mesh::BulkData& recvBulk, const stk::mesh::EntityKey& recvEntityKey, int recvEntityIndex)>;

template <typename TRANSFER, typename INTERPOLATE>
void for_each_node_transfer_pair_run(std::shared_ptr<stk::transfer::TransferBase> transfer,
                                     const stk::mesh::BulkData& sendBulk,
                                     const stk::mesh::BulkData& recvBulk,
                                     NodeLoopFunction& functor)
{
  using SEND        = typename INTERPOLATE::MeshA;
  using RECV        = typename INTERPOLATE::MeshB;
  using RelationVec = typename TRANSFER::EntityProcRelationVec;

  static_assert(std::is_base_of<stk::transfer::spmd::GeometricTransfer<INTERPOLATE>, TRANSFER>::value ,
                "TRANSFER must be a derived class of stk::transfer::spmd/GeometricTransfer.");

  static_assert(std::is_base_of<stk::transfer::spmd::GeometricInterp<SEND, RECV>, INTERPOLATE>::value ,
                "INTERPOLATE must be a derived class of stk::transfer::spmd/GeometricInterp.");

  std::shared_ptr<TRANSFER> castTransfer = std::dynamic_pointer_cast<TRANSFER>(transfer);
  RelationVec& relationVec = castTransfer->get_range_to_domain();

  const std::shared_ptr<SEND> sendMesh = castTransfer->send_mesh();
  const std::shared_ptr<RECV> recvMesh = castTransfer->recv_mesh();

  for(auto& relation : relationVec) {
    const stk::mesh::EntityKey sendEntityKey = relation.second.id();
    const stk::mesh::EntityKey recvEntityKey = relation.first.id();

    functor(sendBulk, sendEntityKey, recvBulk, recvEntityKey);
  }
}

template <typename TRANSFER, typename INTERPOLATE>
void for_each_element_transfer_pair_run(std::shared_ptr<stk::transfer::TransferBase> transfer,
                                        const stk::mesh::BulkData& sendBulk,
                                        const stk::mesh::BulkData& recvBulk,
                                        const ElemLoopFunction& functor)
{
  using SEND        = typename INTERPOLATE::MeshA;
  using RECV        = typename INTERPOLATE::MeshB;
  using RelationVec = typename TRANSFER::EntityProcRelationVec;

  static_assert(std::is_base_of<stk::transfer::spmd::GeometricTransfer<INTERPOLATE>, TRANSFER>::value ,
                "TRANSFER must be a derived class of stk::transfer::spmd/GeometricTransfer.");

  static_assert(std::is_base_of<stk::transfer::spmd::GeometricInterp<SEND, RECV>, INTERPOLATE>::value ,
                "INTERPOLATE must be a derived class of stk::transfer::spmd/GeometricInterp.");

  std::shared_ptr<TRANSFER> castTransfer = std::dynamic_pointer_cast<TRANSFER>(transfer);
  RelationVec& relationVec = castTransfer->get_range_to_domain();

  const std::shared_ptr<SEND> sendMesh = castTransfer->send_mesh();
  const std::shared_ptr<RECV> recvMesh = castTransfer->recv_mesh();

  for(auto& relation : relationVec) {
    const stk::mesh::EntityKey sendEntityKey = relation.second.id();
    const stk::mesh::EntityKey recvEntityKey = relation.first.id().first;
    int recvEntityIndex = relation.first.id().second;

    functor(sendBulk, sendEntityKey, recvBulk, recvEntityKey, recvEntityIndex);
  }
}

stk::mesh::EntityRank get_transfer_rank(const stk::mesh::MetaData& meta,
                                        const std::vector<std::string>& partNames,
                                        const std::string& transferName);

std::vector<std::string> prune_empty_assemblies(const stk::mesh::MetaData& meta,
                                                const std::vector<std::string>& partNames,
                                                const std::string& transferName);

stk::mesh::PartVector get_parts(const stk::mesh::MetaData& meta,
                                const std::vector<std::string>& partNames,
                                const std::string& transferName,
                                const stk::mesh::Part* defaultPart = nullptr);

double mesh_bounding_box_geometric_tolerance(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coordinateField, double scaleFactor);

} // namespace spmd
} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERUTILS_HPP_ */
