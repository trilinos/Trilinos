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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_SIMPLETRANSFER_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_SIMPLETRANSFER_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_transfer_util/spmd/GeometricTransferOptions.hpp>
#include <stk_transfer_util/PointInterpolation.hpp>

#include <ostream>
#include <memory>                                       // for shared_ptr
#include <string>                                       // for string
#include <utility>                                      // for pair
#include <vector>                                       // for vector
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace transfer { namespace spmd { class NodeSendMesh; } } }
namespace stk { namespace transfer { namespace spmd { class NodeRecvMesh; } } }
namespace stk { namespace transfer { namespace spmd { class ElementSendMesh; } } }
namespace stk { namespace transfer { namespace spmd { class ElementRecvMesh; } } }
namespace stk { namespace transfer { namespace spmd { class GeometricTransferDispatchBase; } } }

namespace stk { namespace transfer { class TransferBase; } }

namespace stk {
namespace transfer {
namespace spmd {

class SimpleTransfer {
 public:
  using MasterElementProvider = std::shared_ptr<stk::search::MasterElementProviderInterface>;

  SimpleTransfer(const std::string& transferName)
  : m_transferName(transferName)
  , m_hasParallelMachine(false)
  { }

  SimpleTransfer(const std::string& transferName, stk::ParallelMachine pm)
  : m_transferName(transferName)
  , m_hasParallelMachine(true)
  , m_parallelMachine(pm)
  { }


  virtual ~SimpleTransfer() = default;

  virtual void add_send_field(stk::mesh::FieldBase & field);
  virtual void add_send_field(const stk::transfer::FieldSpec& fieldSpec);
  virtual void add_send_fields(const std::vector<stk::mesh::FieldBase*>& fields);
  virtual void add_send_fields(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);

  virtual void add_recv_field(stk::mesh::FieldBase & field);
  virtual void add_recv_field(const stk::transfer::FieldSpec& fieldSpec);
  virtual void add_recv_fields(const std::vector<stk::mesh::FieldBase*>& fields);
  virtual void add_recv_fields(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);

  virtual void add_send_part_name(const std::string& partName);
  virtual void add_send_part_names(const std::vector<std::string>& partNames);

  virtual void add_recv_part_name(const std::string& partName);
  virtual void add_recv_part_names(const std::vector<std::string>& partNames);

  void setup_patch_recovery_transfer(stk::mesh::BulkData& sendBulk,
                                     stk::mesh::BulkData& recvBulk,
                                     stk::mesh::EntityRank sendRank,
                                     stk::transfer::spmd::RecvMeshType recvType,
                                     MasterElementProvider masterElemProvider,
                                     stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                         stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_patch_recovery_transfer(stk::mesh::BulkData& sendBulk,
                                     stk::mesh::BulkData& recvBulk,
                                     stk::mesh::EntityRank sendRank,
                                     stk::transfer::spmd::RecvMeshType recvType,
                                     MasterElementProvider sendMasterElemProvider,
                                     MasterElementProvider recvMasterElemProvider,
                                     stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                         stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_master_element_transfer(stk::mesh::BulkData& sendBulk,
                                     stk::mesh::BulkData& recvBulk,
                                     stk::mesh::EntityRank sendRank,
                                     stk::transfer::spmd::RecvMeshType recvType,
                                     MasterElementProvider masterElemProvider,
                                     stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                         stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_master_element_transfer(stk::mesh::BulkData& sendBulk,
                                     stk::mesh::BulkData& recvBulk,
                                     stk::mesh::EntityRank sendRank,
                                     stk::transfer::spmd::RecvMeshType recvType,
                                     MasterElementProvider sendMasterElemProvider,
                                     MasterElementProvider recvMasterElemProvider,
                                     stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                         stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_copy_nearest_transfer(stk::mesh::BulkData& sendBulk,
                                   stk::mesh::BulkData& recvBulk,
                                   stk::mesh::EntityRank sendRank,
                                   stk::transfer::spmd::RecvMeshType recvType,
                                   MasterElementProvider masterElemProvider,
                                   stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                       stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_copy_nearest_transfer(stk::mesh::BulkData& sendBulk,
                                   stk::mesh::BulkData& recvBulk,
                                   stk::mesh::EntityRank sendRank,
                                   stk::transfer::spmd::RecvMeshType recvType,
                                   MasterElementProvider sendMasterElemProvider,
                                   MasterElementProvider recvMasterElemProvider,
                                   stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                       stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_sum_nearest_transfer(stk::mesh::BulkData& sendBulk,
                                  stk::mesh::BulkData& recvBulk,
                                  stk::mesh::EntityRank sendRank,
                                  stk::transfer::spmd::RecvMeshType recvType,
                                  MasterElementProvider masterElemProvider,
                                  stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                      stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void setup_sum_nearest_transfer(stk::mesh::BulkData& sendBulk,
                                  stk::mesh::BulkData& recvBulk,
                                  stk::mesh::EntityRank sendRank,
                                  stk::transfer::spmd::RecvMeshType recvType,
                                  MasterElementProvider sendMasterElemProvider,
                                  MasterElementProvider recvMasterElemProvider,
                                  stk::search::ObjectOutsideDomainPolicy extrapolateOption =
                                      stk::search::ObjectOutsideDomainPolicy::IGNORE);

  void initialize(bool overrideParametricSearchWithGeometricClosestNodeToCentroid = false,
                  bool useCentroidForGeometricSearchProximity = true,
                  bool printSearchWarnings = true);

  void apply();

  std::shared_ptr<stk::transfer::TransferBase> get_transfer() const { return m_transfer; }
  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> get_dispatch() const { return m_dispatch; }

  void set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType evalType) { m_defaultPatchRecovery = evalType; }

 protected:
  bool m_committed{false};
  std::string m_transferName{"UN-INITIALIZED TRANSFER NAME"};

  bool m_hasParallelMachine{false};
  stk::ParallelMachine m_parallelMachine{stk::parallel_machine_world()};

  std::vector<stk::transfer::FieldSpec> m_sendFieldSpecs;
  std::vector<stk::transfer::FieldSpec> m_recvFieldSpecs;

  std::vector<std::string> m_sendPartNames;
  std::vector<std::string> m_recvPartNames;

  stk::transfer::PatchRecoveryEvaluationType m_defaultPatchRecovery{stk::transfer::PatchRecoveryEvaluationType::LINEAR_LEAST_SQUARES};

  std::shared_ptr<stk::transfer::TransferBase> m_transfer;
  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> m_dispatch;
  std::shared_ptr<stk::transfer::spmd::GeometricTransferOptions> m_transferOptions;

  stk::mesh::EntityRank recv_type_to_entity_rank(stk::transfer::spmd::RecvMeshType type) const;
  stk::transfer::spmd::SendMeshType entity_rank_to_send_type(stk::mesh::EntityRank rank) const;

  void set_send_mesh_transfer_options(
      stk::mesh::BulkData& sendBulk,
      stk::mesh::EntityRank sendRank,
      stk::search::ObjectOutsideDomainPolicy extrapolateOption);

  void set_send_mesh_transfer_options(
      stk::mesh::BulkData& sendBulk,
      stk::mesh::EntityRank sendRank,
      MasterElementProvider masterElemProvider,
      stk::search::ObjectOutsideDomainPolicy extrapolateOption);

  void set_recv_mesh_transfer_options(stk::mesh::BulkData& recvBulk, stk::transfer::spmd::RecvMeshType recvType);

  void set_recv_mesh_transfer_options(stk::mesh::BulkData& recvBulk,
                                      stk::transfer::spmd::RecvMeshType recvType,
                                      MasterElementProvider masterElemProvider);

  void initialize_meshes();

  void set_transfer_options(bool useClosestNodeToCentroid, bool useCentroidForGeometricProximity, bool printSearchWarnings);

  void execute_search();
};

} // namespace spmd
} // namespace transfer
} // namespace stk


#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_SIMPLETRANSFER_HPP_ */
