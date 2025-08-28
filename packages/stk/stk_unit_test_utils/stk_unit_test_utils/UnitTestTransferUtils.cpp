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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_unit_test_utils/UnitTestTransferUtils.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

namespace impl {
std::shared_ptr<stk::transfer::spmd::ElementSendMesh>
construct_elem_send_mesh(stk::mesh::BulkData& bulk,
                         const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                         double parametricTolerance,
                         double geometricTolerance,
                         const std::vector<std::string>& blocks,
                         stk::search::ObjectOutsideDomainPolicy policy = stk::search::ObjectOutsideDomainPolicy::IGNORE,
                         const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement())
{
  STK_ThrowRequireMsg(blocks.size() > 0, "Input block list is empty");

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = get_parts_list(meta, blocks);

  const stk::mesh::FieldBase* coords = meta.coordinate_field();

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>();

  std::shared_ptr<stk::transfer::InterpolateFieldsInterface> interpolateFields =
      std::make_shared<stk::transfer::MasterElementFieldInterpolator>(bulk, masterElemProvider);

  std::shared_ptr<stk::search::MasterElementParametricCoordsFinder> paramCoordsFinder =
      std::make_shared<stk::search::MasterElementParametricCoordsFinder>(bulk, coords, masterElemProvider, parametricTolerance);

  std::shared_ptr<stk::search::HandleExternalPointInterface> externalPointHandler;

  if(policy == stk::search::ObjectOutsideDomainPolicy::TRUNCATE) {
   externalPointHandler =
       std::make_shared<stk::search::MasterElementExternalPointTruncation>(bulk, coords, masterElemProvider, geometricTolerance);
  } else if(policy == stk::search::ObjectOutsideDomainPolicy::PROJECT) {
    externalPointHandler =
        std::make_shared<stk::search::MasterElementExternalPointProjection>(bulk, coords, masterElemProvider,
                                                                            parametricTolerance, geometricTolerance);
  } else {
    externalPointHandler = std::make_shared<stk::search::ExternalPointNoOpHandler>(parametricTolerance);
  }

  std::shared_ptr<stk::transfer::spmd::ElementSendMesh> mesh =
      std::make_shared<stk::transfer::spmd::ElementSendMesh>(&bulk, coords, fieldSpecs, parts[0]->primary_entity_rank(),
                                                             parts, activeSelector, bulk.parallel(), paramCoordsFinder,
                                                             externalPointHandler, interpolateFields, masterElemProvider);

  mesh->initialize();

  return mesh;
}

std::shared_ptr<stk::transfer::spmd::NodeRecvMesh>
construct_node_recv_mesh(stk::mesh::BulkData& bulk,
                         const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                         double parametricTolerance, double geometricTolerance,
                         const std::vector<std::string>& blocks,
                         const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement())
{
  STK_ThrowRequireMsg(blocks.size() > 0, "Input block list is empty");

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = get_parts_list(meta, blocks);

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::transfer::spmd::NodeRecvMesh> mesh = std::make_shared<stk::transfer::spmd::NodeRecvMesh>(
      &bulk, coords, fieldSpecs, parts, activeSelector, bulk.parallel(), parametricTolerance, geometricTolerance);

  return mesh;
}

std::shared_ptr<stk::transfer::spmd::ElementRecvMesh>
construct_elem_recv_mesh(stk::mesh::BulkData& bulk,
                         const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                         std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator,
                         double parametricTolerance, double geometricTolerance,
                         const std::vector<std::string>& blocks,
                         const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement())
{
  STK_ThrowRequireMsg(blocks.size() > 0, "Input block list is empty");

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = get_parts_list(meta, blocks);

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>();

  std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> mesh =
      std::make_shared<stk::transfer::spmd::ElementRecvMesh>(&bulk, coords, fieldSpecs, parts[0]->primary_entity_rank(),
                                                             parts, activeSelector, bulk.parallel(), pointEvaluator,
                                                             parametricTolerance, geometricTolerance);

  return mesh;
}

}

std::shared_ptr<stk::transfer::spmd::ElementSendMesh>
construct_element_send_mesh(stk::mesh::BulkData& bulk, const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                            double parametricTolerance, double geometricTolerance,
                            stk::search::ObjectOutsideDomainPolicy policy,
                            const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  return impl::construct_elem_send_mesh(bulk, fieldSpecs, parametricTolerance, geometricTolerance,
                                        blocks, policy, activeSelector);
}

std::shared_ptr<stk::transfer::spmd::NodeRecvMesh>
construct_node_recv_mesh(stk::mesh::BulkData& bulk, const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                         double parametricTolerance, double geometricTolerance,
                         const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  return impl::construct_node_recv_mesh(bulk, fieldSpecs, parametricTolerance, geometricTolerance, blocks, activeSelector);
}

std::shared_ptr<stk::transfer::spmd::ElementRecvMesh>
construct_element_centroid_recv_mesh(stk::mesh::BulkData& bulk, const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                                     double parametricTolerance, double geometricTolerance,
                                     const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();

  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::CentroidEvaluator>(bulk, coords);

  return impl::construct_elem_recv_mesh(bulk, fieldSpecs, pointEvaluator,
                                        parametricTolerance, geometricTolerance,
                                        blocks, activeSelector);
}

std::shared_ptr<stk::transfer::spmd::ElementRecvMesh>
construct_element_gauss_point_recv_mesh(stk::mesh::BulkData& bulk, const std::vector<stk::transfer::FieldSpec>& fieldSpecs,
                                        double parametricTolerance, double geometricTolerance,
                                        unsigned integrationOrder, const std::vector<std::string>& blocks,
                                        const stk::mesh::Selector& activeSelector)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>(integrationOrder);

  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::MasterElementGaussPointEvaluator>(bulk, coords, masterElemProvider);

  return impl::construct_elem_recv_mesh(bulk, fieldSpecs, pointEvaluator,
                                        parametricTolerance, geometricTolerance,
                                        blocks, activeSelector);
}



}
}
