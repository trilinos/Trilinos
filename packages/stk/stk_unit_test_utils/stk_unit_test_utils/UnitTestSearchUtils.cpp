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
#include "UnitTestSearchUtils.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

namespace impl {
std::shared_ptr<stk::search::spmd::ElementSendMesh>
construct_send_mesh(stk::mesh::BulkData& bulk, double parametricTolerance,
                    const std::vector<std::string>& blocks,
                    const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement())
{
  STK_ThrowRequireMsg(blocks.size() > 0, "Input block list is empty");

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = get_parts_list(meta, blocks);

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>();
  std::shared_ptr<stk::search::MasterElementParametricCoordsFinder> paramCoordsFinder =
      std::make_shared<stk::search::MasterElementParametricCoordsFinder>(bulk, coords, masterElemProvider, parametricTolerance);
  std::shared_ptr<stk::search::HandleExternalPointInterface> externalPointHandler =
      std::make_shared<stk::search::ExternalPointNoOpHandler>(parametricTolerance);
  std::shared_ptr<stk::search::spmd::ElementSendMesh> mesh =
      std::make_shared<stk::search::spmd::ElementSendMesh>(&bulk, coords, parts[0]->primary_entity_rank(), parts, activeSelector,
                                                           bulk.parallel(), paramCoordsFinder, externalPointHandler, masterElemProvider);

  mesh->initialize();

  return mesh;
}

std::shared_ptr<stk::search::spmd::NodeRecvMesh>
construct_node_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                         const std::vector<std::string>& blocks,
                         const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement())
{
  STK_ThrowRequireMsg(blocks.size() > 0, "Input block list is empty");

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = get_parts_list(meta, blocks);

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::spmd::NodeRecvMesh> mesh = std::make_shared<stk::search::spmd::NodeRecvMesh>(
      &bulk, coords, parts, activeSelector, bulk.parallel(), parametricTolerance, geometricTolerance);

  return mesh;
}

std::shared_ptr<stk::search::spmd::ElementRecvMesh>
construct_elem_recv_mesh(stk::mesh::BulkData& bulk,
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

  std::shared_ptr<stk::search::spmd::ElementRecvMesh> mesh =
      std::make_shared<stk::search::spmd::ElementRecvMesh>(&bulk, coords, parts[0]->primary_entity_rank(), parts,
                                                           activeSelector, bulk.parallel(), pointEvaluator,
                                                           parametricTolerance, geometricTolerance);

  return mesh;
}

}

stk::mesh::PartVector get_parts_list(const stk::mesh::MetaData& meta, const std::vector<std::string>& blocks)
{
  stk::mesh::PartVector parts;
  for(const std::string& block : blocks) {
    stk::mesh::Part* part = meta.get_part(block);
    STK_ThrowRequireMsg(nullptr != part, "Block named: " << block << " does not exist");
    parts.push_back(part);
  }

  return parts;
}

std::vector<std::string> get_all_block_names(const stk::mesh::MetaData& meta)
{
  std::vector<std::string> blocks;

  for(auto part : meta.get_mesh_parts()) {
    if(stk::mesh::is_element_block(*part)) {
      blocks.push_back(part->name());
    }
  }

  return blocks;
}

std::shared_ptr<stk::mesh::BulkData> build_slanted_single_hex_bulk(double slant)
{
  assert(slant > 0);

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);
  std::string meshSpec = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1";
  std::vector<double> coords{ 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,slant, 1,1,slant, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshSpec, coords));

  return bulk;
}

std::shared_ptr<stk::mesh::BulkData> build_single_point_bulk(double x, double y, double z)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);

  std::string meshSpec = "0,1,PARTICLE,1,block_1";
  std::vector<double> coords{ x, y, z };
  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshSpec, coords));

  return bulk;
}

std::shared_ptr<stk::search::spmd::ElementSendMesh>
construct_hex_send_mesh(stk::mesh::BulkData& bulk, double parametricTolerance,
                        const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  return impl::construct_send_mesh(bulk, parametricTolerance, blocks, activeSelector);
}

std::shared_ptr<stk::search::spmd::NodeRecvMesh>
construct_node_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                        const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  return impl::construct_node_recv_mesh(bulk, parametricTolerance, geometricTolerance, blocks, activeSelector);
}

std::shared_ptr<stk::search::spmd::ElementRecvMesh>
construct_element_centroid_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                                     const std::vector<std::string>& blocks, const stk::mesh::Selector& activeSelector)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();

  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::CentroidEvaluator>(bulk, coords);
  return impl::construct_elem_recv_mesh(bulk, pointEvaluator, parametricTolerance, geometricTolerance, blocks, activeSelector);
}

std::shared_ptr<stk::search::spmd::ElementRecvMesh>
construct_hex_gauss_point_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                                    unsigned /*integrationOrder*/, const std::vector<std::string>& blocks,
                                    const stk::mesh::Selector& activeSelector)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>();
  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::MasterElementGaussPointEvaluator>(bulk, coords, masterElemProvider);
  return impl::construct_elem_recv_mesh(bulk, pointEvaluator, parametricTolerance, geometricTolerance, blocks, activeSelector);
}

}
}

