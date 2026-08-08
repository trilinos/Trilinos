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

#include "gtest/gtest.h"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_search/ObjectOutsideDomainPolicy.hpp"

namespace {

TEST(TransferMainSettings, defaultSettings)
{
  stk::transfer_util::TransferMainSettings settings;

  EXPECT_EQ(settings.get_num_input_processors(), 1u); 
  EXPECT_EQ(settings.get_num_output_processors(), 1u); 
  EXPECT_EQ(settings.get_sendMesh_filename(), "");
  EXPECT_EQ(settings.get_recvMesh_filename(), "");
  EXPECT_EQ(settings.get_transfer_fields().size(), 0u);

  auto settingsOODP = settings.get_extrapolate_option();
  std::string settingsOODPConvertedString = stk::search::get_object_outside_domain_policy(settingsOODP);
  std::string settingsOODPString = settings.get_extrapolate_option_string();
  EXPECT_EQ(settingsOODPConvertedString, "IGNORE");
  EXPECT_EQ(settingsOODPString, "IGNORE");

  EXPECT_EQ(settings.get_recv_type_string(), "NODE");

  EXPECT_EQ(settings.get_coord_transf_z_expr(), "z=0");
}

TEST(TransferMainSettings, filenames)
{
  stk::transfer_util::TransferMainSettings settings;
 
  EXPECT_EQ(settings.get_sendMesh_filename(), "");
  EXPECT_EQ(settings.get_recvMesh_filename(), "");
  EXPECT_EQ(settings.get_outputMesh_filename(), "");
  
  settings.set_sendMesh_filename("source.g");

  EXPECT_EQ(settings.get_sendMesh_filename(), "source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "");
  EXPECT_EQ(settings.get_outputMesh_filename(), "");

  settings.set_recvMesh_filename("target.g");

  EXPECT_EQ(settings.get_sendMesh_filename(), "source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "target.g");
  EXPECT_EQ(settings.get_outputMesh_filename(), "");

  settings.set_outputMesh_filename("output.exo");

  EXPECT_EQ(settings.get_sendMesh_filename(), "source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "target.g");
  EXPECT_EQ(settings.get_outputMesh_filename(), "output.exo");

  settings.set_sendMesh_filename("changed_source.g");
  settings.set_recvMesh_filename("changed_target.g");
  settings.set_outputMesh_filename("changed_output.exo");

  EXPECT_EQ(settings.get_sendMesh_filename(), "changed_source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "changed_target.g");
  EXPECT_EQ(settings.get_outputMesh_filename(), "changed_output.exo");
}

TEST(TransferMainSettings, processorCounts)
{
  stk::transfer_util::TransferMainSettings settings;

  EXPECT_EQ(settings.get_num_input_processors(), 1u); 
  EXPECT_EQ(settings.get_num_output_processors(), 1u); 

  settings.set_num_input_processors(8); 
  settings.set_num_output_processors(4); 

  EXPECT_EQ(settings.get_num_input_processors(), 8u); 
  EXPECT_EQ(settings.get_num_output_processors(), 4u); 
}

TEST(TransferMainSettings, fields)
{
  stk::transfer_util::TransferMainSettings settings;
  
  EXPECT_EQ(settings.get_transfer_fields().size(), 0u);
   
  std::vector<std::pair<std::string, std::string>> fieldList = {{"send_field_1", "recv_field_1"},
                                                               {"send_field_2","recv_field_2"}};
  for (auto field : fieldList) {
    settings.set_transfer_field(field);
  }

  std::vector<std::pair<std::string,std::string>> setFields = settings.get_transfer_fields();
  EXPECT_EQ(setFields, fieldList);

  std::string fieldListString = "send_field_1:recv_field_1, send_field_2:recv_field_2";
  EXPECT_EQ(settings.get_field_list_string(), fieldListString);
}

TEST(TransferMainSettings, parts)
{
  stk::transfer_util::TransferMainSettings settings;
  
  EXPECT_EQ(settings.get_transfer_send_parts().size(), 0u);
  EXPECT_EQ(settings.get_transfer_recv_parts().size(), 0u);
   
  std::vector<std::string> sendPartList = {"send_part_1", "send_part_2"};
  std::vector<std::string> recvPartList = {"recv_part_1", "recv_part_2", "recv_part_3"};

  settings.set_transfer_send_parts(sendPartList);
  settings.set_transfer_recv_parts(recvPartList);

  std::vector<std::string> sendParts = settings.get_transfer_send_parts();
  EXPECT_EQ(sendPartList, sendParts);

  std::vector<std::string> recvParts = settings.get_transfer_recv_parts();
  EXPECT_EQ(recvPartList, recvParts);

  std::string partListString = "send_part_1, send_part_2:recv_part_1, recv_part_2, recv_part_3";
  EXPECT_EQ(settings.get_part_list_string(), partListString);
}

TEST(TransferMainSettings, extrapolateOption)
{
  stk::transfer_util::TransferMainSettings settings;

  settings.set_extrapolate_option("IGNORE");
  auto settingsOODPIgnore = settings.get_extrapolate_option();
  std::string settingsOODPIgnoreString = stk::search::get_object_outside_domain_policy(settingsOODPIgnore);
  EXPECT_EQ(settingsOODPIgnoreString, "IGNORE");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "IGNORE");

  settings.set_extrapolate_option("EXTRAPOLATE");
  auto settingsOODPExtrapolate = settings.get_extrapolate_option();
  std::string settingsOODPExtrapolateString = stk::search::get_object_outside_domain_policy(settingsOODPExtrapolate);
  EXPECT_EQ(settingsOODPExtrapolateString, "EXTRAPOLATE");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "EXTRAPOLATE");

  settings.set_extrapolate_option("TRUNCATE");
  auto settingsOODPTruncate = settings.get_extrapolate_option();
  std::string settingsOODPTruncateString = stk::search::get_object_outside_domain_policy(settingsOODPTruncate);
  EXPECT_EQ(settingsOODPTruncateString, "TRUNCATE");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "TRUNCATE");

  settings.set_extrapolate_option("PROJECT");
  auto settingsOODPProject = settings.get_extrapolate_option();
  std::string settingsOODPProjectString = stk::search::get_object_outside_domain_policy(settingsOODPProject);
  EXPECT_EQ(settingsOODPProjectString, "PROJECT");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "PROJECT");

  settings.set_extrapolate_option("ABORT");
  auto settingsOODPAbort = settings.get_extrapolate_option();
  std::string settingsOODPAbortString = stk::search::get_object_outside_domain_policy(settingsOODPAbort);
  EXPECT_EQ(settingsOODPAbortString, "ABORT");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "ABORT");

  settings.set_extrapolate_option("BAD");
  auto settingsOODPBad = settings.get_extrapolate_option();
  std::string settingsOODPBadString = stk::search::get_object_outside_domain_policy(settingsOODPBad);
  EXPECT_EQ(settingsOODPBadString, "UNDEFINED");
  EXPECT_EQ(settings.get_extrapolate_option_string(), "UNDEFINED");

}

TEST(TransferMainSettings, recvType)
{
  stk::transfer_util::TransferMainSettings settings;

  settings.set_recv_type("NODE");
  auto settingsDefaultRecvType = settings.get_recv_type();
  auto settingsDefaultRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsDefaultRecvType, stk::transfer::spmd::RecvMeshType::NODE);
  EXPECT_EQ(settingsDefaultRecvTypeString, "NODE");

  settings.set_recv_type("EDGE_CENTROID");
  auto settingsEdgeCentRecvType = settings.get_recv_type();
  auto settingsEdgeCentRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsEdgeCentRecvType, stk::transfer::spmd::RecvMeshType::EDGE_CENTROID);
  EXPECT_EQ(settingsEdgeCentRecvTypeString, "EDGE_CENTROID");

  settings.set_recv_type("EDGE_GAUSS_POINT");
  auto settingsEdgeGPRecvType = settings.get_recv_type();
  auto settingsEdgeGPRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsEdgeGPRecvType, stk::transfer::spmd::RecvMeshType::EDGE_GAUSS_POINT);
  EXPECT_EQ(settingsEdgeGPRecvTypeString, "EDGE_GAUSS_POINT");

  settings.set_recv_type("FACE_CENTROID");
  auto settingsFaceCentRecvType = settings.get_recv_type();
  auto settingsFaceCentRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsFaceCentRecvType, stk::transfer::spmd::RecvMeshType::FACE_CENTROID);
  EXPECT_EQ(settingsFaceCentRecvTypeString, "FACE_CENTROID");

  settings.set_recv_type("FACE_GAUSS_POINT");
  auto settingsFaceGPRecvType = settings.get_recv_type();
  auto settingsFaceGPRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsFaceGPRecvType, stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT);
  EXPECT_EQ(settingsFaceGPRecvTypeString, "FACE_GAUSS_POINT");

  settings.set_recv_type("ELEMENT_CENTROID");
  auto settingsElemCentRecvType = settings.get_recv_type();
  auto settingsElemCentRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsElemCentRecvType, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID);
  EXPECT_EQ(settingsElemCentRecvTypeString, "ELEMENT_CENTROID");

  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  auto settingsElemGPRecvType = settings.get_recv_type();
  auto settingsElemGPRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsElemGPRecvType, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT);
  EXPECT_EQ(settingsElemGPRecvTypeString, "ELEMENT_GAUSS_POINT");

  settings.set_recv_type("FAKE_CENTROID");
  auto settingsInvalidRecvType = settings.get_recv_type();
  auto settingsInvalidRecvTypeString = settings.get_recv_type_string();
  EXPECT_EQ(settingsInvalidRecvType, stk::transfer::spmd::RecvMeshType::INVALID);
  EXPECT_EQ(settingsInvalidRecvTypeString, "FAKE_CENTROID");
}

TEST(TransferMainSettings, defaultRecvPartRank)
{
  stk::transfer_util::TransferMainSettings settings;

  settings.set_recv_type("NODE");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::ELEMENT_RANK);

  settings.set_recv_type("EDGE_CENTROID");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::EDGE_RANK);

  settings.set_recv_type("EDGE_GAUSS_POINT");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::EDGE_RANK);

  settings.set_recv_type("FACE_CENTROID");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::FACE_RANK);

  settings.set_recv_type("FACE_GAUSS_POINT");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::FACE_RANK);

  settings.set_recv_type("ELEMENT_CENTROID");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::ELEMENT_RANK);

  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::ELEMENT_RANK);

  settings.set_recv_type("FAKE_CENTROID");
  EXPECT_EQ(settings.get_default_recv_part_rank(), stk::topology::INVALID_RANK);
}

TEST(TransferMainSettings, transferType)
{
  stk::transfer_util::TransferMainSettings settings;

  auto settingsDefaultTransferType = settings.get_transfer_type();
  EXPECT_EQ(settingsDefaultTransferType, "INTERP");

  settings.set_transfer_type("COPY");
  EXPECT_EQ(settings.get_transfer_type(), "COPY");

  settings.set_transfer_type("PATCH");
  EXPECT_EQ(settings.get_transfer_type(), "PATCH");

  settings.set_transfer_type("INTERP");
  EXPECT_EQ(settings.get_transfer_type(), "INTERP");

  EXPECT_FALSE(settings.set_transfer_type("FAKE"));
}

TEST(TransferMainSettings, coordTransformZExpression)
{
  stk::transfer_util::TransferMainSettings settings;

  auto settingsDefaultExpr = settings.get_coord_transf_z_expr();
  EXPECT_EQ(settingsDefaultExpr, "z=0");

  settings.set_coord_transf_z_expr("z=x+y");
  EXPECT_EQ(settings.get_coord_transf_z_expr(), "z=x+y");

  settings.set_coord_transf_z_expr("1.5");
  EXPECT_EQ(settings.get_coord_transf_z_expr(), "1.5");
}

}
