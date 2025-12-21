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
}

TEST(TransferMainSettings, filenames)
{
  stk::transfer_util::TransferMainSettings settings;
 
  EXPECT_EQ(settings.get_sendMesh_filename(), "");
  EXPECT_EQ(settings.get_recvMesh_filename(), "");
  
  settings.set_sendMesh_filename("source.g");

  EXPECT_EQ(settings.get_sendMesh_filename(), "source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "");

  settings.set_recvMesh_filename("target.g");

  EXPECT_EQ(settings.get_sendMesh_filename(), "source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "target.g");

  settings.set_sendMesh_filename("changed_source.g");
  settings.set_recvMesh_filename("changed_target.g");

  EXPECT_EQ(settings.get_sendMesh_filename(), "changed_source.g");
  EXPECT_EQ(settings.get_recvMesh_filename(), "changed_target.g");
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

}
