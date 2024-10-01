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

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/NgpFieldBLAS.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include "NgpFieldTestUtils.hpp"
#include <string>
#include <cstdlib>

namespace ngp_field_blas_test {

class NgpFieldBLAS : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldBLAS()
  {
    setup_three_fields_five_hex_three_block_mesh();
  }

  void setup_three_fields_five_hex_three_block_mesh()
  {
    const unsigned numStates = 1;
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stkField1 = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "variableLengthField1", numStates);
    stkField2 = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "variableLengthField2", numStates);
    stkField3 = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "variableLengthField3", numStates);
    stkField4 = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "variableLengthField4", numStates);
    stkField5 = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "variableLengthField5", numStates);

    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);

    const std::vector<double> init1(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField1, block1, numComponent1, init1.data());

    const std::vector<double> init2(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField1, block2, numComponent2, init2.data());

    const std::vector<double> init3(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField2, block1, numComponent1, init3.data());

    const std::vector<double> init4(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField2, block2, numComponent2, init4.data());

    const std::vector<double> init5(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField3, block1, numComponent1, init5.data());

    const std::vector<double> init6(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField3, block2, numComponent2, init6.data());


    const std::vector<double> init7(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField4, block1, numComponent1, init7.data());

    const std::vector<double> init8(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField4, block2, numComponent2, init8.data());

    const std::vector<double> init9(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField5, block1, numComponent1, init9.data());

    const std::vector<double> init10(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField5, block2, numComponent2, init10.data());

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                                 "0,3,HEX_8,9,13,14,15,16,17,18,19,block_2\n"
                                 "0,4,HEX_8,9,20,21,22,23,24,25,26,block_2\n"
                                 "0,5,HEX_8,9,27,28,29,30,31,32,33,block_3";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    EXPECT_FALSE(stkField1->need_sync_to_host());
  }

  const int numComponent1 = 8;
  const int numComponent2 = 3;
  stk::mesh::Field<double>* stkField1 = nullptr;
  stk::mesh::Field<double>* stkField2 = nullptr;
  stk::mesh::Field<double>* stkField3 = nullptr;
  stk::mesh::Field<double>* stkField4 = nullptr;
  stk::mesh::Field<double>* stkField5 = nullptr;
};

class NgpFieldBLASNode : public stk::unit_test_util::MeshFixture
{
public:
  NgpFieldBLASNode()
  {
    setup_three_fields_three_hex_three_block_mesh();
  }

  void setup_three_fields_three_hex_three_block_mesh()
  {
    const unsigned numStates = 1;
    const unsigned bucketCapacity = 2;
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);

    stkField1 = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "variableLengthNodeField1", numStates);
    stkField2 = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "variableLengthNodeField2", numStates);
    stkField3 = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "variableLengthNodeField3", numStates);


    stk::mesh::Part& block1 = get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    get_meta().declare_part_with_topology("block_3", stk::topology::HEX_8);

    const std::vector<double> init1(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField1, block1, numComponent1, init1.data());

    const std::vector<double> init2(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField1, block2, numComponent2, init2.data());

    const std::vector<double> init3(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField2, block1, numComponent1, init3.data());

    const std::vector<double> init4(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField2, block2, numComponent2, init4.data());

    const std::vector<double> init5(numComponent1, -1);
    stk::mesh::put_field_on_mesh(*stkField3, block1, numComponent1, init5.data());

    const std::vector<double> init6(numComponent2, -2);
    stk::mesh::put_field_on_mesh(*stkField3, block2, numComponent2, init6.data());

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2\n"
                                 "0,3,HEX_8,17,18,19,20,21,22,23,24,block_3\n"
                                 "|coordinates:"
                                 // elem 1
                                 "0,0,0, 1,0,0, 1,1,0, 0,1,0, "
                                 "0,0,1, 1,0,1, 1,1,1, 0,1,1, "
                                 // elem 2
                                 "2,0,0, 3,0,0, 3,1,0, 2,1,0, "
                                 "2,0,1, 3,0,1, 3,1,1, 2,1,1, "
                                 // elem 3
                                 "4,0,0, 5,0,0, 5,1,0, 4,1,0, "
                                 "4,0,1, 5,0,1, 5,1,1, 4,1,1";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    EXPECT_FALSE(stkField1->need_sync_to_host());
  }

  const int numComponent1 = 8;
  const int numComponent2 = 3;
  stk::mesh::Field<double>* stkField1 = nullptr;
  stk::mesh::Field<double>* stkField2 = nullptr;
  stk::mesh::Field<double>* stkField3 = nullptr;
};

namespace {
std::vector<double> func1(const double* coords)
{
  double v = coords[0] + 2*coords[1] + 3*coords[2];
  return {v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7};
}

std::vector<double> func2(const double* coords)
{
  double v = coords[0] + 2*coords[1] + 3*coords[2];
  return {v + 8, v + 9, v + 10, v + 11, v + 12, v + 13, v + 14, v + 15};
}

}

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLAS, field_fill_device)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;

  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());

  EXPECT_TRUE(stkField1->need_sync_to_host());

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, selector, myConstantValue);

  const double initialValue = -1;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, initialValue);

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_device_multiple)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;

  std::vector<const stk::mesh::FieldBase*> allFields = {stkField1, stkField2, stkField3, stkField4, stkField5};
  stk::mesh::field_fill(myConstantValue, allFields, stk::ngp::ExecSpace());

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  for (const stk::mesh::FieldBase* field : allFields)
  {
    EXPECT_TRUE(field->need_sync_to_host());

    stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(*field);

    stk::mesh::Selector selector(*field);
    ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField, selector, myConstantValue);

    const double initialValue = -1;
    ngp_field_test_utils::check_field_data_on_host(get_bulk(), *field, selector, initialValue);

    field->sync_to_host();
    ngp_field_test_utils::check_field_data_on_host(get_bulk(), *field, selector, myConstantValue);
  }
}

TEST_F(NgpFieldBLAS, field_fill_selector_device)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, selector, stk::ngp::ExecSpace());

  EXPECT_TRUE(stkField1->need_sync_to_host());

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, selector, myConstantValue);

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_component_selector_device)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  constexpr int component = 1;
  constexpr double myConstantComponentValue = 15.5;
  constexpr double myConstantValue = 55.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, selector, stk::ngp::ExecSpace());
  stk::mesh::field_fill(myConstantComponentValue, *stkField1, component, selector, stk::ngp::ExecSpace());

  EXPECT_TRUE(stkField1->need_sync_to_host());

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<double>& ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, selector, myConstantValue, component, myConstantComponentValue);

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue, component, myConstantComponentValue);
}

TEST_F(NgpFieldBLAS, field_fill_device_component_multiple)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const int component = 1;
  const double myConstantComponentValue = 15.5;
  const double myConstantValue = 55.5;

  std::vector<const stk::mesh::FieldBase*> allFields = {stkField1, stkField2, stkField3, stkField4, stkField5};
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, allFields, selector, stk::ngp::ExecSpace());
  stk::mesh::field_fill(myConstantComponentValue, allFields, component, selector, stk::ngp::ExecSpace());

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  for (const stk::mesh::FieldBase* field : allFields)
  {
    EXPECT_TRUE(field->need_sync_to_host());

    stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(*field);

    ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField, selector, myConstantValue, component, myConstantComponentValue);

    field->sync_to_host();

    ngp_field_test_utils::check_field_data_on_host(get_bulk(), *field, selector, myConstantValue, component, myConstantComponentValue);
  }
}


TEST_F(NgpFieldBLAS, field_fill_host_ngp)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_device());
#else
  EXPECT_TRUE(stkField1->need_sync_to_host());
  stkField1->sync_to_host();
#endif

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

#else /* not-defined STK_USE_DEVICE_MESH */

TEST_F(NgpFieldBLAS, field_fill_host)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());

  stk::mesh::Selector selector(*stkField1);
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_selector_host)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());

  stk::mesh::Selector all(*stkField1);
  stk::mesh::Selector notBlock2 = all - selector;
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_fill_device_with_host_build)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const double myConstantValue = 55.5;
  constexpr bool MarkModOnDevice = true;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{}, MarkModOnDevice);

  EXPECT_TRUE(stkField1->need_sync_to_host());
  stkField1->sync_to_host();

  stk::mesh::Selector selector(*stkField1);
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myConstantValue);
}

#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLAS, field_copy_device)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::ExecSpace());

  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stk::mesh::Selector selector(*stkField1);

  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
  auto ngpField2 = stk::mesh::get_updated_ngp_field<double>(*stkField2);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField2, selector, myConstantValue);

  stkField1->sync_to_host();
  stkField2->sync_to_host();

  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_selector_device)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::ExecSpace());
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, selector, stk::ngp::ExecSpace());

  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stk::mesh::Selector notBlock2(*stkField1);
  notBlock2 -= selector;
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto ngpField1 = stk::mesh::get_updated_ngp_field<double>(*stkField1);
  auto ngpField2 = stk::mesh::get_updated_ngp_field<double>(*stkField2);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField2, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField1, selector, myBlock2Value);
  ngp_field_test_utils::check_field_data_on_device(ngpMesh, ngpField2, selector, myBlock2Value);

  stkField1->sync_to_host();
  stkField2->sync_to_host();

  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myBlock2Value);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_copy_host_ngp)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_device());
  EXPECT_TRUE(stkField2->need_sync_to_device());
#else
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
  stkField1->sync_to_host();
  stkField2->sync_to_host();
#endif


  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

#else /* not-defined STK_USE_DEVICE_MESH */

TEST_F(NgpFieldBLAS, field_copy_host)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stk::mesh::Selector selector(*stkField1);
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_selector_host)
{
  if (get_parallel_size() != 1) { GTEST_SKIP(); }

  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, selector, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stk::mesh::Selector all(*stkField1);
  stk::mesh::Selector notBlock2 = all - selector;
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField2, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myBlock2Value);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField2, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_copy_device_with_host_build)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  const double myConstantValue = 55.5;
  constexpr bool MarkModOnDevice = true;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{}, MarkModOnDevice);
  EXPECT_TRUE(stkField1->need_sync_to_host());

  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{}, MarkModOnDevice);
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stk::mesh::Selector selector(*stkField1);
  auto ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(ngpMesh, *stkField2, selector, myConstantValue);
}

#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_axpbyz_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpbyz(get_bulk(), alpha, *stkField1, beta, *stkField2, *stkField3, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpbyz_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpbyz(get_bulk(), alpha, *stkField1, beta, *stkField2, *stkField3, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField3->need_sync_to_device());
#else
  EXPECT_TRUE(stkField3->need_sync_to_host());
#endif

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_axpbyz_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpbyz(get_bulk(), alpha, *stkField1, beta, *stkField2, *stkField3, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpbyz_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpbyz(get_bulk(), alpha, *stkField1, beta, *stkField2, *stkField3, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_axpby_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpby_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField2->need_sync_to_device());
#else
  EXPECT_TRUE(stkField2->need_sync_to_host());
#endif

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_axpby_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpby_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = alpha * vals1[i] + beta * vals2[i];
    }

    return result;
  };


  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}


#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_axpy_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };


  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpy_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };


  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField2->need_sync_to_device());
#else
  EXPECT_TRUE(stkField2->need_sync_to_host());
#endif

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_axpy_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpy_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_product_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);


  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = vals1[i] * vals2[i];
    }

    return result;
  };


  stk::mesh::field_product(get_bulk(), *stkField1, *stkField2, *stkField3, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_product_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = vals1[i] * vals2[i];
    }

    return result;
  };


  stk::mesh::field_product(get_bulk(), *stkField1, *stkField2, *stkField3, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField3->need_sync_to_device());
#else
  EXPECT_TRUE(stkField3->need_sync_to_host());
#endif

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_product_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = vals1[i] * vals2[i];
    }

    return result;
  };

  stk::mesh::field_product(get_bulk(), *stkField1, *stkField2, *stkField3, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_product_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);


  auto f_expected = [&](const double* coords)
  {
    std::vector<double> vals1 = func1(coords);
    std::vector<double> vals2 = func2(coords);
    std::vector<double> result(vals1.size());
    for (size_t i=0; i < vals1.size(); ++i)
    {
      result[i] = vals1[i] * vals2[i];
    }

    return result;
  };

  stk::mesh::field_product(get_bulk(), *stkField1, *stkField2, *stkField3, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField3->need_sync_to_host());

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLAS, field_axpby)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::field_fill(3.0, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_fill(10.0, *stkField2, stk::ngp::ExecSpace());

  double alpha = 2.0;
  double beta = 5.0;
  stk::mesh::Selector selectRule(*stkField1);

  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, selectRule, stk::ngp::ExecSpace());

  stkField2->sync_to_host();
  stk::mesh::Selector selector(*stkField2);
  constexpr double expectedValue = 56.0;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, expectedValue);
}
#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_scale_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };


  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField1->need_sync_to_host());

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_scale_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };


  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_device());
#else
  EXPECT_TRUE(stkField1->need_sync_to_host());
#endif

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_scale_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField1->need_sync_to_host());

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_scale_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](const double* coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField1->need_sync_to_host());

  stkField1->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

#endif

#ifdef STK_USE_DEVICE_MESH

TEST_F(NgpFieldBLASNode, field_swap_device)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [&](const double* coords)
  {
    return func2(coords);
  };

  auto f2_expected = [&](const double* coords)
  {
    return func1(coords);
  };

  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());

  stkField1->sync_to_host();
  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);
}

TEST_F(NgpFieldBLASNode, field_swap_host)
{
  if (get_parallel_size() > 2) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [&](const double* coords)
  {
    return func2(coords);
  };

  auto f2_expected = [&](const double* coords)
  {
    return func1(coords);
  };


  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::HostExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_device());
  EXPECT_TRUE(stkField2->need_sync_to_device());
#else
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#endif

  stkField1->sync_to_host();
  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);
}

#else

TEST_F(NgpFieldBLASNode, field_swap_exec_space)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [&](const double* coords)
  {
    return func2(coords);
  };

  auto f2_expected = [&](const double* coords)
  {
    return func1(coords);
  };

  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField1->sync_to_host();
  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);
}

TEST_F(NgpFieldBLASNode, field_swap_no_selector)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [&](const double* coords)
  {
    return func2(coords);
  };

  auto f2_expected = [&](const double* coords)
  {
    return func1(coords);
  };

  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, stk::ngp::ExecSpace());
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());

  stkField1->sync_to_host();
  stkField2->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);

}

#endif


}

