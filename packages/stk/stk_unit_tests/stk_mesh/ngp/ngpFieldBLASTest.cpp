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
std::vector<double> func1(stk::mesh::EntityValues<const double> coords)
{
  double v = coords(0_comp) + 2*coords(1_comp) + 3*coords(2_comp);
  return {v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7};
}

std::vector<double> func2(stk::mesh::EntityValues<const double> coords)
{
  double v = coords(0_comp) + 2*coords(1_comp) + 3*coords(2_comp);
  return {v + 8, v + 9, v + 10, v + 11, v + 12, v + 13, v + 14, v + 15};
}

}

TEST_F(NgpFieldBLAS, field_fill_device)
{
  const double myConstantValue = 55.5;

  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
#endif

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_device_multiple)
{
  const double myConstantValue = 55.5;

  std::vector<const stk::mesh::FieldBase*> allFields = {stkField1, stkField2, stkField3, stkField4, stkField5};
  stk::mesh::field_fill(myConstantValue, allFields, stk::ngp::ExecSpace());

  for (const stk::mesh::FieldBase* field : allFields)
  {
#ifdef STK_ENABLE_GPU
    EXPECT_TRUE(field->need_sync_to_host());
#else
    EXPECT_FALSE(field->need_sync_to_host());
#endif

    stk::mesh::Selector selector(*field);
    ngp_field_test_utils::check_field_data_on_host(get_bulk(), *field, selector, myConstantValue);
  }
}

TEST_F(NgpFieldBLAS, field_fill_selector_device)
{
  const double myConstantValue = 55.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, selector, stk::ngp::ExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_component_selector_device)
{
  constexpr int component = 1;
  constexpr double myConstantComponentValue = 15.5;
  constexpr double myConstantValue = 55.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, selector, stk::ngp::ExecSpace());
  stk::mesh::field_fill(myConstantComponentValue, *stkField1, component, selector, stk::ngp::ExecSpace());

#ifdef STK_ENABLE_GPU
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue, component, myConstantComponentValue);
}

TEST_F(NgpFieldBLAS, field_fill_device_component_multiple)
{
  const int component = 1;
  const double myConstantComponentValue = 15.5;
  const double myConstantValue = 55.5;

  std::vector<const stk::mesh::FieldBase*> allFields = {stkField1, stkField2, stkField3, stkField4, stkField5};
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, allFields, selector, stk::ngp::ExecSpace());
  stk::mesh::field_fill(myConstantComponentValue, allFields, component, selector, stk::ngp::ExecSpace());

  for (const stk::mesh::FieldBase* field : allFields)
  {
#ifdef STK_ENABLE_GPU
    EXPECT_TRUE(field->need_sync_to_host());
#else
    EXPECT_FALSE(field->need_sync_to_host());
#endif

    ngp_field_test_utils::check_field_data_on_host(get_bulk(), *field, selector, myConstantValue, component, myConstantComponentValue);
  }
}

TEST_F(NgpFieldBLAS, field_fill_host_ngp)
{
  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});

  EXPECT_TRUE(stkField1->need_sync_to_device());

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_host)
{
  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_fill_selector_host)
{
  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::HostExecSpace{});

  EXPECT_FALSE(stkField1->need_sync_to_host());

  stk::mesh::Selector all(*stkField1);
  stk::mesh::Selector notBlock2 = all - selector;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_fill_device_with_host_build)
{
  const double myConstantValue = 55.5;
  constexpr bool MarkModOnDevice = true;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{}, MarkModOnDevice);

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_device)
{
  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::ExecSpace());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_TRUE(stkField1->need_sync_to_device());
  EXPECT_TRUE(stkField2->need_sync_to_device());
#endif

  stk::mesh::Selector selector(*stkField1);

  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_selector_device)
{
  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::ExecSpace());
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::ExecSpace());
  stk::mesh::field_copy(*stkField1, *stkField2, selector, stk::ngp::ExecSpace());

#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_TRUE(stkField1->need_sync_to_device());
  EXPECT_TRUE(stkField2->need_sync_to_device());
#endif

  stk::mesh::Selector notBlock2(*stkField1);
  notBlock2 -= selector;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myBlock2Value);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_copy_host_ngp)
{
  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_host)
{
  const double myConstantValue = 55.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLAS, field_copy_selector_host)
{
  const double myConstantValue = 55.5;
  const double myBlock2Value = 95.5;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{});
  stk::mesh::Part& block2 = *get_meta().get_part("block_2");
  stk::mesh::Selector selector = block2;
  stk::mesh::field_fill(myBlock2Value, *stkField1, selector, stk::ngp::HostExecSpace{});
  stk::mesh::field_copy(*stkField1, *stkField2, selector, stk::ngp::HostExecSpace{});

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  stk::mesh::Selector all(*stkField1);
  stk::mesh::Selector notBlock2 = all - selector;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, notBlock2, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myBlock2Value);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myBlock2Value);
}

TEST_F(NgpFieldBLAS, field_copy_device_with_host_build)
{
  const double myConstantValue = 55.5;
  constexpr bool MarkModOnDevice = true;
  stk::mesh::field_fill(myConstantValue, *stkField1, stk::ngp::HostExecSpace{}, MarkModOnDevice);

  stk::mesh::field_copy(*stkField1, *stkField2, stk::ngp::HostExecSpace{}, MarkModOnDevice);

  stk::mesh::Selector selector(*stkField1);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField1, selector, myConstantValue);
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, myConstantValue);
}

TEST_F(NgpFieldBLASNode, field_axpbyz_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_FALSE(stkField3->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpbyz_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_TRUE(stkField3->need_sync_to_device());
#endif

  stkField3->sync_to_host();
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpbyz_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_FALSE(stkField3->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpby_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpby_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_TRUE(stkField2->need_sync_to_device());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpby_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  double beta = 5.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpy_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> x = func1(coords);
    std::vector<double> y = func2(coords);
    for (size_t i=0; i < y.size(); ++i)
    {
      y[i] += alpha * x[i];
    }

    return y;
  };

  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpy_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> x = func1(coords);
    std::vector<double> y = func2(coords);
    for (size_t i=0; i < y.size(); ++i)
    {
      y[i] += alpha * x[i];
    }

    return y;
  };

  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, selectRule, stk::ngp::HostExecSpace());

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_TRUE(stkField2->need_sync_to_device());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_axpy_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> x = func1(coords);
    std::vector<double> y = func2(coords);
    for (size_t i=0; i < y.size(); ++i)
    {
      y[i] += alpha * x[i];
    }

    return y;
  };

  stk::mesh::field_axpy(get_bulk(), alpha, *stkField1, *stkField2, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_product_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);


  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_FALSE(stkField3->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_product_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_TRUE(stkField3->need_sync_to_device());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_product_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);


  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
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
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField3->need_sync_to_host());
#else
  EXPECT_FALSE(stkField3->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField3, selectRule, {stkField1, stkField2}, f_expected);
}

TEST_F(NgpFieldBLAS, field_axpby)
{
  stk::mesh::field_fill(3.0, *stkField1, stk::ngp::ExecSpace());
  stk::mesh::field_fill(10.0, *stkField2, stk::ngp::ExecSpace());

  double alpha = 2.0;
  double beta = 5.0;
  stk::mesh::Selector selectRule(*stkField1);

  stk::mesh::field_axpby(get_bulk(), alpha, *stkField1, beta, *stkField2, selectRule, stk::ngp::ExecSpace());

  stk::mesh::Selector selector(*stkField2);
  constexpr double expectedValue = 56.0;
  ngp_field_test_utils::check_field_data_on_host(get_bulk(), *stkField2, selector, expectedValue);
}

TEST_F(NgpFieldBLASNode, field_scale_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, selectRule, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_scale_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, selectRule, stk::ngp::HostExecSpace());

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_TRUE(stkField1->need_sync_to_device());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}

TEST_F(NgpFieldBLASNode, field_scale_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double alpha = 2.0;
  auto f_expected = [&](stk::mesh::EntityValues<const double> coords)
  {
    std::vector<double> result = func1(coords);
    for (size_t i=0; i < result.size(); ++i)
    {
      result[i] = alpha * result[i];
    }

    return result;
  };

  stk::mesh::field_scale(get_bulk(), alpha, *stkField1, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {}, f_expected);
}


TEST_F(NgpFieldBLASNode, field_swap_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func2(coords);
  };

  auto f2_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func1(coords);
  };

  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);
}

TEST_F(NgpFieldBLASNode, field_swap_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func2(coords);
  };

  auto f2_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func1(coords);
  };


  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::HostExecSpace());

#if defined(STK_USE_DEVICE_MESH) && !defined(STK_ENABLE_GPU)
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_TRUE(stkField1->need_sync_to_device());
  EXPECT_TRUE(stkField2->need_sync_to_device());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);
}

TEST_F(NgpFieldBLASNode, field_swap_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField2, selectRule, &func2);

  auto f1_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func2(coords);
  };

  auto f2_expected = [](stk::mesh::EntityValues<const double> coords)
  {
    return func1(coords);
  };

  stk::mesh::field_swap(get_bulk(), *stkField1, *stkField2, stk::ngp::ExecSpace());
#ifdef STK_USE_DEVICE_MESH
  EXPECT_TRUE(stkField1->need_sync_to_host());
  EXPECT_TRUE(stkField2->need_sync_to_host());
#else
  EXPECT_FALSE(stkField1->need_sync_to_host());
  EXPECT_FALSE(stkField2->need_sync_to_host());
#endif

  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField1, selectRule, {stkField2}, f1_expected);
  ngp_field_test_utils::check_field_data_on_host_func(get_bulk(), *stkField2, selectRule, {stkField1}, f2_expected);

}

TEST_F(NgpFieldBLASNode, field_amax_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  stkField1->sync_to_device();

  double hostAMaxVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMaxVal = std::max(hostAMaxVal, myVal(i));
    }
  }

  double globalHostMax = hostAMaxVal;
  stk::all_reduce_max(get_bulk().parallel(), &hostAMaxVal, &globalHostMax, 1u);

  double devAmaxVal = 0.0;
  stk::mesh::field_amax(devAmaxVal, *stkField1, selectRule, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostMax, devAmaxVal);
}

TEST_F(NgpFieldBLASNode, field_amax_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostAMaxVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMaxVal = std::max(hostAMaxVal, myVal(i));
    }
  }

  double globalHostMax = hostAMaxVal;
  stk::all_reduce_max(get_bulk().parallel(), &hostAMaxVal, &globalHostMax, 1u);

  double devAmaxVal = 0.0;
  stk::mesh::field_amax(devAmaxVal, *stkField1, selectRule, stk::ngp::HostExecSpace());

  EXPECT_EQ(globalHostMax, devAmaxVal);
}

TEST_F(NgpFieldBLASNode, field_amax_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostAMaxVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMaxVal = std::max(hostAMaxVal, myVal(i));
    }
  }

  double globalHostMax = hostAMaxVal;
  stk::all_reduce_max(get_bulk().parallel(), &hostAMaxVal, &globalHostMax, 1u);

  double devAmaxVal = 0.0;
  stk::mesh::field_amax(devAmaxVal, *stkField1, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostMax, devAmaxVal);
}

TEST_F(NgpFieldBLASNode, field_amin_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  stkField1->sync_to_device();

  double hostAMinVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMinVal = std::min(hostAMinVal, myVal(i));
    }
  }

  double globalHostMin = hostAMinVal;
  stk::all_reduce_min(get_bulk().parallel(), &hostAMinVal, &globalHostMin, 1u);

  double devAminVal = 0.0;
  stk::mesh::field_amin(devAminVal, *stkField1, selectRule, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostMin, devAminVal);
}

TEST_F(NgpFieldBLASNode, field_amin_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostAMinVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMinVal = std::min(hostAMinVal, myVal(i));
    }
  }

  double globalHostMin = hostAMinVal;
  stk::all_reduce_min(get_bulk().parallel(), &hostAMinVal, &globalHostMin, 1u);

  double devAminVal = 0.0;
  stk::mesh::field_amin(devAminVal, *stkField1, selectRule, stk::ngp::HostExecSpace());

  EXPECT_EQ(globalHostMin, devAminVal);
}

TEST_F(NgpFieldBLASNode, field_amin_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostAMinVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostAMinVal = std::min(hostAMinVal, myVal(i));
    }
  }

  double globalHostMin = hostAMinVal;
  stk::all_reduce_min(get_bulk().parallel(), &hostAMinVal, &globalHostMin, 1u);

  double devAminVal = 0.0;
  stk::mesh::field_amin(devAminVal, *stkField1, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostMin, devAminVal);
}

TEST_F(NgpFieldBLASNode, field_asum_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  stkField1->sync_to_device();

  double hostASumVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostASumVal += myVal(i);
    }
  }

  double globalHostSum = hostASumVal;
  stk::all_reduce_sum(get_bulk().parallel(), &hostASumVal, &globalHostSum, 1u);

  double devAsumVal = 0.0;
  stk::mesh::field_asum(devAsumVal, *stkField1, selectRule, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostSum, devAsumVal);
}

TEST_F(NgpFieldBLASNode, field_asum_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostASumVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostASumVal += myVal(i);
    }
  }

  double globalHostSum = hostASumVal;
  stk::all_reduce_sum(get_bulk().parallel(), &hostASumVal, &globalHostSum, 1u);

  double devAsumVal = 0.0;
  stk::mesh::field_asum(devAsumVal, *stkField1, selectRule, stk::ngp::HostExecSpace());

  EXPECT_EQ(globalHostSum, devAsumVal);
}

TEST_F(NgpFieldBLASNode, field_asum_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostASumVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostASumVal += myVal(i);
    }
  }

  double globalHostSum = hostASumVal;
  stk::all_reduce_sum(get_bulk().parallel(), &hostASumVal, &globalHostSum, 1u);

  double devAsumVal = 0.0;
  stk::mesh::field_asum(devAsumVal, *stkField1, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostSum, devAsumVal);
}

TEST_F(NgpFieldBLASNode, field_nrm2_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  stkField1->sync_to_device();

  double hostNrm2Val = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostNrm2Val += myVal(i) * myVal(i);
    }
  }

  double globalHostNrm2 = hostNrm2Val;
  stk::all_reduce_sum(get_bulk().parallel(), &hostNrm2Val, &globalHostNrm2, 1u);
  globalHostNrm2 = std::sqrt(globalHostNrm2);

  double devNrm2Val = 0.0;
  stk::mesh::field_nrm2(devNrm2Val, *stkField1, selectRule, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostNrm2, devNrm2Val);
}

TEST_F(NgpFieldBLASNode, field_nrm2_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostNrm2Val = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostNrm2Val += myVal(i) * myVal(i);
    }
  }

  double globalHostNrm2 = hostNrm2Val;
  stk::all_reduce_sum(get_bulk().parallel(), &hostNrm2Val, &globalHostNrm2, 1u);
  globalHostNrm2 = std::sqrt(globalHostNrm2);

  double devNrm2Val = 0.0;
  stk::mesh::field_nrm2(devNrm2Val, *stkField1, selectRule, stk::ngp::HostExecSpace());

  EXPECT_EQ(globalHostNrm2, devNrm2Val);
}

TEST_F(NgpFieldBLASNode, field_nrm2_no_selector)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);

  double hostNrm2Val = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  for (auto& n : nodes) {
    auto myVal = stkField1Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal.components()) {
      hostNrm2Val += myVal(i) * myVal(i);
    }
  }

  double globalHostNrm2 = hostNrm2Val;
  stk::all_reduce_sum(get_bulk().parallel(), &hostNrm2Val, &globalHostNrm2, 1u);
  globalHostNrm2 = std::sqrt(globalHostNrm2);

  double devNrm2Val = 0.0;
  stk::mesh::field_nrm2(devNrm2Val, *stkField1, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostNrm2, devNrm2Val);
}

TEST_F(NgpFieldBLASNode, field_dot_device)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func2);

  double hostDotVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  auto stkField2Data = stkField2->data();
  for (auto& n : nodes) {
    auto myVal1 = stkField1Data.entity_values(n);
    auto myVal2 = stkField2Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal1.components()) {
      hostDotVal += myVal1(i) * myVal2(i);
    }
  }

  double globalHostSum = 0.0;
  stk::all_reduce_sum(get_bulk().parallel(), &hostDotVal, &globalHostSum, 1u);

  double devDotVal = 0.0;
  stk::mesh::field_dot(devDotVal, get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::ExecSpace());

  EXPECT_EQ(globalHostSum, devDotVal);
}

TEST_F(NgpFieldBLASNode, field_dot_host)
{
  stk::mesh::Selector selectRule(*stkField1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func1);
  ngp_field_test_utils::set_field_data_on_host(get_bulk(), *stkField1, selectRule, &func2);

  double hostDotVal = 0.0;

  stk::mesh::EntityVector nodes;
  get_bulk().get_entities(stk::topology::NODE_RANK, selectRule, nodes);
  auto stkField1Data = stkField1->data();
  auto stkField2Data = stkField2->data();
  for (auto& n : nodes) {
    auto myVal1 = stkField1Data.entity_values(n);
    auto myVal2 = stkField2Data.entity_values(n);
    for (stk::mesh::ComponentIdx i : myVal1.components()) {
      hostDotVal += myVal1(i) * myVal2(i);
    }
  }

  double globalHostSum = hostDotVal;
  stk::all_reduce_sum(get_bulk().parallel(), &hostDotVal, &globalHostSum, 1u);

  double devDotVal = 0.0;
  stk::mesh::field_dot(devDotVal, get_bulk(), *stkField1, *stkField2, selectRule, stk::ngp::HostExecSpace());

  EXPECT_EQ(globalHostSum, devDotVal);
}

}

