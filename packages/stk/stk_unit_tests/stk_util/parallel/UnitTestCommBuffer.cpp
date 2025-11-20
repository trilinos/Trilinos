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

#include "gtest/gtest.h"
#include <stk_util/stk_config.h>
#include <stk_util/parallel/ParallelComm.hpp>

#if defined(STK_HAS_MPI) && defined(STK_HAVE_STKIO)
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_mesh/base/MetaData.hpp"
#endif

namespace {

class TestCommBuffer : public ::testing::Test
{
public:
  template<typename PACK_LAMBDA>
  void pack(const PACK_LAMBDA& pack_stuff)
  {
    pack_stuff(buf);
    allocate_buf();
    pack_stuff(buf);
  }

  void allocate_buf()
  {
    data.resize(buf.size());
    buf.set_buffer_ptrs(data.data(), data.data(), data.data()+data.size());
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack(const T& gold)
  {
    T val;
    buf.unpack(val);
    EXPECT_NEAR(gold, val, 1.e-6);
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack(const std::vector<T>& gold)
  {
    std::vector<T> val;
    buf.unpack(val);
    ASSERT_EQ(val.size(), gold.size());
    for(unsigned i=0; i<gold.size(); ++i) {
      EXPECT_NEAR(gold[i], val[i], 1.e-6);
    }
  }

  template<typename T>
  typename std::enable_if<!std::is_floating_point<T>::value>::type
  check_unpack(const T& gold)
  {
    T val;
    buf.unpack(val);
    EXPECT_EQ(gold, val);
  }

  template<typename T>
  void check_peek_and_unpack_pair(const T& gold)
  {
    T peekVal;
    buf.peek<T>(peekVal);
    EXPECT_EQ(gold, peekVal);
    T val;
    buf.unpack<T>(val);
    EXPECT_EQ(gold, val);
  }

  stk::CommBuffer buf;
  std::vector<unsigned char> data;
};

TEST_F(TestCommBuffer, pack_unpack_primitives)
{
  char goldc = 's';
  int goldn = 5;
  float goldf = 1.1;
  double goldd = 3.3;
  bool goldb = true;

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldc);
    buffer.pack(goldn);
    buffer.pack(goldf);
    buffer.pack(goldd);
    buffer.pack(goldb);
  });

  buf.reset();

  check_unpack(goldc);
  check_unpack(goldn);
  check_unpack(goldf);
  check_unpack(goldd);
  check_unpack(goldb);
}

TEST_F(TestCommBuffer, pack_unpack_pairs)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpic);
    buffer.pack(goldpdi);
  });

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, skip_pack_pairs_explicit_template_parameter)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  buf.skip<std::pair<int,char>>(1);
  buf.skip<std::pair<unsigned long,int>>(1);
  buf.skip<std::pair<double,int>>(1);
  
  allocate_buf();

  buf.pack<std::pair<int,char>>(goldpic);
  buf.pack<std::pair<unsigned long,int>>(goldpui);
  buf.pack<std::pair<double,int>>(goldpdi);

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpui);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pairs_explicit_template_parameter)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack<std::pair<int,char>>(goldpic);
    buffer.pack<std::pair<unsigned long,int>>(goldpui);
    buffer.pack<std::pair<double,int>>(goldpdi);
  });

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpui);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pairs_explicit_template_parameter_peek_unpack)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpic);
    buffer.pack(goldpui);
    buffer.pack(goldpdi);
  });

  buf.reset();

  check_peek_and_unpack_pair(goldpic);
  check_peek_and_unpack_pair(goldpui);
  check_peek_and_unpack_pair(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pair_string)
{
  std::pair<std::string,int> goldpsi("short string", 5);
  std::pair<double,std::string> goldpds(3.14, "some string");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpsi);
    buffer.pack(goldpds);
  });

  buf.reset();

  check_unpack(goldpsi);
  check_unpack(goldpds);
}

TEST_F(TestCommBuffer, pack_unpack_vector_primitives)
{
  std::vector<char> goldc = {'a', 'b', 'c'};
  std::vector<int> goldn = {5, 6, 7};
  std::vector<float> goldf = {1.1, 1.2, 1.3};
  std::vector<double> goldd = {3.3, 3.4, 3.5};
  std::vector<bool> goldb = {true, false, true};

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldc);
    buffer.pack(goldn);
    buffer.pack(goldf);
    buffer.pack(goldd);
    buffer.pack(goldb);
  });

  buf.reset();

  check_unpack(goldc);
  check_unpack(goldn);
  check_unpack(goldf);
  check_unpack(goldd);
  check_unpack(goldb);
}

TEST_F(TestCommBuffer, pack_unpack_string)
{
  std::string goldstr("stk is awesome");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldstr);
  });

  buf.reset();

  check_unpack(goldstr);
}

TEST_F(TestCommBuffer, pack_unpack_string_explicit_template_parameter)
{
  std::string goldstr("stk is awesome");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack<std::string>(goldstr);
  });

  buf.reset();

  std::string str;
  buf.unpack<std::string>(str);
  EXPECT_EQ(goldstr, str);
}

TEST_F(TestCommBuffer, pack_unpack_vector_string)
{
  std::vector<std::string> goldstr = {"stk", "is", "awesome"};

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldstr);
  });

  buf.reset();

  check_unpack(goldstr);
}

#if defined(STK_HAS_MPI) && defined(STK_HAVE_STKIO)

class TestCommBufferWithMesh : public TestCommBuffer, public stk::unit_test_util::MeshFixtureNoTest
{
public:
  template <typename T, stk::mesh::Layout DataLayout>
  std::pair<stk::mesh::Field<T, DataLayout>*, stk::mesh::Field<T, DataLayout>*>
  build_two_element_mesh_with_vector_fields(const std::vector<T>& initValsPack, const std::vector<T>& initValsUnpack)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    auto& fieldPack = get_meta().declare_field<T, DataLayout>(stk::topology::ELEM_RANK, "elemFieldPack");
    auto& fieldUnpack = get_meta().declare_field<T, DataLayout>(stk::topology::ELEM_RANK, "elemFieldUnpack");
    stk::mesh::put_field_on_mesh(fieldPack, get_meta().universal_part(), initValsPack.size(), initValsPack.data());
    stk::mesh::put_field_on_mesh(fieldUnpack, get_meta().universal_part(), initValsUnpack.size(), initValsUnpack.data());
    stk::io::fill_mesh("generated:1x1x2", get_bulk());
    return std::make_pair(&fieldPack, &fieldUnpack);
  }

  template<typename T, typename EntityValuesType>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack(const std::vector<T>& gold, EntityValuesType& entityValuesUnpack)
  {
    buf.unpack(entityValuesUnpack.pointer(), entityValuesUnpack.num_components(), entityValuesUnpack.component_stride());
    for (stk::mesh::ComponentIdx component : entityValuesUnpack.components()) {
      EXPECT_DOUBLE_EQ(gold[component], entityValuesUnpack(component));
    }
  }

  template<typename T, typename EntityBytesType>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack_bytes(const std::vector<T>& gold, EntityBytesType& entityBytesUnpack)
  {
    buf.unpack_bytes(entityBytesUnpack.pointer(), entityBytesUnpack.num_bytes(), entityBytesUnpack.bytes_per_scalar(),
                     entityBytesUnpack.scalar_byte_stride());

    const std::byte* goldBytes = reinterpret_cast<const std::byte*>(gold.data());
    for (stk::mesh::ByteIdx byte : entityBytesUnpack.bytes()) {
      EXPECT_EQ(goldBytes[byte], entityBytesUnpack(byte));
    }
  }

  template<typename T, typename EntityValuesType>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_peek(const std::vector<T>& gold, EntityValuesType& entityValuesUnpack)
  {
    buf.peek(entityValuesUnpack.pointer(), entityValuesUnpack.num_components(), entityValuesUnpack.component_stride());
    for (stk::mesh::ComponentIdx component : entityValuesUnpack.components()) {
      EXPECT_DOUBLE_EQ(gold[component], entityValuesUnpack(component));
    }
  }

  template<typename T, typename EntityBytesType>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_peek_bytes(const std::vector<T>& gold, EntityBytesType& entityBytesUnpack)
  {
    buf.peek_bytes(entityBytesUnpack.pointer(), entityBytesUnpack.num_bytes(), entityBytesUnpack.bytes_per_scalar(),
                   entityBytesUnpack.scalar_byte_stride());

    const std::byte* goldBytes = reinterpret_cast<const std::byte*>(gold.data());
    for (stk::mesh::ByteIdx byte : entityBytesUnpack.bytes()) {
      EXPECT_EQ(goldBytes[byte], entityBytesUnpack(byte));
    }
  }

};

TEST_F(TestCommBufferWithMesh, pack_layout_right_entity_values)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();

  std::vector<double> initValsPack {1.0, 2.0, 3.0};
  std::vector<double> initValsUnpack(3);
  auto [fieldPack, fieldUnpack] = build_two_element_mesh_with_vector_fields<double, stk::mesh::Layout::Right>(initValsPack,
                                                                                                              initValsUnpack);
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  auto fieldDataPack = fieldPack->data();
  auto entityValuesPack1 = fieldDataPack.entity_values(elem1);
  auto entityValuesPack2 = fieldDataPack.entity_values(elem2);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(entityValuesPack1.pointer(), entityValuesPack1.num_components(), entityValuesPack1.component_stride());
    buffer.pack(entityValuesPack2.pointer(), entityValuesPack2.num_components(), entityValuesPack2.component_stride());
  });

  buf.reset();

  auto fieldDataUnpack = fieldPack->data<stk::mesh::ReadWrite>();
  auto entityValuesUnpack1 = fieldDataUnpack.entity_values(elem1);
  auto entityValuesUnpack2 = fieldDataUnpack.entity_values(elem1);
  check_peek(initValsPack, entityValuesUnpack1);
  check_unpack(initValsPack, entityValuesUnpack1);
  check_peek(initValsPack, entityValuesUnpack2);
  check_unpack(initValsPack, entityValuesUnpack2);
}

TEST_F(TestCommBufferWithMesh, pack_layout_left_entity_values)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();

  std::vector<double> initValsPack {1.0, 2.0, 3.0};
  std::vector<double> initValsUnpack(3);
  auto [fieldPack, fieldUnpack] = build_two_element_mesh_with_vector_fields<double, stk::mesh::Layout::Left>(initValsPack,
                                                                                                             initValsUnpack);
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  auto fieldDataPack = fieldPack->data();
  auto entityValuesPack1 = fieldDataPack.entity_values(elem1);
  auto entityValuesPack2 = fieldDataPack.entity_values(elem2);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(entityValuesPack1.pointer(), entityValuesPack1.num_components(), entityValuesPack1.component_stride());
    buffer.pack(entityValuesPack2.pointer(), entityValuesPack2.num_components(), entityValuesPack2.component_stride());
  });

  buf.reset();

  auto fieldDataUnpack = fieldPack->data<stk::mesh::ReadWrite>();
  auto entityValuesUnpack1 = fieldDataUnpack.entity_values(elem1);
  auto entityValuesUnpack2 = fieldDataUnpack.entity_values(elem1);
  check_peek(initValsPack, entityValuesUnpack1);
  check_unpack(initValsPack, entityValuesUnpack1);
  check_peek(initValsPack, entityValuesUnpack2);
  check_unpack(initValsPack, entityValuesUnpack2);
}

TEST_F(TestCommBufferWithMesh, pack_layout_right_entity_bytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();

  std::vector<double> initValsPack {1.0, 2.0, 3.0};
  std::vector<double> initValsUnpack(3);
  auto [fieldPack, fieldUnpack] = build_two_element_mesh_with_vector_fields<double, stk::mesh::Layout::Right>(initValsPack,
                                                                                                              initValsUnpack);
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  auto fieldBytesPack = fieldPack->data_bytes<const std::byte>();
  auto entityBytesPack1 = fieldBytesPack.entity_bytes<stk::mesh::Layout::Right>(elem1);
  auto entityBytesPack2 = fieldBytesPack.entity_bytes<stk::mesh::Layout::Right>(elem2);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack_bytes(entityBytesPack1.pointer(), entityBytesPack1.num_bytes(), entityBytesPack1.bytes_per_scalar(),
                      entityBytesPack1.scalar_byte_stride());
    buffer.pack_bytes(entityBytesPack2.pointer(), entityBytesPack2.num_bytes(), entityBytesPack2.bytes_per_scalar(),
                      entityBytesPack2.scalar_byte_stride());
  });

  buf.reset();

  auto fieldBytesUnpack = fieldPack->data_bytes<std::byte>();
  auto entityBytesUnpack1 = fieldBytesUnpack.entity_bytes(elem1);
  auto entityBytesUnpack2 = fieldBytesUnpack.entity_bytes(elem2);
  check_peek_bytes(initValsPack, entityBytesUnpack1);
  check_unpack_bytes(initValsPack, entityBytesUnpack1);
  check_peek_bytes(initValsPack, entityBytesUnpack2);
  check_unpack_bytes(initValsPack, entityBytesUnpack2);
}

TEST_F(TestCommBufferWithMesh, pack_layout_left_entity_bytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();

  std::vector<double> initValsPack {1.0, 2.0, 3.0};
  std::vector<double> initValsUnpack(3);
  auto [fieldPack, fieldUnpack] = build_two_element_mesh_with_vector_fields<double, stk::mesh::Layout::Left>(initValsPack,
                                                                                                             initValsUnpack);
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  auto fieldBytesPack = fieldPack->data_bytes<const std::byte>();
  auto entityBytesPack1 = fieldBytesPack.entity_bytes<stk::mesh::Layout::Left>(elem1);
  auto entityBytesPack2 = fieldBytesPack.entity_bytes<stk::mesh::Layout::Left>(elem2);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack_bytes(entityBytesPack1.pointer(), entityBytesPack1.num_bytes(), entityBytesPack1.bytes_per_scalar(),
                      entityBytesPack1.scalar_byte_stride());
    buffer.pack_bytes(entityBytesPack2.pointer(), entityBytesPack2.num_bytes(), entityBytesPack2.bytes_per_scalar(),
                      entityBytesPack2.scalar_byte_stride());
  });

  buf.reset();

  auto fieldBytesUnpack = fieldPack->data_bytes<std::byte>();
  auto entityBytesUnpack1 = fieldBytesUnpack.entity_bytes(elem1);
  auto entityBytesUnpack2 = fieldBytesUnpack.entity_bytes(elem2);
  check_peek_bytes(initValsPack, entityBytesUnpack1);
  check_unpack_bytes(initValsPack, entityBytesUnpack1);
  check_peek_bytes(initValsPack, entityBytesUnpack2);
  check_unpack_bytes(initValsPack, entityBytesUnpack2);
}

#endif // have MPI and STKIO

}

