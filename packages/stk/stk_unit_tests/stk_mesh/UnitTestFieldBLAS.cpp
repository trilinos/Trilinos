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

#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <sstream>
#include <memory>

namespace {

constexpr stk::mesh::Layout LayoutRight = stk::mesh::Layout::Right;
constexpr stk::mesh::Layout LayoutLeft  = stk::mesh::Layout::Left;

//==============================================================================
template <typename T,
          stk::mesh::Layout Layout1 = LayoutRight,
          stk::mesh::Layout Layout2 = LayoutRight,
          stk::mesh::Layout Layout3 = LayoutRight>
struct BLASFixtureScalar
{
  T initial_value1;
  T initial_value2;
  T initial_value3;
  unsigned numEntitiesUniversal;
  unsigned numEntitiesOwned;
  unsigned numEntitiesGlobal;

  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::Field<T, Layout1>* field1;
  stk::mesh::FieldBase* fieldBase1;
  stk::mesh::Field<T, Layout2>* field2;
  stk::mesh::FieldBase* fieldBase2;
  stk::mesh::Field<T, Layout3>* field3;
  stk::mesh::FieldBase* fieldBase3;

  stk::mesh::Part* pPartA;
  stk::mesh::Part* pPartB;
  unsigned numPartAEntitiesOwned;
  unsigned numPartBEntitiesOwned;
  unsigned numPartAEntitiesGlobal;
  unsigned numPartBEntitiesGlobal;

  BLASFixtureScalar(T init1, T init2 = T(), T init3 = T());
  ~BLASFixtureScalar() = default;
};

template <typename T, stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3>
BLASFixtureScalar<T, Layout1, Layout2, Layout3>::BLASFixtureScalar(T init1, T init2, T init3)
{
  initial_value1 = init1;
  initial_value2 = init2;
  initial_value3 = init3;

  MPI_Comm my_comm = MPI_COMM_WORLD;

  const double fractionToPartA = 0.3;
  const double fractionToPartB = 0.3;

  const unsigned int meshSizeX = 8;
  const unsigned int meshSizeY = 4;
  const unsigned int meshSizeZ = 4;

  // Keep small-ish Buckets so that there are many of them, to exercise OpenMP parallelization
  bulk = stk::mesh::MeshBuilder(my_comm).set_spatial_dimension(3)
                                        .set_maximum_bucket_capacity(16).create();
  stk::mesh::MetaData &meta = bulk->mesh_meta_data();

  field1 = &meta.declare_field<T, Layout1>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(*field1, meta.universal_part(), &initial_value1);
  fieldBase1 = field1;

  field2 = &meta.declare_field<T, Layout2>(stk::topology::NODE_RANK, "field2");
  stk::mesh::put_field_on_mesh(*field2, meta.universal_part(), &initial_value2);
  fieldBase2 = field2;

  field3 = &meta.declare_field<T, Layout3>(stk::topology::NODE_RANK, "field3");
  stk::mesh::put_field_on_mesh(*field3, meta.universal_part(), &initial_value3);
  fieldBase3 = field3;

  std::ostringstream osstr;
  osstr << "generated:" << meshSizeX << "x" << meshSizeY << "x" << meshSizeZ;
  stk::io::fill_mesh(osstr.str(), *bulk);

  pPartA = &meta.declare_part("PartA", stk::topology::NODE_RANK);
  pPartB = &meta.declare_part("PartB", stk::topology::NODE_RANK);

  std::vector<stk::mesh::Entity> entities = stk::mesh::get_entities(*bulk, stk::topology::NODE_RANK,
                                                                    meta.locally_owned_part());

  bulk->modification_begin();

  unsigned numToPartA = entities.size() * fractionToPartA;
  for (unsigned i = 0; i < numToPartA; ++i) {
    bulk->change_entity_parts(entities[i], stk::mesh::ConstPartVector{pPartA});
  }

  unsigned numToPartB = entities.size() * fractionToPartB;
  for (unsigned i = numToPartA; i < numToPartA+numToPartB; ++i) {
    bulk->change_entity_parts(entities[i], stk::mesh::ConstPartVector{pPartB});
  }

  bulk->modification_end();


  numEntitiesUniversal = 0;
  unsigned int numPartAEntitiesUniversal = 0;
  numPartAEntitiesOwned = 0;
  unsigned int numPartBEntitiesUniversal = 0;
  numPartBEntitiesOwned = 0;

  unsigned int numPartlessEntities = 0;
  numEntitiesOwned = 0;

  const stk::mesh::BucketVector& nodeBuckets = bulk->get_buckets(stk::topology::NODE_RANK, meta.universal_part());
  for (const stk::mesh::Bucket* bucket : nodeBuckets) {
    const unsigned bucketSize = bucket->size();
    numEntitiesUniversal += bucketSize;

    if (bucket->owned()) {
      numEntitiesOwned += bucketSize;
      if (bucket->member(*pPartA)) {
        numPartAEntitiesOwned += bucketSize;
      }
      if (bucket->member(*pPartB)) {
        numPartBEntitiesOwned += bucketSize;
      }
    }

    if (bucket->member(*pPartA)) {
      numPartAEntitiesUniversal += bucketSize;
    }
    if (bucket->member(*pPartB)) {
      numPartBEntitiesUniversal += bucketSize;
    }
    bool haveABEntities = bucket->member(*pPartA) && bucket->member(*pPartB);
    EXPECT_TRUE(!haveABEntities);
    if (!bucket->member(*pPartA) && !bucket->member(*pPartB)) {
      numPartlessEntities += bucketSize;
    }
  }

  EXPECT_EQ(numEntitiesUniversal, numPartAEntitiesUniversal+numPartBEntitiesUniversal+numPartlessEntities);

  numEntitiesGlobal = 0;
  numPartAEntitiesGlobal = 0;
  numPartBEntitiesGlobal = 0;
  stk::all_reduce_sum(bulk->parallel(), &numEntitiesOwned,      &numEntitiesGlobal,      1u);
  stk::all_reduce_sum(bulk->parallel(), &numPartAEntitiesOwned, &numPartAEntitiesGlobal, 1u);
  stk::all_reduce_sum(bulk->parallel(), &numPartBEntitiesOwned, &numPartBEntitiesGlobal, 1u);

  EXPECT_EQ(numEntitiesGlobal, (meshSizeX+1) * (meshSizeY+1) * (meshSizeZ+1));
}


template <typename T, stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3>
void checkScalarFields(BLASFixtureScalar<T, Layout1, Layout2, Layout3>& fixture,
                       T val1, T val2, T val3, double tol=1.0e-5)
{
  const stk::mesh::Selector selector = stk::mesh::selectField(*fixture.field1) &
                                       fixture.field1->get_mesh().mesh_meta_data().locally_owned_part();
  checkScalarFields(fixture, val1, val2, val3, selector, tol);
}

template <typename T, stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3>
void checkScalarFields(BLASFixtureScalar<T, Layout1, Layout2, Layout3>& fixture,
                       T val1, T val2, T val3, stk::mesh::Selector selector, double tol=1.0e-5)
{
  const stk::mesh::BucketVector& buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(), selector);
  auto field1Data = fixture.field1->template data<>();
  auto field2Data = fixture.field2->template data<>();
  auto field3Data = fixture.field3->template data<>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto field1Values = field1Data.bucket_values(*bucket);
    auto field2Values = field2Data.bucket_values(*bucket);
    auto field3Values = field3Data.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity : bucket->entities()) {
      EXPECT_NEAR(field1Values(entity), val1, tol);
      EXPECT_NEAR(field2Values(entity), val2, tol);
      EXPECT_NEAR(field3Values(entity), val3, tol);
    }
  }
}

template <typename T, stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3>
void checkScalarFields(BLASFixtureScalar<std::complex<T>, Layout1, Layout2, Layout3>& fixture,
                       std::complex<T> val1, std::complex<T> val2, std::complex<T> val3,
                       stk::mesh::Selector selector, double tol=1.0e-5)
{
  const stk::mesh::BucketVector& buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(), selector);
  auto field1Data = fixture.field1->template data<>();
  auto field2Data = fixture.field2->template data<>();
  auto field3Data = fixture.field3->template data<>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto field1Values = field1Data.bucket_values(*bucket);
    auto field2Values = field2Data.bucket_values(*bucket);
    auto field3Values = field3Data.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity : bucket->entities()) {
      EXPECT_LT(std::abs(field1Values(entity) - val1), tol);
      EXPECT_LT(std::abs(field2Values(entity) - val2), tol);
      EXPECT_LT(std::abs(field3Values(entity) - val3), tol);
    }
  }
}


//==============================================================================
template <typename T,
          stk::mesh::Layout Layout1 = LayoutRight,
          stk::mesh::Layout Layout2 = LayoutRight,
          stk::mesh::Layout Layout3 = LayoutRight>
struct BLASFixtureVector
{
  std::array<T, 3> init1;
  std::array<T, 3> init2;
  std::array<T, 3> init3;
  unsigned numEntitiesUniversal;
  unsigned numEntitiesOwned;
  unsigned numEntitiesGlobal;

  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::Field<T, Layout1>* field1;
  stk::mesh::FieldBase* fieldBase1;
  stk::mesh::Field<T, Layout2>* field2;
  stk::mesh::FieldBase* fieldBase2;
  stk::mesh::Field<T, Layout3>* field3;
  stk::mesh::FieldBase* fieldBase3;

  stk::mesh::Part* pPartA;
  stk::mesh::Part* pPartB;
  unsigned numPartAEntitiesOwned;
  unsigned numPartBEntitiesOwned;
  unsigned numPartAEntitiesGlobal;
  unsigned numPartBEntitiesGlobal;

  BLASFixtureVector(const std::array<T, 3>& init1_input, const std::array<T, 3>& init2_input = {T(), T(), T()},
                    const std::array<T, 3>& init3_input = {T(), T(), T()});
  ~BLASFixtureVector() = default;
};

template <typename T, stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3>
BLASFixtureVector<T, Layout1, Layout2, Layout3>::BLASFixtureVector(const std::array<T, 3>& init1_input,
                                                                   const std::array<T, 3>& init2_input,
                                                                   const std::array<T, 3>& init3_input)
{
  init1 = init1_input;
  init2 = init2_input;
  init3 = init3_input;

  MPI_Comm my_comm = MPI_COMM_WORLD;

  const double fractionToPartA = 0.3;
  const double fractionToPartB = 0.3;

  const unsigned int meshSizeX = 8;
  const unsigned int meshSizeY = 4;
  const unsigned int meshSizeZ = 4;

  // Keep small-ish Buckets so that there are many of them, to exercise OpenMP parallelization
  bulk = stk::mesh::MeshBuilder(my_comm).set_spatial_dimension(3)
                                        .set_maximum_bucket_capacity(16).create();
  stk::mesh::MetaData &meta = bulk->mesh_meta_data();

  field1 = &meta.declare_field<T, Layout1>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(*field1, meta.universal_part(), 3, init1.data());
  fieldBase1 = field1;

  field2 = &meta.declare_field<T, Layout2>(stk::topology::NODE_RANK, "field2");
  stk::mesh::put_field_on_mesh(*field2, meta.universal_part(), 3, init2.data());
  fieldBase2 = field2;

  field3 = &meta.declare_field<T, Layout3>(stk::topology::NODE_RANK, "field3");
  stk::mesh::put_field_on_mesh(*field3, meta.universal_part(), 3, init3.data());
  fieldBase3 = field3;

  std::ostringstream osstr;
  osstr << "generated:" << meshSizeX << "x" << meshSizeY << "x" << meshSizeZ;
  stk::io::fill_mesh(osstr.str(), *bulk);

  pPartA = &meta.declare_part("PartA", stk::topology::NODE_RANK);
  pPartB = &meta.declare_part("PartB", stk::topology::NODE_RANK);

  std::vector<stk::mesh::Entity> entities = stk::mesh::get_entities(*bulk, stk::topology::NODE_RANK,
                                                                    meta.locally_owned_part());

  bulk->modification_begin();

  unsigned numToPartA = entities.size() * fractionToPartA;
  for (unsigned i = 0; i < numToPartA; ++i) {
    bulk->change_entity_parts(entities[i], stk::mesh::ConstPartVector{pPartA});
  }

  unsigned numToPartB = entities.size() * fractionToPartB;
  for (unsigned i = numToPartA; i < numToPartA+numToPartB; ++i) {
    bulk->change_entity_parts(entities[i], stk::mesh::ConstPartVector{pPartB});
  }

  bulk->modification_end();


  numEntitiesUniversal = 0;
  unsigned int numPartAEntitiesUniversal = 0;
  numPartAEntitiesOwned = 0;
  unsigned int numPartBEntitiesUniversal = 0;
  numPartBEntitiesOwned = 0;

  unsigned int numPartlessEntities = 0;
  numEntitiesOwned = 0;

  const stk::mesh::BucketVector& nodeBuckets = bulk->get_buckets(stk::topology::NODE_RANK, meta.universal_part());
  for (const stk::mesh::Bucket* bucket : nodeBuckets) {
    const unsigned bucketSize = bucket->size();
    numEntitiesUniversal += bucketSize;

    if (bucket->owned()) {
      numEntitiesOwned += bucketSize;
      if (bucket->member(*pPartA)) {
        numPartAEntitiesOwned += bucketSize;
      }
      if (bucket->member(*pPartB)) {
        numPartBEntitiesOwned += bucketSize;
      }
    }

    if (bucket->member(*pPartA)) {
      numPartAEntitiesUniversal += bucketSize;
    }
    if (bucket->member(*pPartB)) {
      numPartBEntitiesUniversal += bucketSize;
    }
    bool haveABEntities = bucket->member(*pPartA) && bucket->member(*pPartB);
    EXPECT_TRUE(!haveABEntities);
    if (!bucket->member(*pPartA) && !bucket->member(*pPartB)) {
      numPartlessEntities += bucketSize;
    }
  }

  EXPECT_EQ(numEntitiesUniversal, numPartAEntitiesUniversal+numPartBEntitiesUniversal+numPartlessEntities);

  numEntitiesGlobal = 0;
  numPartAEntitiesGlobal = 0;
  numPartBEntitiesGlobal = 0;
  stk::all_reduce_sum(bulk->parallel(), &numEntitiesOwned,      &numEntitiesGlobal,      1u);
  stk::all_reduce_sum(bulk->parallel(), &numPartAEntitiesOwned, &numPartAEntitiesGlobal, 1u);
  stk::all_reduce_sum(bulk->parallel(), &numPartBEntitiesOwned, &numPartBEntitiesGlobal, 1u);

  EXPECT_EQ(numEntitiesGlobal, (meshSizeX+1) * (meshSizeY+1) * (meshSizeZ+1));
}

template <typename T, stk::mesh::Layout Layout>
void checkVectorField(const stk::mesh::Field<T, Layout>& field, const std::array<T, 3>& expectedValue,
                      const double tol=1.5e-3)
{
  const stk::mesh::BucketVector& buckets = field.get_mesh().get_buckets(field.entity_rank(),
                                                                        stk::mesh::selectField(field) &
                                                                        field.get_mesh().mesh_meta_data().locally_owned_part());
  auto fieldData = field.template data<>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto fieldValues = fieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity : bucket->entities()) {
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        EXPECT_NEAR(fieldValues(entity, component), expectedValue[component], tol);
      }
    }
  }
}

template <typename T, stk::mesh::Layout Layout>
void checkVectorField(const stk::mesh::Field<std::complex<T>, Layout>& field,
                      const std::array<std::complex<T>, 3>& expectedValue,
                      const double tol=1.5e-3)
{
  const stk::mesh::BucketVector& buckets = field.get_mesh().get_buckets(field.entity_rank(),
                                                                        stk::mesh::selectField(field) &
                                                                        field.get_mesh().mesh_meta_data().locally_owned_part());
  auto fieldData = field.template data<>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto fieldValues = fieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity : bucket->entities()) {
      for (stk::mesh::ComponentIdx component : fieldValues.components()) {
        EXPECT_LT(std::abs(fieldValues(entity, component) - expectedValue[component]), tol);
      }
    }
  }
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_axpy(T alpha, T initial1, T initial2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(initial1, initial2);

  stk::mesh::field_axpy(alpha, *fixture.field1, *fixture.field2);
  checkScalarFields(fixture, initial1, alpha*initial1+initial2, T());

  stk::mesh::field_axpy(alpha, *fixture.fieldBase1, *fixture.fieldBase2);
  checkScalarFields(fixture, initial1, alpha*initial1*T(2)+initial2, T());
}

TEST(FieldBLAS, scalar_double_axpy)
{
  const double alpha    = 7.11;
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_axpy<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_float_axpy)
{
  const float alpha    = 4.1;
  const float initial1 = 1.2;
  const float initial2 = -3.1;

  test_axpy<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_axpy)
{
  const std::complex<double> alpha    = std::complex<double>(-3.11, 2.00);
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_axpy<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_int_axpy)
{
  const int alpha    = 7;
  const int initial1 = 4;
  const int initial2 = -3;

  test_axpy<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_axpy_selector(T initial1, T initial2, T alpha_1, T alpha_2, T alpha_all)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(initial1, initial2);

  stk::mesh::field_axpy(alpha_1, *fixture.field1, *fixture.field2, stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, initial1, alpha_1*initial1+initial2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, initial1, initial2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, initial1, initial2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_axpy(alpha_2, *fixture.fieldBase1, *fixture.fieldBase2, stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, initial1, alpha_1*initial1+initial2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, initial1, alpha_2*initial1+initial2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, initial1, initial2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_axpy(alpha_all, *fixture.field1, *fixture.field2);
  checkScalarFields(fixture, initial1, (alpha_1+alpha_all)*initial1+initial2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, initial1, (alpha_2+alpha_all)*initial1+initial2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, initial1, alpha_all*initial1+initial2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());
}

TEST(FieldBLAS, scalar_double_selector_axpy)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  const double alpha_1   = 7.11;
  const double alpha_2   = 4.05;
  const double alpha_all = -2.04;

  test_axpy_selector<LayoutRight, LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutRight, LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
}

TEST(FieldBLAS, scalar_float_selector_axpy)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  const float alpha_1   = 7.1;
  const float alpha_2   = 4.7;
  const float alpha_all = -2.3;

  test_axpy_selector<LayoutRight, LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutRight, LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
}

TEST(FieldBLAS, scalar_complex_selector_axpy)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  const std::complex<double> alpha_1   = std::complex<double>(7.11, -42.1);
  const std::complex<double> alpha_2   = std::complex<double>(4.05, 7.22);
  const std::complex<double> alpha_all = std::complex<double>(-2.04, 3.14);

  test_axpy_selector<LayoutRight, LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutRight, LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
}

TEST(FieldBLAS, scalar_int_selector_axpy)
{
  const int initial1 = 4;
  const int initial2 = -3;

  const int alpha_1   = 7;
  const int alpha_2   = 5;
  const int alpha_all = -2;

  test_axpy_selector<LayoutRight, LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutRight, LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutRight>(initial1, initial2, alpha_1, alpha_2, alpha_all);
  test_axpy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2, alpha_1, alpha_2, alpha_all);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_axpy(T alpha, const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_axpy(alpha, *fixture.field1, *fixture.field2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2,
                   {alpha*init1[0] + init2[0], alpha*init1[1] + init2[1], alpha*init1[2] + init2[2]});

  stk::mesh::field_axpy(alpha, *fixture.fieldBase1, *fixture.fieldBase2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2,
                   {alpha*T(2)*init1[0] + init2[0], alpha*T(2)*init1[1] + init2[1], alpha*T(2)*init1[2] + init2[2]});
}

TEST(FieldBLAS, vector_double_axpy)
{
  const double alpha = 4.11;
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_axpy<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_float_axpy)
{
  const float alpha = 4.1;
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_axpy<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_complex_axpy)
{
  const std::complex<double> alpha(4.11, -32.1);
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_axpy<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_int_axpy)
{
  const int alpha = 4;
  std::array<int, 3> init1 = {4, 1, -2};
  std::array<int, 3> init2 = {3, 4, 11};

  test_axpy<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_axpy<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_axpby(T alpha, T init1, T beta, T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_axpby(alpha, *fixture.field1, beta, *fixture.field2);
  checkScalarFields(fixture, init1, alpha*init1 + beta*init2, T());

  stk::mesh::field_axpby(alpha, *fixture.fieldBase1, beta, *fixture.fieldBase2);
  checkScalarFields(fixture, init1, alpha*init1 + alpha*beta*init1 + beta*beta*init2, T());
}

TEST(FieldBLAS, scalar_double_unityBeta_axpby)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;
  const double alpha    = 7.11;
  const double beta     = 1.0;

  test_axpby<LayoutRight, LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, initial1, beta, initial2);
}

TEST(FieldBLAS, scalar_double_axpby)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;
  const double alpha    = 7.11;
  const double beta     = 2.0;

  test_axpby<LayoutRight, LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, initial1, beta, initial2);
}

TEST(FieldBLAS, scalar_float_axpby)
{
  const float initial1 = 1.2;
  const float initial2 = -3.1;
  const float alpha    = 4.1;
  const float beta     = 1.1;

  test_axpby<LayoutRight, LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, initial1, beta, initial2);
}

TEST(FieldBLAS, scalar_complex_axpby)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);
  const std::complex<double> alpha    = std::complex<double>(-3.11, 2.00);
  const std::complex<double> beta     = std::complex<double>(1.0, 0.5);

  test_axpby<LayoutRight, LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, initial1, beta, initial2);
}

TEST(FieldBLAS, scalar_int_axpby)
{
  const int initial1 = 4;
  const int initial2 = -3;
  const int alpha    = 7;
  const int beta     = 2;

  test_axpby<LayoutRight, LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, initial1, beta, initial2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, initial1, beta, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_axpby(T alpha, const std::array<T, 3>& init1, T beta, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_axpby(alpha, *fixture.field1, beta, *fixture.field2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2,
                   {alpha*init1[0] + beta*init2[0], alpha*init1[1] + beta*init2[1], alpha*init1[2] + beta*init2[2]});

  stk::mesh::field_axpby(alpha, *fixture.fieldBase1, beta, *fixture.fieldBase2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2,
                   {alpha*init1[0] + alpha*beta*init1[0] + beta*beta*init2[0],
                    alpha*init1[1] + alpha*beta*init1[1] + beta*beta*init2[1],
                    alpha*init1[2] + alpha*beta*init1[2] + beta*beta*init2[2]});
}

TEST(FieldBLAS, vector_double_unityBeta_axpby)
{
  const double alpha = 4.11;
  const double beta = 1.0;
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_axpby<LayoutRight, LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, init1, beta, init2);
}

TEST(FieldBLAS, vector_double_axpby)
{
  const double alpha = 4.11;
  const double beta = 1.1;
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_axpby<LayoutRight, LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, init1, beta, init2);
}

TEST(FieldBLAS, vector_float_axpby)
{
  const float alpha = 4.1;
  const float beta = 1.1;
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_axpby<LayoutRight, LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, init1, beta, init2);
}

TEST(FieldBLAS, vector_complex_axpby)
{
  const std::complex<double> alpha(4.11, -32.1);
  const std::complex<double> beta(1.1, 1.2);
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_axpby<LayoutRight, LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, init1, beta, init2);
}

TEST(FieldBLAS, vector_int_axpby)
{
  const int alpha = 4;
  const int beta = 2;
  std::array<int, 3> init1 = {4, 1, -2};
  std::array<int, 3> init2 = {3, 4, 11};

  test_axpby<LayoutRight, LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutRight, LayoutLeft >(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutRight>(alpha, init1, beta, init2);
  test_axpby<LayoutLeft,  LayoutLeft >(alpha, init1, beta, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3, typename T>
void test_product(T init1, T init2, T init3)
{
  BLASFixtureScalar<T, Layout1, Layout2, Layout3> fixture(init1, init2, init3);
  checkScalarFields(fixture, init1, init2, init3);

  stk::mesh::field_product(*fixture.field1, *fixture.field2, *fixture.field3);
  checkScalarFields(fixture, init1, init2, init1*init2);

  stk::mesh::field_product(*fixture.field3, *fixture.field1, *fixture.field1);
  checkScalarFields(fixture, init1*init1*init2, init2, init1*init2);
}

TEST(FieldBLAS, scalar_double_product)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;
  const double initial3 = 28.12;

  test_product<LayoutRight, LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
}

TEST(FieldBLAS, scalar_float_product)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;
  const float initial3 = 28.1;

  test_product<LayoutRight, LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
}

TEST(FieldBLAS, scalar_complex_product)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);
  const std::complex<double> initial3 = std::complex<double>(1.28, 3.11);

  test_product<LayoutRight, LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
}

TEST(FieldBLAS, scalar_int_product)
{
  const int initial1 = 4;
  const int initial2 = -3;
  const int initial3 = 1;

  test_product<LayoutRight, LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(initial1, initial2, initial3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(initial1, initial2, initial3);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_product_selector(T init1, T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_product(*fixture.field1, *fixture.field2, *fixture.field2, stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init1*init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_product(*fixture.fieldBase2, *fixture.fieldBase1, *fixture.fieldBase1,
                           stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init1*init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1*init2, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_product(*fixture.field1, *fixture.field2, *fixture.field2,
                           stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());
  stk::mesh::field_product(*fixture.field1, *fixture.field2, *fixture.field2,
                           stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());
  checkScalarFields(fixture, init1, init1*init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1*init2, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, T(pow(init1, 2))*init2, T(),
                    stk::mesh::Selector(*fixture.pPartA).complement() & stk::mesh::Selector(*fixture.pPartB).complement());
}

TEST(FieldBLAS, scalar_double_selector_product)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_product_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_product_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_product)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_product_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_product_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_product)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_product_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_product_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_product)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_product_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_product_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_product_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, stk::mesh::Layout Layout3, typename T>
void test_product(const std::array<T, 3>& init1, const std::array<T, 3>& init2, const std::array<T, 3>& init3)
{
  BLASFixtureVector<T, Layout1, Layout2, Layout3> fixture(init1, init2, init3);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2, init2);
  checkVectorField(*fixture.field3, init3);

  stk::mesh::field_product(*fixture.field1, *fixture.field2, *fixture.field3);
  checkVectorField(*fixture.field1, {init1[0], init1[1], init1[2]});
  checkVectorField(*fixture.field2, {init2[0], init2[1], init2[2]});
  checkVectorField(*fixture.field3, {init1[0]*init2[0], init1[1]*init2[1], init1[2]*init2[2]});

  stk::mesh::field_product(*fixture.fieldBase3, *fixture.fieldBase1, *fixture.fieldBase2);
  checkVectorField(*fixture.field1, {init1[0], init1[1], init1[2]});
  checkVectorField(*fixture.field2, {init1[0]*init1[0]*init2[0], init1[1]*init1[1]*init2[1], init1[2]*init1[2]*init2[2]});
  checkVectorField(*fixture.field3, {init1[0]*init2[0], init1[1]*init2[1], init1[2]*init2[2]});
}

TEST(FieldBLAS, vector_double_product)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};
  std::array<double, 3> init3 {0.24, -7.11, 3.21};

  test_product<LayoutRight, LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(init1, init2, init3);
}

TEST(FieldBLAS, vector_float_product)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};
  std::array<float, 3> init3 {0.2, -7.1, 3.2};

  test_product<LayoutRight, LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(init1, init2, init3);
}

TEST(FieldBLAS, vector_complex_product)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};
  std::array<std::complex<double>, 3> init3 {std::complex<double>(0.24, -1.22),
                                             std::complex<double>(-7.11, 42.1),
                                             std::complex<double>(3.21, 7.11)};

  test_product<LayoutRight, LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(init1, init2, init3);
}

TEST(FieldBLAS, vector_int_product)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};
  std::array<int, 3> init3 {8, -7, 3};

  test_product<LayoutRight, LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutRight, LayoutLeft,  LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutRight, LayoutLeft >(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutRight>(init1, init2, init3);
  test_product<LayoutLeft,  LayoutLeft,  LayoutLeft >(init1, init2, init3);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_copy(T initial1, T initial2)
{
  {
    BLASFixtureScalar<T, Layout1, Layout2> fixture(initial1, initial2);

    stk::mesh::field_copy(*fixture.field1, *fixture.field2);
    checkScalarFields(fixture, initial1, initial1, T(0));

    if (fixture.field1->has_device_data()) {
      EXPECT_EQ(fixture.field1->need_sync_to_device(), false);
    }
    else {
      EXPECT_EQ(fixture.field1->need_sync_to_device(), true);
    }

    EXPECT_EQ(fixture.field2->need_sync_to_device(), true);
  }

  {
    BLASFixtureScalar<T, Layout1, Layout2> fixture(initial1, initial2);

    stk::mesh::field_copy(*fixture.field2, *fixture.field1);
    checkScalarFields(fixture, initial2, initial2, T(0));

    EXPECT_EQ(fixture.field1->need_sync_to_device(), true);
    EXPECT_EQ(fixture.field2->need_sync_to_device(), true);
  }
}

TEST(FieldBLAS, scalar_double_copy)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_copy<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_copy)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_copy<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_copy)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_copy<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_copy)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_copy<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_copy_selector(T init1, T init2)
{
  {
    BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

    stk::mesh::field_copy(*fixture.field1, *fixture.field2, stk::mesh::Selector(*fixture.pPartA));
    checkScalarFields(fixture, init1, init1, T(), stk::mesh::Selector(*fixture.pPartA));
    checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
    checkScalarFields(fixture, init1, init2, T(),
                      stk::mesh::Selector(*fixture.pPartA).complement() &
                      stk::mesh::Selector(*fixture.pPartB).complement());
  }
}

TEST(FieldBLAS, scalar_double_selector_copy)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_copy_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_copy)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_copy_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_copy)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_copy_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_copy)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_copy_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_copy_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_copy(const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  {
    BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
    checkVectorField(*fixture.field1, init1);
    checkVectorField(*fixture.field2, init2);

    stk::mesh::field_copy(*fixture.field1, *fixture.field2);
    checkVectorField(*fixture.field1, init1);
    checkVectorField(*fixture.field2, init1);
  }
  {
    BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
    checkVectorField(*fixture.field1, init1);
    checkVectorField(*fixture.field2, init2);

    stk::mesh::field_copy(*fixture.field2, *fixture.field1);
    checkVectorField(*fixture.field1, init2);
    checkVectorField(*fixture.field2, init2);
  }
}

TEST(FieldBLAS, vector_double_copy)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_copy<LayoutRight, LayoutRight>(init1, init2);
  test_copy<LayoutRight, LayoutLeft >(init1, init2);
  test_copy<LayoutLeft,  LayoutRight>(init1, init2);
  test_copy<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_float_copy)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_copy<LayoutRight, LayoutRight>(init1, init2);
  test_copy<LayoutRight, LayoutLeft >(init1, init2);
  test_copy<LayoutLeft,  LayoutRight>(init1, init2);
  test_copy<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_complex_copy)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_copy<LayoutRight, LayoutRight>(init1, init2);
  test_copy<LayoutRight, LayoutLeft >(init1, init2);
  test_copy<LayoutLeft,  LayoutRight>(init1, init2);
  test_copy<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_int_copy)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_copy<LayoutRight, LayoutRight>(init1, init2);
  test_copy<LayoutRight, LayoutLeft >(init1, init2);
  test_copy<LayoutLeft,  LayoutRight>(init1, init2);
  test_copy<LayoutLeft,  LayoutLeft >(init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_dot(T init1, T init2, double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T field_result = stk::mesh::field_dot(*fixture.field1, *fixture.field2);
  EXPECT_LT(std::abs(field_result - init1 * init2 * T(fixture.numEntitiesGlobal)), tol);

  T fieldBase_result;
  stk::mesh::field_dot(fieldBase_result, *fixture.fieldBase1, *fixture.fieldBase2);
  EXPECT_LT(std::abs(fieldBase_result - init1 * init2 * T(fixture.numEntitiesGlobal)), tol);
}

TEST(FieldBLAS, scalar_double_dot)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_dot<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_dot)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_dot<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_dot)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_dot<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_dot)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_dot<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_dot_selector(T init1, T init2, double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T resultA = stk::mesh::field_dot(*fixture.field1, *fixture.field2, stk::mesh::Selector(*fixture.pPartA));
  EXPECT_LT(std::abs(init1 * init2 *T(fixture.numPartAEntitiesGlobal) - resultA), tol);

  T resultB {};
  stk::mesh::field_dot(resultB, *fixture.fieldBase2, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartB));
  EXPECT_LT(std::abs(init1 * init2 * T(fixture.numPartBEntitiesGlobal) - resultB), tol);

  T resultABc = stk::mesh::field_dot(*fixture.field1, *fixture.field2,
                                     stk::mesh::Selector(*fixture.pPartA).complement() &
                                     stk::mesh::Selector(*fixture.pPartB).complement());
  EXPECT_LT(std::abs(init1 * init2 * T(fixture.numEntitiesGlobal - fixture.numPartAEntitiesGlobal -
                                       fixture.numPartBEntitiesGlobal) - resultABc), tol);
}

TEST(FieldBLAS, scalar_double_selector_dot)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_dot_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_dot)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_dot_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_dot)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_dot_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_dot)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_dot_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_dot_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_dot(const std::array<T, 3>& init1, const std::array<T, 3>& init2, double tol = 1.e-3)
{
  {
    BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
    checkVectorField(*fixture.field1, init1);
    checkVectorField(*fixture.field2, init2);

    const T result = stk::mesh::field_dot(*fixture.field1, *fixture.field2);
    const T expectedResult = (init1[0] * init2[0] + init1[1] * init2[1] + init1[2] * init2[2]) *
                             T(fixture.numEntitiesGlobal);
    EXPECT_LT(std::abs(result - expectedResult), tol);
  }
  {
    BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
    checkVectorField(*fixture.field1, init1);
    checkVectorField(*fixture.field2, init2);

    T result;
    stk::mesh::field_dot(result, *fixture.fieldBase1, *fixture.fieldBase2);
    const T expectedResult = (init1[0] * init2[0] + init1[1] * init2[1] + init1[2] * init2[2]) *
                             T(fixture.numEntitiesGlobal);
    EXPECT_LT(std::abs(result - expectedResult), tol);
  }
}

TEST(FieldBLAS, vector_double_dot)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_dot<LayoutRight, LayoutRight>(init1, init2);
  test_dot<LayoutRight, LayoutLeft >(init1, init2);
  test_dot<LayoutLeft,  LayoutRight>(init1, init2);
  test_dot<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_float_dot)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_dot<LayoutRight, LayoutRight>(init1, init2);
  test_dot<LayoutRight, LayoutLeft >(init1, init2);
  test_dot<LayoutLeft,  LayoutRight>(init1, init2);
  test_dot<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_complex_dot)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_dot<LayoutRight, LayoutRight>(init1, init2);
  test_dot<LayoutRight, LayoutLeft >(init1, init2);
  test_dot<LayoutLeft,  LayoutRight>(init1, init2);
  test_dot<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_int_dot)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_dot<LayoutRight, LayoutRight>(init1, init2);
  test_dot<LayoutRight, LayoutLeft >(init1, init2);
  test_dot<LayoutLeft,  LayoutRight>(init1, init2);
  test_dot<LayoutLeft,  LayoutLeft >(init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_nrm2(T init1, T init2, double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T fieldResult1 = stk::mesh::field_nrm2(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult1 - T(sqrt(std::abs(init1)*std::abs(init1)*double(fixture.numEntitiesGlobal)))), tol);
  T fieldResult2 = stk::mesh::field_nrm2(*fixture.field2);
  EXPECT_LT(std::abs(fieldResult2 - T(sqrt(std::abs(init2)*std::abs(init2)*double(fixture.numEntitiesGlobal)))), tol);

  T fieldBaseResult1;
  stk::mesh::field_nrm2(fieldBaseResult1, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult1 - T(sqrt(std::abs(init1)*std::abs(init1)*double(fixture.numEntitiesGlobal)))), tol);
  T fieldBaseResult2;
  stk::mesh::field_nrm2(fieldBaseResult2, *fixture.fieldBase2);
  EXPECT_LT(std::abs(fieldBaseResult2 - T(sqrt(std::abs(init2)*std::abs(init2)*double(fixture.numEntitiesGlobal)))), tol);
}

TEST(FieldBLAS, scalar_double_nrm2)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_nrm2<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_nrm2)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_nrm2<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_nrm2)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_nrm2<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_nrm2)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_nrm2<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_nrm2_selector(T init1, T init2, double tol = 1.0e-1)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T resultA = stk::mesh::field_nrm2(*fixture.field1, stk::mesh::Selector(*fixture.pPartA));
  EXPECT_LT(std::abs(T(std::abs(init1)*sqrt(T(fixture.numPartAEntitiesGlobal))) - resultA), tol);

  T resultB;
  stk::mesh::field_nrm2(resultB, *fixture.fieldBase2, stk::mesh::Selector(*fixture.pPartB));
  EXPECT_LT(std::abs(T(std::abs(init2)*sqrt(T(fixture.numPartBEntitiesGlobal))) - resultB), tol);

  T resultABc = stk::mesh::field_nrm2(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                      stk::mesh::Selector(*fixture.pPartB).complement());
  EXPECT_LT(std::abs(T(std::abs(init1)*sqrt(T(fixture.numEntitiesGlobal - fixture.numPartAEntitiesGlobal -
                                              fixture.numPartBEntitiesGlobal))) - resultABc), tol);
}

TEST(FieldBLAS, scalar_double_selector_nrm2)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_nrm2_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_nrm2)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_nrm2_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_nrm2)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_nrm2_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_nrm2)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_nrm2_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_nrm2_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_nrm2(const std::array<T, 3>& init1, const std::array<T, 3>& init2, double tol = 1.e-3)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2, init2);

  {
    const T result = stk::mesh::field_nrm2(*fixture.field1);
    const T expectedResult = std::sqrt((std::pow(std::abs(init1[0]), 2) + std::pow(std::abs(init1[1]), 2) +
                                        std::pow(std::abs(init1[2]), 2)) * T(fixture.numEntitiesGlobal));
    EXPECT_LT(std::abs(result - expectedResult), tol);
  }

  {
    T result {};
    stk::mesh::field_nrm2(result, *fixture.fieldBase2);
    const T expectedResult = std::sqrt((std::pow(std::abs(init2[0]), 2) + std::pow(std::abs(init2[1]), 2) +
                                        std::pow(std::abs(init2[2]), 2)) * T(fixture.numEntitiesGlobal));
    EXPECT_LT(std::abs(result - expectedResult), tol);
  }
}

TEST(FieldBLAS, vector_double_nrm2)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_nrm2<LayoutRight, LayoutRight>(init1, init2);
  test_nrm2<LayoutRight, LayoutLeft >(init1, init2);
  test_nrm2<LayoutLeft,  LayoutRight>(init1, init2);
  test_nrm2<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_float_nrm2)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_nrm2<LayoutRight, LayoutRight>(init1, init2);
  test_nrm2<LayoutRight, LayoutLeft >(init1, init2);
  test_nrm2<LayoutLeft,  LayoutRight>(init1, init2);
  test_nrm2<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_complex_nrm2)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_nrm2<LayoutRight, LayoutRight>(init1, init2);
  test_nrm2<LayoutRight, LayoutLeft >(init1, init2);
  test_nrm2<LayoutLeft,  LayoutRight>(init1, init2);
  test_nrm2<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_int_nrm2)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_nrm2<LayoutRight, LayoutRight>(init1, init2);
  test_nrm2<LayoutRight, LayoutLeft >(init1, init2);
  test_nrm2<LayoutLeft,  LayoutRight>(init1, init2);
  test_nrm2<LayoutLeft,  LayoutLeft >(init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_scale(T alpha, T init1)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init1);

  stk::mesh::field_scale(alpha, *fixture.field1);
  stk::mesh::field_scale(alpha, *fixture.fieldBase2);
  checkScalarFields(fixture, alpha*init1, alpha*init1, T{});
}

TEST(FieldBLAS, scalar_double_scale)
{
  const double alpha = 4.27;
  const double initial1 = -3.73;

  test_scale<LayoutRight, LayoutRight>(alpha, initial1);
  test_scale<LayoutRight, LayoutLeft >(alpha, initial1);
  test_scale<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_float_scale)
{
  const float alpha = 4.2;
  const float initial1 = -3.7;

  test_scale<LayoutRight, LayoutRight>(alpha, initial1);
  test_scale<LayoutRight, LayoutLeft >(alpha, initial1);
  test_scale<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_complex_scale)
{
  const std::complex<double> alpha = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial1 = std::complex<double>(-7.21, -1.23);

  test_scale<LayoutRight, LayoutRight>(alpha, initial1);
  test_scale<LayoutRight, LayoutLeft >(alpha, initial1);
  test_scale<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_int_scale)
{
  const int alpha = 4;
  const int initial1 = -3;

  test_scale<LayoutRight, LayoutRight>(alpha, initial1);
  test_scale<LayoutRight, LayoutLeft >(alpha, initial1);
  test_scale<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_scale_selector(T alpha, T init1, T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_scale(alpha, *fixture.field1, stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, alpha*init1, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_scale(alpha, *fixture.fieldBase2, stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, alpha*init1, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, alpha*init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_scale(alpha, *fixture.field2, stk::mesh::Selector(*fixture.pPartA).complement() &
                         stk::mesh::Selector(*fixture.pPartB).complement());
  checkScalarFields(fixture, alpha*init1, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, alpha*init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, alpha*init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());
}

TEST(FieldBLAS, scalar_double_selector_scale)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;
  const double alpha = 2.13;

  test_scale_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_scale)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;
  const float alpha = 7.21;

  test_scale_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_scale)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);
  const std::complex<double> alpha = -4.2;

  test_scale_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_scale)
{
  const int initial1 = 4;
  const int initial2 = -3;
  const int alpha = 2;

  test_scale_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_scale_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_scale(T alpha, const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
  checkVectorField(*fixture.field1, init1);
  checkVectorField(*fixture.field2, init2);

  stk::mesh::field_scale(alpha, *fixture.field1);
  checkVectorField(*fixture.field1, {alpha*init1[0], alpha*init1[1], alpha*init1[2]});

  stk::mesh::field_scale(alpha, *fixture.fieldBase2);
  checkVectorField(*fixture.field2, {alpha*init2[0], alpha*init2[1], alpha*init2[2]});
}

TEST(FieldBLAS, vector_double_scale)
{
  const double alpha = 4.11;
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_scale<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_scale<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_float_scale)
{
  const float alpha = 4.1;
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_scale<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_scale<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_complex_scale)
{
  const std::complex<double> alpha = std::complex<double>(4.11, -32.1);
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_scale<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_scale<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_int_scale)
{
  const int alpha = 4;
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_scale<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_scale<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_scale<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_fill(const T alpha, const T init1)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init1);

  stk::mesh::field_fill(alpha, *fixture.field1);
  stk::mesh::field_fill(alpha, *fixture.fieldBase2);
  checkScalarFields(fixture, alpha, alpha, T());
}

TEST(FieldBLAS, scalar_double_fill)
{
  const double alpha = 4.27;
  const double initial1 = -3.73;

  test_fill<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_float_fill)
{
  const float alpha = 4.2;
  const float initial1 = -3.7;

  test_fill<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_complex_fill)
{
  const std::complex<double> alpha = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial1 = std::complex<double>(-7.21, -1.23);

  test_fill<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_int_fill)
{
  const int alpha = 4;
  const int initial1 = -3;

  test_fill<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, initial1);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_fill_many(const T alpha, const T init1)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init1);

  stk::mesh::field_fill(alpha, {fixture.fieldBase1, fixture.fieldBase2});
  checkScalarFields(fixture, alpha, alpha, T{});
}

TEST(FieldBLAS, scalar_double_many_fill)
{
  const double alpha = 4.27;
  const double initial1 = -3.73;

  test_fill_many<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill_many<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_float_many_fill)
{
  const float alpha = 4.2;
  const float initial1 = -3.7;

  test_fill_many<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill_many<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_complex_many_fill)
{
  const std::complex<double> alpha = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial1 = std::complex<double>(-7.21, -1.23);

  test_fill_many<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill_many<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutLeft >(alpha, initial1);
}

TEST(FieldBLAS, scalar_int_many_fill)
{
  const int alpha = 4;
  const int initial1 = -3;

  test_fill_many<LayoutRight, LayoutRight>(alpha, initial1);
  test_fill_many<LayoutRight, LayoutLeft >(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutRight>(alpha, initial1);
  test_fill_many<LayoutLeft,  LayoutLeft >(alpha, initial1);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_fill_selector(T alpha, T init1, T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_fill(alpha, *fixture.field1, stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, alpha, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_fill(alpha, *fixture.fieldBase2, stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, alpha, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, alpha, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_fill(T(0.0), *fixture.fieldBase1);
  checkScalarFields(fixture, T(0.0), init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, T(0.0), alpha, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, T(0.0), init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_fill(alpha, *fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                        stk::mesh::Selector(*fixture.pPartB).complement());
  checkScalarFields(fixture, T(0.0), init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, T(0.0), alpha, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, alpha, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());
}

TEST(FieldBLAS, scalar_double_selector_fill)
{
  const double alpha = 2.13;
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_fill_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_fill)
{
  const float alpha = 7.21;
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_fill_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_fill)
{
  const std::complex<double> alpha = -4.2;
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_fill_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_fill)
{
  const int alpha = 2;
  const int initial1 = 4;
  const int initial2 = -3;

  test_fill_selector<LayoutRight, LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutRight, LayoutLeft >(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutRight>(alpha, initial1, initial2);
  test_fill_selector<LayoutLeft,  LayoutLeft >(alpha, initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_fill(T alpha, const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_fill(alpha, *fixture.field1);
  stk::mesh::field_fill(alpha, *fixture.fieldBase2);

  checkVectorField(*fixture.field1, {alpha, alpha, alpha});
  checkVectorField(*fixture.field2, {alpha, alpha, alpha});
}

TEST(FieldBLAS, vector_double_fill)
{
  const double alpha = 4.11;
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_fill<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_float_fill)
{
  const float alpha = 4.1;
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_fill<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_complex_fill)
{
  const std::complex<double> alpha = std::complex<double>(4.11, -32.1);
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_fill<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_int_fill)
{
  const int alpha = 4;
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_fill<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_fill_component(const std::array<T, 3>& alpha, const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
  stk::mesh::field_fill_component(alpha.data(), *fixture.field1);
  stk::mesh::field_fill_component(alpha.data(), *fixture.fieldBase2);

  checkVectorField(*fixture.field1, alpha);
  checkVectorField(*fixture.field2, alpha);
}

TEST(FieldBLAS, vector_double_fill_component)
{
  std::array<double, 3> alpha {4.11, 2.11, -3.12};
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_fill_component<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_float_fill_component)
{
  std::array<float, 3> alpha {4.1, 2.1, -3.1};
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_fill_component<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_complex_fill_component)
{
  std::array<std::complex<double>, 3> alpha {std::complex<double>(4.11, -3.1),
                                             std::complex<double>(2.17, -0.25),
                                             std::complex<double>(7.14, -38.1)};
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.13),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.24),
                                             std::complex<double>(11.27, 4.21)};

  test_fill_component<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}

TEST(FieldBLAS, vector_int_fill_component)
{
  std::array<int, 3> alpha {2, -3, 6};
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_fill_component<LayoutRight, LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutRight, LayoutLeft >(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutRight>(alpha, init1, init2);
  test_fill_component<LayoutLeft,  LayoutLeft >(alpha, init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_swap(const T init1, const T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_swap(*fixture.field1, *fixture.fieldBase2);
  checkScalarFields(fixture, init2, init1, T{});
}

TEST(FieldBLAS, scalar_double_swap)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_swap<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_swap)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_swap<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_swap)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_swap<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_swap)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_swap<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_swap_selector(T init1, T init2)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  stk::mesh::field_swap(*fixture.field1, *fixture.field2, stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init2, init1, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_swap(*fixture.fieldBase2, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init2, init1, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init2, init1, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_swap(*fixture.fieldBase2, *fixture.fieldBase1);
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init2, init1, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());

  stk::mesh::field_swap(*fixture.field1, *fixture.field2, stk::mesh::Selector(*fixture.pPartA).complement() &
                        stk::mesh::Selector(*fixture.pPartB).complement());
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartB));
  checkScalarFields(fixture, init1, init2, T(), stk::mesh::Selector(*fixture.pPartA).complement() &
                    stk::mesh::Selector(*fixture.pPartB).complement());
}

TEST(FieldBLAS, swap_selector_double)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_swap_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, swap_selector_float)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_swap_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, swap_selector_complex)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_swap_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, swap_selector_int)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_swap_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_swap_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_swap(const std::array<T, 3>& init1, const std::array<T, 3>& init2)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);
  stk::mesh::field_swap(*fixture.field1, *fixture.fieldBase2);

  checkVectorField(*fixture.field1, init2);
  checkVectorField(*fixture.field2, init1);
}

TEST(FieldBLAS, vector_double_swap)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_swap<LayoutRight, LayoutRight>(init1, init2);
  test_swap<LayoutRight, LayoutLeft >(init1, init2);
  test_swap<LayoutLeft,  LayoutRight>(init1, init2);
  test_swap<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_float_swap)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_swap<LayoutRight, LayoutRight>(init1, init2);
  test_swap<LayoutRight, LayoutLeft >(init1, init2);
  test_swap<LayoutLeft,  LayoutRight>(init1, init2);
  test_swap<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_complex_swap)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_swap<LayoutRight, LayoutRight>(init1, init2);
  test_swap<LayoutRight, LayoutLeft >(init1, init2);
  test_swap<LayoutLeft,  LayoutRight>(init1, init2);
  test_swap<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_int_swap)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_swap<LayoutRight, LayoutRight>(init1, init2);
  test_swap<LayoutRight, LayoutLeft >(init1, init2);
  test_swap<LayoutLeft,  LayoutRight>(init1, init2);
  test_swap<LayoutLeft,  LayoutLeft >(init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_asum(T init1, T init2, double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T fieldResult1 = stk::mesh::field_asum(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult1 - std::abs(init1) * T(fixture.numEntitiesGlobal)), tol);
  T fieldResult2 = stk::mesh::field_asum(*fixture.field2);
  EXPECT_LT(std::abs(fieldResult2 - std::abs(init2) * T(fixture.numEntitiesGlobal)), tol);

  T fieldBaseResult1;
  stk::mesh::field_asum(fieldBaseResult1, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult1 - std::abs(init1) * T(fixture.numEntitiesGlobal)), tol);
  T fieldBaseResult2;
  stk::mesh::field_asum(fieldBaseResult2, *fixture.fieldBase2);
  EXPECT_LT(std::abs(fieldBaseResult2 - std::abs(init2) * T(fixture.numEntitiesGlobal)), tol);
}

TEST(FieldBLAS, scalar_double_asum)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_asum<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_asum)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_asum<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_asum)
{
  const std::complex<double> initial1 = std::complex<double>(4.11, -7.63);
  const std::complex<double> initial2 = std::complex<double>(-7.21, -1.23);

  test_asum<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_asum)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_asum<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_asum_selector(T init1, T init2, double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1, Layout2> fixture(init1, init2);

  T resultA = stk::mesh::field_asum(*fixture.field1, stk::mesh::Selector(*fixture.pPartA));
  EXPECT_LT(std::abs(std::abs(init1) * T(fixture.numPartAEntitiesGlobal) - resultA), tol);

  T resultB;
  stk::mesh::field_asum(resultB, *fixture.fieldBase2, stk::mesh::Selector(*fixture.pPartB));
  EXPECT_LT(std::abs(std::abs(init2) * T(fixture.numPartBEntitiesGlobal) - resultB), tol);

  T resultABc=stk::mesh::field_asum(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                    stk::mesh::Selector(*fixture.pPartB).complement());
  EXPECT_LT(std::abs(std::abs(init1) * T(fixture.numEntitiesGlobal - fixture.numPartAEntitiesGlobal -
                                         fixture.numPartBEntitiesGlobal) - resultABc), tol);
}

TEST(FieldBLAS, scalar_double_selector_asum)
{
  const double initial1 = 4.27;
  const double initial2 = -3.73;

  test_asum_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_float_selector_asum)
{
  const float initial1 = 4.2;
  const float initial2 = -3.7;

  test_asum_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_complex_selector_asum)
{
  const std::complex<double> initial1 = std::complex<double>(4.27, 2.1);
  const std::complex<double> initial2 = std::complex<double>(-3.73, 1.04);

  test_asum_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}

TEST(FieldBLAS, scalar_int_selector_asum)
{
  const int initial1 = 4;
  const int initial2 = -3;

  test_asum_selector<LayoutRight, LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutRight, LayoutLeft >(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutRight>(initial1, initial2);
  test_asum_selector<LayoutLeft,  LayoutLeft >(initial1, initial2);
}


template <stk::mesh::Layout Layout1, stk::mesh::Layout Layout2, typename T>
void test_asum(const std::array<T, 3>& init1, const std::array<T, 3>& init2, double tol = 2.0e-3)
{
  BLASFixtureVector<T, Layout1, Layout2> fixture(init1, init2);

  T result1 = stk::mesh::field_asum(*fixture.field1);

  T result2 {};
  stk::mesh::field_asum(result2, *fixture.fieldBase2);

  EXPECT_LT(std::abs(result1 - (std::abs(init1[0]) + std::abs(init1[1]) + std::abs(init1[2])) *
                     T(fixture.numEntitiesGlobal)), tol);
  EXPECT_LT(std::abs(result2 - (std::abs(init2[0]) + std::abs(init2[1]) + std::abs(init2[2])) *
                     T(fixture.numEntitiesGlobal)), tol);
}

TEST(FieldBLAS, vector_double_asum)
{
  std::array<double, 3> init1 {4.21, 1.23, -2.13};
  std::array<double, 3> init2 {1.32, 4.17, 11.27};

  test_asum<LayoutRight, LayoutRight>(init1, init2);
  test_asum<LayoutRight, LayoutLeft >(init1, init2);
  test_asum<LayoutLeft,  LayoutRight>(init1, init2);
  test_asum<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_float_asum)
{
  std::array<float, 3> init1 {4.2, 1.2, -2.1};
  std::array<float, 3> init2 {1.3, 4.1, 11.2};

  test_asum<LayoutRight, LayoutRight>(init1, init2);
  test_asum<LayoutRight, LayoutLeft >(init1, init2);
  test_asum<LayoutLeft,  LayoutRight>(init1, init2);
  test_asum<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_complex_asum)
{
  std::array<std::complex<double>, 3> init1 {std::complex<double>(4.21, 0.24),
                                             std::complex<double>(1.23, -0.1),
                                             std::complex<double>(-2.13, 4.11)};
  std::array<std::complex<double>, 3> init2 {std::complex<double>(1.32, 23.1),
                                             std::complex<double>(4.17, -0.2),
                                             std::complex<double>(11.27, 4.21)};

  test_asum<LayoutRight, LayoutRight>(init1, init2);
  test_asum<LayoutRight, LayoutLeft >(init1, init2);
  test_asum<LayoutLeft,  LayoutRight>(init1, init2);
  test_asum<LayoutLeft,  LayoutLeft >(init1, init2);
}

TEST(FieldBLAS, vector_int_asum)
{
  std::array<int, 3> init1 {4, 1, -2};
  std::array<int, 3> init2 {3, 4, 11};

  test_asum<LayoutRight, LayoutRight>(init1, init2);
  test_asum<LayoutRight, LayoutLeft >(init1, init2);
  test_asum<LayoutLeft,  LayoutRight>(init1, init2);
  test_asum<LayoutLeft,  LayoutLeft >(init1, init2);
}


//==============================================================================
template <stk::mesh::Layout Layout1, typename T>
void test_amax(T lowVal, T highVal, const double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1> fixture(lowVal);

  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntity(fieldValues.num_entities()/3);
  fieldValues(targetEntity) = highVal*std::abs(lowVal/highVal) + T(highVal - highVal*std::abs(lowVal/highVal))*mpiFrac;

  T fieldResult = stk::mesh::field_amax(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult - std::abs(highVal)), tol);

  T fieldBaseResult {};
  stk::mesh::field_amax(fieldBaseResult, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult - std::abs(highVal)), tol);
}

TEST(FieldBLAS, scalar_double_amax)
{
  const double lowVal  = 2.73;
  const double highVal = -10.27;

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_float_amax)
{
  const float lowVal  = 3.7;
  const float highVal = -10.2;

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_complex_amax)
{
  const std::complex<double> lowVal  = std::complex<double>(-1.11, -2.63);
  const std::complex<double> highVal = std::complex<double>(-100.21, -250.23);

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_int_amax)
{
  const int lowVal  = 2;
  const int highVal = -10;

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}


template <stk::mesh::Layout Layout1, typename T>
void test_amax_empty_selector(T lowVal, T highVal, const double tol = 1.0e-12)
{
  BLASFixtureScalar<T, Layout1> fixture(lowVal);

  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntity(fieldValues.num_entities()/3);
  fieldValues(targetEntity) = highVal*std::abs(lowVal/highVal) + T(highVal - highVal*std::abs(lowVal/highVal))*mpiFrac;

  T fieldResult = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector());
  EXPECT_LT(std::abs(fieldResult - T{}), tol);

  T fieldBaseResult {};
  stk::mesh::field_amax(fieldBaseResult, *fixture.fieldBase1, stk::mesh::Selector());
  EXPECT_LT(std::abs(fieldBaseResult - T{}), tol);
}

TEST(FieldBLAS, scalar_double_amax_empty_selector)
{
  const double lowVal  = 2.73;
  const double highVal = -10.27;

  test_amax_empty_selector<LayoutRight>(lowVal, highVal);
  test_amax_empty_selector<LayoutLeft >(lowVal, highVal);
}

template <stk::mesh::Layout Layout1, typename T>
void test_amax_field_not_on_part(T lowVal, T highVal, const double tol = 1.0e-12)
{
  BLASFixtureScalar<T, Layout1> fixture(lowVal);

  stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  meta.enable_late_fields();

  auto& partialField = meta.declare_field<double, Layout1>(stk::topology::NODE_RANK, "partial_field");
  stk::mesh::put_field_on_mesh(partialField, *fixture.pPartA, &lowVal);

  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(partialField.entity_rank(),
                                                                    *fixture.pPartA & meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = partialField.template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntity(fieldValues.num_entities()/3);
  fieldValues(targetEntity) = highVal*std::abs(lowVal/highVal) + T(highVal - highVal*std::abs(lowVal/highVal))*mpiFrac;

  T fieldResult = stk::mesh::field_amax(partialField, meta.universal_part());
  EXPECT_LT(std::abs(fieldResult - std::abs(highVal)), tol);
}

TEST(FieldBLAS, scalar_double_amax_field_not_on_part)
{
  const double lowVal  = 2.73;
  const double highVal = -10.27;

  test_amax_field_not_on_part<LayoutRight>(lowVal, highVal);
  test_amax_field_not_on_part<LayoutLeft >(lowVal, highVal);
}

template <stk::mesh::Layout Layout1, typename T>
void test_amax_selector(T lowVal, T highValA, T highValAB, T highValABc, const double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1> fixture(lowVal);

  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel()) + 1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));

  {
    const stk::mesh::BucketVector bucketsA = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                       stk::mesh::Selector(*fixture.pPartA) &
                                                                       meta.locally_owned_part());
    stk::mesh::Bucket& bA = *bucketsA[bucketsA.size()/3];

    auto fieldDataA = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesA = fieldDataA.bucket_values(bA);

    stk::mesh::EntityIdx targetEntityA(fieldValuesA.num_entities()/3);
    fieldValuesA(targetEntityA) = highValA*std::abs(lowVal/highValA) +
                                  T(highValA - highValA*std::abs(lowVal/highValA))*mpiFrac;

    T result;
    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(highValA) - result), tol);

    result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowVal) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA) |
                          stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highValA) - result), tol);

    result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                   stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(lowVal) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1);
    EXPECT_LT(std::abs(std::abs(highValA) - result), tol);
  }
  {
    const stk::mesh::BucketVector bucketsB = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                       stk::mesh::Selector(*fixture.pPartB) &
                                                                       meta.locally_owned_part());
    stk::mesh::Bucket& bB = *bucketsB[bucketsB.size()/3];

    auto fieldDataB = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesB = fieldDataB.bucket_values(bB);

    stk::mesh::EntityIdx targetEntityB(fieldValuesB.num_entities()/3);
    fieldValuesB(targetEntityB) = highValAB*std::abs(lowVal/highValAB) +
                                  T(highValAB - highValAB*std::abs(lowVal/highValAB))*mpiFrac;

    T result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(highValA) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highValAB) - result), tol);

    result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartA) |
                                   stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highValAB) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA).complement() &
                          stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(lowVal) - result), tol);
    EXPECT_LT(std::abs(std::abs(highValAB) - stk::mesh::field_amax(*fixture.field1)), tol);
  }
  {
    const stk::mesh::BucketVector bucketsABc = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                         stk::mesh::Selector(*fixture.pPartA).complement() &
                                                                         stk::mesh::Selector(*fixture.pPartB).complement() &
                                                                         meta.locally_owned_part());
    stk::mesh::Bucket& bABc = *bucketsABc[bucketsABc.size()/3];

    auto fieldDataABc = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesABc = fieldDataABc.bucket_values(bABc);

    stk::mesh::EntityIdx targetEntityABc(fieldValuesABc.num_entities()/3);
    fieldValuesABc(targetEntityABc) = highValABc*std::abs(lowVal/highValABc) +
                                      T(highValABc - highValABc*std::abs(lowVal/highValABc))*mpiFrac;

    T result;
    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(highValA) - result), tol);

    result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highValAB) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA) |
                          stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highValAB) - result), tol);

    result = stk::mesh::field_amax(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                   stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(highValABc) - result), tol);

    stk::mesh::field_amax(result, *fixture.fieldBase1);
    EXPECT_LT(std::abs(std::abs(highValABc) - result), tol);
  }
}

TEST(FieldBLAS, scalar_double_selector_amax)
{
  const double lowVal     = 1.27;
  const double highValA   = -3.73;
  const double highValAB  = -4.43;
  const double highValABc = -5.03;

  test_amax_selector<LayoutRight>(lowVal, highValA, highValAB, highValABc);
  test_amax_selector<LayoutLeft >(lowVal, highValA, highValAB, highValABc);
}

TEST(FieldBLAS, scalar_float_selector_amax)
{
  const float lowVal     = 1.3;
  const float highValA   = -3.7;
  const float highValAB  = -4.4;
  const float highValABc = -5.0;

  test_amax_selector<LayoutRight>(lowVal, highValA, highValAB, highValABc);
  test_amax_selector<LayoutLeft >(lowVal, highValA, highValAB, highValABc);
}

TEST(FieldBLAS, scalar_complex_selector_amax)
{
  const std::complex<double> lowVal     = std::complex<double>(0.51, 0.32);
  const std::complex<double> highValA   = std::complex<double>(-3.73, 4.04);
  const std::complex<double> highValAB  = std::complex<double>(4.95, -5.12);
  const std::complex<double> highValABc = std::complex<double>(-6.03, 6.11);

  test_amax_selector<LayoutRight>(lowVal, highValA, highValAB, highValABc);
  test_amax_selector<LayoutLeft >(lowVal, highValA, highValAB, highValABc);
}

TEST(FieldBLAS, scalar_int_selector_amax)
{
  const int lowVal     = 1;
  const int highValA   = -3;
  const int highValAB  = -4;
  const int highValABc = -5;

  test_amax_selector<LayoutRight>(lowVal, highValA, highValAB, highValABc);
  test_amax_selector<LayoutLeft >(lowVal, highValA, highValAB, highValABc);
}


template <stk::mesh::Layout Layout1, typename T>
void test_amax(const std::array<T, 3>& lowVal, const std::array<T, 3>& highVal, const double tol = 1.0e-3)
{
  BLASFixtureVector<T, Layout1> fixture(lowVal);

  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntity(fieldValues.num_entities()/3);
  for (stk::mesh::ScalarIdx scalar : fieldValues.scalars()) {
    fieldValues(targetEntity, scalar) = highVal[scalar]*std::abs(lowVal[scalar]/highVal[scalar]) +
        T(highVal[scalar] - highVal[scalar]*std::abs(lowVal[scalar]/highVal[scalar]))*mpiFrac;
  }

  T maxHighVal = std::max({std::abs(highVal[0]), std::abs(highVal[1]), std::abs(highVal[2])});

  T fieldResult = stk::mesh::field_amax(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult - maxHighVal), tol);

  T fieldBaseResult {};
  stk::mesh::field_amax(fieldBaseResult, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult - maxHighVal), tol);
}

TEST(FieldBLAS, vector_double_amax)
{
  std::array<double, 3> lowVal  {2.73, 1.23, 0.49};
  std::array<double, 3> highVal {-10.27, -9.66, 8.52};

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_float_amax)
{
  std::array<float, 3> lowVal  {2.7, 1.2, 0.4};
  std::array<float, 3> highVal {-9.6, -10.2, 8.5};

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_complex_amax)
{
  std::array<std::complex<double>, 3> lowVal  = {std::complex<double>(-1.11, -2.63),
                                                 std::complex<double>(-2.22, -1.65),
                                                 std::complex<double>(-1.58, -1.27)};
  std::array<std::complex<double>, 3> highVal = {std::complex<double>(-100.21, -250.23),
                                                 std::complex<double>(-110.54, -260.28),
                                                 std::complex<double>(-150.61, -331.87)};

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_int_amax)
{
  std::array<int, 3> lowVal  {2, 3, 5};
  std::array<int, 3> highVal {-10, -11, -28};

  test_amax<LayoutRight>(lowVal, highVal);
  test_amax<LayoutLeft >(lowVal, highVal);
}

//==============================================================================
template <stk::mesh::Layout Layout1, typename T>
void test_amin(T lowVal, T highVal, const double tol = 1.0e-3)
{
  BLASFixtureScalar<T, Layout1> fixture(highVal);
  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntityIdx(fieldValues.num_entities()/3);
  fieldValues(targetEntityIdx) = lowVal*std::abs(highVal/lowVal) + T(lowVal - lowVal*std::abs(highVal/lowVal))*mpiFrac;

  T fieldResult = stk::mesh::field_amin(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult - std::abs(lowVal)), tol);

  T fieldBaseResult {};
  stk::mesh::field_amin(fieldBaseResult, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult - std::abs(lowVal)), tol);
}

TEST(FieldBLAS, scalar_double_amin)
{
  const double lowVal  = 2.73;
  const double highVal = -10.27;

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_float_amin)
{
  const float lowVal  = 3.7;
  const float highVal = -10.2;

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_complex_amin)
{
  const std::complex<double> lowVal  = std::complex<double>(-1.11, -2.63);
  const std::complex<double> highVal = std::complex<double>(-100.21, -250.23);

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, scalar_int_amin)
{
  const int lowVal  = 2;
  const int highVal = -10;

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

template <stk::mesh::Layout Layout1, typename T>
void test_amin_empty_selector(T lowVal, T highVal, const double tol = 1.0e-12)
{
  BLASFixtureScalar<T, Layout1> fixture(highVal);
  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntityIdx(fieldValues.num_entities()/3);
  fieldValues(targetEntityIdx) = lowVal*std::abs(highVal/lowVal) + T(lowVal - lowVal*std::abs(highVal/lowVal))*mpiFrac;

  T fieldResult = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector());
  EXPECT_LT(std::abs(fieldResult - std::numeric_limits<T>::max()), tol);

  T fieldBaseResult {};
  stk::mesh::field_amin(fieldBaseResult, *fixture.fieldBase1, stk::mesh::Selector());
  EXPECT_LT(std::abs(fieldBaseResult - std::numeric_limits<T>::max()), tol);
}

TEST(FieldBLAS, scalar_double_amin_empty_selector)
{
  const double lowVal  = 2.73;
  const double highVal = -10.27;

  test_amin_empty_selector<LayoutRight>(lowVal, highVal);
  test_amin_empty_selector<LayoutLeft >(lowVal, highVal);
}


template <stk::mesh::Layout Layout1, typename T>
void test_amin_selector(T highVal, T lowValA, T lowValAB, T lowValABc, const double tol = 1.0e-1)
{
  BLASFixtureScalar<T, Layout1> fixture(highVal);
  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));

  {
    const stk::mesh::BucketVector bucketsA = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                       stk::mesh::Selector(*fixture.pPartA) &
                                                                       meta.locally_owned_part());
    stk::mesh::Bucket& bA = *bucketsA[bucketsA.size()/3];

    auto fieldDataA = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesA = fieldDataA.bucket_values(bA);

    stk::mesh::EntityIdx targetEntityA(fieldValuesA.num_entities()/3);
    fieldValuesA(targetEntityA) = lowValA*std::abs(highVal/lowValA) +
                                  T(lowValA - lowValA*std::abs(highVal/lowValA))*mpiFrac;

    T result;
    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(lowValA) - result), tol);

    result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(highVal) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA) |
                          stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowValA) - result), tol);

    result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                   stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(highVal) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1);
    EXPECT_LT(std::abs(std::abs(lowValA) - result), tol);
  }
  {
    const stk::mesh::BucketVector bucketsB = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                       stk::mesh::Selector(*fixture.pPartB) &
                                                                       meta.locally_owned_part());
    stk::mesh::Bucket& bB = *bucketsB[bucketsB.size()/3];

    auto fieldDataB = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesB = fieldDataB.bucket_values(bB);

    stk::mesh::EntityIdx targetEntityB(fieldValuesB.num_entities()/3);
    fieldValuesB(targetEntityB) = lowValAB*std::abs(highVal/lowValAB) +
                                  T(lowValAB - lowValAB*std::abs(highVal/lowValAB))*mpiFrac;

    T result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(lowValA) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowValAB) - result), tol);

    result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartA) |
                                   stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowValAB) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA).complement() &
                          stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(highVal) - result), tol);
    EXPECT_LT(std::abs(std::abs(lowValAB) - stk::mesh::field_amin(*fixture.field1)), tol);
  }
  {
    const stk::mesh::BucketVector bucketsABc = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                         stk::mesh::Selector(*fixture.pPartA).complement() &
                                                                         stk::mesh::Selector(*fixture.pPartB).complement() &
                                                                         meta.locally_owned_part());
    stk::mesh::Bucket& bABc = *bucketsABc[bucketsABc.size()/3];

    auto fieldDataABc = fixture.field1->template data<stk::mesh::Unsynchronized>();
    auto fieldValuesABc = fieldDataABc.bucket_values(bABc);

    stk::mesh::EntityIdx targetEntityABc(fieldValuesABc.num_entities()/3);
    fieldValuesABc(targetEntityABc) = lowValABc*std::abs(highVal/lowValABc) +
                                      T(lowValABc - lowValABc*std::abs(highVal/lowValABc))*mpiFrac;

    T result;
    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA));
    EXPECT_LT(std::abs(std::abs(lowValA) - result), tol);

    result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowValAB) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1, stk::mesh::Selector(*fixture.pPartA) |
                          stk::mesh::Selector(*fixture.pPartB));
    EXPECT_LT(std::abs(std::abs(lowValAB) - result), tol);

    result = stk::mesh::field_amin(*fixture.field1, stk::mesh::Selector(*fixture.pPartA).complement() &
                                   stk::mesh::Selector(*fixture.pPartB).complement());
    EXPECT_LT(std::abs(std::abs(lowValABc) - result), tol);

    stk::mesh::field_amin(result, *fixture.fieldBase1);
    EXPECT_LT(std::abs(std::abs(lowValABc) - result), tol);
  }
}

TEST(FieldBLAS, scalar_double_selector_amin)
{
  const double highVal   = -6.27;
  const double lowValA   = 5.73;
  const double lowValAB  = 4.43;
  const double lowValABc = 3.03;

  test_amin_selector<LayoutRight>(highVal, lowValA, lowValAB, lowValABc);
  test_amin_selector<LayoutLeft >(highVal, lowValA, lowValAB, lowValABc);
}

TEST(FieldBLAS, scalar_float_selector_amin)
{
  const float highVal   = -6.3;
  const float lowValA   = 5.7;
  const float lowValAB  = 4.4;
  const float lowValABc = 3.1;

  test_amin_selector<LayoutRight>(highVal, lowValA, lowValAB, lowValABc);
  test_amin_selector<LayoutLeft >(highVal, lowValA, lowValAB, lowValABc);
}

TEST(FieldBLAS, scalar_complex_selector_amin)
{
  const std::complex<double> highVal   = std::complex<double>(7.51, -8.32);
  const std::complex<double> lowValA   = std::complex<double>(-6.73, 6.04);
  const std::complex<double> lowValAB  = std::complex<double>(4.95, -5.12);
  const std::complex<double> lowValABc = std::complex<double>(-4.03, 4.11);

  test_amin_selector<LayoutRight>(highVal, lowValA, lowValAB, lowValABc);
  test_amin_selector<LayoutLeft >(highVal, lowValA, lowValAB, lowValABc);
}

TEST(FieldBLAS, scalar_int_selector_amin)
{
  const int highVal   = -6;
  const int lowValA   = 5;
  const int lowValAB  = 4;
  const int lowValABc = 3;

  test_amin_selector<LayoutRight>(highVal, lowValA, lowValAB, lowValABc);
  test_amin_selector<LayoutLeft >(highVal, lowValA, lowValAB, lowValABc);
}


template <stk::mesh::Layout Layout1, typename T>
void test_amin(const std::array<T, 3>& lowVal, const std::array<T, 3>& highVal, const double tol = 1.0e-3)
{
  BLASFixtureVector<T, Layout1> fixture(highVal);

  const stk::mesh::MetaData& meta = fixture.bulk->mesh_meta_data();
  const stk::mesh::BucketVector buckets = fixture.bulk->get_buckets(fixture.field1->entity_rank(),
                                                                    stk::mesh::selectField(*fixture.field1) &
                                                                    meta.locally_owned_part());

  double mpiFrac = double(stk::parallel_machine_rank(fixture.bulk->parallel())+1) /
                   double(stk::parallel_machine_size(fixture.bulk->parallel()));
  stk::mesh::Bucket& b = *buckets[buckets.size()/3];

  auto fieldData = fixture.field1->template data<stk::mesh::Unsynchronized>();
  auto fieldValues = fieldData.bucket_values(b);

  stk::mesh::EntityIdx targetEntity(fieldValues.num_entities()/3);
  for (stk::mesh::ScalarIdx scalar : fieldValues.scalars()) {
    fieldValues(targetEntity, scalar) = lowVal[scalar]*std::abs(highVal[scalar]/lowVal[scalar]) +
        T(lowVal[scalar] - lowVal[scalar]*std::abs(highVal[scalar]/lowVal[scalar]))*mpiFrac;
  }

  T minLowVal = std::min({std::abs(lowVal[0]), std::abs(lowVal[1]), std::abs(lowVal[2])});

  T fieldResult = stk::mesh::field_amin(*fixture.field1);
  EXPECT_LT(std::abs(fieldResult - minLowVal), tol);

  T fieldBaseResult {};
  stk::mesh::field_amin(fieldBaseResult, *fixture.fieldBase1);
  EXPECT_LT(std::abs(fieldBaseResult - minLowVal), tol);
}

TEST(FieldBLAS, vector_double_amin)
{
  std::array<double, 3> lowVal  {2.73, 1.23, 0.49};
  std::array<double, 3> highVal {-10.27, -9.66, 8.52};

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_float_amin)
{
  std::array<float, 3> lowVal  {2.7, 1.2, 0.4};
  std::array<float, 3> highVal {-9.6, -10.2, 8.5};

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_complex_amin)
{
  std::array<std::complex<double>, 3> lowVal  = {std::complex<double>(-1.11, -2.63),
                                                 std::complex<double>(-2.22, -1.65),
                                                 std::complex<double>(-1.58, -1.27)};
  std::array<std::complex<double>, 3> highVal = {std::complex<double>(-100.21, -250.23),
                                                 std::complex<double>(-110.54, -260.28),
                                                 std::complex<double>(-150.61, -331.87)};

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}

TEST(FieldBLAS, vector_int_amin)
{
  std::array<int, 3> lowVal  {2, 3, 5};
  std::array<int, 3> highVal {-10, -11, -28};

  test_amin<LayoutRight>(lowVal, highVal);
  test_amin<LayoutLeft >(lowVal, highVal);
}
}
