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
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include "FieldDataAccessFixture.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>

namespace {

#ifndef STK_UNIFIED_MEMORY

class LegacyFieldDataBucketAccess : public FieldDataAccessFixture {};

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, host_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      int* bucketValues = stk::mesh::field_data(field, *bucket);
      bucketValues[entity] = ++value;
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      const int* constBucketValues = stk::mesh::field_data(field, *bucket);
      EXPECT_EQ(constBucketValues[entity], ++value);
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      int* bucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
      bucketValuesBase[entity] = ++value*10;
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      const int* constBucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
      EXPECT_EQ(constBucketValuesBase[entity], ++value*10);
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, host_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const unsigned numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned component = 0; component < numComponents; ++component) {
        int* bucketValues = stk::mesh::field_data(field, *bucket);
        bucketValues[entity*numComponents + component] = ++value;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned component = 0; component < numComponents; ++component) {
        const int* constBucketValues = stk::mesh::field_data(field, *bucket);
        EXPECT_EQ(constBucketValues[entity*numComponents + component], ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned component = 0; component < numComponents; ++component) {
        int* bucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
        bucketValuesBase[entity*numComponents + component] = ++value*10;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned component = 0; component < numComponents; ++component) {
        const int* constBucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
        EXPECT_EQ(constBucketValuesBase[entity*numComponents + component], ++value*10);
      }
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, host_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const unsigned numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        int* bucketValues = stk::mesh::field_data(field, *bucket);
        bucketValues[entity*numCopies + copy] = ++value;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        const int* constBucketValues = stk::mesh::field_data(field, *bucket);
        EXPECT_EQ(constBucketValues[entity*numCopies + copy], ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        int* bucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
        bucketValuesBase[entity*numCopies + copy] = ++value*10;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        const int* constBucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
        EXPECT_EQ(constBucketValuesBase[entity*numCopies + copy], ++value*10);
      }
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, host_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const unsigned numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());
  const unsigned numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        for (unsigned component = 0; component < numComponents; ++component) {
          int* bucketValues = stk::mesh::field_data(field, *bucket);
          bucketValues[entity*numComponents*numCopies + copy*numComponents + component] = ++value;
        }
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        for (unsigned component = 0; component < numComponents; ++component) {
          const int* constBucketValues = stk::mesh::field_data(field, *bucket);
          EXPECT_EQ(constBucketValues[entity*numComponents*numCopies + copy*numComponents + component], ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        for (unsigned component = 0; component < numComponents; ++component) {
          int* bucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
          bucketValuesBase[entity*numComponents*numCopies + copy*numComponents + component] = ++value*10;
        }
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (unsigned entity = 0; entity < bucket->size(); ++entity) {
      for (unsigned copy = 0; copy < numCopies; ++copy) {
        for (unsigned component = 0; component < numComponents; ++component) {
          const int* constBucketValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, *bucket));
          EXPECT_EQ(constBucketValuesBase[entity*numComponents*numCopies + copy*numComponents + component], ++value*10);
        }
      }
    }
  }
}


//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, device_scalar) {
 // We don't have a legacy Bucket-based device-side data access method
 // that is equivalent to the new API.
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, device_multiComponent) {
 // We don't have a legacy Bucket-based device-side data access method
 // that is equivalent to the new API.
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, device_multiCopy) {
 // We don't have a legacy Bucket-based device-side data access method
 // that is equivalent to the new API.
}

//==============================================================================
TEST_F(LegacyFieldDataBucketAccess, device_multiCopy_multiComponent) {
 // We don't have a legacy Bucket-based device-side data access method
 // that is equivalent to the new API.
}

#endif

} //namespace <anonymous>
