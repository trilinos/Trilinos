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

#ifndef STK_MESH_FIXTURES_TESTHEXFIXTURE_HPP
#define STK_MESH_FIXTURES_TESTHEXFIXTURE_HPP

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <vector>
#include <string>

namespace stk {
namespace mesh {
namespace fixtures {

class TestHexFixture : public ::ngp_testing::Test
{
protected:
  TestHexFixture()
  : m_bulk(stk::mesh::MeshBuilder(MPI_COMM_WORLD).create()),
    m_meta(m_bulk->mesh_meta_data()),
    m_hexFixture(nullptr)
  {
    m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");
  }

  virtual ~TestHexFixture() { delete m_hexFixture; }

  void setup_mesh(size_t nx, size_t ny, size_t nz,
                  const std::vector<std::string>& entityRankNames = std::vector<std::string>())
  {
    STK_ThrowRequireMsg(m_hexFixture == nullptr, "TestHexFixture::setup_mesh may only be called once.");
    m_meta.initialize(3, entityRankNames, m_coord_field->name());
    m_hexFixture = new HexFixture(m_meta, *m_bulk, nx, ny, nz, 1, 1);
    m_hexFixture->m_meta.commit();
    m_hexFixture->generate_mesh();
  }

  MPI_Comm get_comm() const { return m_bulk->parallel(); }

  int get_parallel_rank() const { return m_bulk->parallel_rank(); }

  int get_parallel_size() const { return m_bulk->parallel_size(); }

  stk::mesh::MetaData& get_meta() { return m_meta; }
  const stk::mesh::MetaData& get_meta() const { return m_meta; }

  stk::mesh::BulkData& get_bulk() { return *m_bulk; }
  const stk::mesh::BulkData& get_bulk() const { return *m_bulk; }

  HexFixture::CoordFieldType& get_coord_field() { return *m_coord_field; }

private:
  std::shared_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData& m_meta;
  HexFixture::CoordFieldType* m_coord_field;
  HexFixture* m_hexFixture;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
TestHexFixture : public ::ngp_testing::Test
{
protected:
  TestHexFixture()
  : m_bulk(stk::mesh::MeshBuilder(MPI_COMM_WORLD).create()),
    m_meta(m_bulk->mesh_meta_data()),
    m_hexFixture(nullptr)
  {
    m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");
  }

  virtual ~TestHexFixture() { delete m_hexFixture; }

  void setup_mesh(size_t nx, size_t ny, size_t nz,
                  const std::vector<std::string>& entityRankNames = std::vector<std::string>())
  {
    STK_ThrowRequireMsg(m_hexFixture == nullptr, "TestHexFixture::setup_mesh may only be called once.");
    m_meta.initialize(3, entityRankNames, m_coord_field->name());
    m_hexFixture = new stk::mesh::fixtures::HexFixture(m_meta, *m_bulk, nx, ny, nz, 1, 1);
    m_hexFixture->m_meta.commit();
    m_hexFixture->generate_mesh();
  }

  MPI_Comm get_comm() const { return m_bulk->parallel(); }

  int get_parallel_rank() const { return m_bulk->parallel_rank(); }

  int get_parallel_size() const { return m_bulk->parallel_size(); }

  stk::mesh::MetaData& get_meta() { return m_meta; }
  const stk::mesh::MetaData& get_meta() const { return m_meta; }

  stk::mesh::BulkData& get_bulk() { return *m_bulk; }
  const stk::mesh::BulkData& get_bulk() const { return *m_bulk; }

  stk::mesh::Field<double>& get_coord_field() { return *m_coord_field; }

private:
  std::shared_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData& m_meta;
  stk::mesh::Field<double>* m_coord_field;
  stk::mesh::fixtures::HexFixture* m_hexFixture;
};

} // namespace simple_fields

}}}

#endif

