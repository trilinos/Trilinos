// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_UNIT_SINGLE_ELEMENT_FIXTURES_H_
#define AKRI_UNIT_SINGLE_ELEMENT_FIXTURES_H_

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData

#include "Akri_AuxMetaData.hpp"
#include "Akri_FieldRef.hpp"

namespace krino {

class SimpleStkFixture
{
public:
  SimpleStkFixture(unsigned dimension, MPI_Comm comm = MPI_COMM_WORLD) {
    bulk = stk::mesh::MeshBuilder(comm)
               .set_spatial_dimension(dimension)
               .create();
    
    meta = bulk->mesh_meta_data_ptr();
    AuxMetaData::create(*meta);
  }
  void commit() { meta->commit(); }
  void write_results(const std::string & filename) { write_results(filename, *bulk); }
  static void write_results(const std::string & filename, stk::mesh::BulkData & mesh, const bool use64bitIds = true);
  const stk::mesh::MetaData & meta_data() const { return *meta; }
  stk::mesh::MetaData & meta_data() { return *meta; }
  stk::mesh::BulkData & bulk_data() { return *bulk; }

private:
  std::shared_ptr<stk::mesh::MetaData> meta;
  std::unique_ptr<stk::mesh::BulkData> bulk;
};

class SimpleStkFixture2d : public SimpleStkFixture
{
public:
  SimpleStkFixture2d(MPI_Comm comm = MPI_COMM_WORLD) : SimpleStkFixture(2) {}
};

class SimpleStkFixture3d : public SimpleStkFixture
{
public:
  SimpleStkFixture3d(MPI_Comm comm = MPI_COMM_WORLD) : SimpleStkFixture(3) {}
};

// Fixture to create a single element stk mesh of the given topology.
class SingleElementFixture
{
public:
  SingleElementFixture(const stk::topology & topology);

  void generate_mesh();

  stk::topology my_topology;
  SimpleStkFixture stk_fixture;
  stk::mesh::Part * block_part;
  FieldRef coord_field;
  FieldRef scalar_field;
  stk::mesh::Entity my_elem;
};


}

#endif /* AKRI_UNIT_SINGLE_ELEMENT_FIXTURES_H_ */
