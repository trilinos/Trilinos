// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_PhaseInput.hpp>
#include <Akri_Surface_Manager.hpp>

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <memory>

namespace krino {

class SurfaceManagerFixture : public ::testing::Test
{
protected:
  SurfaceManagerFixture()
    : bulk(stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create()),
      meta(bulk->mesh_meta_data()) {}

  std::unique_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData & meta;
};

TEST_F(SurfaceManagerFixture, SetSurfaceNames_KeepsParsedIdentifiers)
{
  SurfaceNames parsedNames;
  const Surface_Identifier ls1 = parsedNames.get_identifier("ls1");
  const Surface_Identifier ls2 = parsedNames.get_identifier("ls2");

  Surface_Manager & manager = Surface_Manager::get(meta);
  manager.set_surface_names(parsedNames);

  EXPECT_EQ(ls1, manager.get_identifier("ls1"));
  EXPECT_EQ(ls2, manager.get_identifier("LS2"));
  EXPECT_EQ("LS2", manager.get_name(ls2));
  EXPECT_EQ(Surface_Identifier(2), manager.get_identifier("ls3")); // new surfaces follow the parsed ones
}

TEST_F(SurfaceManagerFixture, SetSurfaceNames_AfterSurfaceIsNamed_Throws)
{
  Surface_Manager & manager = Surface_Manager::get(meta);
  manager.get_identifier("ls1");
  EXPECT_ANY_THROW(manager.set_surface_names(SurfaceNames()));
}

TEST_F(SurfaceManagerFixture, SeparateMeshes_HaveSeparateSurfaceManagers)
{
  auto otherBulk = stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create();
  Surface_Manager & manager = Surface_Manager::get(meta);
  Surface_Manager & otherManager = Surface_Manager::get(otherBulk->mesh_meta_data());
  EXPECT_NE(&manager, &otherManager);

  manager.get_identifier("ls1");
  EXPECT_EQ(Surface_Identifier(0), otherManager.get_identifier("other"));
}

} // namespace krino
