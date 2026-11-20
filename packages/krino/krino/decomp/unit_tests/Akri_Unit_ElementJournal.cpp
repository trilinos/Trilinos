// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_UnitCubitJournalUtils.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace krino
{

namespace
{

class RightTetJournalFixture : public DecompositionFixture<RightTet, LSPerInterfacePolicy, 1>
{
public:
  RightTetJournalFixture()
  {
    if (stk::parallel_machine_size(mComm) > 1) return;
    build_mesh_from_spec(meshSpec);
    setup_ls_fields();
  }

  AnalyticSurfaceInterfaceGeometry build_geometry(const Surface & surface)
  {
    return AnalyticSurfaceInterfaceGeometry({levelset_fields()[0].identifier},
        {&surface},
        get_aux_meta().active_part(),
        cdfem_support(),
        Phase_Support::get(mMesh.mesh_meta_data()));
  }
};

} // namespace

static int count_block_names_starting_with(const std::string & journal, const std::string & start)
{
  std::istringstream lines(journal);
  std::string line;
  int count = 0;
  while (std::getline(lines, line))
  {
    const size_t pos = line.find(" name \"");
    if (line.rfind("block ", 0) == 0 && pos != std::string::npos && line.compare(pos+7, start.size(), start) == 0)
      ++count;
  }
  return count;
}

TEST_F(RightTetJournalFixture, sphereThroughTheRegion_drawnInItsNamedGroup)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Sphere sphere(stk::math::Vector3d(0., 0., 0.), 0.5);
  const AnalyticSurfaceInterfaceGeometry geometry = build_geometry(sphere);

  CubitJournal journal;
  geometry.add_surfaces_to_journal(journal,
      BoundingBox(stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., 1., 1.)),
      [](const Surface_Identifier) { return std::string("the_sphere"); });
  const std::string str = journal.str();

  EXPECT_NE(std::string::npos, str.find("name \"the_sphere_tri\""));
  EXPECT_FALSE(parse_cubit_journal(str).tris.empty());
}

TEST_F(RightTetJournalFixture, facetedSurfaceThroughTheRegion_itsFacetDrawn)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  MeshFacetedSurface surface(mComm, 3, {{-1., -1., 0.3}, {3., -1., 0.3}, {-1., 3., 0.3}}, {{{0, 1, 2}}});
  const BoundingBox region(stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., 1., 1.));
  surface.prepare_to_compute(0., region, 0.);
  const AnalyticSurfaceInterfaceGeometry geometry = build_geometry(surface);

  CubitJournal journal;
  geometry.add_surfaces_to_journal(journal, region, [](const Surface_Identifier) { return std::string("the_facets"); });

  EXPECT_EQ(1u, parse_cubit_journal(journal.str()).tris.size());
}

TEST_F(RightTetJournalFixture, sphereAwayFromTheRegion_nothingDrawn)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Sphere sphere(stk::math::Vector3d(5., 5., 5.), 0.5);
  const AnalyticSurfaceInterfaceGeometry geometry = build_geometry(sphere);

  CubitJournal journal;
  geometry.add_surfaces_to_journal(journal,
      BoundingBox(stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., 1., 1.)),
      [](const Surface_Identifier) { return std::string("the_sphere"); });

  EXPECT_TRUE(parse_cubit_journal(journal.str()).nodes.empty());
}


TEST_F(RightTetJournalFixture, tetCutBySphere_journalHasParentLeavesOfBothPhasesAndTheSphere)
{
  if (stk::parallel_machine_size(mComm) > 1) return;
  const Sphere sphere(stk::math::Vector3d(0., 0., 0.), 0.5);
  const AnalyticSurfaceInterfaceGeometry geometry = build_geometry(sphere);
  CDMesh::decompose_mesh(mMesh, geometry);
  const CDMesh & decomposedMesh = CDMesh::get_decomposed_mesh(mMesh);

  stk::mesh::EntityVector elements;
  stk::mesh::get_entities(mMesh, stk::topology::ELEMENT_RANK, elements);
  stk::mesh::EntityId parentId = 0;
  for (auto element : elements)
    if (nullptr != decomposedMesh.find_mesh_element(mMesh.identifier(element)))
      parentId = mMesh.identifier(element);
  ASSERT_NE(0u, parentId);
  const std::string prefix = "element_" + std::to_string(parentId) + "_";

  const std::string fileBaseName = "tetCutBySphereJournal";
  decomposedMesh.write_elements_journal(fileBaseName, {parentId}, geometry);
  const std::string fileName = fileBaseName + ".0.jou";
  std::ostringstream contents;
  contents << std::ifstream(fileName).rdbuf();
  std::remove(fileName.c_str());
  const std::string str = contents.str();
  const ParsedCubitJournal parsed = parse_cubit_journal(str);

  EXPECT_EQ(1, count_block_names_starting_with(str, prefix + "parent_tet"));
  EXPECT_GE(parsed.tets.size(), 3u);
  EXPECT_EQ(3, count_block_names_starting_with(str, prefix)); // the parent and both phases
  EXPECT_EQ(1, count_block_names_starting_with(str, "surface_"));
  EXPECT_FALSE(parsed.tris.empty());
}

} // namespace krino
