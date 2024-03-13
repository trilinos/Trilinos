// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourElement.hpp>
#include <gtest/gtest.h>

#include <Akri_Unit_Single_Element_Fixtures.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>

namespace krino {

TEST(SingleElementFixture, LS_Element_Tet4)
{
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  stk::mesh::MetaData & meta = test_fixture.stk_fixture.meta_data();
  stk::mesh::BulkData & bulk = test_fixture.stk_fixture.bulk_data();

  const unsigned ndim = meta.spatial_dimension();
  const stk::mesh::FieldBase & coords_field = test_fixture.coord_field.field();
  const stk::mesh::FieldBase & scalar_field = test_fixture.scalar_field.field();

  stk::mesh::Entity elem = test_fixture.my_elem;
  const stk::mesh::Entity* elem_nodes = bulk.begin_nodes(elem);
  const unsigned npe = bulk.num_nodes(elem);

  // Set coordinates
  const double coords[4][3] = {{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,4.}};
  for (unsigned node=0; node<npe; ++node)
  {
    double * coords_data = static_cast<double*>(stk::mesh::field_data(coords_field, elem_nodes[node]));
    for (unsigned dim=0; dim<ndim; ++dim)
    {
      coords_data[dim] = coords[node][dim];
    }
  }

  // Set isovariable
  const double isovar[4] = {-1., -1., -1., 1.};
  for (unsigned node=0; node<npe; ++node)
  {
    double * scalar_data = static_cast<double*>(stk::mesh::field_data(scalar_field, elem_nodes[node]));
    *scalar_data = isovar[node];
  }

  //
  // Create facets on 0 isosurface
  //

  const double isoval = 0.0;
  krino::ContourElement ls_elem( bulk, elem, coords_field, scalar_field, isoval );
  const double length_scale = 1.0; // Used for snapping facets to vertices of element when distance is small compared to length_scale
  ls_elem.compute_subelement_decomposition(length_scale);

  Faceted_Surface<Facet3d> faceted_surface;
  ls_elem.build_subelement_facets(faceted_surface);
  const auto & facets = faceted_surface.get_facets_3d();

  ASSERT_EQ(1u, facets.size());

  const Facet3d & facet = facets[0];

  EXPECT_EQ(2.0, facet.facet_vertex(0)[2]);
  EXPECT_EQ(2.0, facet.facet_vertex(1)[2]);
  EXPECT_EQ(2.0, facet.facet_vertex(2)[2]);

  const stk::math::Vector3d normal = facet.facet_normal();
  const double area = facet.facet_area();

  EXPECT_EQ(0.0, normal[0]);
  EXPECT_EQ(0.0, normal[1]);
  EXPECT_EQ(1.0, normal[2]);

  EXPECT_EQ(0.125, area);
}

template <typename MESHSPEC>
class ContourElementFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  ContourElementFixture()
  : myLs(LevelSet::build(mMesh.mesh_meta_data(), "LS", sierra::Diag::sierraTimer()))
  {
    if(stk::parallel_machine_size(mComm) == 1)
    {
      myLs.setup(); // registers field and sets field refs on the object
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    }
  }
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;
  void assign_nodal_level_set(const std::vector<double> & nodalLevelset)
  {
    const FieldRef distanceField = myLs.get_distance_field();
    for (size_t n=0; n<nodalLevelset.size(); ++n)
    {
      stk::mesh::Entity node = this->get_assigned_node_for_index(n);
      ASSERT_TRUE(mMesh.is_valid(node));
      double * dist = field_data<double>(distanceField, node);
      ASSERT_TRUE(dist != nullptr);
      *dist = nodalLevelset[n];
    }
  }
  void build_facets_in_element(const stk::mesh::Entity elem, FacetedSurfaceBase & facetedSurface)
  {
    facetedSurface.clear();
    ContourElement contourElement(mMesh, elem, myLs.get_coordinates_field(), myLs.get_distance_field(), 0.);
    contourElement.compute_subelement_decomposition(avgElemSize);
    contourElement.build_subelement_facets(facetedSurface);
  }
  void expect_num_facets_in_element(const int dim, const stk::mesh::Entity elem, const size_t goldNumFacets)
  {
    std::unique_ptr<FacetedSurfaceBase> facetedSurface = FacetedSurfaceBase::build(dim);
    build_facets_in_element(elem, *facetedSurface);
    EXPECT_EQ(goldNumFacets, facetedSurface->size());
  }

protected:
  const double eps{1.e-9};
  const double avgElemSize{1.};
  MESHSPEC meshSpec;
  LevelSet & myLs;
};

typedef ContourElementFixture<TwoTri306090> ContourElementTwoTri306090;

TEST_F(ContourElementTwoTri306090, twoElementsWhereSnappedSignPatternIsDifferentThanActualSignPattern_buildExpectedFacets)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    assign_nodal_level_set({{-eps, 1., -eps, -1.}});
    const std::vector<stk::mesh::Entity> & elems = get_owned_elements();
    ASSERT_EQ(2u, elems.size());

    expect_num_facets_in_element(2, elems[0], 0);
    expect_num_facets_in_element(2, elems[1], 1);
  }
}

typedef ContourElementFixture<TwoQuads> ContourElementTwoQuads;

TEST_F(ContourElementTwoQuads, twoQuadElementsWithSmallVariationAcrossElementAndSnapping_buildExpectedFacets)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const double epsCloseToSnapTol = 0.9e-4;
    assign_nodal_level_set({-epsCloseToSnapTol, -1.5*epsCloseToSnapTol, -1.5*epsCloseToSnapTol, -epsCloseToSnapTol, 1., 1.});
    const std::vector<stk::mesh::Entity> & elems = get_owned_elements();
    ASSERT_EQ(2u, elems.size());

    expect_num_facets_in_element(2, elems[0], 1);
    expect_num_facets_in_element(2, elems[1], 0);
  }
}

}
