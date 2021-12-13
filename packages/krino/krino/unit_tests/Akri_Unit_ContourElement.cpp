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

  Faceted_Surface faceted_surface("tmp");
  ls_elem.build_subelement_facets(faceted_surface);
  const auto & facets = faceted_surface.get_facets();

  ASSERT_EQ(1u, facets.size());

  Facet & facet = *facets[0];

  EXPECT_EQ(2.0, facet.facet_vertex(0)[2]);
  EXPECT_EQ(2.0, facet.facet_vertex(1)[2]);
  EXPECT_EQ(2.0, facet.facet_vertex(2)[2]);

  const Vector3d normal = facet.facet_normal();
  const double area = facet.facet_area();

  EXPECT_EQ(0.0, normal[0]);
  EXPECT_EQ(0.0, normal[1]);
  EXPECT_EQ(1.0, normal[2]);

  EXPECT_EQ(0.125, area);
}

}
