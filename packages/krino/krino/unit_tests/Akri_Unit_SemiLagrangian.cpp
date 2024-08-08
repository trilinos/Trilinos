#include <Akri_AnalyticSurf.hpp>
#include <Akri_Composite_Surface.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_SemiLagrangian.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

class SemiLagrangianTriElements : public StkMeshTriFixture
{
public:
  SemiLagrangianTriElements() {}

  void build_quad_split_4tri()
  {
    myDistanceField = mMesh.mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, "distance");
    stk::mesh::put_field_on_mesh(myDistanceField.field(), mMesh.mesh_meta_data().universal_part(), 1, 1, nullptr);
    QuadSplit4Tri meshSpec;
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, mBuilder.get_processor_distribution_for_num_elements(meshSpec.allElementConn.size()));
  }

  bool point_has_matching_Y_value(const stk::math::Vector3d & pt, const std::vector<double> & goldYValues)
  {
    const double absTol = 1.e-6;
    for(double goldYValue : goldYValues)
      if (std::abs(goldYValue - pt[1]) < absTol)
        return true;
    return false;
  }

  bool facet_points_all_have_matching_Y_value(const Facet2d & facet, const std::vector<double> & goldFacetYValues)
  {
    for (int pt=0; pt<2; ++pt)
      if(!point_has_matching_Y_value(facet.facet_vertex(pt), goldFacetYValues))
        return false;
    return true;
  }

  void test_facets(const FacetedSurfaceBase & facets, const std::vector<double> & goldFacetYValues)
  {
    ASSERT_FALSE(facets.get_facets_2d().empty());
    for (auto facet : facets.get_facets_2d())
      EXPECT_TRUE(facet_points_all_have_matching_Y_value(facet, goldFacetYValues)) << "Non matching facet " << facet;
  }

  void initialize_planar_surfaces(Composite_Surface & initSurfaces)
  {
    initSurfaces.add(new Plane(stk::math::Vector3d(0,-1,0), 0.2));
    initSurfaces.add(new Plane(stk::math::Vector3d(0,1,0), 0.2));
  }

  void set_interface_velocity(const std::vector<std::string> & interfaceVelocity)
  {
    initialize_expression_vector(interfaceVelocity, myInterfaceVelocity);
  }

protected:
  FieldRef myDistanceField;
  double avgEdgeLength{1.0};
  std::vector<String_Function_Expression> myInterfaceVelocity;
};

TEST_F(SemiLagrangianTriElements, initialTwoParallelPlanesIntersectingElementsMultipleTimes_adaptivelyContourElementsAndThenAdvect_exactlyRecoverPlanes)
{
  build_quad_split_4tri();

  Composite_Surface initSurfaces("init");
  initialize_planar_surfaces(initSurfaces);

  const stk::mesh::Selector activeFieldSelector = mMesh.mesh_meta_data().universal_part();
  BoundingBox nodeBBox = krino::compute_nodal_bbox(mMesh, activeFieldSelector, get_coordinates_field());
  constexpr double zeroNarrowBandSize = 0.;
  initSurfaces.prepare_to_compute(0.0, nodeBBox, zeroNarrowBandSize);

  // initialize
  const double time0 = 0.;
  std::unique_ptr<FacetedSurfaceBase> facetsOrig = FacetedSurfaceBase::build(2);
  compute_nodal_surface_distance(mMesh, get_coordinates_field(), myDistanceField, initSurfaces, time0, zeroNarrowBandSize);
  build_initial_adaptive_facets_after_nodal_distance_is_initialized_from_initial_surfaces(mMesh, activeFieldSelector, time0, get_coordinates_field(), myDistanceField, avgEdgeLength, initSurfaces, *facetsOrig);
  test_facets(*facetsOrig, {-0.2, 0.2});

  set_interface_velocity({"0.", "1."});

  const BoundingBox paddedNodeBBox = compute_padded_node_bounding_box_for_semilagrangian(mMesh, activeFieldSelector, 0., 0.5, get_coordinates_field(), myInterfaceVelocity, *facetsOrig);
  facetsOrig->prepare_to_compute(paddedNodeBBox, 0.);

  // single step advection (still surprised that this is unstable in the circle advection, molenkamp, test (not here))
  std::unique_ptr<FacetedSurfaceBase> facetsEnd = FacetedSurfaceBase::build(2);
  calc_single_step_semilagrangian_nodal_distance_and_build_facets(mMesh, activeFieldSelector, 0.0, 0.5, get_coordinates_field(), myDistanceField, myInterfaceVelocity, zeroNarrowBandSize, avgEdgeLength, *facetsOrig, *facetsEnd);
  test_facets(*facetsEnd, {0.3});

  // predistor-corrector (seems super stable and accurate, but expensive)
  std::unique_ptr<FacetedSurfaceBase> facetsPred = FacetedSurfaceBase::build(2);
  predict_semilagrangian_nodal_distance_and_build_facets(mMesh, activeFieldSelector, 0.0, 0.5, get_coordinates_field(), myDistanceField, myInterfaceVelocity, zeroNarrowBandSize, avgEdgeLength, *facetsOrig, *facetsPred);

  facetsPred->prepare_to_compute(paddedNodeBBox, 0.);
  correct_semilagrangian_nodal_distance_and_build_facets(mMesh, activeFieldSelector, 0.0, 0.5, get_coordinates_field(), myDistanceField, myInterfaceVelocity, zeroNarrowBandSize, avgEdgeLength, *facetsOrig, *facetsPred, *facetsEnd);
  test_facets(*facetsEnd, {0.3});
}
}

