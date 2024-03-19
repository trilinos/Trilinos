// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_Facet.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Snap.hpp>
#include <Akri_Unit_BoundingBoxMesh.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <gtest/gtest.h>
#include <stk_io/IossBridge.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <memory>

namespace krino {

static void expect_good_mesh(const stk::mesh::BulkData & mesh, const CDMesh & cdmesh)
{
  EXPECT_TRUE(check_induced_parts(mesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_TRUE(check_face_and_edge_relations(mesh));
  EXPECT_TRUE(check_shared_entity_nodes(mesh));
  EXPECT_TRUE(cdmesh.check_element_side_parts());
}

template <class MESH_FIXTURE>
class AnalyticDecompositionFixture : public ::testing::Test
{
public:
  AnalyticDecompositionFixture()
  : fixture(), cdfemSupport(CDFEM_Support::get(fixture.meta_data()))
  {
    AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());
    auto & vec_type = fixture.meta_data().spatial_dimension() == 2 ? FieldType::VECTOR_2D : FieldType::VECTOR_3D;
    coord_field = aux_meta.register_field("coordinates", vec_type, stk::topology::NODE_RANK, 1u, 1u, fixture.meta_data().universal_part());
    cdfemSupport.set_coords_field(coord_field);
    cdfemSupport.add_edge_interpolation_field(coord_field);

    cdfemSupport.set_prolongation_model(INTERPOLATION);
  }

  void setup_phase_support(const stk::mesh::PartVector & blocks)
  {
    PhaseTag p, n;
    p.add(surfaceIdentifier, 1);
    n.add(surfaceIdentifier, -1);
    PhaseVec named_phases{{"A", p}, {"B", n}};

    Phase_Support & phase_support = Phase_Support::get(fixture.meta_data());
    Block_Surface_Connectivity block_surface_info;
    phase_support.set_input_block_surface_connectivity(block_surface_info);
    phase_support.register_blocks_for_level_set(surfaceIdentifier, blocks);
    std::vector<std::tuple<stk::mesh::PartVector, 
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(blocks,interface_name_gen,named_phases));
    phase_support.decompose_blocks(ls_sets);
  }

  stk::mesh::Part & declare_input_block(const std::string & name, const stk::topology topo)
  {
    auto & block_part = fixture.meta_data().declare_part_with_topology(name, topo);
    stk::io::put_io_part_attribute(block_part);
    return block_part;
  }

  void decompose_mesh(const InterfaceGeometry & interfaceGeometry)
  {
    NodeToCapturedDomainsMap nodesToCapturedDomains;
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
    {
      const double minIntPtWeightForEstimatingCutQuality = cdfemSupport.get_snapper().get_edge_tolerance();
      nodesToCapturedDomains = snap_as_much_as_possible_while_maintaining_quality(krino_mesh->stk_bulk(),
          krino_mesh->get_active_part(),
          cdfemSupport.get_snap_fields(),
          interfaceGeometry,
          cdfemSupport.get_global_ids_are_parallel_consistent(),
          cdfemSupport.get_snapping_sharp_feature_angle_in_degrees(),
          minIntPtWeightForEstimatingCutQuality,
          cdfemSupport.get_max_edge_snap());
    }
    interfaceGeometry.prepare_to_decompose_elements(krino_mesh->stk_bulk(), nodesToCapturedDomains);

    krino_mesh->generate_nonconformal_elements();
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      krino_mesh->snap_nearby_intersections_to_nodes(interfaceGeometry, nodesToCapturedDomains);
    krino_mesh->set_phase_of_uncut_elements(interfaceGeometry);
    krino_mesh->triangulate(interfaceGeometry);
    krino_mesh->decompose(interfaceGeometry);
    krino_mesh->stash_field_data(-1);
    krino_mesh->modify_mesh();
    krino_mesh->prolongation();
  }

  void commit()
  {
    krino_mesh = std::make_unique<CDMesh>(fixture.bulk_data());
  }

  void write_results(const std::string & filename)
  {
    stk::io::write_mesh(filename, fixture.bulk_data());
  }

  void test_build_good_mesh(const bool doSnapping,
      const InterfaceGeometry & interfaceGeometry,
      const typename BoundingBoxMesh::BoundingBoxType & domain,
      const double meshSize,
      const std::string & filename = "")
  {
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int parallel_size = stk::parallel_machine_size(pm);

    stk::mesh::MetaData & meta = fixture.meta_data();
    AuxMetaData & aux_meta = AuxMetaData::get(meta);

    if (parallel_size > 2) return;

    if (doSnapping)
      cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
    else
      cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_NODE);

    auto & block1_part = aux_meta.get_part("block_1");
    setup_phase_support({&block1_part});

    fixture.set_domain(domain, meshSize);
    fixture.populate_mesh();

    stk::mesh::BulkData & mesh = fixture.bulk_data();
    stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(),&aux_meta.active_part()});

    commit();

    try
    {
      decompose_mesh(interfaceGeometry);
    }
    catch (const std::exception & exception)
    {
      std::cout << "Decomposing mesh failed with exception:\n";
      std::cout << exception.what() << "\n";
      ASSERT_TRUE(false);
    }

    expect_good_mesh(mesh, *krino_mesh);

    std::cout << log.get_log() << std::endl;
    if (!filename.empty())
      write_results(filename);
  }

  MESH_FIXTURE fixture;
  FieldRef coord_field;
  CDFEM_Support & cdfemSupport;
  std::unique_ptr<CDMesh> krino_mesh;
  LogRedirecter log;
  Surface_Identifier surfaceIdentifier{0};
};

template <class MESH_FIXTURE>
class SphereDecompositionFixture : public AnalyticDecompositionFixture<MESH_FIXTURE>
{
public:
  SphereDecompositionFixture()
  {
    mySphereGeometry.reset(new AnalyticSurfaceInterfaceGeometry({this->surfaceIdentifier}, {&mySphere}, AuxMetaData::get(this->fixture.meta_data()).active_part(), this->cdfemSupport, Phase_Support::get(this->fixture.meta_data())));
  }
protected:
  const InterfaceGeometry & get_interface_geometry() const { return *mySphereGeometry; }
  typename BoundingBoxMesh::BoundingBoxType get_domain() const
  {
    const stk::math::Vector3d extents = (2 == this->fixture.meta_data().spatial_dimension()) ? stk::math::Vector3d(1.,1.,0.) : stk::math::Vector3d(1.,1.,1.);
    typename BoundingBoxMesh::BoundingBoxType domain(-0.5*extents, 0.5*extents);
    return domain;
  }
  double get_mesh_size() const { return 1./6.; }
  Sphere mySphere{stk::math::Vector3d::ZERO, 0.35};
  std::unique_ptr<AnalyticSurfaceInterfaceGeometry> mySphereGeometry;
};

typedef SphereDecompositionFixture<BoundingBoxMeshTri3> CDMeshSphereTestsBboxMesh2D;
TEST_F(CDMeshSphereTestsBboxMesh2D, Sphere_SnapMesh)
{
  const bool doSnap = true;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "2DSnapSphere.e");
}

TEST_F(CDMeshSphereTestsBboxMesh2D, Sphere_CutMesh)
{
  const bool doSnap = false;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "2DCutSphere.e");
}

typedef SphereDecompositionFixture<BoundingBoxMeshTet4> CDMeshSphereTestsBboxMesh3D;
TEST_F(CDMeshSphereTestsBboxMesh3D, Sphere_SnapMesh)
{
  const bool doSnap = true;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "3DSnapSphere.e");
}

TEST_F(CDMeshSphereTestsBboxMesh3D, Sphere_CutMesh)
{
  const bool doSnap = false;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "3DCutSphere.e");
}

template <class MESH_FIXTURE>
class CubeDecompositionFixture : public AnalyticDecompositionFixture<MESH_FIXTURE>
{
public:
  CubeDecompositionFixture()
  {
    const double x = 0.25;
    const std::array<stk::math::Vector3d,8> cubeVerts =
      {{
        {-x,-x,-x}, {+x,-x,-x}, {+x,+x,-x}, {-x,+x,-x},
        {-x,-x,+x}, {+x,-x,+x}, {+x,+x,+x}, {-x,+x,+x}
      }};

    const std::array<std::array<int,3>,12> facetsVerts
      {{
        {{0,4,7}}, {{0,7,3}},
        {{1,2,6}}, {{1,6,5}},
        {{0,1,5}}, {{0,5,4}},
        {{2,3,7}}, {{2,7,6}},
        {{0,3,2}}, {{0,2,1}},
        {{4,5,6}}, {{4,6,7}}
      }};
    for (auto && facetVerts : facetsVerts)
      myCube.emplace_back_3d( cubeVerts[facetVerts[0]], cubeVerts[facetVerts[1]], cubeVerts[facetVerts[2]] );

    myCubeGeometry.reset(new AnalyticSurfaceInterfaceGeometry({this->surfaceIdentifier}, {&myCube}, AuxMetaData::get(this->fixture.meta_data()).active_part(), this->cdfemSupport, Phase_Support::get(this->fixture.meta_data())));
  }
protected:
  const InterfaceGeometry & get_interface_geometry() const { return *myCubeGeometry; }
  typename BoundingBoxMesh::BoundingBoxType get_domain() const
  {
    const stk::math::Vector3d extents = (2 == this->fixture.meta_data().spatial_dimension()) ? stk::math::Vector3d(1.,1.,0.) : stk::math::Vector3d(1.,1.,1.);
    typename BoundingBoxMesh::BoundingBoxType domain(-0.5*extents, 0.5*extents);
    return domain;
  }
  double get_mesh_size() const { return 1./6.; }
  Faceted_Surface<Facet3d> myCube;
  std::unique_ptr<AnalyticSurfaceInterfaceGeometry> myCubeGeometry;
};

typedef CubeDecompositionFixture<BoundingBoxMeshTri3> CDMeshCubeTestsBboxMesh2D;
TEST_F(CDMeshCubeTestsBboxMesh2D, Cube_SnapMesh)
{
  const bool doSnap = true;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "2DSnapCube.e");
}

typedef CubeDecompositionFixture<BoundingBoxMeshTet4> CDMeshCubeTestsBboxMesh3D;
TEST_F(CDMeshCubeTestsBboxMesh3D, Sphere_SnapMesh)
{
  const bool doSnap = true;
  test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size()); //test_build_good_mesh(doSnap, get_interface_geometry(), get_domain(), get_mesh_size(), "3DSnapCube.e");
}

} // namespace krino

