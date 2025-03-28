#include <Akri_AuxMetaData.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Smooth.hpp>
#include <Akri_UnitTestUtils.hpp>

namespace krino {

void move_node(const CoordinatesFieldRef coordsField, const stk::mesh::Entity node, const stk::math::Vector3d & newNodeLocation)
{
  double * coords = field_data<double>(coordsField, node);
    for (unsigned d = 0; d < coordsField.dim(); ++d)
      coords[d] = newNodeLocation[d];
}

void move_node(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::math::Vector3d & newNodeLocation)
{
  const CoordinatesFieldRef coordsField(mesh.mesh_meta_data().coordinate_field(), mesh.mesh_meta_data().spatial_dimension());
  return move_node(coordsField, node, newNodeLocation);
}

stk::math::Vector3d get_node_location(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  const CoordinatesFieldRef coordsField(mesh.mesh_meta_data().coordinate_field(), mesh.mesh_meta_data().spatial_dimension());
  return get_vector_field(mesh, coordsField, node, coordsField.dim());
}

template <typename MESHSPEC>
class SmoothingFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  SmoothingFixture()
  {
  }

  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::write_mesh;

  double get_element_quality(const stk::mesh::Entity elem)
  {
    return element_quality(mMesh, mMesh.mesh_meta_data().spatial_dimension(), mMesh.mesh_meta_data().coordinate_field(), elem);
  }
protected:
  MESHSPEC meshSpec;
};

class CubeOf12TetsSmoothing : public SmoothingFixture<CubeOf12Tets>
{
public:
  CubeOf12TetsSmoothing()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1,1,1,1,1});
  }
  stk::mesh::Entity get_center_node() { return get_assigned_node_for_index(8); }
  stk::math::Vector3d get_original_center_node_location() { return stk::math::Vector3d::ZERO; }
  stk::math::Vector3d get_center_node_location() { return get_node_location(mMesh, get_center_node()); }
  void displace_center_node(const stk::math::Vector3d & disp)
  {
    const stk::math::Vector3d newNodeLocation = get_original_center_node_location() + disp;
    move_node(mMesh, get_center_node(), newNodeLocation);
  }
  CoordinatesFieldRef get_coordinates_field() { const CoordinatesFieldRef coordsField(mMesh.mesh_meta_data().coordinate_field(), mMesh.mesh_meta_data().spatial_dimension()); return coordsField; }

  double get_mesh_size() { return 1.0; }
  double get_location_tolerance() { return 1.e-3; }
};

TEST_F(CubeOf12TetsSmoothing, perturbedNode_smoothNodeUsingODT_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    improve_quality_by_ODT_smoothing(mMesh, get_coordinates_field(), elemSelector);

    expect_near_absolute(get_center_node_location(), get_original_center_node_location(), get_location_tolerance());
  }
}

TEST_F(CubeOf12TetsSmoothing, perturbedNode_smoothNodeUsingQuality_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    improve_quality_by_optimized_mean_ratio_smoothing(mMesh, get_coordinates_field(), elemSelector);

    expect_near_absolute(get_center_node_location(), get_original_center_node_location(), get_location_tolerance());
  }
}

}


