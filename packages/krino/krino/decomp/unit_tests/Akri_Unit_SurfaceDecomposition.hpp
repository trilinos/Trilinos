// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_SURFACEDECOMPOSITION_HPP_
#define KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_SURFACEDECOMPOSITION_HPP_

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <Akri_AllReduce.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMeshRespectingMeshBasedGeometry.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

// The mean ratios below which DecompositionMeasures counts the elements.
constexpr std::array<double, 3> lowMeanRatios{{0.05, 0.1, 0.2}};

struct DecompositionMeasures
{
  double negativeVolume;
  double interfaceArea;
  double worstQuality; // mean ratio
  double worstScaledJacobian;
  std::array<size_t, lowMeanRatios.size()> numElementsBelowLowMeanRatios;
  size_t numElements;
};

// The number of the selected elements whose mean ratio is below each of lowMeanRatios.
inline std::array<size_t, lowMeanRatios.size()> count_elements_below_low_mean_ratios(
    const stk::mesh::BulkData & mesh, const stk::mesh::Selector & ownedElements)
{
  const MeanRatioQualityMetric meanRatio;
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::array<size_t, lowMeanRatios.size()> counts{};
  std::vector<stk::math::Vector3d> nodeCoords;
  for (auto * bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedElements))
  {
    for (auto element : *bucket)
    {
      fill_element_node_coordinates(mesh, element, coordsField, nodeCoords);
      const double quality = meanRatio.get_element_quality_metric(dim, nodeCoords);
      for (size_t i = 0; i < lowMeanRatios.size(); ++i)
        if (quality < lowMeanRatios[i]) ++counts[i];
    }
  }
  return all_reduce_sum_array(mesh.parallel(), counts);
}

// A bounding box tet4 or tet10 mesh, as krino's flat-walled BCC generated meshes are (or tri3 or
// tri6 mesh in 2D, a flat-walled triangular lattice), with one surface whose
// negative side is one phase and positive side another, set up for decomposition with snapping, or
// only cutting, and with the sharp features cut before the surface or not (by default, they are).
class SurfaceDecomposition
{
public:
  SurfaceDecomposition(const BoundingBoxMesh::BoundingBoxType & domain,
      const double meshSize,
      const Edge_Degeneracy_Handling edgeDegeneracyHandling =
          SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE,
      const std::optional<bool> cutSharpFeaturesBeforeSurfaces = std::nullopt,
      const stk::topology elemTopology = stk::topology::TETRAHEDRON_4,
      const std::optional<double> edgeTol = std::nullopt)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    myLSFields = LSPerInterfacePolicy::setup_levelsets_on_all_blocks(meta, 1);

    CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
    cdfemSupport.set_coords_field(meta.coordinate_field());
    cdfemSupport.set_cdfem_edge_degeneracy_handling(edgeDegeneracyHandling);
    if (cutSharpFeaturesBeforeSurfaces)
      cdfemSupport.set_cut_sharp_features_before_surfaces(*cutSharpFeaturesBeforeSurfaces);
    if (edgeTol) cdfemSupport.set_cdfem_edge_tol(*edgeTol);
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.register_parent_node_ids_field();
    cdfemSupport.finalize_fields();

    myBBoxMesh.set_mesh_structure_type(elemTopology.base() == stk::topology::TETRAHEDRON_4 ?
            FLAT_WALLED_BCC_BOUNDING_BOX_MESH :
            FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH);
    populate_bounding_box_mesh_and_activate(
        myBBoxMesh, domain.get_min(), domain.get_max(), meshSize);
  }

  DecompositionMeasures decompose_and_measure(const Surface & surface)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    AuxMetaData & auxMeta = AuxMetaData::get(meta);
    const Phase_Support & phaseSupport = Phase_Support::get(meta);
    const Surface_Identifier surfaceIdentifier = myLSFields[0].identifier;

    const AnalyticSurfaceInterfaceGeometry geometry({surfaceIdentifier},
        {&surface},
        auxMeta.active_part(),
        CDFEM_Support::get(meta),
        phaseSupport);
    CDMesh::decompose_mesh(mesh, geometry);
    return measure(surfaceIdentifier);
  }

  // Smooths the decomposed mesh, moving only the nodes off every boundary (select_nodes_on_any_boundary),
  // and measures it again.
  DecompositionMeasures smooth_keeping_boundary_nodes_fixed_and_measure()
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const CoordinatesFieldRef coordsField(meta.coordinate_field(), meta.spatial_dimension());
    smooth_mesh_keeping_boundary_nodes_fixed(myBBoxMesh.bulk_data(), coordsField, AuxMetaData::get(meta).active_part());
    return measure(myLSFields[0].identifier);
  }

  // Smooths the decomposed mesh, the nodes on the topological features of the mesh, and with a feature
  // angle also on its sharp features, moving within their surfaces or along their curves
  // (smooth_mesh_respecting_mesh_based_geometry), and measures it again.
  DecompositionMeasures smooth_respecting_mesh_based_geometry_and_measure(const std::optional<double> cosFeatureAngle = std::nullopt)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const CoordinatesFieldRef coordsField(meta.coordinate_field(), meta.spatial_dimension());
    smooth_mesh_respecting_mesh_based_geometry(
        myBBoxMesh.bulk_data(), coordsField, AuxMetaData::get(meta).active_part(), cosFeatureAngle);
    return measure(myLSFields[0].identifier);
  }

  DecompositionMeasures measure(const Surface_Identifier surfaceIdentifier) const
  {
    const stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    const stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const AuxMetaData & auxMeta = AuxMetaData::get(meta);
    const Phase_Support & phaseSupport = Phase_Support::get(meta);

    const stk::mesh::Selector ownedActive = meta.locally_owned_part() & auxMeta.active_part();
    DecompositionMeasures measures;
    measures.negativeVolume = compute_volume_or_surface_area(mesh,
        stk::topology::ELEMENT_RANK,
        phaseSupport.get_negative_levelset_block_selector(surfaceIdentifier) & ownedActive);
    all_reduce_sum(mesh.parallel(), measures.negativeVolume);
    measures.interfaceArea = compute_volume_or_surface_area(mesh,
        meta.side_rank(),
        phaseSupport.get_negative_levelset_interface_selector(surfaceIdentifier) & ownedActive);
    all_reduce_sum(mesh.parallel(), measures.interfaceArea);
    measures.worstQuality = compute_mesh_quality(
        mesh, auxMeta.active_part(), meta.coordinate_field(), MeanRatioQualityMetric());
    measures.worstScaledJacobian = compute_mesh_quality(
        mesh, auxMeta.active_part(), meta.coordinate_field(), ScaledJacobianQualityMetric());
    measures.numElementsBelowLowMeanRatios = count_elements_below_low_mean_ratios(mesh, ownedActive);
    std::vector<size_t> counts;
    const stk::mesh::Selector activeSelector = auxMeta.active_part();
    stk::mesh::comm_mesh_counts(mesh, counts, &activeSelector);
    measures.numElements = counts[stk::topology::ELEMENT_RANK];
    return measures;
  }

  // The sharp features captured only by snapping to them (CDFEM_Support).
  void capture_sharp_features_only_during_snapping()
  {
    CDFEM_Support::get(myBBoxMesh.meta_data()).set_capture_sharp_features_only_during_snapping(true);
  }

  // Whether CDFEM would decompose the mesh again for this surface, as an application asks before
  // each step.
  bool decomposition_needs_update(const Surface & surface)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const AnalyticSurfaceInterfaceGeometry geometry({myLSFields[0].identifier},
        {&surface},
        AuxMetaData::get(meta).active_part(),
        CDFEM_Support::get(meta),
        Phase_Support::get(meta));
    return CDMesh::decomposition_needs_update(myBBoxMesh.bulk_data(), geometry);
  }

  // Writes the active elements, composed into one file, as file fileIndex of the series fileBaseName
  // (fileBaseName, then fileBaseName-s0002, ...), at the given time.  On more than one processor
  // this needs parallel IO.
  void write_active_mesh(const std::string & fileBaseName, const int fileIndex, const double time) const
  {
    output_composed_mesh_with_fields(myBBoxMesh.bulk_data(),
        AuxMetaData::get(myBBoxMesh.meta_data()).active_part(),
        create_file_name(fileBaseName, fileIndex),
        1,
        time);
  }

private:
  BoundingBoxMesh myBBoxMesh;
  std::vector<LS_Field> myLSFields;
};

} // namespace krino

#endif /* KRINO_KRINO_DECOMP_UNIT_TESTS_AKRI_UNIT_SURFACEDECOMPOSITION_HPP_ */
