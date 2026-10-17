// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_CurvatureErrorEstimator.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_ErrorEstimateUtils.hpp>

#include <Akri_MeshHelpers.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_NodeMarker.hpp>
#include <Akri_RefinementSupport.hpp>
#include <algorithm>
#include <numbers>

namespace krino {

namespace {

// The interface normal, for the surface, of each selected element that has a marked node; other entries are left
// unchanged.
void fill_interface_normals_of_elements_with_marked_nodes(const stk::mesh::BulkData& mesh,
    const InterfaceGeometry & interfaceGeometry,
    const Surface_Identifier surfId,
    const stk::mesh::Selector & elemSelector,
    const FieldRef nodeMarkerField,
    std::vector<stk::math::Vector3d> & interfaceNormals)
{
  for (auto && bucket : mesh.buckets(stk::topology::ELEMENT_RANK))
    if (elemSelector(bucket))
      for (auto && elem : *bucket)
        if (element_has_marked_node(mesh, elem, nodeMarkerField))
          interfaceNormals[elem.local_offset()] = interfaceGeometry.compute_interface_normal(mesh, surfId, elem);
}

}

std::unique_ptr<CurvatureErrorEstimator> build_curvature_error_estimator(const RefinementSupport & refinementSupport, const InterfaceGeometry & interfaceGeometry, const stk::diag::Timer & parentTimer)
{
  if (refinementSupport.get_interface_refinement_chordal_error_tolerance() > 0.)
  {
    const double chordalErrorTol = refinementSupport.get_interface_refinement_chordal_error_tolerance();
    krinolog << "Refining interface elements that have a chordal error estimate > " << chordalErrorTol << ".\n";
    return std::make_unique<FacetChordalErrorFromGradientJump>(interfaceGeometry, chordalErrorTol, parentTimer);
  }
  else if (refinementSupport.get_interface_refinement_curvature_tolerance() > 0.)
  {
    if (dynamic_cast<const AnalyticSurfaceInterfaceGeometry *>(&interfaceGeometry) == nullptr)
    {
      std::ostringstream err;
      err << "interface_refinement_curvature_tolerance is only supported for bounding surfaces.\n"
          << "With level sets, use either interface_refinement_curvature_angle_tolerance "
          << "or interface_refinement_chordal_error_tolerance.\n";
      throw std::runtime_error(err.str());
    }
    const double curvatureTol = refinementSupport.get_interface_refinement_curvature_tolerance();
    krinolog << "Refining interface elements that have estimated curvature time element size > " << curvatureTol << " (possibly expensive).\n";
    return std::make_unique<CurvatureTimesElementSizeError>(interfaceGeometry, curvatureTol, parentTimer);
  }
  const double surfaceAngleTolDeg = refinementSupport.get_interface_refinement_curvature_angle_tolerance();
  krinolog << "Refining interface elements that have estimated angle between interface normals > " << surfaceAngleTolDeg << "(deg).\n";
  const double surfaceAngleTolRadians = surfaceAngleTolDeg * std::numbers::pi / 180.;
  return std::make_unique<FacetAngleCurvatureError>(interfaceGeometry, surfaceAngleTolRadians, parentTimer);
}

double CurvatureTimesElementSizeError::estimate_element_error(const stk::mesh::BulkData& mesh, const stk::mesh::Entity elem) const
{
  stk::diag::TimeBlock timer__(myTimer);
  return myInterfaceGeometry.estimate_curvature_times_element_size(mesh, elem);
};

void FacetAngleCurvatureError::precompute(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const FieldRef nodeMarkerField)
{
  STK_ThrowRequireMsg(mesh.is_automatic_aura_on() || mesh.parallel_size() == 1, "Method requires automatic aura.");
  stk::mesh::communicate_field_data(mesh, {&nodeMarkerField.field()}); // need node marker on ghost nodes
  stk::diag::TimeBlock timer__(myTimer);

  std::vector<stk::math::Vector3d> interfaceNormals = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, stk::math::Vector3d::ZERO);
  myErrorByElemOffset.assign(interfaceNormals.size(), 0.);

  for (auto & surfId : myInterfaceGeometry.get_surface_identifiers())
  {
    std::fill(interfaceNormals.begin(), interfaceNormals.end(), stk::math::Vector3d::ZERO);
    fill_interface_normals_of_elements_with_marked_nodes(mesh, myInterfaceGeometry, surfId, elemSelector, nodeMarkerField, interfaceNormals);

    accumulate_facet_angle_errors_from_interface_normals(mesh, elemSelector, interfaceNormals, myErrorByElemOffset);
  }
}

double FacetAngleCurvatureError::estimate_element_error(const stk::mesh::BulkData&, const stk::mesh::Entity elem) const
{
  return myErrorByElemOffset[elem.local_offset()];
};

void FacetChordalErrorFromGradientJump::precompute(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const FieldRef nodeMarkerField)
{
  STK_ThrowRequireMsg(mesh.is_automatic_aura_on() || mesh.parallel_size() == 1, "Method requires automatic aura.");
  stk::mesh::communicate_field_data(mesh, {&nodeMarkerField.field()}); // need node marker on ghost nodes
  stk::diag::TimeBlock timer__(myTimer);

  std::vector<stk::math::Vector3d> interfaceNormals = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, stk::math::Vector3d::ZERO);
  myErrorByElemOffset.assign(interfaceNormals.size(), 0.);

  for (auto & surfId : myInterfaceGeometry.get_surface_identifiers())
  {
    std::fill(interfaceNormals.begin(), interfaceNormals.end(), stk::math::Vector3d::ZERO);
    fill_interface_normals_of_elements_with_marked_nodes(mesh, myInterfaceGeometry, surfId, elemSelector, nodeMarkerField, interfaceNormals);

    accumulate_facet_chordal_errors_from_interface_normals(mesh, coordsField, elemSelector, interfaceNormals, myErrorByElemOffset);
  }
}

double FacetChordalErrorFromGradientJump::estimate_element_error(const stk::mesh::BulkData&, const stk::mesh::Entity elem) const
{
  return myErrorByElemOffset[elem.local_offset()];
};

}


