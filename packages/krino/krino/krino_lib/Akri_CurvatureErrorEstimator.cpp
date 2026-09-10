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

#include <Akri_MeshHelpers.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_NodeMarker.hpp>
#include <Akri_RefinementSupport.hpp>
#include <numbers>

namespace krino {

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
  myErrorByElemOffset = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, 0.);
  std::vector<stk::mesh::Entity> sideNbrs;

  for (auto & surfId : myInterfaceGeometry.get_surface_identifiers())
  {
    for (auto && bucket : mesh.buckets(stk::topology::ELEMENT_RANK))
      if (elemSelector(bucket))
        for (auto && elem : *bucket)
          if (element_has_marked_node(mesh, elem, nodeMarkerField))
            interfaceNormals[elem.local_offset()] = myInterfaceGeometry.compute_interface_normal(mesh, surfId, elem);

    for (auto * bucket : mesh.buckets(stk::topology::ELEMENT_RANK))
    {
      if (bucket->owned() && elemSelector(bucket))
      {
        for (auto && elem : *bucket)
        {
          const stk::math::Vector3d & elemNormal = interfaceNormals[elem.local_offset()];
          if (!elemNormal.zero_length())
          {
            double nbrSum = 0.;
            unsigned nbrCount = 0;
            fill_selected_side_attached_elements(mesh, elemSelector, elem, sideNbrs);
            for (auto & nbr : sideNbrs)
            {
              const stk::math::Vector3d & nbrNormal = interfaceNormals[nbr.local_offset()];
              if (!nbrNormal.zero_length())
              {
                nbrSum += std::acos(Dot(elemNormal,nbrNormal));
                nbrCount += 1;
              }
            }
            // average face contribution over neighbors
            myErrorByElemOffset[elem.local_offset()] = nbrSum/nbrCount;
          }
        }
      }
    }
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
  /*
  For two adjacent equilateral triangular facets whose vertices lie on a sphere, let
      d = n1 · n2
  be the dot product of the two outward-facing unit normals.  The normalized
  facet chordal error, err/h, where err is the maximum deviation
  between the spherical surface and a facet and h is the
  facet edge length, can be written directly in terms of d as
      err/h = (sqrt(5 - 3*d) - sqrt(1 + d)) / sqrt(12 * (1 - d))
  The nearly coplanar, small-facet limit corresponds to d -> 1.  In that case,
      err/h ≈ sqrt((1 - d) / 6)
  so the chordal error relative to edge length decreases like the square root
  of the normal mismatch 1 - d.

  For two adjacent equal-length line segments whose endpoints lie on a circle,
  the normalized chordal error, err/h, where epsilon is the maximum
  sagitta deviation between the circular arc and one straight segment and h is
  the segment chord length, can be written directly in terms of d as
      err/h = (1 - sqrt((1 + d) / 2)) / sqrt(2 * (1 - d))
  The nearly collinear, small-segment limit corresponds to d -> 1.  In that
      err/h ≈ sqrt((1 - d) / 32)
  */

  STK_ThrowRequireMsg(mesh.is_automatic_aura_on() || mesh.parallel_size() == 1, "Method requires automatic aura.");
  stk::mesh::communicate_field_data(mesh, {&nodeMarkerField.field()}); // need node marker on ghost nodes
  stk::diag::TimeBlock timer__(myTimer);
  std::vector<stk::math::Vector3d> elemNodesCoords;

  std::vector<stk::math::Vector3d> interfaceNormals = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, stk::math::Vector3d::ZERO);
  myErrorByElemOffset = create_vector_indexable_by_entity_offset(mesh, stk::topology::ELEMENT_RANK, 0.);
  std::vector<stk::mesh::Entity> sideNbrs;
  const double coeff = mesh.mesh_meta_data().spatial_dimension() == 2 ? std::sqrt(1./32.) : std::sqrt(1./6.);

  for (auto & surfId : myInterfaceGeometry.get_surface_identifiers())
  {
    for (auto * bucket : mesh.buckets(stk::topology::ELEMENT_RANK))
      if (elemSelector(bucket))
        for (auto && elem : *bucket)
          if (element_has_marked_node(mesh, elem, nodeMarkerField))
            interfaceNormals[elem.local_offset()] = myInterfaceGeometry.compute_interface_normal(mesh, surfId, elem);

    for (auto && bucket : mesh.buckets(stk::topology::ELEMENT_RANK))
    {
      if (bucket->owned() && elemSelector(bucket))
      {
        for (auto && elem : *bucket)
        {
          const stk::math::Vector3d & elemNormal = interfaceNormals[elem.local_offset()];
          if (!elemNormal.zero_length())
          {
            fill_element_node_coordinates(mesh, elem, coordsField, elemNodesCoords);
            const double elemSize = compute_simplex_RMS_edge_lengths(elemNodesCoords);
            double nbrSum = 0.;
            unsigned nbrCount = 0;

            fill_selected_side_attached_elements(mesh, elemSelector, elem, sideNbrs);
            for (auto & nbr : sideNbrs)
            {
              const stk::math::Vector3d & nbrNormal = interfaceNormals[nbr.local_offset()];
              if (!nbrNormal.zero_length())
              {
                nbrSum += std::sqrt((1.-Dot(elemNormal,nbrNormal)));
                nbrCount += 1;
              }
            }
            // average face contribution over neighbors
            myErrorByElemOffset[elem.local_offset()] = elemSize*coeff*nbrSum/nbrCount;
          }
        }
      }
    }
  }
}

double FacetChordalErrorFromGradientJump::estimate_element_error(const stk::mesh::BulkData&, const stk::mesh::Entity elem) const
{
  return myErrorByElemOffset[elem.local_offset()];
};

}


