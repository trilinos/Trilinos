// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ErrorEstimateUtils.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <cmath>

namespace krino {

// For each owned selected element with an interface normal, adds the mean angle between its normal and the normals of
// its selected side neighbors that have one.
void accumulate_facet_angle_errors_from_interface_normals(const stk::mesh::BulkData& mesh,
    const stk::mesh::Selector & elemSelector,
    const std::vector<stk::math::Vector3d> & interfaceNormals,
    std::vector<double> & errorByElemOffset)
{
  std::vector<stk::mesh::Entity> sideNbrs;
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
          // average face contribution over neighbors, summed over surfaces
          if (nbrCount > 0)
            errorByElemOffset[elem.local_offset()] += nbrSum/nbrCount;
        }
      }
    }
  }
}

// For each owned selected element with an interface normal, adds the chordal error estimated from the jump between its
// normal and the normals of its selected side neighbors that have one, scaled by the element size.
void accumulate_facet_chordal_errors_from_interface_normals(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const std::vector<stk::math::Vector3d> & interfaceNormals,
    std::vector<double> & errorByElemOffset)
{
  std::vector<stk::mesh::Entity> sideNbrs;
  std::vector<stk::math::Vector3d> elemNodesCoords;
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
  const double coeff = mesh.mesh_meta_data().spatial_dimension() == 2 ? std::sqrt(1./32.) : std::sqrt(1./6.);
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
          // average face contribution over neighbors, summed over surfaces
          if (nbrCount > 0)
            errorByElemOffset[elem.local_offset()] += elemSize*coeff*nbrSum/nbrCount;
        }
      }
    }
  }
}

} // namespace krino
