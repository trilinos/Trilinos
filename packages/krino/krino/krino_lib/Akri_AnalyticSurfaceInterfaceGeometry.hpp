// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_AnalyticSurfaceInterfaceGeometry_h
#define Akri_AnalyticSurfaceInterfaceGeometry_h

#include <Akri_DetermineElementSign.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_Surface.hpp>
#include <Akri_Surface_Identifier.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace krino {

class CDFEM_Support;
class Phase_Support;
class FieldRef;

class SurfaceElementCutter : public ElementCutter
{
public:
  SurfaceElementCutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::vector<const Surface *> & surfaces,
    const std::vector<int8_t> & elementSigns,
    const double edgeTol);
  virtual ~SurfaceElementCutter() {}

  virtual bool might_have_interior_or_face_intersections() const override { return false; }
  virtual void fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const override {}
  virtual std::vector<InterfaceID> get_sorted_cutting_interfaces() const override;
  virtual std::vector<int> get_interface_signs_based_on_crossings(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const override;
  virtual void fill_tetrahedron_face_interior_intersections(const std::array<stk::math::Vector3d,3> & faceNodes,
    const InterfaceID & interface1,
    const InterfaceID & interface2,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const override {}
  virtual std::string visualize(const stk::mesh::BulkData & mesh) const override { std::string empty; return empty; }
  virtual int interface_sign_for_uncrossed_element(const InterfaceID interface, const std::vector<stk::math::Vector3d> & elemNodesCoords) const override;
  virtual std::pair<int, double> interface_edge_crossing_sign_and_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const override;
  virtual int get_starting_phase_for_cutting_surfaces() const override { return 0; }

  const std::vector<int8_t> & get_element_signs() const { return myElementSigns; }

private:
  const Surface & get_surface(const InterfaceID interface) const;
  stk::math::Vector3d parametric_to_global_coordinates(const stk::math::Vector3d & pCoords) const;

  const MasterElement & myMasterElem;
  std::vector<stk::math::Vector3d> myElementNodeCoords;
  const std::vector<const Surface*> & mySurfaces;
  std::vector<int8_t> myElementSigns;
  double myEdgeCrossingTol;
};

class AnalyticSurfaceInterfaceGeometry : public InterfaceGeometry {

public:
  AnalyticSurfaceInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport);
  AnalyticSurfaceInterfaceGeometry(const std::vector<Surface_Identifier> & surfaceIdentifiers,
    const std::vector<const Surface*> & surfaces,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport);

  virtual ~AnalyticSurfaceInterfaceGeometry() {}

  void add_surface(const Surface_Identifier surfaceIdentifier, const Surface & surface);

  virtual bool might_have_interior_or_face_intersections() const override { return true; }

  virtual void prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh) const override;
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;

  virtual std::vector<stk::mesh::Entity> get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const override;
  virtual void fill_elements_that_intersect_distance_interval(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval) const override;

  virtual bool snapped_elements_may_have_new_intersections() const override { return false; }

  virtual std::vector<IntersectionPoint> get_edge_intersection_points(const stk::mesh::BulkData & mesh,
      const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;
  virtual void append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const override;

  virtual void store_phase_for_uncut_elements(const stk::mesh::BulkData & mesh) const override {}
  virtual void store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & independentSnapInfos,
      const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;
  virtual const ElementToDomainMap & get_phase_for_uncut_elements() const override { return myUncutElementPhases; }

  virtual std::unique_ptr<ElementCutter> build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const override;

  virtual PhaseTag get_starting_phase(const ElementCutter * cutter) const override;

  const std::vector<Surface_Identifier> & get_surface_identifiers() const override { return mySurfaceIdentifiers; }
  FieldRef get_coordinates_field(const stk::mesh::BulkData & mesh) const;

protected:
  const stk::mesh::Part & get_active_part() const { return myActivePart; }
  const CDFEM_Support & get_cdfem_support() const { return myCdfemSupport; }
  const Phase_Support & get_phase_support() const { return myPhaseSupport; }
  void set_elements_to_intersect_and_prepare_to_compute_with_surfaces(const stk::mesh::BulkData & mesh,
      const std::vector<stk::mesh::Entity> & elementsToIntersect) const;
  void set_element_signs(const stk::mesh::BulkData & mesh,
      const std::vector<stk::mesh::Selector> & perSurfaceElementSelector) const;
  stk::mesh::Selector get_mesh_parent_element_selector() const;
  std::vector<stk::mesh::Entity> get_mesh_parent_elements(const stk::mesh::BulkData & mesh) const;

private:
  const Surface & get_surface_with_identifer(const Surface_Identifier surfaceIdentifier) const;

  std::vector<const Surface*> mySurfaces;
  const stk::mesh::Part & myActivePart;
  const CDFEM_Support & myCdfemSupport;
  const Phase_Support & myPhaseSupport;
  std::vector<Surface_Identifier> mySurfaceIdentifiers;
  double myEdgeCrossingTol;
  mutable ElementToSignsMap myElementsToSigns;
  mutable ElementToDomainMap myUncutElementPhases;
  mutable std::vector<stk::mesh::Entity> myElementsToIntersect;
};

} // namespace krino

#endif
