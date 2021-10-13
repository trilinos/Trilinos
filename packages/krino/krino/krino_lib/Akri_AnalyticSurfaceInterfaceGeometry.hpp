// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_AnalyticSurfaceInterfaceGeometry_h
#define Akri_AnalyticSurfaceInterfaceGeometry_h

#include <Akri_InterfaceID.hpp>
#include <Akri_Surface.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"

namespace krino {

class CDFEM_Support;
class Phase_Support;

class SurfaceElementCutter : public ElementCutter
{
public:
  SurfaceElementCutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const Surface & surface);
  virtual ~SurfaceElementCutter() {}

  virtual bool might_have_interior_or_face_intersections() const override { return false; }
  virtual void fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const override {}
  virtual std::vector<InterfaceID> get_sorted_cutting_interfaces() const override
    { std::vector<InterfaceID> interfaces; interfaces.push_back(InterfaceID(0,0)); return interfaces; }
  virtual std::vector<int> get_interface_signs_based_on_crossings(const std::vector<Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const override
    { std::vector<int> interfaceSigns(1,0); return interfaceSigns; }
  virtual void fill_tetrahedron_face_interior_intersections(const std::array<Vector3d,3> & faceNodes,
    const InterfaceID & interface1,
    const InterfaceID & interface2,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const override {}
  virtual std::string visualize(const stk::mesh::BulkData & mesh) const override { std::string empty; return empty; }
  virtual bool have_crossing(const InterfaceID interface, const Segment3d & edge) const override;
  virtual double interface_crossing_position(const InterfaceID interface, const Segment3d & edge) const override;
  virtual int sign_at_position(const InterfaceID interface, const Vector3d & paramCoords) const override;
  virtual int get_starting_phase_for_cutting_surfaces() const override { return 0; }

private:
  Vector3d parametric_to_global_coordinates(const Vector3d & pCoords) const;

  const MasterElement & myMasterElem;
  std::vector<Vector3d> myElementNodeCoords;
  const Surface & mySurface;
};

class AnalyticSurfaceInterfaceGeometry : public InterfaceGeometry {

public:
  AnalyticSurfaceInterfaceGeometry(const Surface & surface,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
: mySurface(surface),
  myActivePart(activePart),
  myCdfemSupport(cdfemSupport),
  myPhaseSupport(phaseSupport) {}

  virtual ~AnalyticSurfaceInterfaceGeometry() {}

  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;
  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override;

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
      const NodeToCapturedDomainsMap & nodesToSnappedDomains) const override {}
  virtual const ElementToDomainMap & get_phase_for_uncut_elements() const override { return myUncutElementPhases; }

  virtual std::unique_ptr<ElementCutter> build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const override;

  virtual PhaseTag get_starting_phase(const ElementCutter * cutter) const override;

private:
  const Surface & mySurface;
  const stk::mesh::Part & myActivePart;
  const CDFEM_Support & myCdfemSupport;
  const Phase_Support & myPhaseSupport;
  ElementToDomainMap myUncutElementPhases;
  mutable std::vector<stk::mesh::Entity> myElementsToIntersect;
};

} // namespace krino

#endif
