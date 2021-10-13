// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_INTERFACEGEOMETRY_H_
#define AKRI_INTERFACEGEOMETRY_H_
#include <Akri_Element_Cutter.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_PhaseTag.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <string>

namespace krino {

class SnapInfo;

class ElementCutter
{
public:
  virtual ~ElementCutter() {}
  virtual std::vector<InterfaceID> get_sorted_cutting_interfaces() const = 0;
  virtual std::vector<int> get_interface_signs_based_on_crossings(const std::vector<Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const = 0;
  virtual void fill_tetrahedron_face_interior_intersections(const std::array<Vector3d,3> & faceNodes,
    const InterfaceID & interface1,
    const InterfaceID & interface2,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const = 0;
  virtual bool might_have_interior_or_face_intersections() const = 0;
  virtual void fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const = 0;
  virtual std::string visualize(const stk::mesh::BulkData & mesh) const = 0;
  virtual bool have_crossing(const InterfaceID interface, const Segment3d & edge) const = 0;
  virtual double interface_crossing_position(const InterfaceID interface, const Segment3d & edge) const = 0;
  virtual int sign_at_position(const InterfaceID interface, const Vector3d & paramCoords) const = 0;
  virtual int get_starting_phase_for_cutting_surfaces() const = 0;
};

class InterfaceGeometry {
public:
  InterfaceGeometry() {}

  virtual ~InterfaceGeometry() {}
  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh, const NodeToCapturedDomainsMap & nodesToCapturedDomains) const = 0;
  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const = 0;

  virtual std::vector<IntersectionPoint> get_edge_intersection_points(const stk::mesh::BulkData & mesh,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const = 0;

  virtual void append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const = 0;

  // FIXME: Temporary methods
  virtual void store_phase_for_uncut_elements(const stk::mesh::BulkData & mesh) const = 0;
  virtual void store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & independentSnapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const = 0;

  virtual const ElementToDomainMap & get_phase_for_uncut_elements() const = 0;

  virtual std::unique_ptr<ElementCutter> build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const = 0;

  virtual PhaseTag get_starting_phase(const ElementCutter * cutter) const = 0;
};
}

#endif // AKRI_INTERFACEGEOMETRY_H_
