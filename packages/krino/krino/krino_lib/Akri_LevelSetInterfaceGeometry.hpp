// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_LEVELSETINTERFACEGEOMETRY_H_
#define KRINO_INCLUDE_AKRI_LEVELSETINTERFACEGEOMETRY_H_
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_ParentsToChildMapper.hpp>
#include <stk_mesh/base/Part.hpp>
#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"

namespace krino {

class CDFEM_Support;
class Phase_Support;

class LevelSetElementCutter : public ElementCutter
{
public:
  LevelSetElementCutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const ParentEdgeMap & parentEdges,
    const Phase_Support & phaseSupport,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker);
  virtual ~LevelSetElementCutter() {}

  virtual bool might_have_interior_or_face_intersections() const override
    { return myElementInterfaceCutter->get_num_cutting_surfaces() > 1; }
  virtual void fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const override
    { myElementInterfaceCutter->fill_interior_intersections(intersectionPointFilter, intersections); }
  virtual std::vector<InterfaceID> get_sorted_cutting_interfaces() const override
    { std::vector<InterfaceID> interfaces; myElementInterfaceCutter->fill_interfaces_with_cutting_surface(interfaces); return interfaces; }
  virtual std::vector<int> get_interface_signs_based_on_crossings(const std::vector<Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const override;
  virtual void fill_tetrahedron_face_interior_intersections(const std::array<Vector3d,3> & faceNodes,
    const InterfaceID & interface1,
    const InterfaceID & interface2,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const override
    { myElementInterfaceCutter->fill_tetrahedron_face_interior_intersections(faceNodes, interface1, interface2, intersectionPointFilter, intersections); }
  virtual std::string visualize(const stk::mesh::BulkData & mesh) const override;
  virtual bool have_crossing(const InterfaceID interface, const Segment3d & edge) const override
    { return myElementInterfaceCutter->have_crossing(interface, edge); }
  virtual double interface_crossing_position(const InterfaceID interface, const Segment3d & edge) const override
    { return myElementInterfaceCutter->interface_crossing_position(interface, edge); }
  virtual int sign_at_position(const InterfaceID interface, const Vector3d & paramCoords) const override
    { return myElementInterfaceCutter->sign_at_position(interface, paramCoords); }
  virtual int get_starting_phase_for_cutting_surfaces() const override
    { return myElementInterfaceCutter->get_starting_phase_for_cutting_surfaces(); }

  PhaseTag get_starting_phase(const CDFEM_Support & cdfemSupport, const Phase_Support & phaseSupport) const;
  void update_edge_crossings(const unsigned iEdge, const std::vector<std::vector<double>> & nodesIsovar);
  void add_interfaces_with_uncaptured_intersection_within_element(const std::vector<Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings) const
    { return myElementInterfaceCutter->add_interfaces_with_uncaptured_intersection_within_element(elemNodesCoords, elemNodesSnappedDomains, interfacesWithUncapturedCrossings); }

private:
  std::vector<const CDFEM_Parent_Edge *> myParentEdges;
  std::vector<bool> myParentEdgesAreOrientedSameAsElementEdges;
  std::unique_ptr<Element_Cutter> myElementInterfaceCutter;
};

class LevelSetInterfaceGeometry : public InterfaceGeometry {

public:
  LevelSetInterfaceGeometry(const stk::mesh::MetaData & meta);
  LevelSetInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
: myActivePart(activePart),
  myCdfemSupport(cdfemSupport),
  myPhaseSupport(phaseSupport) {}
  virtual ~LevelSetInterfaceGeometry() {}

  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;

  virtual std::vector<IntersectionPoint> get_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const override;

  // FIXME: Temporary methods
  virtual void store_phase_for_uncut_elements(const stk::mesh::BulkData & mesh) const override;
  virtual void store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & independentSnapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;

  virtual const ElementToDomainMap & get_phase_for_uncut_elements() const override { return myUncutElementPhases; }

  virtual std::unique_ptr<ElementCutter> build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const override;

  virtual PhaseTag get_starting_phase(const ElementCutter * cutter) const override;

private:
  const stk::mesh::Part & myActivePart;
  const CDFEM_Support & myCdfemSupport;
  const Phase_Support & myPhaseSupport;
  mutable ParentEdgeMap myParentEdges;
  mutable ParentsToChildMapper myParentsToChildMapper;
  mutable ElementToDomainMap myUncutElementPhases;
};

}


#endif /* KRINO_INCLUDE_AKRI_LEVELSETINTERFACEGEOMETRY_H_ */
