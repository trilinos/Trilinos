// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ELEMENT_CUTTER_H_
#define KRINO_INCLUDE_AKRI_ELEMENT_CUTTER_H_
#include <Akri_InterfaceID.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_Segment.hpp>
#include <stk_topology/topology.hpp>
#include <array>
#include <map>
#include <memory>
#include <functional>

namespace krino {

class Cutting_Surface;
class CDFEM_Parent_Edge;
struct ElementIntersection;

typedef std::function<bool(const std::vector<int> &)> ElementIntersectionPointFilter;
typedef std::map<InterfaceID, std::shared_ptr<Cutting_Surface>> InterfaceToSurface;

class Element_Cutter
{
public:
  Element_Cutter(const MasterElement & masterElement) : myTopology(masterElement.get_topology()) {}
  virtual ~Element_Cutter() {}
  double interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const;
  int sign_at_position(const InterfaceID interface, const stk::math::Vector3d & p_coords) const;
  virtual bool have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const = 0;
  virtual int get_ls_per_interface_phase_at_location(const stk::math::Vector3d & pCoords) const = 0;
  virtual void add_interfaces_with_uncaptured_intersection_within_element(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings) const = 0;

  int get_num_cutting_surfaces() const {return cutting_surfaces.size(); }
  bool have_cutting_surface(const InterfaceID interface) const;
  void fill_interfaces_with_cutting_surface(std::vector<InterfaceID> & interfaces) const;
  Cutting_Surface * get_cutting_surface(const InterfaceID interface) const;

  virtual bool is_one_ls_per_phase() const = 0;
  void fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const;
  void fill_interior_intersections(std::vector<ElementIntersection> & intersections) const;
  void append_tetrahedron_face_interior_intersections(const std::array<stk::math::Vector3d,3> & faceNodes, const InterfaceID & interface1, const InterfaceID & interface2, const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const;
  static std::string visualize_cutting_surfaces(const stk::topology & topology, const InterfaceToSurface & cuttingSurfaces);
  virtual std::string visualize() const { return visualize_cutting_surfaces(myTopology, cutting_surfaces); }

  int get_starting_phase_for_cutting_surfaces() const;

  struct Edge_Crossing {
    Edge_Crossing(unsigned e, double p, int s) : edge(e), pos(p), sign(s) {}
    unsigned edge;
    double pos;
    int sign;
  };

  static constexpr double theSnapToNodeTol{1.e-12};

protected:
  stk::topology myTopology;
  InterfaceToSurface cutting_surfaces;

private:
  virtual bool intersection_point_is_real(const stk::math::Vector3d & intersection, const std::vector<int> & sortedDomains) const = 0;
  void fill_triangle_interior_intersection_points(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const;
  void fill_tetrahedron_interior_intersection_points(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const;
};

class LS_Per_Interface_Cutter : public Element_Cutter
{
public:
  LS_Per_Interface_Cutter(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker);
  virtual bool is_one_ls_per_phase() const override { return false; }
  virtual bool have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const override;
  virtual int get_ls_per_interface_phase_at_location(const stk::math::Vector3d & pCoords) const override;
  virtual void add_interfaces_with_uncaptured_intersection_within_element(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings) const override;
private:
  virtual bool intersection_point_is_real(const stk::math::Vector3d & /*intersection*/, const std::vector<int> & /*sortedDomains*/) const override { return true; }
  static std::vector<Edge_Crossing> build_cut_edges(const InterfaceID & interface,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & parentEdgeIsOrientedSameAsElementEdge);
};

class One_LS_Per_Phase_Cutter : public Element_Cutter
{
public:
  One_LS_Per_Phase_Cutter(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker);
  virtual bool is_one_ls_per_phase() const override { return true; }
  virtual bool have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const override;
  virtual int get_ls_per_interface_phase_at_location(const stk::math::Vector3d & pCoords) const override;
  virtual void add_interfaces_with_uncaptured_intersection_within_element(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings) const override;
  virtual std::string visualize() const override { return visualize_cutting_surfaces(myTopology, all_cutting_surfaces); }

private:
  virtual bool intersection_point_is_real(const stk::math::Vector3d & intersection, const std::vector<int> & sortedDomains) const override;
  static std::vector<Edge_Crossing> build_cut_edges(const InterfaceID & interface,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & parentEdgeIsOrientedSameAsElementEdge);
  static std::shared_ptr<Cutting_Surface> attempt_to_build_cutting_surface(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker,
    const InterfaceID& interface);

  InterfaceToSurface all_cutting_surfaces;
};

std::unique_ptr<Element_Cutter> create_element_cutter(const bool oneLSPerPhase,
    const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker);

}

#endif /* KRINO_INCLUDE_AKRI_ELEMENT_CUTTER_H_ */
