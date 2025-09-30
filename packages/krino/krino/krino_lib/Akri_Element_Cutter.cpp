// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_Cutting_Surface.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Plane_Intersections.hpp>
#include "Akri_MasterElementDeterminer.hpp"

namespace krino {

bool all_interfaces_have_single_level_set(const std::vector<InterfaceID> &interfaces)
{
  for (auto && interface : interfaces)
    if (interface.first_ls() != interface.second_ls())
      return false;
  return true;
}

void fill_sorted_domains(const bool isOneLSPerPhase,const InterfaceID & interface0, const InterfaceID & interface1, std::vector<int> & sortedDomains)
{
  sortedDomains.clear();
  if (isOneLSPerPhase)
  {

    sortedDomains = {interface0.first_ls(), interface0.second_ls(), interface1.first_ls(), interface1.second_ls()};
    stk::util::sort_and_unique(sortedDomains);
  }
  else
  {
    STK_ThrowAssert(all_interfaces_have_single_level_set({interface0, interface1}));
    sortedDomains = {interface0.first_ls(), interface1.first_ls()};
  }
}

void fill_sorted_domains(const bool isOneLSPerPhase,const InterfaceID & interface0, const InterfaceID & interface1, const InterfaceID & interface2, std::vector<int> & sortedDomains)
{
  sortedDomains.clear();
  if (isOneLSPerPhase)
  {

    sortedDomains =  {interface0.first_ls(), interface0.second_ls(), interface1.first_ls(), interface1.second_ls(), interface2.first_ls(), interface2.second_ls()};
    stk::util::sort_and_unique(sortedDomains);
  }
  else
  {
    STK_ThrowAssert(all_interfaces_have_single_level_set({interface0, interface1, interface2}));
    sortedDomains = {interface0.first_ls(), interface1.first_ls(), interface2.first_ls()};
  }
}

bool sorted_domains_form_triple_point(const bool is_one_interface_per_phase, const std::vector<int> & sortedDomains)
{
  return !is_one_interface_per_phase || sortedDomains.size() == 3;
}

bool sorted_domains_form_quadruple_point(const bool is_one_interface_per_phase, const std::vector<int> & sortedDomains)
{
  return !is_one_interface_per_phase || sortedDomains.size() == 4;
}

void fill_tetrahedron_face_surface_intersection_points(const int iFace, const Cutting_Surface * surf0, const Cutting_Surface * surf1, const std::vector<int> & sortedDomains, std::vector<ElementIntersection> & intersections)
{
  stk::math::Vector3d intersectionPoint;
  if (find_intersection_of_two_planes_and_side_of_tet(iFace, surf0->get_plane(), surf1->get_plane(), intersectionPoint))
    intersections.emplace_back(intersectionPoint, sortedDomains);
}

void fill_tetrahedron_interior_surface_intersection_points(const Cutting_Surface * surf0, const Cutting_Surface * surf1, const Cutting_Surface * surf2, const std::vector<int> & sortedDomains, std::vector<ElementIntersection> & intersections)
{
  stk::math::Vector3d intersectionPoint;
  if (find_intersection_of_three_planes_within_tet(surf0->get_plane(), surf1->get_plane(), surf2->get_plane(), intersectionPoint))
    intersections.emplace_back(intersectionPoint, sortedDomains);
}

std::string visualize_element(const stk::topology & topology, std::array<int,4> & usedGeomIds)
{
  const auto & masterElem = MasterElementDeterminer::getMasterElement(topology);
  const double * elemNodeParamCoords = masterElem.nodal_parametric_coordinates();

  const unsigned numNodes = topology.num_nodes();
  const int dim = topology.dimension();

  std::ostringstream os;
  for (unsigned n=0; n<numNodes; ++n)
  {
    stk::math::Vector3d nodeCoord(&elemNodeParamCoords[n*dim],dim);
    os << "Create vertex " << nodeCoord[0] << " " << nodeCoord[1] << " " << nodeCoord[2] << std::endl;
    usedGeomIds[0]++;
    os << "Create node " << usedGeomIds[0] << " vertex " << usedGeomIds[0] << std::endl;
  }
  const int nodeIds = usedGeomIds[0]-numNodes+1;
  if (numNodes == 3)
    os << "Create tri node " << nodeIds << " " << nodeIds+1 << " " << nodeIds+2 << std::endl;
  else
    os << "Create tet node " << nodeIds << " " << nodeIds+1 << " " << nodeIds+2 << " " << nodeIds+3 << std::endl;
  return os.str();
}

std::string Element_Cutter::visualize_cutting_surfaces(const stk::topology & topology, const InterfaceToSurface & cuttingSurfaces)
{
  std::array<int,4> usedGeomIds{0,0,0,0};
  std::string viz = visualize_element(topology, usedGeomIds);
  for (auto && cuttingSurf : cuttingSurfaces)
  {
    std::ostringstream os;
    os << cuttingSurf.first;
    const std::string interfaceDesc = "on interface " + os.str();
    viz += cuttingSurf.second->print_visualization(usedGeomIds, interfaceDesc);
  }
  return viz;
}


void Element_Cutter::fill_triangle_interior_intersection_points(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const
{
  std::vector<InterfaceID> interfacesWithCuttingSurface;
  fill_interfaces_with_cutting_surface(interfacesWithCuttingSurface);

  std::vector<int> sortedDomains;
  std::set<std::vector<int>> alreadyFound;
  stk::math::Vector3d intersectionPoint;

  intersections.clear();
  for (size_t index0=0; index0<interfacesWithCuttingSurface.size(); ++index0)
  {
    for (size_t index1=index0+1; index1<interfacesWithCuttingSurface.size(); ++index1)
    {
      fill_sorted_domains(is_one_ls_per_phase(), interfacesWithCuttingSurface[index0], interfacesWithCuttingSurface[index1], sortedDomains);
      if (alreadyFound.find(sortedDomains) == alreadyFound.end() &&
          intersectionPointFilter(sortedDomains) &&
          sorted_domains_form_triple_point(is_one_ls_per_phase(), sortedDomains))
      {
        alreadyFound.insert(sortedDomains);
        const stk::math::Plane3d & plane0 = cutting_surfaces.at(interfacesWithCuttingSurface[index0])->get_plane();
        const stk::math::Plane3d & plane1 = cutting_surfaces.at(interfacesWithCuttingSurface[index1])->get_plane();
        if (find_intersection_of_two_2D_planes_within_tri(plane0, plane1, intersectionPoint) &&
            intersection_point_is_real(intersectionPoint, sortedDomains))
        {
          intersections.emplace_back(intersectionPoint, sortedDomains);
        }
      }
    }
  }
}

void Element_Cutter::append_tetrahedron_face_interior_intersections(const std::array<stk::math::Vector3d,3> & faceNodes,
    const InterfaceID & interface1,
    const InterfaceID & interface2,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const
{
  std::vector<int> sortedDomains;
  fill_sorted_domains(is_one_ls_per_phase(), interface1, interface2, sortedDomains);
  if (sorted_domains_form_triple_point(is_one_ls_per_phase(), sortedDomains) &&
      intersectionPointFilter(sortedDomains))
  {
    const stk::math::Plane3d facePlane{faceNodes[0], faceNodes[1], faceNodes[2]};
    stk::math::Vector3d intersectionPoint;

    const stk::math::Plane3d & plane0 = cutting_surfaces.at(interface1)->get_plane();
    const stk::math::Plane3d & plane1 = cutting_surfaces.at(interface2)->get_plane();
    if (find_intersection_of_three_planes(plane0, plane1, facePlane, intersectionPoint) &&
        intersection_point_is_real(intersectionPoint, sortedDomains))
    {
      const stk::math::Vector3d faceCoords = triangle_parametric_coordinates_of_projected_point(faceNodes, intersectionPoint);
      if (within_tri_bounds(faceCoords))
        intersections.emplace_back(faceCoords, sortedDomains);
    }
  }
}

void Element_Cutter::fill_tetrahedron_interior_intersection_points(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const
{
  std::vector<InterfaceID> interfacesWithCuttingSurface;
  fill_interfaces_with_cutting_surface(interfacesWithCuttingSurface);

  std::vector<int> sortedDomains;
  std::set<std::vector<int>> alreadyFound;
  stk::math::Vector3d intersectionPoint;

  intersections.clear();
  for (size_t index0=0; index0<interfacesWithCuttingSurface.size(); ++index0)
  {
    for (size_t index1=index0+1; index1<interfacesWithCuttingSurface.size(); ++index1)
    {
      for (size_t index2=index1+1; index2<interfacesWithCuttingSurface.size(); ++index2)
      {
        fill_sorted_domains(is_one_ls_per_phase(), interfacesWithCuttingSurface[index0], interfacesWithCuttingSurface[index1], interfacesWithCuttingSurface[index2], sortedDomains);
        if (alreadyFound.find(sortedDomains) == alreadyFound.end() &&
            intersectionPointFilter(sortedDomains) &&
            sorted_domains_form_quadruple_point(is_one_ls_per_phase(), sortedDomains))
        {
          alreadyFound.insert(sortedDomains);
          const stk::math::Plane3d & plane0 = cutting_surfaces.at(interfacesWithCuttingSurface[index0])->get_plane();
          const stk::math::Plane3d & plane1 = cutting_surfaces.at(interfacesWithCuttingSurface[index1])->get_plane();
          const stk::math::Plane3d & plane2 = cutting_surfaces.at(interfacesWithCuttingSurface[index2])->get_plane();
          if (find_intersection_of_three_planes_within_tet(plane0, plane1, plane2, intersectionPoint) &&
              intersection_point_is_real(intersectionPoint, sortedDomains))
          {
            intersections.emplace_back(intersectionPoint, sortedDomains);
          }
        }
      }
    }
  }
}

static ElementIntersectionPointFilter
keep_all_intersecion_points_filter()
{
  auto filter =
  [](const std::vector<int> & /*intersectionPointSortedDomains*/)
  {
    return true;
  };
  return filter;
}

void Element_Cutter::fill_interior_intersections(std::vector<ElementIntersection> & intersections) const
{
  fill_interior_intersections(keep_all_intersecion_points_filter(), intersections);
}

void Element_Cutter::fill_interior_intersections(const ElementIntersectionPointFilter & intersectionPointFilter, std::vector<ElementIntersection> & intersections) const
{
  if (myTopology.base() == stk::topology::TETRAHEDRON_4)
    fill_tetrahedron_interior_intersection_points(intersectionPointFilter, intersections);
  else if (myTopology.base() == stk::topology::TRIANGLE_3_2D || myTopology.base() == stk::topology::TRIANGLE_3)
    fill_triangle_interior_intersection_points(intersectionPointFilter, intersections);
  else
    STK_ThrowErrorMsg("Unexepected topology " << myTopology.name());
}

std::unique_ptr<Element_Cutter> create_element_cutter(const bool oneLSPerPhase,
    const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
{
  std::unique_ptr<Element_Cutter> cutter;
  if(oneLSPerPhase)
    cutter.reset(new One_LS_Per_Phase_Cutter(masterElem, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges, intersectingPlanesDiagonalPicker));
  else
    cutter.reset(new LS_Per_Interface_Cutter(masterElem, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges, intersectingPlanesDiagonalPicker));

  return cutter;
}


int
Element_Cutter::sign_at_position(const InterfaceID interface, const stk::math::Vector3d & p_coords) const
{
  const auto iter = cutting_surfaces.find(interface);
  if(iter == cutting_surfaces.end())
  {
    // Right now, One_LS_Per_Phase_Cutter does not build up the cutting surfaces for the single_phase case.  So
    // we need a "default" answer.
    return +1;
  }
  const Cutting_Surface & interface_surface = *(iter->second);

  return interface_surface.sign_at_position(p_coords);
}

double
Element_Cutter::interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  STK_ThrowErrorMsgIf(cutting_surfaces.find(interface) == cutting_surfaces.end(),
      "No cutting surface found for interface " << interface);
  const Cutting_Surface & interface_surface = *(cutting_surfaces.find(interface)->second);
  return interface_surface.interface_crossing_position(edgeNodeCoords);
}

void
Element_Cutter::fill_interfaces_with_cutting_surface(std::vector<InterfaceID> & interfaces) const
{
  interfaces.clear();
  for (auto && entry : cutting_surfaces)
    interfaces.push_back(entry.first);
}

int
Element_Cutter::get_starting_phase_for_cutting_surfaces() const
{
  if (cutting_surfaces.empty())
    return -1;
  return cutting_surfaces.begin()->first.first_ls();
}

bool
Element_Cutter::have_cutting_surface(const InterfaceID interface) const
{
  return cutting_surfaces.find(interface) != cutting_surfaces.end();
}

Cutting_Surface *
Element_Cutter::get_cutting_surface(const InterfaceID interface) const
{
  const auto iter = cutting_surfaces.find(interface);
  if (iter != cutting_surfaces.end())
    return iter->second.get();
  return nullptr;
}

static std::vector< std::array<unsigned,2> > get_element_edge_node_ordinals(const stk::topology baseTopology)
{
  const unsigned numEdges = baseTopology.num_edges();
  std::vector< std::array<unsigned,2> > elemEdgeNodeOrdinals(numEdges);
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * node_ordinals = get_edge_node_ordinals(baseTopology, i);
    elemEdgeNodeOrdinals[i][0] = node_ordinals[0];
    elemEdgeNodeOrdinals[i][1] = node_ordinals[1];
  }
  return elemEdgeNodeOrdinals;
}

static std::vector<stk::math::Vector3d>
get_parametric_node_coords_on_element_nodes_and_cut_edges(const MasterElement & masterElem,
    const std::vector<Element_Cutter::Edge_Crossing> & cutEdges,
    const std::vector< std::array<unsigned,2> > elemEdgeNodeOrdinals)
{
  const stk::topology baseTopology = masterElem.get_topology().base();
  const double * elemNodeParamCoords = masterElem.nodal_parametric_coordinates();

  const unsigned numNodes = baseTopology.num_nodes();
  const unsigned numEdges = baseTopology.num_edges();
  const int dim = baseTopology.dimension();

  std::vector<stk::math::Vector3d> nodeParamCoords(numNodes+numEdges);

  for (unsigned n=0; n<numNodes; ++n)
    nodeParamCoords[n] = stk::math::Vector3d(&elemNodeParamCoords[n*dim],dim);

  for (auto && cut_edge : cutEdges)
  {
    const unsigned edge = cut_edge.edge;
    const double pos = cut_edge.pos;

    nodeParamCoords[numNodes+edge] = (1.-pos)*nodeParamCoords[elemEdgeNodeOrdinals[edge][0]] + pos*nodeParamCoords[elemEdgeNodeOrdinals[edge][1]];
  }

  return nodeParamCoords;
}

static std::vector<int>
get_node_signs_on_element_nodes(const stk::topology baseTopology,
    const std::vector<Element_Cutter::Edge_Crossing> & cutEdges,
    const std::vector< std::array<unsigned,2> > elemEdgeNodeOrdinals)
{
  const unsigned numNodes = baseTopology.num_nodes();
  std::vector<int> nodeSigns(numNodes,-2);

  unsigned edge_error = 0;
  for (auto && cut_edge : cutEdges)
  {
    const unsigned edge = cut_edge.edge;
    const unsigned node0_ord = elemEdgeNodeOrdinals[edge][0];
    const unsigned node1_ord = elemEdgeNodeOrdinals[edge][1];

    int node0_sign = -cut_edge.sign;
    int node1_sign = cut_edge.sign;

    if (cut_edge.pos < Element_Cutter::theSnapToNodeTol)
      node0_sign = 0;
    else if (cut_edge.pos > 1.-Element_Cutter::theSnapToNodeTol)
      node1_sign = 0;

    if ((nodeSigns[node0_ord]*node0_sign == -1) ||
        (nodeSigns[node1_ord]*node1_sign == -1))
      edge_error += 1<<edge;

    if (nodeSigns[node0_ord] != 0)
      nodeSigns[node0_ord] = node0_sign;
    if (nodeSigns[node1_ord] != 0)
      nodeSigns[node1_ord] = node1_sign;
  }

  if (edge_error > 0)
  {
    std::ostringstream err_msg;
    err_msg << "Edge crossings are not consistent: " << std::endl;
    for(auto && cut_edge : cutEdges)
    {
      const unsigned edge = cut_edge.edge;
      const double pos = cut_edge.pos;
      const int sign = cut_edge.sign;
      const unsigned node0_ord = elemEdgeNodeOrdinals[edge][0];
      const unsigned node1_ord = elemEdgeNodeOrdinals[edge][1];
      err_msg << "  Edge: " << edge << " error = " << ((edge_error & (1<<edge)) > 0) << " pos = " << pos << " sign = " << sign << " node0_ord = " << node0_ord << " node1_ord = " << node1_ord << std::endl;
    }
    krinolog << err_msg.str();
    throw std::runtime_error(err_msg.str());
  }

  return nodeSigns;
}

static std::shared_ptr<Cutting_Surface>
build_triangle_cutting_surface(const std::vector<int> & nodeSigns, const std::vector<stk::math::Vector3d> & nodeParamCoords)
{
  const int case_id = (nodeSigns[0]+1) +
                      (nodeSigns[1]+1)*3 +
                      (nodeSigns[2]+1)*9;

  static const unsigned case_permutations[] =
  { 0, 0, 0, 2, 0, 0, 2, 1, 1, 1, //  0-9
    1, 2, 2, 0, 2, 2, 1, 1, 1, 1, //  10-19
    2, 0, 0, 2, 0, 0, 0 };        //  20-26

  stk::topology topo = stk::topology::TRIANGLE_6_2D;
  std::vector<unsigned> permute(6);
  topo.permutation_node_ordinals(case_permutations[case_id], permute.begin());

  const unsigned i0 = permute[0];
  const unsigned i1 = permute[1];
  const unsigned i2 = permute[2];
  const unsigned i3 = permute[3];
  const unsigned i5 = permute[5];

  const int permute_case_id =
      (nodeSigns[i0]+1) +
      (nodeSigns[i1]+1)*3 +
      (nodeSigns[i2]+1)*9;

  switch (permute_case_id)
  {
    case 1:  // ls[0]==0 && ls[1]<0 && ls[2]<0
    case 25: // ls[0]==0 && ls[1]>0 && ls[2]>0
    {
      STK_ThrowAssert((nodeSigns[i0] == 0 && nodeSigns[i1] < 0 && nodeSigns[i2] < 0) || (nodeSigns[i0] == 0 && nodeSigns[i1] > 0 && nodeSigns[i2] > 0));
      return (permute_case_id==1) ?
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i0], (nodeParamCoords[i0] + nodeParamCoords[i1] - nodeParamCoords[i2]))) :
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i0], (nodeParamCoords[i0] + nodeParamCoords[i2] - nodeParamCoords[i1])));
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0
    case 24: // ls[0]<0 && ls[1]>0 && ls[2]>0)
    {
      STK_ThrowAssert((nodeSigns[i0] > 0 && nodeSigns[i1] < 0 && nodeSigns[i2] < 0) || (nodeSigns[i0] < 0 && nodeSigns[i1] > 0 && nodeSigns[i2] > 0));
      return (permute_case_id==2) ?
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i5], nodeParamCoords[i3])) :
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i3], nodeParamCoords[i5]));
    }
    break;

    case 5:  // ls[0]>0 && ls[1]==0 && ls[2]<0
    case 21: // ls[0]<0 && ls[1]==0 && ls[2]>0
    {
      STK_ThrowAssert((nodeSigns[i0] > 0 && nodeSigns[i1] == 0 && nodeSigns[i2] < 0) || (nodeSigns[i0] < 0 && nodeSigns[i1] == 0 && nodeSigns[i2] > 0));
      return (permute_case_id==5) ?
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i5], nodeParamCoords[i1])) :
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i1], nodeParamCoords[i5]));
    }
    break;

    case 4:  // ls[0]==0 && ls[1]==0 && ls[2]<0
    case 22: // ls[0]==0 && ls[1]==0 && ls[2]>0
    {
      STK_ThrowAssert((nodeSigns[i0] == 0 && nodeSigns[i1] == 0 && nodeSigns[i2] < 0) || (nodeSigns[i0] == 0 && nodeSigns[i1] == 0 && nodeSigns[i2] > 0));
      return (permute_case_id==4) ?
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i0], nodeParamCoords[i1])) :
          (std::make_shared<Plane_Cutting_Surface_2D>(nodeParamCoords[i1], nodeParamCoords[i0]));
    }
    break;

    default:
    {
      krinolog << "Case id " << case_id << " " << permute_case_id << stk::diag::dendl;
      ThrowRuntimeError("Subelement decomposition error.");
    }
  }

  ThrowRuntimeError("Subelement decomposition error.");
}

static std::shared_ptr<Cutting_Surface>
build_tetrahedral_cutting_surface(const std::vector<int> & nodeSigns,
    const std::vector<stk::math::Vector3d> & nodeParamCoords,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
{
  const int case_id = (nodeSigns[0]+1) +
                      (nodeSigns[1]+1)*3 +
                      (nodeSigns[2]+1)*9 +
                      (nodeSigns[3]+1)*27;

  static const unsigned case_permutations[] =
    { 0, 0, 0, 1, 0, 0, 1, 5, 0, 2, //  0-9
      2, 6, 1, 0, 0, 1, 1,10, 2, 2, //  10-19
      2,11, 2, 4, 1, 8, 4, 4, 3, 3, //  20-29
      4, 3, 3, 9, 5, 7, 7, 6, 6, 9, //  30-39
      0, 9, 9, 6, 7, 7, 7, 9,11, 3, //  40-49
      4, 3, 3, 4, 4, 8, 3, 4, 4,11, //  50-59
      4, 2, 2,10, 8, 1,10, 0, 1, 6, //  60-69
      2, 2, 7, 5, 1, 0, 0, 1, 0, 0, //  70-79
      0 };

  stk::topology full_topology = stk::topology::TETRAHEDRON_10;
  std::vector<unsigned> permute_node_ordinals(10);
  full_topology.permutation_node_ordinals(case_permutations[case_id], permute_node_ordinals.begin());

  const unsigned i0 = permute_node_ordinals[0];
  const unsigned i1 = permute_node_ordinals[1];
  const unsigned i2 = permute_node_ordinals[2];
  const unsigned i3 = permute_node_ordinals[3];
  const unsigned i4 = permute_node_ordinals[4];
  const unsigned i5 = permute_node_ordinals[5];
  const unsigned i6 = permute_node_ordinals[6];
  const unsigned i7 = permute_node_ordinals[7];
  const unsigned i8 = permute_node_ordinals[8];

  const int permute_case_id = (nodeSigns[i0]+1) +
                              (nodeSigns[i1]+1)*3 +
                              (nodeSigns[i2]+1)*9 +
                              (nodeSigns[i3]+1)*27;

  switch (permute_case_id)
  {
    case 1:  // ls[0]=0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 79: // ls[0]=0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x0p13 = nodeParamCoords[i0] + nodeParamCoords[i3]-nodeParamCoords[i1];
      const stk::math::Vector3d x0p12 = nodeParamCoords[i0] + nodeParamCoords[i2]-nodeParamCoords[i1];
      return (permute_case_id==1) ?
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i0],x0p13,x0p12)) :
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i0],x0p12,x0p13));
    }
    break;

    case 4:  // ls[0]=0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 76: // ls[0]=0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    {
      const stk::math::Vector3d x0p23 = nodeParamCoords[i0] + nodeParamCoords[i3]-nodeParamCoords[i2];
      return (permute_case_id==4) ?
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i0],nodeParamCoords[i1],x0p23)) :
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i1],nodeParamCoords[i0],x0p23));
    }
    break;

    case 2:  // ls[0]>0 && ls[1]<0 && ls[2]<0 && ls[3]<0
    case 78: // ls[0]<0 && ls[1]>0 && ls[2]>0 && ls[3]>0
    {
      STK_ThrowAssert((nodeSigns[i0] > 0 && nodeSigns[i1] < 0 && nodeSigns[i2] < 0 && nodeSigns[i3] < 0) || (nodeSigns[i0] < 0 && nodeSigns[i1] > 0 && nodeSigns[i2] > 0 && nodeSigns[i3] > 0));

      return (permute_case_id==2) ?
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i4],nodeParamCoords[i7],nodeParamCoords[i6])) :
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i4],nodeParamCoords[i6],nodeParamCoords[i7]));
    }
    break;

    case 5:  // ls[0]>0 && ls[1]=0 && ls[2]<0 && ls[3]<0
    case 75: // ls[0]<0 && ls[1]=0 && ls[2]>0 && ls[3]>0
    {
      return (permute_case_id==5) ?
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i1],nodeParamCoords[i7],nodeParamCoords[i6])) :
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i1],nodeParamCoords[i6],nodeParamCoords[i7]));
    }
    break;

    case 8:  // ls[0]>0 && ls[1]>0 && ls[2]<0 && ls[3]<0
    {
      const bool face_diag = intersectingPlanesDiagonalPicker({{i0,i1,i2,i3}});
      return (face_diag) ?
          (std::make_shared<Intersecting_Planes_Cutting_Surface>(nodeParamCoords[i8],nodeParamCoords[i7],nodeParamCoords[i6],nodeParamCoords[i5])) :
          (std::make_shared<Intersecting_Planes_Cutting_Surface>(nodeParamCoords[i7],nodeParamCoords[i6],nodeParamCoords[i5],nodeParamCoords[i8]));
    }
    break;

    case 13: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    case 67: // ls[0]=0 && ls[1]=0 && ls[2]=0 && ls[3]>0
    {
      return (permute_case_id==13) ?
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i0],nodeParamCoords[i2],nodeParamCoords[i1])) :
          (std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i0],nodeParamCoords[i1],nodeParamCoords[i2]));
    }
    break;

    case 14: // ls[0]>0 && ls[1]=0 && ls[2]=0 && ls[3]<0
    {
      STK_ThrowAssert(nodeSigns[i0] > 0 && nodeSigns[i1] == 0 && nodeSigns[i2] == 0 && nodeSigns[i3] < 0);
      return std::make_shared<Plane_Cutting_Surface>(nodeParamCoords[i2],nodeParamCoords[i1],nodeParamCoords[i7]);
    }
    break;

    default:
    {
      krinolog << "Case id " << case_id << " " << permute_case_id << stk::diag::dendl;
      ThrowRuntimeError("Subelement decomposition error.");
    }
  }

  ThrowRuntimeError("Subelement decomposition error.");
}

static std::shared_ptr<Cutting_Surface>
build_cutting_surface(const MasterElement & masterElem,
    const std::vector<Element_Cutter::Edge_Crossing> & cutEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
{
  const stk::topology baseTopology = masterElem.get_topology().base();
  const std::vector< std::array<unsigned,2> > elemEdgeNodeOrdinals = get_element_edge_node_ordinals(baseTopology);

  const auto nodeParamCoords = get_parametric_node_coords_on_element_nodes_and_cut_edges(masterElem, cutEdges, elemEdgeNodeOrdinals);
  const auto nodeSigns = get_node_signs_on_element_nodes(baseTopology, cutEdges, elemEdgeNodeOrdinals);

  if (baseTopology == stk::topology::TRIANGLE_3_2D || baseTopology == stk::topology::TRIANGLE_3)
  {
    STK_ThrowRequire(cutEdges.size() == 2);
    return build_triangle_cutting_surface(nodeSigns, nodeParamCoords);
  }

  STK_ThrowRequireMsg(baseTopology == stk::topology::TETRAHEDRON_4, "Unsupported base topology: " << baseTopology.name());
  STK_ThrowRequire(cutEdges.size() == 3 || cutEdges.size() == 4);
  return build_tetrahedral_cutting_surface(nodeSigns, nodeParamCoords, intersectingPlanesDiagonalPicker);
}

bool node_captures_intersection_point_domains(const std::vector<int> & nodeDomains, const std::vector<int> & intersectionPointDomains)
{
  return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(nodeDomains, intersectionPointDomains);
}

template <class NODEDOMAINS>
bool any_node_captures_intersection_point_domains(const NODEDOMAINS & nodesSnappedDomains, const std::vector<int> & intersectionPointDomains)
{
  if (nodesSnappedDomains.empty())
    return false;
  for (auto && nodeSnappedDomains : nodesSnappedDomains)
    if (node_captures_intersection_point_domains(*nodeSnappedDomains, intersectionPointDomains))
      return true;
  return false;
}

bool intersection_is_already_captured(const std::vector<const std::vector<int>*> & elementNodesSnappedDomains, const std::vector<unsigned> & nodes, const std::vector<int> & intersectionPointDomains)
{
  if (elementNodesSnappedDomains.empty())
    return false;
  for (auto && node : nodes)
    if (node_captures_intersection_point_domains(*elementNodesSnappedDomains[node], intersectionPointDomains))
      return true;
  return false;
}

LS_Per_Interface_Cutter::LS_Per_Interface_Cutter(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
  : Element_Cutter(masterElem)
{
  for(auto && interface : get_interfaces_present(parentEdges))
  {
    const std::vector<Element_Cutter::Edge_Crossing> cutEdges = build_cut_edges(interface, parentEdges, areParentEdgesOrientedSameAsElementEdges);
    STK_ThrowRequire(!cutEdges.empty());

    cutting_surfaces[interface] = build_cutting_surface(masterElem, cutEdges, intersectingPlanesDiagonalPicker);
  }
}

bool
LS_Per_Interface_Cutter::have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const auto iter = cutting_surfaces.find(interface);
  STK_ThrowRequire(iter != cutting_surfaces.end());
  const Cutting_Surface & interface_surface = *(iter->second);
  return interface_surface.sign_at_position(edgeNodeCoords[0]) == -interface_surface.sign_at_position(edgeNodeCoords[1]);
}

int LS_Per_Interface_Cutter::get_ls_per_interface_phase_at_location(const stk::math::Vector3d & /*pCoords*/) const
{
  STK_ThrowRequireMsg(false, "Improper usage of LS_Per_Interface_Cutter.");
  return -1;
}

void LS_Per_Interface_Cutter::add_interfaces_with_uncaptured_intersection_within_element(const std::vector<stk::math::Vector3d> & /*elemNodesCoords*/,
    const std::vector<const std::vector<int> *> & /*elemNodesSnappedDomains*/,
    std::set<InterfaceID> & /*interfacesWithUncapturedCrossings*/) const
{
  // For LS_Per_Interface_Cutter, internal intersections only happen if there are edge intersections so nothing to do here.
}

std::vector<Element_Cutter::Edge_Crossing>
LS_Per_Interface_Cutter::build_cut_edges(const InterfaceID & interface,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & parentEdgeIsOrientedSameAsElementEdge)
{
  std::vector<Edge_Crossing> cutEdges;

  for(unsigned i=0; i < parentEdges.size(); ++i)
  {
    const CDFEM_Parent_Edge * parent_edge = parentEdges[i];
    if(parent_edge)
    {
      const int crossing_sign = parent_edge->get_crossing_sign(interface);
      const double crossing_pos = parent_edge->get_crossing_position(interface);

      if( crossing_pos >= 0. )
      {
        const double pos = parentEdgeIsOrientedSameAsElementEdge[i] ? crossing_pos : (1 -crossing_pos);
        const int sign = parentEdgeIsOrientedSameAsElementEdge[i] ? crossing_sign : (-crossing_sign);
        cutEdges.push_back(Edge_Crossing(i,pos,sign));
        if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "Crossing of edge " << i << " at position " << pos << "\n";
      }
    }
  }
  return cutEdges;
}

std::vector<Element_Cutter::Edge_Crossing>
One_LS_Per_Phase_Cutter::build_cut_edges(const InterfaceID & interface,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & parentEdgeIsOrientedSameAsElementEdge)
{
  std::vector<Edge_Crossing> cutEdges;
  const bool isOneLSPerPhase = true;

  for(unsigned i=0; i < parentEdges.size(); ++i)
  {
    const CDFEM_Parent_Edge * parent_edge = parentEdges[i];
    if(parent_edge)
    {
      double crossing_pos;
      int crossing_sign;
      bool is_fake;
      std::tie(crossing_pos, crossing_sign, is_fake) = parent_edge->get_crossing_position_and_sign(isOneLSPerPhase, interface);

      if( crossing_pos >= 0. )
      {
        const double pos = parentEdgeIsOrientedSameAsElementEdge[i] ? crossing_pos : (1 -crossing_pos);
        const int sign = parentEdgeIsOrientedSameAsElementEdge[i] ? crossing_sign : (-crossing_sign);
        cutEdges.push_back(Edge_Crossing(i,pos,sign));
        if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "Crossing of edge " << i << " for interface " << interface << " at position " << pos << " fake = " << std::boolalpha << is_fake << "\n";
      }
    }
  }
  return cutEdges;
}

std::shared_ptr<Cutting_Surface> One_LS_Per_Phase_Cutter::attempt_to_build_cutting_surface(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker,
    const InterfaceID& interface)
{
  std::shared_ptr<Cutting_Surface> cuttingSurface;
  try
  {
    const std::vector<Edge_Crossing> cutEdges = build_cut_edges(interface, parentEdges, areParentEdgesOrientedSameAsElementEdges);
    STK_ThrowRequire(!cutEdges.empty());
    cuttingSurface = build_cutting_surface(masterElem, cutEdges, intersectingPlanesDiagonalPicker);
  }
  catch (const std::exception & err)
  {
    std::stringstream err_msg;
    err_msg << "Error while cutting element with interface " << interface << std::endl;
    err_msg << err.what() << std::endl;
    throw std::runtime_error(err_msg.str());
  }
  return cuttingSurface;
}

static std::set<int> determine_optimal_phases_at_location(const stk::math::Vector3d & location,
    const InterfaceToSurface & allSurfaces)
{
  // optimal phases are phases those that are on the "right" side of at least one surface,
  // and never on the "wrong" side of an interface
  std::set<int> toPhases;
  std::set<int> fromPhases;

  for (auto && surface : allSurfaces)
  {
    const int sign = surface.second->sign_at_position(location);
    const int toPhase = (sign < 0) ? surface.first.first_ls() : surface.first.second_ls();
    const int fromPhase = (sign < 0) ? surface.first.second_ls() : surface.first.first_ls();
    toPhases.insert(toPhase);
    fromPhases.insert(fromPhase);
  }

  for (int phase : fromPhases)
    toPhases.erase(phase);

  return toPhases;
}

bool intersection_point_is_real(const stk::math::Vector3d & intersectionPoint,
    const InterfaceToSurface & allSurfaces,
    const std::vector<int> & sortedDomains)
{
  for (auto && surface : allSurfaces)
  {
    if (!surface.second->on_surface(intersectionPoint, Element_Cutter::theSnapToNodeTol))
    {
      const int sign = surface.second->sign_at_position(intersectionPoint);
      const int toPhase = (sign < 0) ? surface.first.first_ls() : surface.first.second_ls();
      const int fromPhase = (sign < 0) ? surface.first.second_ls() : surface.first.first_ls();
      const bool toPhaseIsInDomains = std::binary_search(sortedDomains.begin(), sortedDomains.end(), toPhase);
      const bool fromPhaseIsInDomains = std::binary_search(sortedDomains.begin(), sortedDomains.end(), fromPhase);
      if (!toPhaseIsInDomains && fromPhaseIsInDomains)
        return false;
    }
  }
  return true;
}

static void add_triple_point_interfaces(const bool isOneLSPerPhase, const std::vector<int> & triplePointDomains, std::set<InterfaceID> & interfaces)
{
  STK_ThrowRequire(triplePointDomains.size() == 3);
  if (isOneLSPerPhase)
  {
    interfaces.insert(InterfaceID(triplePointDomains[0],triplePointDomains[1]));
    interfaces.insert(InterfaceID(triplePointDomains[1],triplePointDomains[2]));
    interfaces.insert(InterfaceID(triplePointDomains[0],triplePointDomains[2]));
  }
  else
  {
    interfaces.insert(InterfaceID(triplePointDomains[0],triplePointDomains[0]));
    interfaces.insert(InterfaceID(triplePointDomains[1],triplePointDomains[1]));
    interfaces.insert(InterfaceID(triplePointDomains[2],triplePointDomains[2]));
  }
}

static void add_interfaces_that_have_uncaptured_intersection_within_sides_of_tet(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface & allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  stk::math::Vector3d intersectionPoint;
  std::vector<int> sortedDomains;
  constexpr bool oneLSPerPhase(true);

  stk::topology topology = stk::topology::TETRAHEDRON_4;
  std::vector<unsigned> faceNodes(3);

  for (auto && surface1 : allSurfaces)
  {
    for (auto && surface2 : allSurfaces)
    {
      if (surface1.first < surface2.first &&
          (interfacesWithUncapturedCrossings.count(surface1.first) == 0 || interfacesWithUncapturedCrossings.count(surface2.first)))
      {
        fill_sorted_domains(oneLSPerPhase, surface1.first, surface2.first, sortedDomains);

        if (sorted_domains_form_triple_point(oneLSPerPhase, sortedDomains))
        {
          for (int iFace=0; iFace<4; ++iFace)
          {
            topology.face_node_ordinals(iFace, faceNodes);
            const std::array<stk::math::Vector3d,3> faceNodeCoords{elemNodesCoords[faceNodes[0]], elemNodesCoords[faceNodes[1]], elemNodesCoords[faceNodes[2]]};
            const stk::math::Plane3d facePlane{elemNodesCoords[faceNodes[0]], elemNodesCoords[faceNodes[1]], elemNodesCoords[faceNodes[2]]};
            if (!intersection_is_already_captured(nodesSnappedDomains, faceNodes, sortedDomains) &&
                find_intersection_of_three_planes(facePlane, surface1.second->get_plane(), surface2.second->get_plane(), intersectionPoint) &&
                intersection_point_is_real(intersectionPoint, allSurfaces, sortedDomains))
            {
              const stk::math::Vector3d faceCoords = triangle_parametric_coordinates_of_projected_point(faceNodeCoords, intersectionPoint);
              if (within_tri_bounds(faceCoords))
                add_triple_point_interfaces(oneLSPerPhase, sortedDomains, interfacesWithUncapturedCrossings);
            }
          }
        }
      }
    }
  }
}

static void add_interfaces_that_have_uncaptured_intersection_within_sides_of_tet(const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface & allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  stk::math::Vector3d intersectionPoint;
  std::vector<int> sortedDomains;
  constexpr bool oneLSPerPhase(true);

  stk::topology topology = stk::topology::TETRAHEDRON_4;
  std::vector<unsigned> faceNodes(3);

  for (auto && surface1 : allSurfaces)
  {
    for (auto && surface2 : allSurfaces)
    {
      if (surface1.first < surface2.first &&
          (interfacesWithUncapturedCrossings.count(surface1.first) == 0 || interfacesWithUncapturedCrossings.count(surface2.first)))
      {
        fill_sorted_domains(oneLSPerPhase, surface1.first, surface2.first, sortedDomains);
        if (sorted_domains_form_triple_point(oneLSPerPhase, sortedDomains))
        {
          for (int iFace=0; iFace<4; ++iFace)
          {
            topology.face_node_ordinals(iFace, faceNodes);
            if (!intersection_is_already_captured(nodesSnappedDomains, faceNodes, sortedDomains) &&
                find_intersection_of_two_planes_and_side_of_tet(iFace, surface1.second->get_plane(), surface2.second->get_plane(), intersectionPoint) &&
                intersection_point_is_real(intersectionPoint, allSurfaces, sortedDomains))
              add_triple_point_interfaces(oneLSPerPhase, sortedDomains, interfacesWithUncapturedCrossings);
          }
        }
      }
    }
  }
}

static void add_interfaces_that_have_uncaptured_intersection_within_tri(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface & allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  stk::math::Vector3d intersectionPoint;
  std::vector<int> sortedDomains;
  constexpr bool oneLSPerPhase(true);
  const std::vector<unsigned> triNodes = {0,1,2};

  for (auto && surface1 : allSurfaces)
  {
    for (auto && surface2 : allSurfaces)
    {
      if (surface1.first < surface2.first &&
          (interfacesWithUncapturedCrossings.count(surface1.first) == 0 || interfacesWithUncapturedCrossings.count(surface2.first)))
      {
        fill_sorted_domains(oneLSPerPhase, surface1.first, surface2.first, sortedDomains);
        if (sorted_domains_form_triple_point(oneLSPerPhase, sortedDomains) &&
            !intersection_is_already_captured(nodesSnappedDomains, triNodes, sortedDomains) &&
            find_intersection_of_two_2D_planes_within_tri(surface1.second->get_plane(), surface2.second->get_plane(), intersectionPoint) && // ASSUME THAT PROVIDED TRIANGLE CONTAINED WITHIN PARENT TRIANGLE
            intersection_point_is_real(intersectionPoint, allSurfaces, sortedDomains))
        {
          const std::array<stk::math::Vector3d,3> triNodeCoords{elemNodesCoords[0], elemNodesCoords[1], elemNodesCoords[2]};
          const stk::math::Vector3d triCoords = triangle_parametric_coordinates_of_projected_point(triNodeCoords, intersectionPoint);
          if (within_tri_bounds(triCoords))
            add_triple_point_interfaces(oneLSPerPhase, sortedDomains, interfacesWithUncapturedCrossings);
        }
      }
    }
  }
}

static void add_interfaces_that_have_uncaptured_intersection_within_tri(const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface& allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  stk::math::Vector3d intersectionPoint;
  std::vector<int> sortedDomains;
  constexpr bool oneLSPerPhase(true);
  const std::vector<unsigned> triNodes = {0,1,2};

  for (auto && surface1 : allSurfaces)
  {
    for (auto && surface2 : allSurfaces)
    {
      if (surface1.first < surface2.first &&
          (interfacesWithUncapturedCrossings.count(surface1.first) == 0 || interfacesWithUncapturedCrossings.count(surface2.first)))
      {
        fill_sorted_domains(oneLSPerPhase, surface1.first, surface2.first, sortedDomains);
        if (sorted_domains_form_triple_point(oneLSPerPhase, sortedDomains) &&
            !intersection_is_already_captured(nodesSnappedDomains, triNodes, sortedDomains) &&
            find_intersection_of_two_2D_planes_within_tri(surface1.second->get_plane(), surface2.second->get_plane(), intersectionPoint) &&
            intersection_point_is_real(intersectionPoint, allSurfaces, sortedDomains))
          add_triple_point_interfaces(oneLSPerPhase, sortedDomains, interfacesWithUncapturedCrossings);
      }
    }
  }
}

static void add_interfaces_that_have_uncaptured_intersection_within_element(const stk::topology & baseTopology,
    const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface & allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  if (baseTopology == stk::topology::TETRAHEDRON_4)
    add_interfaces_that_have_uncaptured_intersection_within_sides_of_tet(nodesSnappedDomains, allSurfaces, interfacesWithUncapturedCrossings);
  else if (baseTopology == stk::topology::TRIANGLE_3 || baseTopology == stk::topology::TRIANGLE_3_2D)
    add_interfaces_that_have_uncaptured_intersection_within_tri(nodesSnappedDomains, allSurfaces, interfacesWithUncapturedCrossings);
  else
    STK_ThrowRequireMsg(false, "Unsupported topology " << baseTopology);
}

static void add_interfaces_that_have_uncaptured_intersection_within_element(const stk::topology & baseTopology,
    const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int>*> & nodesSnappedDomains,
    const InterfaceToSurface & allSurfaces,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings)
{
  if (baseTopology == stk::topology::TETRAHEDRON_4)
    add_interfaces_that_have_uncaptured_intersection_within_sides_of_tet(elemNodesCoords, nodesSnappedDomains, allSurfaces, interfacesWithUncapturedCrossings);
  else if (baseTopology == stk::topology::TRIANGLE_3 || baseTopology == stk::topology::TRIANGLE_3_2D)
    add_interfaces_that_have_uncaptured_intersection_within_tri(elemNodesCoords, nodesSnappedDomains, allSurfaces, interfacesWithUncapturedCrossings);
  else
    STK_ThrowRequireMsg(false, "Unsupported topology " << baseTopology);
}

void One_LS_Per_Phase_Cutter::add_interfaces_with_uncaptured_intersection_within_element(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    std::set<InterfaceID> & interfacesWithUncapturedCrossings) const
{
  // Does not consider uncaptured edge crossings, these should be considered already
  add_interfaces_that_have_uncaptured_intersection_within_element(myTopology.base(), elemNodesCoords, elemNodesSnappedDomains, all_cutting_surfaces, interfacesWithUncapturedCrossings);
}

bool
One_LS_Per_Phase_Cutter::have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const auto iter = cutting_surfaces.find(interface);
  STK_ThrowRequire(iter != cutting_surfaces.end());
  const Cutting_Surface & interface_surface = *(iter->second);
  if (interface_surface.sign_at_position(edgeNodeCoords[0]) == -interface_surface.sign_at_position(edgeNodeCoords[1]))
  {
    const double loc = interface_crossing_position(interface, edgeNodeCoords);
    const stk::math::Vector3d intersectionPoint = (1.-loc)*edgeNodeCoords[0] + loc*edgeNodeCoords[1];
    return krino::intersection_point_is_real(intersectionPoint, all_cutting_surfaces, {interface.first_ls(), interface.second_ls()});
  }
  return false;
}

int One_LS_Per_Phase_Cutter::get_ls_per_interface_phase_at_location(const stk::math::Vector3d & pCoords) const
{
  const std::set<int> optimalPhases = determine_optimal_phases_at_location(pCoords, cutting_surfaces);
  STK_ThrowRequireMsg(optimalPhases.size()==1, "Unexpected phase configuration with " << optimalPhases.size() << " optimal phases when evaluated phase at " << pCoords << "\n" << visualize());
  return *optimalPhases.begin();
}

bool One_LS_Per_Phase_Cutter::intersection_point_is_real(const stk::math::Vector3d & intersectionPoint,
    const std::vector<int> & sortedDomains) const
{
  return krino::intersection_point_is_real(intersectionPoint, all_cutting_surfaces, sortedDomains);
}

One_LS_Per_Phase_Cutter::One_LS_Per_Phase_Cutter(const MasterElement & masterElem,
    const std::vector<const CDFEM_Parent_Edge *> & parentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
  : Element_Cutter(masterElem)
{
  std::set<int> edgePhases = get_phases_present_on_edges_and_interior(parentEdges);
  for (int phase1 : edgePhases)
  {
    for (int phase2 : edgePhases)
    {
      if (phase2 > phase1)
      {
        const InterfaceID interface(phase1, phase2);
        all_cutting_surfaces[interface] = attempt_to_build_cutting_surface(masterElem, parentEdges, areParentEdgesOrientedSameAsElementEdges, intersectingPlanesDiagonalPicker, interface);
      }
    }
  }

  std::set<InterfaceID> interfacesWithRealCrossings;
  for (auto && interface : get_interfaces_present(parentEdges))
    interfacesWithRealCrossings.insert(interface);

  const std::vector<const std::vector<int>*> emptyNodesSnappedDomains;
  add_interfaces_that_have_uncaptured_intersection_within_element(myTopology.base(), emptyNodesSnappedDomains, all_cutting_surfaces, interfacesWithRealCrossings);

  for (auto && interface : interfacesWithRealCrossings)
    cutting_surfaces[interface] = all_cutting_surfaces.at(interface);
}

}
