/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SIERRA_MESH_ARRAY_TOPOLOGIES_HPP
#define STK_SIERRA_MESH_ARRAY_TOPOLOGIES_HPP

#include <string>

namespace sierra {
namespace mesh {

/** TopologyEnum can be used for specifying topology at run-time. */
enum TopologyEnum {
  InvalidTopologyEnum = -1,
  Tet4Enum            =  0,
  Hex8Enum            =  1,
  NodeEnum            =  2
};

/** Structs allow for describing each topology as a different compile-time type. */
struct InvalidTopology { static const int value = InvalidTopologyEnum; };
struct Node
{
	static const int value = NodeEnum;
};

struct Tet4
{
  static const int value = Tet4Enum;
  static const int num_sides = 4;
  static const int nodes_per_side = 3;
};

static const int Tet4_side_nodes[Tet4::num_sides][Tet4::nodes_per_side] =
{
  { 0, 1, 3 }, //side 1
  { 1, 2, 3 }, //side 2
  { 0, 3, 2 }, //side 3
  { 0, 2, 1 }, //side 4
};

struct Hex8
{
  static const int value = Hex8Enum;
  static const int num_sides = 6;
  static const int nodes_per_side = 4;
};

static const int Hex8_side_nodes[Hex8::num_sides][Hex8::nodes_per_side] =
{
  { 0, 1, 5, 4 }, //side 1
  { 1, 2, 6, 5 }, //side 2
  { 2, 3, 7, 6 }, //side 3
  { 0, 4, 7, 3 }, //side 4
  { 0, 3, 2, 1 }, //side 5
  { 4, 5, 6, 7 }  //side 6
};

inline
std::string topology_name(int topology)
{
  std::string name("unknown");
  switch(topology)
  {
  case Node::value: name = "node"; break;
  case Tet4::value: name = "tetra"; break;
  case Hex8::value: name = "hex"; break;
  default: break;
  }
  return name;
}

template<typename T>
struct num_nodes {
  static const int value = 0;
};

template<> struct num_nodes<Tet4> { static const int value = 4; };
template<> struct num_nodes<Hex8> { static const int value = 8; };

//runtime num-nodes for a given topology:
inline
int num_topology_nodes(int topo)
{
  switch(topo) {
   case Tet4::value: return num_nodes<Tet4>::value; break;
   case Hex8::value: return num_nodes<Hex8>::value; break;
   default: return 0;
  }
  return 0;
}

} // mesh
} // sierra

#endif

