/*--------------------------------------------------------------------*/
/*    Copyright 2001 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Tet10.h>
#include <Ioss_ElementVariableType.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Tet10 : public ElementVariableType
  {
  public:
    static void factory() {static St_Tet10 registerThis;}

  protected:
    St_Tet10()
      : ElementVariableType("tetra10", 10) {}
  };
}

// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 10;
    static const int nedge     =  6;
    static const int nedgenode =  3;
    static const int nface     =  4;
    static const int nfacenode =  6;
    static const int nfaceedge =  3;
    static int edge_node_order[nedge][nedgenode];
    static int face_node_order[nface][nfacenode];
    static int face_edge_order[nface][nfaceedge];
    static int nodes_per_face[nface+1];
    static int edges_per_face[nface+1];
  };
  
  // Edge numbers are zero-based [0..number_edges)
  int Constants::edge_node_order[nedge][nedgenode] = // [edge][edge_node]
    { {0,1,4}, {1,2,5}, {2,0,6}, {0,3,7}, {1,3,8}, {2,3,9} };
  
  // Face numbers are zero-based [0..number_faces)
  int Constants::face_node_order[nface][nfacenode] = // [face][face_node]
    { {0,1,3,4,8,7}, {1,2,3,5,9,8},
      {0,3,2,7,9,6}, {0,2,1,6,5,4} };
  
  int Constants::face_edge_order[nface][nfaceedge] = // [face][face_edge]
    { {0,4,3}, {1,5,4}, {3,5,2}, {2,1,0} };

  // face 0 returns number of nodes for all faces if homogenous
  //        returns -1 if faces have differing topology
  int Constants::nodes_per_face[nface+1] =
    {6,6,6,6,6};
  
  // face 0 returns number of edges for all faces if homogenous
  //        returns -1 if faces have differing topology
  int Constants::edges_per_face[nface+1] =
    {3,3,3,3,3};
}

void Ioss::Tet10::factory()
{
  static Ioss::Tet10 registerThis;
  Ioss::St_Tet10::factory();
}

Ioss::Tet10::Tet10()
  : Ioss::ElementTopology("tetra10", "Tetrahedron_10")
{
  Ioss::ElementTopology::alias("tetra10", "tet10");
  Ioss::ElementTopology::alias("tetra10", "Solid_Tet_10_3D");}

Ioss::Tet10::~Tet10() {}

int Ioss::Tet10::parametric_dimension()           const {return  3;}
int Ioss::Tet10::spatial_dimension()           const {return  3;}
int Ioss::Tet10::order()               const {return  2;}

int Ioss::Tet10::number_corner_nodes() const {return 4;}
int Ioss::Tet10::number_nodes()        const {return Constants::nnode;}
int Ioss::Tet10::number_edges()        const {return Constants::nedge;}
int Ioss::Tet10::number_faces()        const {return Constants::nface;}

int Ioss::Tet10::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Tet10::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return  Constants::nodes_per_face[face];
}

int Ioss::Tet10::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::edges_per_face[face];
}

Ioss::IntVector Ioss::Tet10::edge_connectivity(int edge_number) const
{
  assert(edge_number > 0 && edge_number <= Constants::nedge);
  Ioss::IntVector connectivity(Constants::nedgenode);

  for (int i=0; i < Constants::nedgenode; i++)
    connectivity[i] = Constants::edge_node_order[edge_number-1][i];

  return connectivity;
}

Ioss::IntVector Ioss::Tet10::face_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= number_faces());
  Ioss::IntVector connectivity(Constants::nodes_per_face[face_number]);

  for (int i=0; i < Constants::nodes_per_face[face_number]; i++)
    connectivity[i] = Constants::face_node_order[face_number-1][i];

  return connectivity;
}

Ioss::IntVector Ioss::Tet10::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Tet10::face_type(int face_number) const
{
  // face_number == 0 returns topology for all faces if
  // all faces are the same topology; otherwise, returns NULL
  // face_number is 1-based.

  assert(face_number >= 0 && face_number <= number_faces());
//  return Ioss::ElementTopology::factory("triface6");
  return Ioss::ElementTopology::factory("tri6");
}

Ioss::ElementTopology* Ioss::Tet10::edge_type(int edge_number) const
{
  // edge_number == 0 returns topology for all edges if
  // all edges are the same topology; otherwise, returns NULL
  // edge_number is 1-based.

  assert(edge_number >= 0 && edge_number <= number_edges());
  return Ioss::ElementTopology::factory("edge3");
}

Ioss::IntVector Ioss::Tet10::face_edge_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= Constants::nface);

  int nface_edge = number_edges_face(face_number);
  Ioss::IntVector fcon(nface_edge);

  for (int i=0; i < nface_edge; i++)
    fcon[i] = Constants::face_edge_order[face_number-1][i];

  return fcon;
}
