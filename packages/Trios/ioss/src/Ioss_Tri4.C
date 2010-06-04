/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Tri4.h>
#include <Ioss_ElementVariableType.h>
#include <string>

#include <assert.h>

// ========================================================================
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Tri4 : public ElementVariableType
  {
  public:
    static void factory() {static St_Tri4 registerThis;}
  protected:
    St_Tri4()
      : ElementVariableType("tri4", 4) {}
  };
}
//------------------------------------------------------------------------
namespace {
  struct Constants {
    static const int nnode = 4;
    static const int nedge = 3;
    static const int nedgenode = 2;
    static const int nface = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
    static int edge_node_order[nedge][nedgenode];
  };
  // Edge numbers are zero-based [0..number_edges)
  int Constants::edge_node_order[nedge][nedgenode] = // [edge][edge_node]
    { {0,1}, {1,2}, {2,0} };
}

void Ioss::Tri4::factory()
{
  static Ioss::Tri4 registerThis;
  Ioss::St_Tri4::factory();
}

Ioss::Tri4::Tri4()
  : Ioss::ElementTopology("tri4", "Triangle_4")
{
  Ioss::ElementTopology::alias("tri4", "triangle4");
  Ioss::ElementTopology::alias("tri4", "Solid_Tri_4_2D");
  Ioss::ElementTopology::alias("tri4", "Face_Tri_4_3D");
  Ioss::ElementTopology::alias("tri4", "triface4");
}

Ioss::Tri4::~Tri4() {}

int Ioss::Tri4::parametric_dimension()           const {return  2;}
int Ioss::Tri4::spatial_dimension()           const {return  2;}
int Ioss::Tri4::order()               const {return  1;}

int Ioss::Tri4::number_corner_nodes() const {return     3;}
int Ioss::Tri4::number_nodes()        const {return Constants::nnode;}
int Ioss::Tri4::number_edges()        const {return Constants::nedge;}
int Ioss::Tri4::number_faces()        const {return Constants::nface;}

int Ioss::Tri4::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}
int Ioss::Tri4::number_nodes_face(int /* face */) const {return  Constants::nfacenode;}
int Ioss::Tri4::number_edges_face(int /* face */) const {return  Constants::nfaceedge;}

Ioss::IntVector Ioss::Tri4::edge_connectivity(int edge_number) const
{
  assert(edge_number > 0 && edge_number <= number_edges());
  Ioss::IntVector connectivity(Constants::nedgenode);

  for (int i=0; i < Constants::nedgenode; i++)
    connectivity[i] = Constants::edge_node_order[edge_number-1][i];

  return connectivity;
}

Ioss::IntVector Ioss::Tri4::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Tri4::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());

  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;

  return connectivity;
}

Ioss::ElementTopology* Ioss::Tri4::face_type(int /* face_number */) const
{
  return (Ioss::ElementTopology*)NULL;
}

Ioss::ElementTopology* Ioss::Tri4::edge_type(int edge_number) const
{
  assert(edge_number >= 0 && edge_number <= number_edges());
  return Ioss::ElementTopology::factory("edge2d2");
}
