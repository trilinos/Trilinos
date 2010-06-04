/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Quad8.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
#include <Ioss_ElementVariableType.h>
namespace Ioss {
  class St_Quad8 : public ElementVariableType
  {
  public:
    static void factory() {static St_Quad8 registerThis;}

  protected:
    St_Quad8()
      : ElementVariableType("quad8", 8) {}
  };
}
// ========================================================================
Ioss::Quad8 Ioss::Quad8::instance_;

namespace {
  struct Constants {
    static const int nnode     = 8;
    static const int nedge     = 4;
    static const int nedgenode = 3;
    static const int nface     = 0;
    static int edge_node_order[nedge][nedgenode];
  };

  // Edge numbers are zero-based [0..number_edges)
  int Constants::edge_node_order[nedge][nedgenode] = // [edge][edge_node]
    { {0,1,4}, {1,2,5}, {2,3,6}, {3,0,7} };
}

void Ioss::Quad8::factory()
{
  static Ioss::Quad8 registerThis;
  Ioss::St_Quad8::factory();
}

Ioss::Quad8::Quad8()
  : Ioss::ElementTopology("quad8", "Quadrilateral_8")
{
  Ioss::ElementTopology::alias("quad8", "Solid_Quad_8_2D");
  Ioss::ElementTopology::alias("quad8", "Face_Quad_8_3D");
  Ioss::ElementTopology::alias("quad8", "quadface8");
}

Ioss::Quad8::~Quad8() {}

int Ioss::Quad8::parametric_dimension()           const {return  2;}
int Ioss::Quad8::spatial_dimension()           const {return  2;}
int Ioss::Quad8::order()               const {return  2;}

int Ioss::Quad8::number_corner_nodes() const {return     4;}
int Ioss::Quad8::number_nodes()        const {return Constants::nnode;}
int Ioss::Quad8::number_edges()        const {return Constants::nedge;}
int Ioss::Quad8::number_faces()        const {return Constants::nface;}

int Ioss::Quad8::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Quad8::number_nodes_face(int /* face */) const
{
  return 0;
}
int Ioss::Quad8::number_edges_face(int /* face */) const
{
  return 0;
}

Ioss::IntVector Ioss::Quad8::edge_connectivity(int edge_number) const
{
  assert(edge_number > 0 && edge_number <= number_edges());
  Ioss::IntVector connectivity(Constants::nedgenode);
  assert(edge_number > 0 && edge_number <= Constants::nedge);

  for (int i=0; i < Constants::nedgenode; i++)
    connectivity[i] = Constants::edge_node_order[edge_number-1][i];

  return connectivity;
}

Ioss::IntVector Ioss::Quad8::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Quad8::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Quad8::face_type(int /* face_number */) const
{
  return (Ioss::ElementTopology*)NULL;
}

Ioss::ElementTopology* Ioss::Quad8::edge_type(int edge_number) const
{
  assert(edge_number >= 0 && edge_number <= number_edges());
  return Ioss::ElementTopology::factory("edge2d3");
}
