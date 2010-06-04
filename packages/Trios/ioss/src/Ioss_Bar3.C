/*--------------------------------------------------------------------*/
/*    Copyright 2004, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Bar3.h>
#include <Ioss_ElementVariableType.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Bar3 : public ElementVariableType
  {
  public:
    static void factory() {static St_Bar3 registerThis;}

  protected:
    St_Bar3()
      : ElementVariableType("bar3", 3) {}
  };
}
// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 3;
    static const int nedge     = 1;
    static const int nedgenode = 3;
    static const int nface     = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
  };
}

void Ioss::Bar3::factory()
{
  static Ioss::Bar3 registerThis;
  Ioss::St_Bar3::factory();
}

Ioss::Bar3::Bar3()
  : Ioss::ElementTopology("bar3", "Beam_3")
{
  Ioss::ElementTopology::alias("bar3", "Rod_3_3D");
  Ioss::ElementTopology::alias("bar3", "rod3");
  Ioss::ElementTopology::alias("bar3", "rod3d3");
  Ioss::ElementTopology::alias("bar3", "truss3");
  Ioss::ElementTopology::alias("bar3", "beam3");
  Ioss::ElementTopology::alias("bar3", "Rod_3_2D");
  Ioss::ElementTopology::alias("bar3", "rod2d3");
}

Ioss::Bar3::~Bar3() {}

int Ioss::Bar3::parametric_dimension()           const {return  1;}
int Ioss::Bar3::spatial_dimension()           const {return  3;}
int Ioss::Bar3::order()               const {return  2;}

int Ioss::Bar3::number_corner_nodes() const {return 2;}
int Ioss::Bar3::number_nodes()        const {return Constants::nnode;}
int Ioss::Bar3::number_edges()        const {return Constants::nedge;}
int Ioss::Bar3::number_faces()        const {return Constants::nface;}

int Ioss::Bar3::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Bar3::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

int Ioss::Bar3::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfaceedge;
}

Ioss::IntVector Ioss::Bar3::edge_connectivity(int /* edge_number */) const
{
  Ioss::IntVector connectivity(Constants::nedgenode);
  connectivity[0] = 0;
  connectivity[1] = 1;
  connectivity[2] = 2;
  return connectivity;
}

Ioss::IntVector Ioss::Bar3::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Bar3::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Bar3::face_type(int /* face_number */) const
{ return (Ioss::ElementTopology*)NULL; }

Ioss::ElementTopology* Ioss::Bar3::edge_type(int /* edge_number */) const
{ return Ioss::ElementTopology::factory("edge3"); }
