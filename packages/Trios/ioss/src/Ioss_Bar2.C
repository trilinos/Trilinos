/*--------------------------------------------------------------------*/
/*    Copyright 2004, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Bar2.h>
#include <Ioss_ElementVariableType.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Bar2 : public ElementVariableType
  {
  public:
    static void factory() {static St_Bar2 registerThis;}

  protected:
    St_Bar2()
      : ElementVariableType("bar2", 2) {}
  };
}
// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 2;
    static const int nedge     = 1;
    static const int nedgenode = 2;
    static const int nface     = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
  };
}

void Ioss::Bar2::factory()
{
  static Ioss::Bar2 registerThis;
  Ioss::St_Bar2::factory();
}

Ioss::Bar2::Bar2()
  : Ioss::ElementTopology("bar2", "Beam_2")
{
  Ioss::ElementTopology::alias("bar2", "Rod_2_3D");
  Ioss::ElementTopology::alias("bar2", "rod2");
  Ioss::ElementTopology::alias("bar2", "rod");
  Ioss::ElementTopology::alias("bar2", "beam2");
  Ioss::ElementTopology::alias("bar2", "bar");
  Ioss::ElementTopology::alias("bar2", "truss");
  Ioss::ElementTopology::alias("bar2", "truss2");
  Ioss::ElementTopology::alias("bar2", "beam");
  Ioss::ElementTopology::alias("bar2", "rod3d2");
  Ioss::ElementTopology::alias("bar2", "Rod_2_2D");
  Ioss::ElementTopology::alias("bar2", "rod2d2");
  Ioss::ElementTopology::alias("bar2", "beam-r2");
}


Ioss::Bar2::~Bar2() {}

int Ioss::Bar2::parametric_dimension() const {return  1;}
int Ioss::Bar2::spatial_dimension()    const {return  3;}
int Ioss::Bar2::order()                const {return  1;}

int Ioss::Bar2::number_corner_nodes()  const {return 2;}
int Ioss::Bar2::number_nodes()         const {return Constants::nnode;}
int Ioss::Bar2::number_edges()         const {return Constants::nedge;}
int Ioss::Bar2::number_faces()         const {return Constants::nface;}

int Ioss::Bar2::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Bar2::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

int Ioss::Bar2::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfaceedge;
}

Ioss::IntVector Ioss::Bar2::edge_connectivity(int /* edge_number */) const
{
  Ioss::IntVector connectivity(Constants::nedgenode);
  connectivity[0] = 0;
  connectivity[1] = 1;
  return connectivity;
}

Ioss::IntVector Ioss::Bar2::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Bar2::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Bar2::face_type(int /* face_number */) const
{ return (Ioss::ElementTopology*)NULL; }

Ioss::ElementTopology* Ioss::Bar2::edge_type(int /* edge_number */) const
{ return Ioss::ElementTopology::factory("edge2"); }
