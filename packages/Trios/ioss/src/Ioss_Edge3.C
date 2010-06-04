/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Edge3.h>
#include <Ioss_ElementVariableType.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Edge3 : public ElementVariableType
  {
  public:
    static void factory() {static St_Edge3 registerThis;}

  protected:
    St_Edge3()
      : ElementVariableType("edge3", 3) {}
  };
}
// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 3;
    static const int nedge     = 0;
    static const int nedgenode = 3;
    static const int nface     = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
  };
}

void Ioss::Edge3::factory()
{
  static Ioss::Edge3 registerThis;
  Ioss::St_Edge3::factory();
}

Ioss::Edge3::Edge3()
  : Ioss::ElementTopology("edge3", "Line_3")
{
  Ioss::ElementTopology::alias("edge3", "edge3d3");
}

Ioss::Edge3::~Edge3() {}

int Ioss::Edge3::parametric_dimension()           const {return  1;}
int Ioss::Edge3::spatial_dimension()           const {return  3;}
int Ioss::Edge3::order()               const {return  2;}

int Ioss::Edge3::number_corner_nodes() const {return     2;}
int Ioss::Edge3::number_nodes()        const {return Constants::nnode;}
int Ioss::Edge3::number_edges()        const {return Constants::nedge;}
int Ioss::Edge3::number_faces()        const {return Constants::nface;}

int Ioss::Edge3::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Edge3::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

int Ioss::Edge3::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfaceedge;
}

Ioss::IntVector Ioss::Edge3::edge_connectivity(int /* edge_number */) const
{
  Ioss::IntVector connectivity(Constants::nedgenode);
  connectivity[0] = 0;
  connectivity[1] = 1;
  connectivity[2] = 2;
  return connectivity;
}

Ioss::IntVector Ioss::Edge3::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Edge3::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Edge3::face_type(int /* face_number */) const
{ return (Ioss::ElementTopology*)NULL; }

Ioss::ElementTopology* Ioss::Edge3::edge_type(int /* edge_number */) const
{ return Ioss::ElementTopology::factory("edge3"); }
