/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Unknown.h>
#include <Ioss_ElementVariableType.h>
#include <string>

#include <assert.h>

// ========================================================================
namespace Ioss {
  class St_Unknown : public ElementVariableType
  {
  public:
    static void factory();

  protected:
    St_Unknown()
      : ElementVariableType("unknown", 0) {}
  };
}
void Ioss::St_Unknown::factory()
{ static Ioss::St_Unknown registerThis; }

// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 0;
    static const int nedge     = 0;
    static const int nedgenode = 0;
    static const int nface     = 0;
    static const int nfacenode = 0;
  };
}

void Ioss::Unknown::factory()
{
  static Ioss::Unknown registerThis;
  Ioss::St_Unknown::factory();
}

Ioss::Unknown::Unknown()
  : Ioss::ElementTopology("unknown", "unknown") {}

Ioss::Unknown::~Unknown() {}

int Ioss::Unknown::parametric_dimension()           const {return  0;}
int Ioss::Unknown::spatial_dimension()           const {return  0;}
int Ioss::Unknown::order()               const {return  0;}

int Ioss::Unknown::number_corner_nodes() const {return number_nodes();}
int Ioss::Unknown::number_nodes()        const {return Constants::nnode;}
int Ioss::Unknown::number_edges()        const {return Constants::nedge;}
int Ioss::Unknown::number_faces()        const {return Constants::nface;}

int Ioss::Unknown::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Unknown::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

int Ioss::Unknown::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

Ioss::IntVector Ioss::Unknown::edge_connectivity(int edge_number) const
{
  Ioss::IntVector connectivity;
  assert(edge_number >= 0 && edge_number <= Constants::nedge);
  return connectivity;
}

Ioss::IntVector Ioss::Unknown::face_connectivity(int face_number) const
{
  assert(face_number >= 0 && face_number <= number_faces());
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Unknown::element_connectivity() const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Unknown::face_type(int face_number) const
{
  // face_number == 0 returns topology for all faces if
  // all faces are the same topology; otherwise, returns NULL
  // face_number is 1-based.

  assert(face_number >= 0 && face_number <= number_faces());
  return Ioss::ElementTopology::factory("unknown");
}

Ioss::ElementTopology* Ioss::Unknown::edge_type(int edge_number) const
{
  // edge_number == 0 returns topology for all edges if
  // all edges are the same topology; otherwise, returns NULL
  // edge_number is 1-based.

  assert(edge_number >= 0 && edge_number <= number_edges());
  return Ioss::ElementTopology::factory("unknown");
}
