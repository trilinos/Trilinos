/*--------------------------------------------------------------------*/
/*    Copyright 2004 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_ShellLine2D3.h>
#include <Ioss_ElementVariableType.h>
#include <string>

#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_ShellLine2D3 : public ElementVariableType
  {
  public:
    static void factory() {static St_ShellLine2D3 registerThis;}

  protected:
    St_ShellLine2D3()
      : ElementVariableType("shellline2d3", 3) {}
  };
}
// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 3;
    static const int nedge     = 2;
    static const int nedgenode = 3;
    static const int nface     = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
  };
}

void Ioss::ShellLine2D3::factory()
{
  static Ioss::ShellLine2D3 registerThis;
  Ioss::St_ShellLine2D3::factory();
}

Ioss::ShellLine2D3::ShellLine2D3()
  : Ioss::ElementTopology("shellline2d3", "ShellLine_3")
{
  Ioss::ElementTopology::alias("shellline2d3", "Shell_Line_3_2D");
}

Ioss::ShellLine2D3::~ShellLine2D3() {}

int Ioss::ShellLine2D3::parametric_dimension()           const {return  1;}
int Ioss::ShellLine2D3::spatial_dimension()           const {return  2;}
int Ioss::ShellLine2D3::order()               const {return  2;}

int Ioss::ShellLine2D3::number_corner_nodes() const {return 2;}
int Ioss::ShellLine2D3::number_nodes()        const {return Constants::nnode;}
int Ioss::ShellLine2D3::number_edges()        const {return Constants::nedge;}
int Ioss::ShellLine2D3::number_faces()        const {return Constants::nface;}

int Ioss::ShellLine2D3::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::ShellLine2D3::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfacenode;
}

int Ioss::ShellLine2D3::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfaceedge;
}

Ioss::IntVector Ioss::ShellLine2D3::edge_connectivity(int edge_number) const
{
  Ioss::IntVector connectivity(Constants::nedgenode);
  if (edge_number == 1) {
    connectivity[0] = 0;
    connectivity[1] = 1;
    connectivity[2] = 2;
  } else {
    connectivity[0] = 1;
    connectivity[1] = 0;
    connectivity[2] = 2;
  }
  return connectivity;
}

Ioss::IntVector Ioss::ShellLine2D3::face_connectivity(int /* face_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::ShellLine2D3::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::ShellLine2D3::face_type(int /* face_number */) const
{ return (Ioss::ElementTopology*)NULL; }

Ioss::ElementTopology* Ioss::ShellLine2D3::edge_type(int /* edge_number */) const
{ return Ioss::ElementTopology::factory("edge2d3"); }
