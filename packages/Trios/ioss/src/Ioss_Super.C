/*--------------------------------------------------------------------*/
/*    Copyright 2000, 2002 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Super.h>
#include <Ioss_ElementVariableType.h>
#include <cstring>
#include <cstdlib>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Super : public ElementVariableType
  {
  public:
    St_Super(const std::string &my_name, int node_count)
      : ElementVariableType(my_name, node_count) {}
  };
}

void Ioss::Super::factory() {}

// ========================================================================
Ioss::Super::Super(const std::string &my_name, int node_count)
  : Ioss::ElementTopology(my_name, "Unknown"), nodeCount(node_count),
    storageType(new St_Super(my_name, node_count))
{}

Ioss::Super::~Super()
{
  delete storageType;
}

void Ioss::Super::make_super(const std::string &type)
{
  // Decode name to determine number of nodes...
  std::string node_count_str = type.substr(5);
  int node_count = std::atoi(node_count_str.c_str());
  new Ioss::Super(type, node_count);
}

int Ioss::Super::parametric_dimension() const {return  3;}
int Ioss::Super::spatial_dimension()    const {return  3;}
int Ioss::Super::order()                const {return  1;}

int Ioss::Super::number_corner_nodes() const {return nodeCount;}
int Ioss::Super::number_nodes()        const {return nodeCount;}
int Ioss::Super::number_edges()        const {return 0;}
int Ioss::Super::number_faces()        const {return 0;}

int Ioss::Super::number_nodes_edge(int /* edge */) const {return  0;}

int Ioss::Super::number_nodes_face(int face) const {return 0;}
int Ioss::Super::number_edges_face(int face) const {return 0;}

Ioss::IntVector Ioss::Super::edge_connectivity(int edge_number) const
{
  Ioss::IntVector connectivity(0);
  return connectivity;
}

Ioss::IntVector Ioss::Super::face_connectivity(int face_number) const
{
  Ioss::IntVector connectivity(0);
  return connectivity;
}

Ioss::IntVector Ioss::Super::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Super::face_type(int face_number) const
{
  return Ioss::ElementTopology::factory("unknown");
}

Ioss::ElementTopology* Ioss::Super::edge_type(int edge_number) const
{
  return Ioss::ElementTopology::factory("unknown");
}

Ioss::IntVector Ioss::Super::face_edge_connectivity(int face_number) const
{
  Ioss::IntVector fcon(0);
  return fcon;
}
