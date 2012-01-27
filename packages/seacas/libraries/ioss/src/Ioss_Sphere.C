// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_Sphere.h>
#include <Ioss_ElementVariableType.h>
#include <string>
#include <assert.h>

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_Sphere : public ElementVariableType
  {
  public:
    static void factory() {static St_Sphere registerThis;}

  protected:
    St_Sphere()
      : ElementVariableType("sphere", 1) {}
  };
}
// ========================================================================
namespace {
  struct Constants {
    static const int nnode     = 1;
    static const int nedge     = 0;
    static const int nedgenode = 0;
    static const int nface     = 0;
    static const int nfacenode = 0;
    static const int nfaceedge = 0;
  };
}
void Ioss::Sphere::factory()
{
  static Ioss::Sphere registerThis;
  Ioss::St_Sphere::factory();
}

Ioss::Sphere::Sphere()
  : Ioss::ElementTopology("sphere", "Particle")
{
  Ioss::ElementTopology::alias("sphere", "sphere1");
  Ioss::ElementTopology::alias("sphere", "particle");
  Ioss::ElementTopology::alias("sphere", "particles");
  Ioss::ElementTopology::alias("sphere", "sphere-mass");
  Ioss::ElementTopology::alias("sphere", "Particle_1_3D");
  Ioss::ElementTopology::alias("sphere", "Particle_1_2D");
  Ioss::ElementTopology::alias("sphere", "circle");  
  Ioss::ElementTopology::alias("sphere", "circle1");
}

Ioss::Sphere::~Sphere() {}

int Ioss::Sphere::parametric_dimension()           const {return  0;}
int Ioss::Sphere::spatial_dimension()           const {return  3;}
int Ioss::Sphere::order()               const {return  1;}

int Ioss::Sphere::number_corner_nodes() const {return number_nodes();}
int Ioss::Sphere::number_nodes()        const {return Constants::nnode;}
int Ioss::Sphere::number_edges()        const {return Constants::nedge;}
int Ioss::Sphere::number_faces()        const {return Constants::nface;}

int Ioss::Sphere::number_nodes_edge(int /* edge */) const {return  Constants::nedgenode;}

int Ioss::Sphere::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return  Constants::nfacenode;
}

int Ioss::Sphere::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nfaceedge;
}

Ioss::IntVector Ioss::Sphere::edge_connectivity(int /* edge_number */) const
{
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Sphere::face_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= number_faces());
  Ioss::IntVector connectivity;
  return connectivity;
}

Ioss::IntVector Ioss::Sphere::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i=0; i < number_nodes(); i++)
    connectivity[i] = i;
  return connectivity;
}

Ioss::ElementTopology* Ioss::Sphere::face_type(int face_number) const
{
  // face_number == 0 returns topology for all faces if
  // all faces are the same topology; otherwise, returns NULL
  // face_number is 1-based.

  assert(face_number >= 0 && face_number <= number_faces());
  return (Ioss::ElementTopology*)NULL;
}

Ioss::ElementTopology* Ioss::Sphere::edge_type(int edge_number) const
{
  // edge_number == 0 returns topology for all edges if
  // all edges are the same topology; otherwise, returns NULL
  // edge_number is 1-based.

  assert(edge_number >= 0 && edge_number <= number_edges());
  return NULL;
}
