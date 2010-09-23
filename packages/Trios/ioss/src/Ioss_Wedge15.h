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

// -*- Mode: c++ -*-
#ifndef IOSS_Ioss_Wedge15_h
#define IOSS_Ioss_Wedge15_h

#include <Ioss_CodeTypes.h>
#include <Ioss_ElementTopology.h>
#include <string>

// STL Includes
#include <vector>

namespace Ioss {
  class Wedge15 : public Ioss::ElementTopology {

  public:
    static void factory();
    ~Wedge15();

    int spatial_dimension()           const;
    int parametric_dimension()       const;
    bool is_element()                 const {return true;}
    int order()               const;

    int number_corner_nodes() const;
    int number_nodes()        const;
    int number_edges()        const;
    int number_faces()        const;

    int number_nodes_edge(int edge=0)   const;
    int number_nodes_face(int face=0)   const;
    int number_edges_face(int face=0)   const;

    bool faces_similar()      const;

    Ioss::IntVector edge_connectivity(int edge_number) const;
    Ioss::IntVector face_connectivity(int face_number) const;
    Ioss::IntVector element_connectivity()             const;

    Ioss::IntVector face_edge_connectivity(int face_number) const;

    Ioss::ElementTopology* face_type(int face_number = 0) const;
    Ioss::ElementTopology* edge_type(int edge_number = 0) const;

  protected:
    Wedge15();

  private:
    Wedge15(const Wedge15&); // Do not implement
  };
}
#endif
