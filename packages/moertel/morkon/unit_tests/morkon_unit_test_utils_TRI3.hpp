/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef MRK_UNIT_TEST_UTILS_HPP
#define MRK_UNIT_TEST_UTILS_HPP

#include "mrk_default_kokkos_device_type.hpp"
#include <mrk_api_classes.hpp>

namespace morkon_exp {

struct Mrk_InterfaceFixtureBase
{
  typedef Morkon_Manager<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                        manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>     interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                    interface_3d_ptr;
  typedef InterfaceBase::SideEnum                           interface_side_t;
};

struct Mrk_2x2_TriangleInterfaceFixtureBase : public Mrk_InterfaceFixtureBase
{
  static const size_t NumNodes = 8;
  static const global_idx_t node_gids[NumNodes];
  static const interface_side_t node_side[NumNodes];

  static const size_t NumFaces = 4;
  static const size_t NodesPerFace = TopoConsts<MRK_TRI3>::NODES_PER_FACE;
  static const global_idx_t face_gids[NumFaces];
  static const interface_side_t face_sides[NumFaces];
  static const global_idx_t face_node_gids[NumFaces][NodesPerFace];

  void connect_faces_to_nodes(interface_3d_ptr interface)
  {
    for (size_t face_i = 0; face_i < NumFaces; ++face_i)
    {
      interface->hsa_add_face(face_sides[face_i], face_gids[face_i], NodesPerFace, face_node_gids[face_i]);
    }
  }
};

struct Mrk_2x2_aligned_TriangleInterfaceFixture : public Mrk_2x2_TriangleInterfaceFixtureBase
{
  static const double node_coords[NumNodes][3];

  Mrk_2x2_aligned_TriangleInterfaceFixture(interface_3d_ptr interface)
  {
    for (size_t node_i = 0; node_i < NumNodes; ++node_i)
    {
      interface->hsa_add_node(node_side[node_i], node_gids[node_i], node_coords[node_i]);
    }
    connect_faces_to_nodes(interface);
  }
};

struct Mrk_2x2_offset_TriangleInterfaceFixture : public Mrk_2x2_TriangleInterfaceFixtureBase
{
  static const double node_coords[NumNodes][3];

  Mrk_2x2_offset_TriangleInterfaceFixture(interface_3d_ptr interface)
  {
    for (size_t node_i = 0; node_i < NumNodes; ++node_i)
    {
      interface->hsa_add_node(node_side[node_i], node_gids[node_i], node_coords[node_i]);
    }
    connect_faces_to_nodes(interface);
  }
};

}

#endif
