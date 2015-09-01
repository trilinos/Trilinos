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

#include "morkon_unit_test_utils_TRI3.hpp"

using namespace morkon_exp;

const global_idx_t Mrk_2x2_TriangleInterfaceFixtureBase::node_gids[NumNodes] = {1,2,3,4,5,6,7,8};
const Mrk_InterfaceFixtureBase::interface_side_t Mrk_2x2_TriangleInterfaceFixtureBase::node_side[NumNodes] =
                        { InterfaceBase::NON_MORTAR_SIDE, InterfaceBase::NON_MORTAR_SIDE,
                          InterfaceBase::NON_MORTAR_SIDE, InterfaceBase::NON_MORTAR_SIDE,
                          InterfaceBase::MORTAR_SIDE, InterfaceBase::MORTAR_SIDE,
                          InterfaceBase::MORTAR_SIDE, InterfaceBase::MORTAR_SIDE };

const global_idx_t Mrk_2x2_TriangleInterfaceFixtureBase::face_gids[NumFaces] = {1,2,3,4};
const Mrk_InterfaceFixtureBase::interface_side_t Mrk_2x2_TriangleInterfaceFixtureBase::face_sides[NumFaces] =
                        { InterfaceBase::NON_MORTAR_SIDE, InterfaceBase::NON_MORTAR_SIDE,
                          InterfaceBase::MORTAR_SIDE, InterfaceBase::MORTAR_SIDE};
const global_idx_t Mrk_2x2_TriangleInterfaceFixtureBase::face_node_gids[NumFaces][NodesPerFace] =
                                                                            { {1,2,3}, {1,3,4}, {5,6,7}, {5,7,8}};



const double Mrk_2x2_aligned_TriangleInterfaceFixture::node_coords[NumNodes][3] =
                                                            { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                                                              {0,0,0.1}, {0,1,0.1}, {1,1,0.1}, {1,0,0.1} };

const double Mrk_2x4_offset_TriangleInterfaceFixture::node_coords[NumNodes][3] =
                                            { {0,0,0},         {1,0,0},         {1,1,0},         {0,1,0},
                                              {0.6, 0.6, 0.1}, {0.6, 1.6, 0.1}, {1.6, 1.6, 0.1}, {1.6, 0.6, 0.1} };
