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

#ifndef MORKON_EXP_INTERFACE_HOST_SIDE_ADAPTER_H
#define MORKON_EXP_INTERFACE_HOST_SIDE_ADAPTER_H

namespace morkon_exp {

template <unsigned int DIM>
struct Interface_HostSideAdapter
{
  struct NodeInfo
  {
    global_idx_t                   m_id;
    InterfaceBase::SideEnum      m_side;
    double[DIM]                m_coords;
  };

  struct FaceInfo
  {
    global_idx_t                     m_id;
    InterfaceBase::SideEnum        m_side;
    std::vector<global_idx_t>     m_nodes;
  };

  typedef typename std::map<typename global_idx_t, NodeInfo>   node_map_type;
  typedef typename std::map<typename global_idx_t, FaceInfo> face_map_type;

  node_map_type   m_nodes;
  face_map_type m_faces;
};


}
#endif

