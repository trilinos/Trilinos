// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef TextMesh_hpp
#define TextMesh_hpp

#include "stk_util/stk_config.h"
#include <string>
#include <vector>
namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace unit_test_util
{

std::string get_full_text_mesh_desc(const std::string& textMeshConnectivityDesc, const std::vector<double>& coordVec);
void setup_text_mesh(stk::mesh::BulkData& bulkData, const std::string& meshDesc);


inline std::string line_2_1d_mesh_desc() {
  return "0,1,LINE_2_1D,1,2,block_1"
         "|dimension: 1"
         "|coordinates: 0,1";
}

inline std::string line_3_1d_mesh_desc() {
  return "0,1,LINE_3_1D,1,2,3,block_1"
         "|dimension: 1"
         "|coordinates: 0,1,0.5";
}

inline std::string beam_2_mesh_desc() {
  return "0,1,BEAM_2,1,2,block_1"
         "|coordinates: 0,0,0, 1,0,0";
}

inline std::string beam_3_mesh_desc() {
  return "0,1,BEAM_3,1,2,3,block_1"
         "|coordinates: 0,0,0, 1,0,0, 0.5,0,0";
}

inline std::string shell_line_2_mesh_desc() {
  return "0,1,SHELL_LINE_2,1,2,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0";
}

inline std::string shell_line_3_mesh_desc() {
  return "0,1,SHELL_LINE_3,1,2,3,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 0.5,0";
}

inline std::string tri_3_2d_mesh_desc() {
  return "0,1,TRI_3_2D,1,2,3,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 0,1";
}

inline std::string tri_6_2d_mesh_desc() {
  return "0,1,TRI_6_2D,1,2,3,4,5,6,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 0,1, 0.5,0, 0.5,0.5, 0,0.5";
}

inline std::string quad_4_2d_mesh_desc() {
  return "0,1,QUAD_4_2D,1,2,3,4,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 1,1, 0,1";
}

inline std::string quad_8_2d_mesh_desc() {
  return "0,1,QUAD_8_2D,1,2,3,4,5,6,7,8,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 1,1, 0,1, 0.5,0, 1,0.5, 0.5,1, 0,0.5";
}

inline std::string quad_9_2d_mesh_desc() {
  return "0,1,QUAD_9_2D,1,2,3,4,5,6,7,8,9,block_1"
         "|dimension: 2"
         "|coordinates: 0,0, 1,0, 1,1, 0,1, 0.5,0, 1,0.5, 0.5,1, 0,0.5, 0.5,0.5";
}

inline std::string shell_tri_3_mesh_desc() {
  return "0,1,SHELL_TRI_3,1,2,3,block_1"
         "|coordinates: 0,0,0, 1,0,0, 0,1,0";
}

inline std::string shell_tri_6_mesh_desc() {
  return "0,1,SHELL_TRI_6,1,2,3,4,5,6,block_1"
         "|coordinates: 0,0,0, 1,0,0, 0,1,0, 0.5,0,0, 0.5,0.5,0, 0,0.5,0";
}

inline std::string shell_quad_4_mesh_desc() {
  return "0,1,SHELL_QUAD_4,1,2,3,4,block_1"
         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0";
}

inline std::string shell_quad_8_mesh_desc() {
  return "0,1,SHELL_QUAD_8,1,2,3,4,5,6,7,8,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0 ";
}

inline std::string shell_quad_9_mesh_desc() {
  return "0,1,SHELL_QUAD_9,1,2,3,4,5,6,7,8,9,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0, 0.5,0.5,0.0";
}

inline std::string tet_4_mesh_desc() {
  return "0,1,TET_4,1,2,3,4,block_1"
         "|coordinates: 0,0,0, 1,0,0, 0,1,0, 0,0,1";
}

inline std::string tet_10_mesh_desc() {
  return "0,1,TET_10,1,2,3,4,5,6,7,8,9,10,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, "
                       "0.5,0.0,0.0, 0.5,0.5,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5, 0.5,0.0,0.5, 0.0,0.5,0.5";
}

inline std::string pyr_5_mesh_desc() {
  return "0,1,PYRAMID_5,1,2,3,4,5,block_1"
         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1";
}

inline std::string pyr_13_mesh_desc() {
  return "0,1,PYRAMID_13,1,2,3,4,5,6,7,8,9,10,11,12,13,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0, "
                       "0.0,0.0,0.5, 0.5,0.0,0.5, 0.5,0.5,0.5, 0.0,0.5,0.5  ";
}

inline std::string wedge_6_mesh_desc() {
  return "0,1,WEDGE_6,1,2,3,4,5,6,block_1"
         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,0,1, 1,0,1, 1,1,1";
}

inline std::string wedge_15_mesh_desc() {
  return "0,1,WEDGE_15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, "
                       "0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,0.5,0.0, "
                       "0.0,0.0,0.5, 1.0,0.0,0.5, 1.0,1.0,0.5, "
                       "0.5,0.0,1.0, 1.0,0.5,1.0, 0.5,0.5,1.0  ";
}

inline std::string hex_8_mesh_desc() {
  return "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";
}

inline std::string hex_20_mesh_desc() {
  return "0,1,HEX_20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, "
                       "0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, 0.0,1.0,1.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0, "
                       "0.0,0.0,0.5, 1.0,0.0,0.5, 1.0,1.0,0.5, 0.0,1.0,0.5, "
                       "0.5,0.0,1.0, 1.0,0.5,1.0, 0.5,1.0,1.0, 0.0,0.5,1.0  ";
}

inline std::string hex_27_mesh_desc() {
  return "0,1,HEX_27,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, "
                       "0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, 0.0,1.0,1.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0, "
                       "0.0,0.0,0.5, 1.0,0.0,0.5, 1.0,1.0,0.5, 0.0,1.0,0.5, "
                       "0.5,0.0,1.0, 1.0,0.5,1.0, 0.5,1.0,1.0, 0.0,0.5,1.0, "
                       "0.5,0.5,0.5, "
                       "0.5,0.5,0.0, 0.5,0.5,1.0, 0.0,0.5,0.5, 1.0,0.5,0.5, 0.5,0.0,0.5, 0.5,1.0,0.5";
}

inline std::string line_2_edge_mesh_desc() {
  return shell_tri_3_mesh_desc() + "|sideset:name=surface_1; data=1,3";
}

inline std::string line_3_edge_mesh_desc() {
  return shell_tri_6_mesh_desc() + "|sideset:name=surface_1; data=1,3";
}

inline std::string tri_3_face_mesh_desc() {
  return tet_4_mesh_desc() + "|sideset:name=surface_1; data=1,1";
}

inline std::string tri_6_face_mesh_desc() {
  return tet_10_mesh_desc() + "|sideset:name=surface_1; data=1,1";
}

inline std::string quad_4_face_mesh_desc() {
  return hex_8_mesh_desc() +  "|sideset:name=surface_1; data=1,1";
}

inline std::string quad_8_face_mesh_desc() {
  return hex_20_mesh_desc() +  "|sideset:name=surface_1; data=1,1";
}

inline std::string quad_9_face_mesh_desc() {
  return hex_27_mesh_desc() +  "|sideset:name=surface_1; data=1,1";
}


} // namespace unit_test_util
} // namespace stk

#endif
