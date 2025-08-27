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
// 

#ifndef stk_mesh_CreateFaces_hpp
#define stk_mesh_CreateFaces_hpp

#include <stk_util/stk_config.h>

namespace stk {
  namespace mesh {

    class BulkData;
    class Selector;
    class Part;
    
    /** Create faces for all elements in "element_selector" and attach them to
         * existing elements.
         *
         * If connect_faces_to_edges is true, connect pre-existing edges to faces
         *
         * This is a parallel collective function (it should be called on all
         * processors at the same time
         *
         */
    void create_faces(  BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges );

    /** Create faces for all elements in "element_selector" and attach them to
     * existing elements.
     *
     * This is a parallel collective function (it should be called on all
     * processors at the same time
     *
     */
    void create_faces(  BulkData & mesh, const Selector & element_selector );

    /** Create faces for all elements in the mesh and attach them to
         * existing elements.
         *
         * If connect_faces_to_edges is true, connect pre-existing edges to faces
         *
         * This is a parallel collective function (it should be called on all
         * processors at the same time
         *
         */
    void create_faces( BulkData & mesh, bool connect_faces_to_edges );

    /** Create faces for all elements in the mesh and attach them to
     * existing elements.
     *
     * This is a parallel collective function (it should be called on all
     * processors at the same time
     *
     */
    void create_faces( BulkData & mesh );

    enum class FaceCreationBehavior {
        CREATE_FACES_FACE_CREATION_CLASSIC = 42,
        CREATE_FACES_FACE_CREATION_CURRENT = 73
    };

  }
}

#endif
