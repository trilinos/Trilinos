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

#ifndef stk_mesh_FieldState_hpp
#define stk_mesh_FieldState_hpp

namespace stk {
namespace mesh {

  /** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Enumeration of states for multi-state
 *          \ref stk::mesh::Field "fields".
 *
 *  A field may be declared to have field data for multiple states.
 *  -  Field states <b>StateNone</b>, <b>StateNew</b>, or <b>StateNP1</b>
 *     refer to the current or newest of a field.
 *  -  Field states <b>StateOld</b> or <b>StateN</b> refer to
 *     the previous state of a field with two or more states.
 *   - The remaining field states <b>StateNM1</b>, <b>StateNM2</b>,
 *     <b>StateNM3</b>, <b>StateNM4</b> refer to prior states
 *     N-1, N-2, N-3, and N-4 accordingly.
 *
 * See Field.hpp for more.
 */
enum FieldState {
  StateNone    = 0,  ///< \brief State of a field with one state
  StateNew     = 0,  ///< \brief Newest state of a field with two or more states
  StateNP1     = 0,  ///< \brief Newest state of a field with two or more states
  StateOld     = 1,  ///< \brief Previous state of a field with two states
  StateN       = 1,  ///< \brief Previous state of a field with three+ states
  StateNM1     = 2,  ///< \brief Previous-1 state of a field with three+ states
  StateNM2     = 3,  ///< \brief Previous-2 state of a field with four+ states
  StateNM3     = 4,  ///< \brief Previous-3 state of a field with five+ states
  StateNM4     = 5,  ///< \brief Previous-4 state of a field with six states
  StateInvalid = 6
};

/** \brief Maximum number of states that a \ref stk::mesh::Field "field"
 *         can have.
 */
enum { MaximumFieldStates = 6 };

/** \} */

} //namespace mesh
} //namespace stk

#endif //stk_mesh_FieldState_hpp
