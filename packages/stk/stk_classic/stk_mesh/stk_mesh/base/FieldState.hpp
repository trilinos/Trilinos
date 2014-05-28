/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


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
  StateNone = 0,  ///< \brief State of a field with one state
  StateNew  = 0,  ///< \brief Newest state of a field with two states
  StateNP1  = 0,  ///< \brief Newest state of a field with three+ states
  StateOld  = 1,  ///< \brief Previous state of a field with two states
  StateN    = 1,  ///< \brief Previous state of a field with three+ states
  StateNM1  = 2,  ///< \brief Previous-1 state of a field with three+ states
  StateNM2  = 3,  ///< \brief Previous-2 state of a field with four+ states
  StateNM3  = 4,  ///< \brief Previous-3 state of a field with five+ states
  StateNM4  = 5   ///< \brief Previous-4 state of a field with six states
};

/** \brief Maximum number of states that a \ref stk::mesh::Field "field"
 *         can have.
 */
enum { MaximumFieldStates = 6 };

/** \} */

} //namespace mesh
} //namespace stk

#endif //stk_mesh_FieldState_hpp
