/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Field_hpp
#define stk_mesh_Field_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/FieldTraits.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \ingroup stk_mesh_module
 *  \brief  Field with defined data type and multi-dimensions (if any)
 */
template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class Field : public FieldBase {
public:

  /** \brief  Query this field for a given field state. */
  Field & field_of_state( FieldState input_state ) const {
    return static_cast<Field &>( * FieldBase::field_state(input_state) );
  }

private:

#ifndef DOXYGEN_COMPILE

  ~Field();
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

#endif /* DOXYGEN_COMPILE */
};


} // namespace mesh
} // namespace stk


#endif /* stk_mesh_Field_hpp */

