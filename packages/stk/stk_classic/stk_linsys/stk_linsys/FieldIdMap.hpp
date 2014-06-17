/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_linsys_FieldIdMap_hpp
#define stk_linsys_FieldIdMap_hpp

#include <map>
#include <stk_mesh/base/Field.hpp>

namespace stk_classic {
namespace linsys {

/** Mappings from stk_classic::mesh::Field objects to integer ids used by fei objects.
 */
typedef std::map<const stk_classic::mesh::FieldBase*,int> FieldIdMap;

}//namespace linsys
}//namespace stk_classic

#endif
