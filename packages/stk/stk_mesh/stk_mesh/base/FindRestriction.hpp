/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_base_FindRestriction_hpp
#define stk_mesh_base_FindRestriction_hpp

#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase, etc
#include <stk_mesh/base/Types.hpp>      // for EntityRank, PartVector
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {

//Given a field and a vector of parts, determine whether the field has a restriction
//for the part vector. (Common usage is to provide the part-vector from a bucket; i.e.,
//determine whether the field should be allocated in the bucket.)
const FieldBase::Restriction& find_restriction(const FieldBase& field,
                                               EntityRank erank,
                                               const PartVector& parts);

const FieldBase::Restriction& find_and_check_restriction(const FieldBase& field,
                                                         EntityRank erank,
                                                         const PartVector& parts);

const FieldBase::Restriction& find_restriction(const FieldBase& field,
                                               EntityRank erank,
                                               const Part & part);

} // namespace mesh
} // namespace stk

#endif // stk_mesh_base_FindRestriction_hpp
