/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef unit_test_UnitTestModificationEndWrapper_hpp
#define unit_test_UnitTestModificationEndWrapper_hpp

#include <stk_mesh/base/BulkData.hpp>

// must be in same namespace as BulkData to make friending w/out exposure possible

namespace stk_classic {
namespace mesh {

class UnitTestModificationEndWrapper
{
 public:
  static bool wrap(stk_classic::mesh::BulkData& mesh, bool generate_aura);
};

} // namespace mesh

namespace unit_test {

bool modification_end_wrapper(stk_classic::mesh::BulkData& mesh, bool generate_aura);

} // namespace unit_test

} // namespace stk_classic

#endif
