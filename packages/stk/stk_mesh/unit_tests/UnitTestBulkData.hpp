/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef unit_tests_UnitTestBulkData_hpp
#define unit_tests_UnitTestBulkData_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class UnitTestBulkData {
public:
  static void testBulkData( ParallelMachine );
  static void testCreateMore( ParallelMachine );
  static void testChangeOwner_nodes( ParallelMachine );
  static void testChangeOwner_loop( ParallelMachine );
  static void testChangeOwner_box( ParallelMachine );
  static void testChangeParts( ParallelMachine );
  static void testChangeParts_loop( ParallelMachine );
  static void testDestroy_nodes( ParallelMachine );
  static void testDestroy_loop( ParallelMachine );

  /** Generate a box mesh which is globally ( ngx X ngy X ngz )
   *  elements where:
   *    ngx = root_box[0][1] - root_box[0][0] ;
   *    ngy = root_box[1][1] - root_box[1][0] ;
   *    ngz = root_box[2][1] - root_box[2][0] ;
   *
   *  The box is partitioned via recursive coordinate bisection
   *  and the resulting local box are given in 'local_box'.
   */
  static void generate_boxes( BulkData  & mesh ,
                              const bool  generate_aura ,
                              const int   root_box[][2] ,
                                    int   local_box[][2] );

  static bool modification_end( BulkData & , bool );
};

}
}

#endif // unit_tests_UnitTestBulkData_hpp
