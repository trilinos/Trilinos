/*------------------------------------------------------------------------*/
/*  Copyright 2013 Sandia Corporation.                                    */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_IO_DATABASEPURPOSE_HPP
#define STK_IO_DATABASEPURPOSE_HPP
namespace stk {
  namespace io {
    enum DatabasePurpose {
      PURPOSE_UNKNOWN=0,
      WRITE_RESULTS=1,
      WRITE_RESTART=2,
      READ_MESH=4,
      READ_RESTART=8
    };
  }
}
#endif
