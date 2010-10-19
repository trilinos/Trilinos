/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef Ioi_Events_h
#define Ioi_Events_h

#include <Ioss_DBUsage.h>

  //: Enumerations used in Publish/Subscribe classes
namespace Ioi {

enum EventState {
  DEFINE_MODEL   = 0,
  DEFINE_FIELDS  = 1,
  FIELDS_DEFINED = 2,
  DUMP_FIELDS    = 3,
  FIELDS_DUMPED  = 4,
  RESTORE_FIELDS = 5,
  FIELDS_RESTORED= 6
};

enum EventInterest {
  WRITE_RESTART   =  Ioss::WRITE_RESTART,
  READ_RESTART    =  Ioss::READ_RESTART,
  WRITE_RESULTS   =  Ioss::WRITE_RESULTS,
  READ_MODEL      =  Ioss::READ_MODEL,
  WRITE_HISTORY   =  Ioss::WRITE_HISTORY,
  WRITE_HEARTBEAT =  Ioss::WRITE_HEARTBEAT
};

}//namespace Ioi
#endif // Ioi_Events_h
