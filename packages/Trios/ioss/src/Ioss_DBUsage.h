/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_DBUsage_h
#define SIERRA_Ioss_DBUsage_h

namespace Ioss {

enum DatabaseUsage {
  WRITE_RESTART   =  1,
  READ_RESTART    =  2,
  WRITE_RESULTS   =  4,
  READ_MODEL      =  8,
  WRITE_HISTORY   = 16,
  WRITE_HEARTBEAT = 32
};

inline bool is_input_event(Ioss::DatabaseUsage db_usage) {
  return db_usage == Ioss::READ_MODEL || db_usage == Ioss::READ_RESTART;
}

}
#endif // SIERRA_Ioss_DBUsage_h
