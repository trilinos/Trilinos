/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef SIERRA_Ioss_State_h
#define SIERRA_Ioss_State_h
namespace Ioss {
  enum State {
    STATE_INVALID = -1,
    STATE_UNKNOWN,
    STATE_READONLY,
    STATE_CLOSED,
    STATE_DEFINE_MODEL,
    STATE_MODEL,
    STATE_DEFINE_TRANSIENT,
    STATE_TRANSIENT,
    STATE_LAST_ENTRY
  };
}
#endif
