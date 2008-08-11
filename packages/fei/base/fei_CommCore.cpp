/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_CommCore.hpp"

namespace fei {

//----------------------------------------------------------------------------
fei::CommCore::CommCore()
  :
#ifndef FEI_SER
    mpiReqs_(),
    mpiStatuses_(),
#endif
    tmpIntData_(),
    tag_(713)
{
}

fei::CommCore::~CommCore()
{
}

fei::CommCore* fei::CommCore::create()
{
  static fei::CommCore fei_comm_core;
  return(&fei_comm_core);
}

int fei::CommCore::get_tag()
{
  return(tag_);
}

void fei::CommCore::increment_tag()
{
  ++tag_;
}

} //namespace fei

