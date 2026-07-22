// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// ************************************************************************
// 
// C++ wrappers for Zoltan Timer library.  
//
// ************************************************************************

#ifndef ZOLTAN_TIMER_CPP_H_
#define ZOLTAN_TIMER_CPP_H_

#include "zoltan_timer.h"
#include <stdio.h>

#ifdef TFLOP
#include <string.h>
#else
#include <string>
#endif

class Zoltan_Timer_Object {

public:

  // Constructor
  Zoltan_Timer_Object(int flag = ZOLTAN_TIME_WALL) {
    // Assumption: MPI has been initialized prior to this call.
    ZTStruct = Zoltan_Timer_Create(flag);
  }

  // Copy constructor
  Zoltan_Timer_Object(const Zoltan_Timer_Object &zt)
  {
  this->ZTStruct = Zoltan_Timer_Copy(zt.ZTStruct);
  }

  // Copy operator
  Zoltan_Timer_Object & operator= (const Zoltan_Timer_Object &zt)
  {
    Zoltan_Timer_Copy_To(&(this->ZTStruct), zt.ZTStruct);
    return *this;
  }

  // Destructor
  ~Zoltan_Timer_Object() {
    Zoltan_Timer_Destroy(&ZTStruct);
  }

  int Init(const int &use_barrier, const std::string & name) {
    return Zoltan_Timer_Init(this->ZTStruct, use_barrier, name.c_str());
  }

  int Reset(const int &ts_idx, const int &use_barrier, 
            const std::string & name) {
    return Zoltan_Timer_Reset(this->ZTStruct, ts_idx, use_barrier, name.c_str());
  }

  int Start(const int &ts_idx, const MPI_Comm &comm) {
    return Zoltan_Timer_Start(this->ZTStruct, ts_idx, comm,
                             (char *) __FILE__, __LINE__);
  }

  int Stop(const int &ts_idx, const MPI_Comm &comm) {
    return Zoltan_Timer_Stop(this->ZTStruct, ts_idx, comm, 
                             (char *) __FILE__, __LINE__);
  }

  int Print(const int &ts_idx, const int &proc, 
            const MPI_Comm &comm, FILE *os) const {
    return Zoltan_Timer_Print(this->ZTStruct, ts_idx, proc, comm, os);
  }

  int PrintAll(const int &proc, const MPI_Comm &comm, FILE *os) const {
    return Zoltan_Timer_PrintAll(this->ZTStruct, proc, comm, os);
  }

private:

  struct Zoltan_Timer *ZTStruct;
};

#endif
