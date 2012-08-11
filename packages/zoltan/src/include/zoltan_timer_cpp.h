/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
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
