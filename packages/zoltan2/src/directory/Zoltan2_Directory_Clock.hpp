/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan2 Directory for Load-balancing, Partitioning, Ordering and Coloring
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

#ifndef ZOLTAN2_DIRECTORY_CLOCK_H_
#define ZOLTAN2_DIRECTORY_CLOCK_H_

// a temporary class to get some timing information out - to be deleted
#include <time.h>
class Zoltan2_Directory_Clock {
  public:
  Zoltan2_Directory_Clock(std::string name_,
    Teuchos::RCP<const Teuchos::Comm<int> > comm_) : comm(comm_), name(name_),
    bCompleted(false) {
    comm->barrier();
    startTime = getTime();
  }
  ~Zoltan2_Directory_Clock() {
    complete();
  }
  void complete() {
    /*
    if(!bCompleted) {
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        comm->barrier();
        if(proc == comm->getRank()) {
          if(proc == 0) {
            std::cout << "CLOCK: " << name << " ";
          }
          double deltaTime = getTime() - startTime;
          std::cout << " " << comm->getRank() << ": " << deltaTime;
          if(proc == comm->getSize()-1) {
            std::cout << std::endl;
          }
        }
      }
      comm->barrier();
    }
    */
    bCompleted = true;
  }
  private:
    double getTime() { return static_cast<double>(clock()) /
      static_cast<double>(CLOCKS_PER_SEC); }
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    std::string name;
    double startTime;
    bool bCompleted;
};

#endif
