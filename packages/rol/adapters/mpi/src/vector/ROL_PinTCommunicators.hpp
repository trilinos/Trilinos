// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTCOMMUNICATORS_H
#define ROL_PINTCOMMUNICATORS_H

#include <ostream>

#include "mpi.h"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

/** 
 * This class handles PinT communication. Its
 * based on the approach taken in X-Braid.
 */
class PinTCommunicators {
public:
  PinTCommunicators(MPI_Comm parent,int spatialProcs, bool freeParent=false,const std::string & label="Default")
    : freeParent_(freeParent), label_(label)
  { 
    parentComm_ = parent; 

    int myGlobalRank = -1;
    int globalProcs = -1;
    MPI_Comm_size(parentComm_, &globalProcs);
    MPI_Comm_rank(parentComm_, &myGlobalRank);

    // make sure they divide evenly (for sanity!)
    assert(globalProcs % spatialProcs == 0);

    int space = myGlobalRank / spatialProcs;
    int time  = myGlobalRank % spatialProcs;

    // this decomposition comes from X-Braid 
    MPI_Comm_split(parentComm_,space,myGlobalRank,&spaceComm_); 
    MPI_Comm_split(parentComm_, time,myGlobalRank, &timeComm_); 

    MPI_Comm_size(timeComm_, &timeSize_);
    MPI_Comm_rank(timeComm_, &timeRank_);

    MPI_Comm_size(spaceComm_, &spaceSize_);
    MPI_Comm_rank(spaceComm_, &spaceRank_);
  }

  /**
   * Deleted the copy constructor because the destructor is non-trivial. This could
   * be fixed with a more complex copy constructor implementation and more logic in
   * the destructor. However, following the pointer and reference semantics
   * seems more straightforward.
   */
  PinTCommunicators(const PinTCommunicators &) = delete;

  // cleanup
  ~PinTCommunicators() {
    MPI_Comm_free(&spaceComm_);  
    MPI_Comm_free(&timeComm_); 
    if(freeParent_)
      MPI_Comm_free(&parentComm_); 
   }

   MPI_Comm getParentCommunicator() const { return parentComm_; }
   MPI_Comm getSpaceCommunicator() const { return spaceComm_; }
   MPI_Comm getTimeCommunicator() const { return timeComm_; }

   int getTimeRank() const { return timeRank_; }
   int getTimeSize() const { return timeSize_; }

   int getSpaceRank() const { return spaceRank_; }
   int getSpaceSize() const { return spaceSize_; }

   // support for coarsening
   ROL::Ptr<PinTCommunicators> buildCoarseCommunicators() const
   {
      // this would be easy to generalize by just adjusting the subdivide
      // variable

      int subdivide = 2;

      int myGlobalRank = -1;
      int globalProcs = -1;
      MPI_Comm_size(parentComm_, &globalProcs);
      MPI_Comm_rank(parentComm_, &myGlobalRank);

      if(globalProcs  % 2 != 0) 
        throw std::logic_error("Building coarse communicators requires a power of two.");

      int halfCount = globalProcs / subdivide;
      
      // split the communicator in half, using only the lower half
      MPI_Comm halfComm;
      if(myGlobalRank<halfCount) {
        MPI_Comm_split(parentComm_,0,myGlobalRank,&halfComm); 
      }
      else {
        MPI_Comm_split(parentComm_,MPI_UNDEFINED,myGlobalRank,&halfComm); 
        assert(MPI_COMM_NULL==halfComm);

        return ROL::nullPtr;
      }

      // build a new pint communicator that works on half the parent processors
      return ROL::makePtr<PinTCommunicators>(halfComm,getSpaceSize(),true,label_+"->Coarse");
   }

private:

  bool freeParent_;    // for when we coarsen
  MPI_Comm parentComm_;

  MPI_Comm spaceComm_;
  MPI_Comm timeComm_;
  int timeSize_;
  int timeRank_;
  int spaceSize_;
  int spaceRank_;

  std::string label_;
};

} // namespace ROL

#endif
