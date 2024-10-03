// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTCOMMUNICATIONUTILITIES_HPP
#define ROL_PINTCOMMUNICATIONUTILITIES_HPP

#include <vector>
#include <mpi.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_TimeStamp.hpp"
#include "ROL_PinTCommunicators.hpp"
#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_PinTVector.hpp"

namespace ROL {
namespace PinT { // a utilities namespace for PinT

//! Used for serialization of a vector of time stamps so it can be moved using MPI
template <typename RealT>
void serializeTimeStamps(const std::vector<ROL::TimeStamp<RealT>> & stamps,std::vector<double> & serialized,int startIndex) 
{
  serialized.clear();
    
  // note that initial step ignores the false step
  for(size_t i=startIndex;i<stamps.size();i++) { 
    const ROL::TimeStamp<RealT> & ts = stamps[i];

    // only thing that works now
    assert(ts.t.size()==2);

    serialized.push_back(ts.t.at(0));
    serialized.push_back(ts.t.at(1));
  }

  // sanity check of post condition
  assert(serialized.size()==(stamps.size()-startIndex)*2); 
}

//! Used for serialization of a vector of time stamps so it can be moved using MPI
template <typename RealT>
void deserializeTimeStamps(const std::vector<double> & serialized,std::vector<ROL::TimeStamp<RealT>> & stamps,int startIndex) 
{
  stamps.clear();

  // sanity check
  assert(serialized.size() % 2 == 0);

  // allocate space
  stamps.resize(serialized.size() / 2 - startIndex);
    
  // note that initial step ignores the false step
  for(size_t i=startIndex,ind=2*startIndex;i<stamps.size();i++,ind+=2) { 
    ROL::TimeStamp<RealT> & ts = stamps[i];

    ts.t.resize(2);
    ts.t[0] = serialized[ind];
    ts.t[1] = serialized[ind+1];
  }

  // sanity check of post condition
  assert(serialized.size()==(stamps.size()-startIndex)*2); 
}

/**
 * \brief This method distributes time stamps towards a coarse subset of processors.
 *
 * Currently the implementation sends all the time stamps (equally distributed in sequence) to 
 * half of the coarse processors.
 *
 * \param[in] fineDistrStamps The time stamps distributed on the fine level processors
 * \param[out] coarseDistrStamps The time stamps distributed on the coarse level processors
 * \param[in] communicators PinT communicators used on the fine level
 * \param[in] startIndex Index to start the fine level distributed time stamps.
 *
 * \note I use the phrase export because its going from fine to coarse and is called from the
 *       fine level. That is the communicators that are being used are the fine PinTCommunicators.
 */ 
template <typename RealT>
bool exportToCoarseDistribution_TimeStamps(const std::vector<ROL::TimeStamp<RealT>> & fineDistrStamps,
                                           std::vector<ROL::TimeStamp<RealT>> & coarseDistrStamps,
                                           const ROL::PinTCommunicators & communicators,
                                           int startIndex)
{
  MPI_Comm comm = communicators.getTimeCommunicator();
  int myRank    = communicators.getTimeRank();
  int procSize  = communicators.getTimeSize();

  std::vector<RealT> serialized;
  serializeTimeStamps(fineDistrStamps,serialized,startIndex);

  // send to your lower processor
  if(myRank % 2 ==0) {
    MPI_Send(const_cast<RealT*>(&serialized[0]),int(serialized.size()),MPI_DOUBLE,myRank/2,0,comm);
  }
  else {
    MPI_Send(const_cast<RealT*>(&serialized[0]),int(serialized.size()),MPI_DOUBLE,(myRank-1)/2,0,comm);
  }

  if(myRank < procSize / 2) {
    std::vector<double> newSerialized;

    int count = 0;
    MPI_Status status;
    std::vector<double> buffer(2*serialized.size());

    // recieve first block
    MPI_Recv(&buffer[0],int(buffer.size()),MPI_DOUBLE,2*myRank,0,comm,&status);
    MPI_Get_count(&status,MPI_DOUBLE,&count);
    newSerialized.insert(newSerialized.end(),buffer.begin(),buffer.begin()+count); 

    // recieve second block
    MPI_Recv(&buffer[0],int(buffer.size()),MPI_DOUBLE,2*myRank+1,0,comm,&status);
    MPI_Get_count(&status,MPI_DOUBLE,&count);
    newSerialized.insert(newSerialized.end(),buffer.begin(),buffer.begin()+count); 

    deserializeTimeStamps(newSerialized,coarseDistrStamps,0);

    return true;
  }

  return false;
}

// next two functions communicate from the fine processor mesh to the coarse processor mesh

/**
 * \brief Send a block of vectors to the corresponding processor in the coarse processor grid.
 *
 * Currently the implmentation will send two processors of vectors in the fine grid to one processor
 * of vectors in the coarse grid
 *
 * \param[in] fineVectors Fine vectors to send to the coarse grid.
 * \param[in] vectorComm Object that encapsulates vector communication
 * \param[in] communicators Object that seperates spatial parallelism from temporal parallelism
 */
template <typename RealT>
void sendToCoarseDistribution_Vector(const std::vector<ROL::Ptr<ROL::Vector<RealT>>> & fineVectors,
                                     const ROL::PinTVectorCommunication<RealT> & vectorComm,
                                     const ROL::PinTCommunicators & communicators)
{
  MPI_Comm comm = communicators.getTimeCommunicator();
  int myRank    = communicators.getTimeRank();

  int targetRank = -1;
  if(myRank % 2 == 0) {
    targetRank = myRank/2;
  } else {
    targetRank = (myRank-1)/2;
  }

  // send to your lower processor
  for(size_t i=0;i<fineVectors.size();i++) {
    auto & vec = *fineVectors[i];
    // std::cout << "    sending " << myRank << "/" << communicators.getTimeSize() << " " << fineVectors[i] << std::endl;
    vectorComm.send(comm,targetRank,vec,i);
  }
}

/**
 * \brief Recieve a block of vectors from the corresponding processor in the fine processor grid.
 *
 * Currently the implmentation will recieve from two processors vectors in the fine grid to one processor
 * of vectors in the coarse grid
 *
 * \param[in] coarseVectors Destination for any coarse vectors
 * \param[in] vectorComm Object that encapsulates vector communication
 * \param[in] communicators Object that seperates spatial parallelism from temporal parallelism
 * \param[in] size_a Number of vectors to recieve from the first processor in the fine grid.
 * \param[in] size_b Number of vectors to recieve from the second processor in the fine grid.
 *
 * \note This asserts the coarseVectors size is "size_a + size_b"
 */
template <typename RealT>
void recvFromFineDistribution_Vector(std::vector<ROL::Ptr<ROL::Vector<RealT>>> & coarseVectors,
                                     const ROL::PinTVectorCommunication<RealT> & vectorComm,
                                     const ROL::PinTCommunicators & communicators,
                                     int size_a,
                                     int size_b)
{
  MPI_Comm comm = communicators.getTimeCommunicator();
  int myRank    = communicators.getTimeRank();
  int procSize  = communicators.getTimeSize();

  assert(int(coarseVectors.size())==size_a+size_b); // this is the expected size

  // only if this is in the right set
  if(myRank < procSize / 2) {
    // recieve from previously sent processors
    for(int i=0;i<std::max(size_a,size_b);i++) {
      if(i<size_a)
        vectorComm.recv(comm,2*myRank,*coarseVectors[i],i);
      if(i<size_b)
        vectorComm.recv(comm,2*myRank+1,*coarseVectors[size_a+i],i);
    }
  }
  else {
    assert(false);

    return;
  }
}

// next two functions communicate from the coarse processor mesh to the fine processor mesh

/**
 * \brief Send a block of vectors to the corresponding processor in the coarse processor grid.
 *
 * Currently the implmentation will send two processors of vectors in the fine grid to one processor
 * of vectors in the coarse grid
 *
 * \param[in] coarseVectors_a Coarse vectors to send to the fine grid processors a.
 * \param[in] coarseVectors_a Coarse vectors to send to the fine grid processors b.
 * \param[in] vectorComm Object that encapsulates vector communication
 * \param[in] communicators Object that seperates spatial parallelism from temporal parallelism
 *
 * \note This asserts the coarseVectors size is "size_a + size_b"
 */
template <typename RealT>
void sendToFineDistribution_Vector(const std::vector<ROL::Ptr<ROL::Vector<RealT>>> & coarseVectors_a,
                                   const std::vector<ROL::Ptr<ROL::Vector<RealT>>> & coarseVectors_b,
                                   const ROL::PinTVectorCommunication<RealT> & vectorComm,
                                   const ROL::PinTCommunicators & communicators)
{
  MPI_Comm comm = communicators.getTimeCommunicator();
  int myRank    = communicators.getTimeRank();

  // send to your lower processor
  for(size_t i=0;i<coarseVectors_a.size();i++) {
    vectorComm.send(comm,2*myRank,*coarseVectors_a[i],i);
  }

  for(size_t i=0;i<coarseVectors_b.size();i++) {
    vectorComm.send(comm,2*myRank+1,*coarseVectors_b[i],i);
  }
}

/**
 * \brief Recieve a block of vectors from the corresponding processor in the fine processor grid.
 *
 * Currently the implmentation will recieve from two processors vectors in the fine grid to one processor
 * of vectors in the coarse grid
 *
 * \param[in] fineVectors Destination for any fine vectors
 * \param[in] vectorComm Object that encapsulates vector communication
 * \param[in] communicators Object that seperates spatial parallelism from temporal parallelism
 */
template <typename RealT>
void recvFromCoarseDistribution_Vector(std::vector<ROL::Ptr<ROL::Vector<RealT>>> & fineVectors,
                                     const ROL::PinTVectorCommunication<RealT> & vectorComm,
                                     const ROL::PinTCommunicators & communicators)
{
  MPI_Comm comm = communicators.getTimeCommunicator();
  int myRank    = communicators.getTimeRank();

  int targetRank = -1;
  if(myRank % 2 == 0) {
    targetRank = myRank/2;
  } else {
    targetRank = (myRank-1)/2;
  }

  // send to your lower processor
  for(int i=0;i<int(fineVectors.size());i++) {
    vectorComm.recv(comm,targetRank,*fineVectors[i],i);
  }
}

} // namespace PinT
} // namespace ROL 

#endif 
