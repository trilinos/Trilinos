// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <mpi.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_PinTCommunicationUtilities.hpp"
#include "ROL_PinTVectorCommunication_StdVector.hpp"
#include "ROL_PinTVector.hpp"

typedef double Real;

void run_TimeStamp_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int numSteps);
void run_VectorExportToCoarse_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int numSteps);
void run_VectorExportToFine_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int numSteps);

int main(int argc, char* argv[]) 
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int myRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // split the comm so if three processors are used you can run a 3, 2 and 1 processor
  // case all at once
 
  int numSteps = 4;
  int errorFlag  = 0;

  try {

    (*outStream) << "Running TimeStamp test" << std::endl;
    run_TimeStamp_test(MPI_COMM_WORLD, outStream,numSteps);

    // (*outStream) << "Running VectorExportToCoarse test" << std::endl;
    // run_VectorExportToCoarse_test(MPI_COMM_WORLD, outStream,numSteps);

    // (*outStream) << "Running VectorExportToFine test" << std::endl;
    // run_VectorExportToFine_test(MPI_COMM_WORLD, outStream,numSteps);
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  int errors = std::abs(errorFlag);
  MPI_Allreduce(&errors,&errorFlag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if (errorFlag != 0 && myRank==0)
    std::cout << "End Result: TEST FAILED\n";
  else if(myRank==0)
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}

void printSerializedTimeStamps(const std::vector<ROL::TimeStamp<Real>> & timeStamps,
                               int myRank,
                               std::ostream & outStream)
{
  std::vector<double> serialized;
  ROL::PinT::serializeTimeStamps(timeStamps,serialized,0);
  
  std::stringstream ss;
  ss << "Processor " << myRank << ": ";
  for(auto t : serialized) 
    ss << t << " ";
  outStream << ss.str() << std::endl;
}

void run_TimeStamp_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int local_Nt)
{
  int numRanks = -1;
  int myRank = -1;

  std::stringstream ss;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  double totaltime = 1.75; // we will do 1.75 periods, this way the final result is non-trivial
  auto  Nt  = numRanks*local_Nt; // use a power of two
  double dt = totaltime/Nt;
  double timeOffset  = 0.0;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  // compute the time offset for each processor
  {
    double myFinalTime = dt*local_Nt;
    MPI_Exscan(&myFinalTime,&timeOffset,1,MPI_DOUBLE,MPI_SUM,communicators->getTimeCommunicator());
  }

  std::vector<ROL::TimeStamp<Real>> timeStamps(local_Nt);
  for( size_t k=0; k< static_cast<size_t>(local_Nt); ++k ) {
    timeStamps.at(k).t.resize(2);
    timeStamps.at(k).t.at(0) = k*dt+timeOffset;
    timeStamps.at(k).t.at(1) = (k+1)*dt+timeOffset;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printSerializedTimeStamps(timeStamps,myRank,*outStream);

  // build the coarse communicators
  auto coarseComms = communicators->buildCoarseCommunicators();

  // check the coarse communicators are allocated correctly
  if(myRank<2) {
    if(ROL::is_nullPtr(coarseComms)) {
      ss << "Coarse communicator is unexpectely null on rank " << myRank << std::endl;
      throw std::logic_error(ss.str());
    }  
  }
  else {
    if(not ROL::is_nullPtr(coarseComms)) {
      ss << "Coarse communicator is unexpectely NOT null on rank " << myRank << std::endl;
      throw std::logic_error(ss.str());
    }  
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<ROL::TimeStamp<Real>> coarseStamps;
  bool coarseStampsExist = ROL::PinT::exportToCoarseDistribution_TimeStamps(timeStamps,coarseStamps,*communicators,0);

  // get the coarse communicators
  if(not ROL::is_nullPtr(coarseComms)) {
    if(not coarseStampsExist) {
      ss << "Incorrect distribution of coarse stamps on rank " << myRank << std::endl;
      throw std::logic_error(ss.str());
    }

    if(static_cast<long long>(coarseStamps.size())!=local_Nt*2) {
      ss << "Incorrect number of coarse stamps on rank " << myRank << std::endl;
      throw std::logic_error(ss.str());
    }

    printSerializedTimeStamps(coarseStamps,coarseComms->getTimeRank(),*outStream);
  }
  else {
    if(coarseStampsExist) {
      ss << "Incorrect distribution of coarse stamps on rank " << myRank << std::endl;
      throw std::logic_error(ss.str());
    }
  }
}

void run_VectorExportToCoarse_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int local_Nt)
{
  typedef ROL::Ptr<ROL::Vector<Real>> PtrVector;
  typedef ROL::Ptr<ROL::PinTVector<Real>> PtrPinTVector;
  
  int numRanks = -1;
  int myRank = -1;

  std::stringstream ss;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  int spatialProcs = 1;
  ROL::Ptr<ROL::PinTCommunicators> fineComm = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spatialProcs);
  ROL::Ptr<ROL::PinTCommunicators> crseComm = fineComm->buildCoarseCommunicators();
  ROL::Ptr<const ROL::PinTVectorCommunication<Real>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<Real>>();

  // allocate fine vector
  PtrPinTVector fine_vec;
  int spaceSize = 20;
  {
    std::vector<Real> data(spaceSize,0.0);
    PtrVector vec = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtrFromRef(data));

    std::vector<int> stencil = {-1,0};
    fine_vec = ROL::makePtr<ROL::PinTVector<Real>>(fineComm,vectorComm,vec,local_Nt*numRanks,stencil);
  }

  // allocate coarse vector (if relevant)
  PtrPinTVector crse_vec;
  if(not ROL::is_nullPtr(crseComm)) {
    std::vector<Real> data(spaceSize,0.0);
    PtrVector vec = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtrFromRef(data));

    std::vector<int> stencil = {-1,0};
    crse_vec = ROL::makePtr<ROL::PinTVector<Real>>(crseComm,vectorComm,vec,local_Nt*numRanks,stencil);
  }

  // allocate and fill a fine vector
  for(int i=-1;i<fine_vec->numOwnedSteps();i++) {
    auto v = fine_vec->getVectorPtr(i);
    v->setScalar(Real(myRank)+(i+1.0)/(spaceSize+1.0)); // EXPECTED VALUE
  }

  MPI_Barrier(fineComm->getTimeCommunicator());

/*
  for(int i=0;i<4;i++) {
    if(myRank==i) {
      std::cout << "FINE PROCESS " << myRank << std::endl;
      std::cout << "===============================" << std::endl;
      fine_vec->print(std::cout);
      std::cout << "===============================" << std::endl;
      std::cout.flush();
    }

    MPI_Barrier(fineComm->getTimeCommunicator());
  }

  MPI_Barrier(fineComm->getTimeCommunicator());
*/

  if(myRank==0)
    std::cout << std::endl << std::endl;

  // build a STL vector with the fine vectors
  int startIndex = (myRank % 2==0) ? -1 : 0;
  std::vector<ROL::Ptr<ROL::Vector<Real>>> fineVectors(fine_vec->numOwnedSteps()-startIndex);
  for(int i=startIndex;i<fine_vec->numOwnedSteps();i++)
    fineVectors[i-startIndex] = fine_vec->getVectorPtr(i); 
  
  // send data
  ROL::PinT::sendToCoarseDistribution_Vector(fineVectors,*vectorComm,*fineComm);

  // recv data
  if(not ROL::is_nullPtr(crseComm)) {
    // build a STL vector with the fine vectors
    int fineSize = fine_vec->numOwnedSteps();
    std::vector<ROL::Ptr<ROL::Vector<Real>>> crseVectors(2*fineSize+1);
    for(int i=-1;i<crse_vec->numOwnedSteps();i++)
      crseVectors[i+1] = crse_vec->getVectorPtr(i); 
   
    ROL::PinT::recvFromFineDistribution_Vector(crseVectors,*vectorComm,*fineComm,fineSize+1,fineSize);

    /*
    for(int i=0;i<2;i++) {
      if(myRank==i) {
        std::cout << "COARSE PROCESS " << myRank << std::endl;
        std::cout << "===============================" << std::endl;
        crse_vec->print(std::cout);
        std::cout << "===============================" << std::endl;
        std::cout.flush();
      }

      MPI_Barrier(crseComm->getTimeCommunicator());
    }
    MPI_Barrier(crseComm->getTimeCommunicator());
    */

    int numFineSteps = fine_vec->numOwnedSteps();
    for(int i=-1;i<crse_vec->numOwnedSteps();i++) {
      Real norm = crse_vec->getVectorPtr(i)->norm();

      // the expected value is a function of the rank that sent it
      // and the spatial index of the vector, see line labeled "EXPECTED VALUE"
      Real expectedValue = 0.0;
      if(i<numFineSteps) 
        // the first numFineSteps+1 (not iteration start) expect this
        expectedValue = (Real(2.0*myRank)+(i+1.0)/(spaceSize+1.0));
      else
        // the last numFineSteps expect this
        expectedValue = (Real(2.0*myRank+1)+(i-numFineSteps+1.0)/(spaceSize+1.0));

      Real expectedNorm  = std::sqrt(spaceSize*expectedValue*expectedValue);

      if(std::fabs(norm-expectedNorm)/expectedNorm > 1e-15) {
        ss << "Expected norm (" << expectedNorm << ") does not match computed norm (" << norm 
           << ") on processor " << myRank << " vector " << i;

        throw std::logic_error(ss.str());
      }
    }
  }
}

void run_VectorExportToFine_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int local_Nt)
{
  typedef ROL::Ptr<ROL::Vector<Real>> PtrVector;
  typedef ROL::Ptr<ROL::PinTVector<Real>> PtrPinTVector;
  
  int numRanks = -1;
  int myRank = -1;

  std::stringstream ss;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  int spatialProcs = 1;
  ROL::Ptr<ROL::PinTCommunicators> fineComm = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spatialProcs);
  ROL::Ptr<ROL::PinTCommunicators> crseComm = fineComm->buildCoarseCommunicators();
  ROL::Ptr<const ROL::PinTVectorCommunication<Real>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<Real>>();

  // allocate fine vector
  PtrPinTVector fine_vec;
  int spaceSize = 20;
  {
    std::vector<Real> data(spaceSize,0.0);
    PtrVector vec = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtrFromRef(data));

    std::vector<int> stencil = {-1,0};
    fine_vec = ROL::makePtr<ROL::PinTVector<Real>>(fineComm,vectorComm,vec,local_Nt*numRanks,stencil);
  }

  // allocate coarse vector
  PtrPinTVector crse_vec;
  if(not ROL::is_nullPtr(crseComm)) {
    std::vector<Real> data(spaceSize,0.0);
    PtrVector vec = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtrFromRef(data));

    std::vector<int> stencil = {-1,0};
    crse_vec = ROL::makePtr<ROL::PinTVector<Real>>(crseComm,vectorComm,vec,2*local_Nt*numRanks/2,stencil);

    // allocate and fill a fine vector
    for(int i=-1;i<crse_vec->numOwnedSteps();i++) {
      auto v = crse_vec->getVectorPtr(i);
      v->setScalar(Real(myRank)+(i+1.0)/(crse_vec->numOwnedSteps()+1.0)); // EXPECTED VALUE
    }
  }

  MPI_Barrier(fineComm->getTimeCommunicator());

  // recv data
  if(not ROL::is_nullPtr(crseComm)) {

    // build a STL vector with the fine vectors
    std::vector<ROL::Ptr<ROL::Vector<Real>>> crseVectors_a(local_Nt+1,ROL::nullPtr);
    std::vector<ROL::Ptr<ROL::Vector<Real>>> crseVectors_b(local_Nt+1,ROL::nullPtr);
    for(int i=-1;i<local_Nt;i++)
      crseVectors_a[i+1] = crse_vec->getVectorPtr(i); 

    for(int i=local_Nt-1;i<crse_vec->numOwnedSteps();i++)
      crseVectors_b[i-local_Nt+1] = crse_vec->getVectorPtr(i); 

    // send data
    ROL::PinT::sendToFineDistribution_Vector(crseVectors_a,crseVectors_b,*vectorComm,*fineComm);
  }

  // build a STL vector with the fine vectors
  std::vector<ROL::Ptr<ROL::Vector<Real>>> fineVectors(fine_vec->numOwnedVectors());
  for(int i=-1;i<fine_vec->numOwnedSteps();i++)
    fineVectors[i+1] = fine_vec->getVectorPtr(i); 

  ROL::PinT::recvFromCoarseDistribution_Vector(fineVectors,*vectorComm,*fineComm);

  for(int i=-1;i<fine_vec->numOwnedSteps();i++) {
    Real norm = fine_vec->getVectorPtr(i)->norm();

    int coarseRank = (myRank % 2 ==0) ? myRank/2 : (myRank-1)/2;
    int ts_i = (myRank % 2 ==0) ? i : i+local_Nt;

    // the expected value is a function of the rank that sent it
    // and the spatial index of the vector, see line labeled "EXPECTED VALUE"
    Real expectedValue = 0.0;

    expectedValue = Real(coarseRank)+(ts_i+1.0)/(8+1.0); // EXPECTED VALUE

    Real expectedNorm  = std::sqrt(spaceSize*expectedValue*expectedValue);

    if(std::fabs(norm-expectedNorm)/expectedNorm > 1e-15) {
      ss << "Expected norm (" << expectedNorm << ") does not match computed norm (" << norm 
         << ") on processor " << myRank << " vector " << i;

      throw std::logic_error(ss.str());
    }
  }
}
