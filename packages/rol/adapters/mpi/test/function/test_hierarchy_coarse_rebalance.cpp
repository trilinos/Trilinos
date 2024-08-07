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
#include "ROL_PinTConstraint.hpp"
#include "ROL_PinTHierarchy.hpp"

typedef double Real;

#define ROL_TEST_ASSERT(a,print) { \
    if(not a) if(print) {\
      std::stringstream ss;\
      ss << "\nFAILURE: ******************************************************" << std::endl;\
      ss <<   "FAILURE: * Test \"" << #a << "\" failed: line " << __LINE__ << std::endl;\
      ss <<   "FAILURE: ******************************************************" << std::endl;\
      throw std::logic_error(ss.str());\
    }\
  }

#define ROL_TEST_EQUALITY(a,b,print) { \
    if(not (a==b)) {\
      std::stringstream ss;\
      ss << "\nFAILURE: ******************************************************" << std::endl;\
      ss <<   "FAILURE: * Test (" << #a << " != " << #b << ") failed: line " << __LINE__ << std::endl;\
      ss <<   "FAILURE: *      (" <<  a << " != " <<  b << ")" << std::endl;\
      ss <<   "FAILURE: ******************************************************" << std::endl;\
      throw std::logic_error(ss.str());\
    }\
    else if(print) { \
      std::stringstream ss;\
      ss << "\nPASS: ******************************************************" << std::endl;\
      ss <<   "PASS: * Test (" << #a << " == " << #b << ") passed: line " << __LINE__ << std::endl;\
      ss <<   "PASS: *      (" <<  a << " == " <<  b << ")" << std::endl;\
      ss <<   "PASS: ******************************************************" << std::endl;\
      *outStream << ss.str() << std::endl;\
    } \
  }

#define ROL_TEST_FLOAT_EQUALITY(a,b,tol,print) { \
    if(not (std::fabs(a-b)<=tol)) {\
      std::stringstream ss;\
      ss << "\nFAILURE: ******************************************************" << std::endl;\
      ss <<   "FAILURE: * Test (" << #a << " != " << #b << ") failed: line " << __LINE__ << std::endl;\
      ss <<   "FAILURE: *      (" << "abs(" << a << " - " <<  b << ") > " << tol << std::endl;\
      ss <<   "FAILURE: ******************************************************" << std::endl;\
      throw std::logic_error(ss.str());\
    }\
    else if(print) { \
      std::stringstream ss;\
      ss << "\nPASS: ******************************************************" << std::endl;\
      ss <<   "PASS: * Test (" << #a << " == " << #b << ") passed: line " << __LINE__ << std::endl;\
      ss <<   "PASS: *      (" << "abs(" << a << " - " <<  b << ") <= " << tol << std::endl;\
      ss <<   "PASS: ******************************************************" << std::endl;\
      *outStream << ss.str() << std::endl;\
    } \
  }

void testHiearchyCoarseTimeStamps(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);

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

  int errorFlag  = 0;

  try {
    testHiearchyCoarseTimeStamps(MPI_COMM_WORLD, outStream);
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

std::string printTimeStamps(const std::string & prefix, const std::vector<ROL::TimeStamp<Real>> & stamps)
{
  std::stringstream ss;


  ss << std::setprecision(3);
  ss << std::fixed;
  ss << prefix;
  ss << "Stamps:  ";
  for(size_t i=0;i<stamps.size();i++) {
    ss << std::setw(6) << stamps[i].t[0] << " ";
  }
  ss << std::endl;
  ss << prefix;
  ss << "         ";
  for(size_t i=0;i<stamps.size();i++) {
    ss << std::setw(6) << stamps[i].t[1] << " ";
  }
  ss << std::endl;

  ss << prefix << "-------------------------------------------------------------"
                  "-------------------------------------------------------------";

  return ss.str();
}

void testHiearchyCoarseTimeStamps(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  std::stringstream ss;

  int spatialProcs = 1;
  int spaceDim = 2;
  int numSteps = 16+1; // power of two plus one
  double dt = 0.1;

  // build communicators
  ROL::Ptr<ROL::PinTCommunicators> pintComm = ROL::makePtr<ROL::PinTCommunicators>(comm,spatialProcs);
  ROL::Ptr<const ROL::PinTVectorCommunication<Real>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<Real>>();

  std::string rank;
  {
    std::stringstream rank_ss;
    rank_ss << "P" << pintComm->getTimeRank() << ". ";
    rank = rank_ss.str();
  }

  // build local vector
  std::vector<Real> data(spaceDim,0.0);
  ROL::Ptr<ROL::Vector<Real>> localVector = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtrFromRef(data));

  // build fine pint vector
  ROL::Ptr<ROL::Vector<Real>> simVec_0 = ROL::buildStatePinTVector(pintComm,vectorComm,numSteps,localVector);
  ROL::PinTVector<Real> & pintSimVec_0 = dynamic_cast<ROL::PinTVector<Real>&>(*simVec_0);

  // Build the hiearchy
  //////////////////////////////////////////////////////////////////////////////////////
 
  // do a scan to get the starting time
  int localNt = pintSimVec_0.numOwnedSteps();
  double myFinalTime = dt*localNt;
  double timeOffset  = 0.0;
  MPI_Exscan(&myFinalTime,&timeOffset,1,MPI_DOUBLE,MPI_SUM,pintComm->getTimeCommunicator());

  *outStream << rank << "Setting the time steps: " << timeOffset << " " << localNt << "/" << numSteps << std::endl; 
 
  // build time stamps
  auto timeStamps_ptr = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>(localNt);
  auto & timeStamps = *timeStamps_ptr;
  for(int k=0; k<localNt; ++k ) {
    timeStamps[k].t.resize(2);
    timeStamps[k].t.at(0) = k*dt + timeOffset;
    timeStamps[k].t.at(1) = (k+1)*dt + timeOffset;
  }

  bool rebalance = true;
  ROL::PinTHierarchy<Real> hierarchy(timeStamps_ptr);
  hierarchy.setMaxLevels(3);
  hierarchy.buildLevels(pintComm,vectorComm,rebalance);

  // Check the time stamps by level
  //////////////////////////////////////////////////////////////////////////////////////
  
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_0 = hierarchy.getTimeStampsByLevel(0);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_1 = hierarchy.getTimeStampsByLevel(1);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_2 = hierarchy.getTimeStampsByLevel(2);

  if(stamps_0!=ROL::nullPtr) *outStream << printTimeStamps(rank + "STAMP 0: ",*stamps_0) << std::endl;;
  if(stamps_1!=ROL::nullPtr) *outStream << printTimeStamps(rank + "STAMP 1: ",*stamps_1) << std::endl;;
  if(stamps_2!=ROL::nullPtr) *outStream << printTimeStamps(rank + "STAMP 2: ",*stamps_2) << std::endl;;

  // send everyone's end points to the root: we are assuming 4 processors
  double endpoint = stamps_0->at(localNt-1).t[1];
  MPI_Send(&endpoint,1,MPI_DOUBLE,0,pintComm->getTimeRank(),pintComm->getTimeCommunicator());

  // error check on the root processors
  if(pintComm->getTimeRank()==0) {
 
    std::vector<double> endpoints;
    for(int p=0;p<pintComm->getTimeSize();p++) {
      double buffer = -1234.5; 
      MPI_Status status;
      MPI_Recv(&buffer,1,MPI_DOUBLE,p,p,pintComm->getTimeCommunicator(),&status);
      
      // *outStream << "ENDPOINT from proccessor " << p << " is " << buffer << std::endl;
      endpoints.push_back(buffer); 
    }

    Real tol = 1e-14;

    ROL_TEST_EQUALITY(static_cast<long long>(stamps_0->size()),localNt,false);
    ROL_TEST_EQUALITY(static_cast<long long>(stamps_1->size()),localNt,false);
    ROL_TEST_EQUALITY(static_cast<long long>(stamps_2->size()),localNt,false);

    ROL_TEST_FLOAT_EQUALITY(stamps_0->at(0).t[0],stamps_1->at(0).t[0],tol,false);
    ROL_TEST_FLOAT_EQUALITY(stamps_0->at(0).t[0],stamps_2->at(0).t[0],tol,false);
  
    // protect the first time step
    ROL_TEST_FLOAT_EQUALITY(stamps_0->at(0).t[1],stamps_1->at(0).t[1],tol,false);
    ROL_TEST_FLOAT_EQUALITY(stamps_0->at(0).t[1],stamps_2->at(0).t[1],tol,false);

    for(int i=0;i<localNt-1;i++) {
      // make sure time stamps are tied together
      ROL_TEST_FLOAT_EQUALITY(stamps_0->at(i).t[1],stamps_0->at(i+1).t[0],tol,false);
      ROL_TEST_FLOAT_EQUALITY(stamps_1->at(i).t[1],stamps_1->at(i+1).t[0],tol,false);
      ROL_TEST_FLOAT_EQUALITY(stamps_2->at(i).t[1],stamps_2->at(i+1).t[0],tol,false);
    }

    ROL_TEST_EQUALITY(stamps_0->at(localNt-1).t[1],endpoints[0],false);
    ROL_TEST_EQUALITY(stamps_1->at(localNt-1).t[1],endpoints[1],false);
    ROL_TEST_EQUALITY(stamps_2->at(localNt-1).t[1],endpoints[3],false);
  }
}
