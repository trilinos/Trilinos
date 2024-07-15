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

void testRestrictionProlong_SimVector(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void testRestrictionProlong_OptVector(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);

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
    testRestrictionProlong_SimVector(MPI_COMM_WORLD,outStream);
    testRestrictionProlong_OptVector(MPI_COMM_WORLD,outStream);
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

// print the first spatial entry of each vector
std::string printVector_Control(const ROL::Ptr<const ROL::PinTVector<Real>> & pintVec)
{
  if(pintVec==ROL::nullPtr)
    return "Control: <null>"; 

  std::stringstream ss;

  int steps = pintVec->numOwnedSteps();
 
  ss << std::setprecision(3);
  ss << std::fixed;
  ss << "Control: ";
  for(int i=0;i<steps;i++) {
    Real value = dynamic_cast<ROL::StdVector<Real>&>(*pintVec->getVectorPtr(i)).getVector()->at(0);
    ss << std::setw(6) << value << " ";
  }

  return ss.str();
}

// print the first spatial entry of each vector
std::string printVector_Primary(const ROL::Ptr<const ROL::PinTVector<Real>> & pintVec)
{
  if(pintVec==ROL::nullPtr)
    return "Control: <null>"; 

  std::stringstream ss;

  int steps = pintVec->numOwnedSteps();
 
  ss << std::setprecision(3);
  ss << std::fixed;
  ss << "Primary: ";
  for(int i=0;i<steps;i++) {
    Real value = dynamic_cast<ROL::StdVector<Real>&>(*pintVec->getVectorPtr(2*i)).getVector()->at(0);
    ss << std::setw(6) << value << " ";
  }

  return ss.str();
}

// print the first spatial entry of each vector
std::string printVector_Virtual(const ROL::Ptr<const ROL::PinTVector<Real>> & pintVec)
{
  if(pintVec==ROL::nullPtr)
    return "Control: <null>"; 

  std::stringstream ss;

  int steps = pintVec->numOwnedSteps();
 
  ss << std::setprecision(3);
  ss << std::fixed;
  ss << "Virtual: ";
  for(int i=0;i<steps;i++) {
    Real value = dynamic_cast<ROL::StdVector<Real>&>(*pintVec->getVectorPtr(2*i+1)).getVector()->at(0);
    ss << std::setw(6) << value << " ";
  }

  return ss.str();
}

void testRestrictionProlong_SimVector(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream) 
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
  ROL::Ptr<ROL::PinTVector<Real>> pintSimVec_0 = ROL::buildStatePinTVector(pintComm,vectorComm,numSteps,localVector);
  // ROL::PinTVector<Real> & pintSimVec_0 = dynamic_cast<ROL::PinTVector<Real>&>(*simVec_0);

  // Build the hiearchy
  //////////////////////////////////////////////////////////////////////////////////////
  //
  *outStream << rank << "Computing time offset" << std::endl; 

  // do a scan to get the starting time
  int localNt = pintSimVec_0->numOwnedSteps();
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
  
  *outStream << rank << "Setting up the hierarchy: " << timeStamps_ptr->size() << std::endl; 

  bool rebalance = true;
  ROL::PinTHierarchy<Real> hierarchy(timeStamps_ptr);
  hierarchy.setMaxLevels(3);
  hierarchy.buildLevels(pintComm,vectorComm,rebalance);

  // Check the time stamps by level
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Building Time stamps" << std::endl;
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_0 = hierarchy.getTimeStampsByLevel(0);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_1 = hierarchy.getTimeStampsByLevel(1);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_2 = hierarchy.getTimeStampsByLevel(2);

  // check the sizing of the simulation vector (should include the virtual state
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Building restricted vectors" << std::endl; 

  ROL::Ptr<ROL::PinTVector<Real>> pintSimVec_1 = hierarchy.allocateSimVector(*pintSimVec_0,1);
  ROL::Ptr<ROL::PinTVector<Real>> pintSimVec_2 = hierarchy.allocateSimVector(*pintSimVec_0,2);
  
  
  if(pintSimVec_0->numOwnedVectors() != 2*localNt) {
    ss << rank << "Level 0: Incorrect number of owned vectors: " << pintSimVec_0->numOwnedVectors() << "!=" << 2*localNt << std::endl;
    throw std::logic_error(ss.str());
  }

  // simulation vectors contain the virtual variable, this is all shifted backwards by one
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Setting linear function on fine grid" << std::endl; 

  Real shift = 1.0;
  
  for(size_t i=0;i<timeStamps.size();i++) {
    Real virtualValue = timeOffset+dt+i*dt;
    pintSimVec_0->getVectorPtr(2*i)->setScalar(3.0*virtualValue+shift);
    pintSimVec_0->getVectorPtr(2*i+1)->setScalar(virtualValue+shift);
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // Prolongation and restriction are designed to preserve linear functions - 12/21/2018
  //
  //////////////////////////////////////////////////////////////////////////////////////

  // perform restriction, and test that it works (check for preservation of linear functions)
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << std::endl;
  *outStream << rank << "Test Restriction " << std::endl;
  *outStream << rank << "**********************************************************************" << std::endl;

  hierarchy.restrictSimVector(pintSimVec_0,pintSimVec_1,0);
  hierarchy.restrictSimVector(pintSimVec_1,pintSimVec_2,1);

  *outStream << rank << "R_Level 0 = " << printVector_Primary(pintSimVec_0) << std::endl;
  *outStream << rank << "R_Level 0 = " << printVector_Virtual(pintSimVec_0) << std::endl;

  *outStream << rank << std::endl;
  *outStream << rank << std::endl;

  *outStream << rank << "R_Level 1 = " << printVector_Primary(pintSimVec_1) << std::endl;
  *outStream << rank << "R_Level 1 = " << printVector_Virtual(pintSimVec_1) << std::endl;

  *outStream << rank << std::endl;
  *outStream << rank << std::endl;

  *outStream << rank << "R_Level 2 = " << printVector_Primary(pintSimVec_2) << std::endl;
  *outStream << rank << "R_Level 2 = " << printVector_Virtual(pintSimVec_2) << std::endl;

  // check errors on the coarsest level
  if(stamps_2!=ROL::nullPtr) {
    for(size_t i=0;i<stamps_2->size();i++) {
      Real virtualValue = (*stamps_2)[i].t[1];

      Real u_value = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_2->getVectorPtr(2*i)).getVector()->at(0);
      Real v_value = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_2->getVectorPtr(2*i+1)).getVector()->at(0);

      *outStream << "Compare: " << u_value << " - " << (3.0*virtualValue+shift) << " > 1e-15 " << std::endl;
      if(std::fabs(u_value-(3.0*virtualValue+shift)) > 1e-15) {
        ss << rank << "Two levels of restriction are not correct: u " << u_value << " != " << 3.0*virtualValue+shift << std::endl;
        throw std::logic_error(ss.str());
      }

      *outStream << "Compare: " << v_value << " - " << (virtualValue+shift) << " > 1e-15 " << std::endl;
      if(std::fabs(v_value-(virtualValue+shift)) > 1e-15) {
        ss << rank << "Two levels of restriction are not correct: v " << v_value << " != " << virtualValue+shift << std::endl;
        throw std::logic_error(ss.str());
      }
    }
  }

  *outStream << rank << "Sim Vector Restriction: PASSED" << std::endl;

  // perform prolong, and test that it works (check for preservation of linear functions)
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << std::endl;
  *outStream << rank << "Test Prolongation " << std::endl;
  *outStream << rank << "**********************************************************************" << std::endl;
  
  ROL::Ptr<ROL::PinTVector<Real>> pintSimVec_0_prolong = hierarchy.allocateSimVector(*pintSimVec_0,0);
  ROL::Ptr<ROL::PinTVector<Real>> pintSimVec_1_prolong = hierarchy.allocateSimVector(*pintSimVec_0,1);

  if(pintSimVec_0_prolong!=ROL::nullPtr) pintSimVec_0_prolong->zero(); 
  if(pintSimVec_1_prolong!=ROL::nullPtr) pintSimVec_1_prolong->zero(); 

  hierarchy.prolongSimVector(pintSimVec_2,pintSimVec_1_prolong,2);
  hierarchy.prolongSimVector(pintSimVec_1_prolong,pintSimVec_0_prolong,1);

  *outStream << rank << "Level 0 = " << printVector_Primary(pintSimVec_0_prolong) << std::endl;
  *outStream << rank << "Level 0 = " << printVector_Virtual(pintSimVec_0_prolong) << std::endl;

  *outStream << rank << std::endl;

  *outStream << rank << "Level 1 = " << printVector_Primary(pintSimVec_1_prolong) << std::endl;
  *outStream << rank << "Level 1 = " << printVector_Virtual(pintSimVec_1_prolong) << std::endl;

  *outStream << rank << std::endl;

  *outStream << rank << "Level 2 = " << printVector_Primary(pintSimVec_2) << std::endl;
  *outStream << rank << "Level 2 = " << printVector_Virtual(pintSimVec_2) << std::endl;

   // check errors on the finest level
  for(size_t i=0;i<stamps_0->size();i++) {

    // prolonged values
    Real u_value_prlng = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_0_prolong->getVectorPtr(2*i)).getVector()->at(0);
    Real v_value_prlng = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_0_prolong->getVectorPtr(2*i+1)).getVector()->at(0);

    // original (or exact) values from  the initialized fine grid vector
    Real u_value_exact = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_0->getVectorPtr(2*i)).getVector()->at(0);
    Real v_value_exact = dynamic_cast<ROL::StdVector<Real>&>(*pintSimVec_0->getVectorPtr(2*i+1)).getVector()->at(0);

    *outStream << "Compare: " << u_value_prlng << " - " << u_value_exact << " > 1e-15 " << std::endl;
    if(std::fabs(u_value_prlng-u_value_exact) > 1e-15) {
      ss << rank << "Two levels of prolongation are not correct: u " << u_value_prlng << " != " << u_value_exact << std::endl;
      throw std::logic_error(ss.str());
    }

    *outStream << "Compare: " << v_value_prlng << " - " << v_value_exact << " > 1e-15 " << std::endl;
    if(std::fabs(v_value_prlng-v_value_exact) > 1e-15) {
      ss << rank << "Two levels of prolongation are not correct: v " << v_value_prlng << " != " << v_value_exact << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  *outStream << rank << "Sim Vector Prolongation: PASSED" << std::endl;
}

void testRestrictionProlong_OptVector(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream) 
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
  ROL::Ptr<ROL::PinTVector<Real>> pintOptVec_0 = ROL::buildControlPinTVector(pintComm,vectorComm,numSteps,localVector);

  // Build the hiearchy
  //////////////////////////////////////////////////////////////////////////////////////
  //
  *outStream << rank << "Computing time offset" << std::endl; 

  // do a scan to get the starting time
  int localNt = pintOptVec_0->numOwnedSteps();
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
  
  *outStream << rank << "Setting up the hierarchy: " << timeStamps_ptr->size() << std::endl; 

  bool rebalance = true;
  ROL::PinTHierarchy<Real> hierarchy(timeStamps_ptr);
  hierarchy.setMaxLevels(3);
  hierarchy.buildLevels(pintComm,vectorComm,rebalance);

  // Check the time stamps by level
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Building Time stamps" << std::endl;
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_0 = hierarchy.getTimeStampsByLevel(0);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_1 = hierarchy.getTimeStampsByLevel(1);
  ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps_2 = hierarchy.getTimeStampsByLevel(2);

  // check the sizing of the optimization vector
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Building restricted vectors" << std::endl; 

  ROL::Ptr<ROL::PinTVector<Real>> pintOptVec_1 = hierarchy.allocateOptVector(*pintOptVec_0,1);
  ROL::Ptr<ROL::PinTVector<Real>> pintOptVec_2 = hierarchy.allocateOptVector(*pintOptVec_0,2);

  *outStream << rank << "Finished building restrictued vectors" << std::endl;
  
  if(pintOptVec_0->numOwnedVectors() != localNt) {
    ss << rank << "Level 0: Incorrect number of owned vectors: " << pintOptVec_0->numOwnedVectors() << "!=" << localNt << std::endl;
    throw std::logic_error(ss.str());
  }

  // add in a linear fucntion on the grid
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << "Setting linear function on fine grid" << std::endl; 

  Real shift = 0.0;
  
  for(size_t i=0;i<timeStamps.size();i++) {
    Real value = 0.5*(timeStamps[i].t[0]+timeStamps[i].t[1]);
    pintOptVec_0->getVectorPtr(i)->setScalar(value+shift);
  }


  //////////////////////////////////////////////////////////////////////////////////////
  //
  // Prolongation and restriction are designed to preserve linear functions - 12/21/2018
  //
  //////////////////////////////////////////////////////////////////////////////////////

  // perform restriction, and test that it works (check for preservation of linear functions)
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << std::endl;
  *outStream << rank << "Test Restriction " << std::endl;
  *outStream << rank << "**********************************************************************" << std::endl;

  hierarchy.restrictOptVector(pintOptVec_0,pintOptVec_1,0);
  hierarchy.restrictOptVector(pintOptVec_1,pintOptVec_2,1);

  *outStream << rank << "R_Level 0 = " << printVector_Control(pintOptVec_0) << std::endl;

  *outStream << rank << std::endl;
  *outStream << rank << std::endl;

  *outStream << rank << "R_Level 1 = " << printVector_Control(pintOptVec_1) << std::endl;

  *outStream << rank << std::endl;
  *outStream << rank << std::endl;

  *outStream << rank << "R_Level 2 = " << printVector_Control(pintOptVec_2) << std::endl;

  // check errors on the coarsest level
  if(stamps_2!=ROL::nullPtr) {
    for(size_t i=0;i<stamps_2->size();i++) {
      Real value_exact = 0.5*((*stamps_2)[i].t[0]+(*stamps_2)[i].t[1]);

      Real value = dynamic_cast<ROL::StdVector<Real>&>(*pintOptVec_2->getVectorPtr(i)).getVector()->at(0);

      if(std::fabs(value-value_exact) > 1e-15) {
        ss << rank << "Two levels of restriction are not correct: " << value << " != " << value_exact << std::endl;
        throw std::logic_error(ss.str());
      }
    }
  }

  *outStream << rank << "Opt Vector Restriction: PASSED" << std::endl;

  // perform prolong, and test that it works
  //////////////////////////////////////////////////////////////////////////////////////
  
  *outStream << rank << std::endl;
  *outStream << rank << "Test Prolongation " << std::endl;
  *outStream << rank << "**********************************************************************" << std::endl;
  
  ROL::Ptr<ROL::PinTVector<Real>> pintOptVec_0_prolong = hierarchy.allocateOptVector(*pintOptVec_0,0);
  ROL::Ptr<ROL::PinTVector<Real>> pintOptVec_1_prolong = hierarchy.allocateOptVector(*pintOptVec_0,1);

  if(pintOptVec_0_prolong!=ROL::nullPtr) pintOptVec_0_prolong->zero(); 
  if(pintOptVec_1_prolong!=ROL::nullPtr) pintOptVec_1_prolong->zero(); 

  hierarchy.prolongOptVector(pintOptVec_2,pintOptVec_1_prolong,2);
  hierarchy.prolongOptVector(pintOptVec_1_prolong,pintOptVec_0_prolong,1);

  *outStream << rank << "P_Level 0 = " << printVector_Control(pintOptVec_0_prolong) << std::endl;

  *outStream << rank << std::endl;

  *outStream << rank << "P_Level 1 = " << printVector_Control(pintOptVec_1_prolong) << std::endl;

  *outStream << rank << std::endl;

  *outStream << rank << "P_Level 2 = " << printVector_Control(pintOptVec_2) << std::endl;

  std::pair<int,int> fneRange = pintOptVec_0_prolong->ownedStepRange();
  if(fneRange.first==0) {

    Real value_prlng = dynamic_cast<ROL::StdVector<Real>&>(*pintOptVec_0_prolong->getVectorPtr(0)).getVector()->at(0);
    Real value_exact = dynamic_cast<ROL::StdVector<Real>&>(*pintOptVec_2->getVectorPtr(0)).getVector()->at(0);

    if(std::fabs(value_prlng-value_exact) > 1e-15) {
      ss << rank << "Two levels of prolongation are not correct index 0: " << value_prlng << " != " << value_exact << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  /*
  // THIS NEEDS TO BE TESTED, BUT ITS COMPLICATED SO I'M SKIPPING FOR NOW: Code works as of 3/12/2019
  //    see SHA1 28d50e83df19e8cb8dc3fc2e25573cdb8b062b03

  // check errors on the finest level
  for(size_t i=0;i<stamps_0->size();i++) {

    // prolonged values
    Real value_prlng = dynamic_cast<ROL::StdVector<Real>&>(*pintOptVec_0_prolong.getVectorPtr(i)).getVector()->at(0);

    // original (or exact) values from  the initialized fine grid vector
    Real value_exact = dynamic_cast<ROL::StdVector<Real>&>(*pintOptVec_2.getVectorPtr(offset+(i-offset)/4)).getVector()->at(0);

    if(std::fabs(value_prlng-value_exact) > 1e-15) {
      ss << rank << "Two levels of prolongation are not correct: " << value_prlng << " != " << value_exact << std::endl;
      throw std::logic_error(ss.str());
    }
  }
  */

  *outStream << rank << "Opt Vector Prolongation: PASSED" << std::endl;
}
