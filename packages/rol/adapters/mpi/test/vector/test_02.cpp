// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_PinTVector.hpp"
#include "ROL_PinTVectorCommunication_StdVector.hpp"

typedef double RealT;

int main(int argc, char* argv[]) 
{

//  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  std::string procStr = std::to_string(myRank) + "/" + std::to_string(numRanks) + ": ";

  assert(numRanks==4);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  try {

    int spatialProcs = 1;
    ROL::Ptr<ROL::PinTCommunicators> pintComm = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spatialProcs);
 
    ROL::Ptr<ROL::PinTCommunicators> coarseComm = pintComm->buildCoarseCommunicators();

    std::stringstream ss;
    if(myRank>=2) {
      if(not ROL::is_nullPtr(coarseComm)) {
        ss << "Rank " << myRank << " is included in the coarse communicator! This is wrong." << std::endl;
        throw std::logic_error(ss.str());
      }
    }
    else {
      if(ROL::is_nullPtr(coarseComm)) {
        ss << "Rank " << myRank << " is not included in the coarse communicator! This is wrong." << std::endl;
        throw std::logic_error(ss.str());
      }
     
      if(coarseComm->getSpaceSize()!=1) {
        ss << "Rank " << myRank << " has space size not equal to 1" << std::endl;
        throw std::logic_error(ss.str());
      }

      if(coarseComm->getTimeSize()!=2) {
        ss << "Rank " << myRank << " has time size not equal to 2" << std::endl;
        throw std::logic_error(ss.str());
      }

      if(coarseComm->getTimeRank()!=myRank) {
        ss << "Rank " << myRank << " has time rank not equal to itself" << std::endl;
        throw std::logic_error(ss.str());
      }
    }
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  int errors = std::abs(errorFlag);
  MPI_Allreduce(&errors,&errorFlag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
