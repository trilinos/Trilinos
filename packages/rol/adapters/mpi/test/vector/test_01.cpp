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

  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;

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
  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  std::string procStr = std::to_string(myRank) + "/" + std::to_string(numRanks) + ": ";

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  try {

    int spatialProcs = 1;
    ROL::Ptr<ROL::PinTCommunicators> pintComm = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spatialProcs);
 
    // if(numRanks!=3) {
    //   throw std::logic_error("Three processors are required to run this test!");
    // }
    
    *outStream << "Testing checkVector" << std::endl; 

    ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<RealT>>();
    
    {
      // allocate state vector
      std::vector<RealT> x_data(2); x_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); x_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
      std::vector<RealT> y_data(2); y_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); y_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
      std::vector<RealT> z_data(2); z_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); z_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
   
      PtrVector x_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(x_data));
      PtrVector y_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(y_data));
      PtrVector z_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(z_data));
  
      std::vector<int> stencil = {-1,0};
  
      ROL::Ptr<ROL::PinTVector<RealT>> x_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,x_vec,3*numRanks,-1,1,2);
      ROL::Ptr<ROL::PinTVector<RealT>> y_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,y_vec,3*numRanks,-1,1,2);
      ROL::Ptr<ROL::PinTVector<RealT>> z_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,z_vec,3*numRanks,-1,1,2);
  
      /*
      *outStream << "X = " << std::endl;
      x_pint->print(*outStream);
      *outStream << "\n\nY = " << std::endl;
      y_pint->print(*outStream);
      *outStream << "\n\nZ = " << std::endl;
      z_pint->print(*outStream);
      */

      if(x_pint->numOwnedVectors()!=6) {
        std::stringstream ss;
        ss << procStr << "Number owned vectors is " << x_pint->numOwnedVectors() << " is not 6!";
        throw std::logic_error("Rank " + ss.str());
      }
  
      if(x_pint->numOwnedSteps()!=3) {
        std::stringstream ss;
        ss << procStr << "Number owned steps is " << x_pint->numOwnedSteps() << " is not 3!";
        throw std::logic_error("Rank " + ss.str());
      }
  
      std::vector<RealT> consistency = x_pint->checkVector(*y_pint,*z_pint,(myRank==0 ? true : false),*outStream);
      ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );
      if (checkvec.norm() > std::sqrt(errtol)) {
        errorFlag++;
        std::stringstream ss;
        ss << procStr << "Failed check vector!";
        throw std::logic_error("Rank " + ss.str());
      }
    }

    // test boundary exchange (left to right)
    /////////////////////////////////////////////////////////////////////////////////////////////

    *outStream << "Testing boundary exchange" << std::endl;

    {
      int replicate = 2;
      std::vector<int> stencil = {-1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,-1,1,replicate);

      TEUCHOS_ASSERT(  ROL::is_nullPtr(p_pint->getVectorPtr(-1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 2)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 3)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 4)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 5)) );
      TEUCHOS_ASSERT(  ROL::is_nullPtr(p_pint->getVectorPtr( 6)) );

      p_pint->getVectorPtr( 0)->scale((myRank+1)*100+0);
      p_pint->getVectorPtr( 1)->scale((myRank+1)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank+1)*100+2);
      p_pint->getVectorPtr( 3)->scale((myRank+1)*100+3);
      p_pint->getVectorPtr( 4)->scale((myRank+1)*100+4);
      p_pint->getVectorPtr( 5)->scale((myRank+1)*100+5);

      p_pint->boundaryExchangeLeftToRight();

      if(myRank!=0) { // no left boundary exchange to check
        for(int i=0;i<2;i++) {
          const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(i)).getVector();

          for(auto v : p_std) {
            bool correct = (v== (myRank)*100+4+i); 
            if(not correct) { 
              std::stringstream ss;
              ss << procStr << "Checking of left boundary exchange failed: expected " << myRank*100+5 << " found " << v << std::endl;
              throw std::logic_error("Rank " + ss.str());
            }
          }
        }
      } // end if myRank

      p_pint->boundaryExchangeRightToLeft();

      if(myRank!=2) { // no right boundary exchange to check
        for(int i=0;i<2;i++) {
          const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(i)).getVector();

          for(auto v : p_std) {
            bool correct = (v== (myRank+2)*100+i); 
            if(not correct) { 
              std::stringstream ss;
              ss << procStr << "Checking of right to left boundary exchange failed: expected " << (myRank+2)*100+0 << " found " << v << std::endl;
              throw std::logic_error("Rank " + ss.str());
            }
          }
        }
      } // end if myRank

    } // end check left

    *outStream << "Passed boundary exchange" << std::endl;

#if 0
    // test boundary exchange (insert into buffer from vector)
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    *outStream << "Testing boundary exchange (left stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {-1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(-1)) );  // backwards time is owned for this stencil
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 2)) );

      p_pint->getVectorPtr(-1)->scale((myRank+1)*100-1); // backwards time is owned for this stencil
      p_pint->getVectorPtr( 0)->scale((myRank+1)*100+0);
      p_pint->getVectorPtr( 1)->scale((myRank+1)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank+1)*100+2);

      p_pint->boundaryExchange();

      if(myRank!=0) { // no left boundary exchange to check
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector();

        for(auto v : p_std) {
          bool correct = (v== (myRank)*100+2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left boundary exchange failed: expected " << myRank*100+2 << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left

    *outStream << "Passed left boundary exchange" << std::endl;

    *outStream << "Testing boundary exchange (right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(2)) );

      p_pint->getVectorPtr( 0)->scale((myRank+1)*100+-2);
      p_pint->getVectorPtr( 1)->scale((myRank+1)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank+1)*100+2);
      p_pint->getVectorPtr( 3)->scale((myRank+1)*100+3);

      p_pint->boundaryExchange();

      if(myRank!=2) { // no right boundary exchange to check
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector();

        for(auto v : p_std) {
          bool correct = (v== (myRank+2)*100-2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of right boundary exchange failed: expected " << (myRank+2)*100-2 << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check right

    *outStream << "Passed right boundary exchange" << std::endl;

    *outStream << "Testing boundary exchange (left,right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {-1,0,1};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(2)) );

      p_pint->getVectorPtr(-1)->scale((myRank+1)*100-3);
      p_pint->getVectorPtr( 0)->scale((myRank+1)*100-2); // this vector is actually sent (because of the stencil)
      p_pint->getVectorPtr( 1)->scale((myRank+1)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank+1)*100+2); // this vector is actually sent (because of the stencil)
      p_pint->getVectorPtr( 3)->scale((myRank+1)*100+3);

      p_pint->boundaryExchange();

      if(myRank!=0) { // no left boundary exchange to check
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector();

        for(auto v : p_std) {
          bool correct = (v== (myRank)*100+2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << myRank*100+2 << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

      if(myRank!=2) { // no right boundary exchange to check
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector();

        for(auto v : p_std) {
          bool correct = (v== (myRank+2)*100-2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << (myRank+2)*100-2 << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left/rigth

    *outStream << "Passed left/right boundary exchange" << std::endl;

    // test boundary exchange (sum into vector from buffer)
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    *outStream << "Testing sum boundary exchange (left stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {-1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(-1)) );  // backwards time is owned for this stencil
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 2)) );
      TEUCHOS_ASSERT(  ROL::is_nullPtr(p_pint->getVectorPtr( 3)) );

      p_pint->getRemoteBufferPtr(-1)->scale((myRank)*100-5);

      p_pint->getVectorPtr(-1)->scale((myRank)*100-1); 
      p_pint->getVectorPtr( 0)->scale((myRank)*100+0);
      p_pint->getVectorPtr( 1)->scale((myRank)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank)*100+2);

      p_pint->boundaryExchangeSumInto();

      if(myRank!=2) { 
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector();

        for(auto v : p_std) {
          double correct_value = ((myRank)*100+2) + ((myRank+1)*100-5); 
          bool correct = (v== correct_value);
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left boundary exchange failed: expected " << correct_value << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left

    *outStream << "Passed left sum boundary exchange" << std::endl;

    *outStream << "Testing sum boundary exchange (right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(  ROL::is_nullPtr(p_pint->getVectorPtr(-1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 2)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 3)) );

      p_pint->getRemoteBufferPtr(3)->scale((myRank)*100+5); 

      p_pint->getVectorPtr( 0)->scale((myRank)*100+-2);
      p_pint->getVectorPtr( 1)->scale((myRank)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank)*100+2);
      p_pint->getVectorPtr( 3)->scale((myRank)*100+3);

      p_pint->boundaryExchangeSumInto();

      if(myRank!=0) { 
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector();

        for(auto v : p_std) {
          double correct_value = ((myRank)*100-2) + ((myRank-1)*100+5); 
          bool correct = (v== correct_value);
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of right boundary exchange failed: expected " << correct_value << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check right

    *outStream << "Passed right sum boundary exchange" << std::endl;

    *outStream << "Testing sum boundary exchange (left,right stencil)" << std::endl;
 
    // comments are correct in this block: ECC 1/25/2018
    {
      std::vector<int> stencil = {-1,0,1};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr(-1)) ); // this comes from the stencil
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 0)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 1)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 2)) );
      TEUCHOS_ASSERT( !ROL::is_nullPtr(p_pint->getVectorPtr( 3)) ); // this comes from the stencil

      p_pint->getRemoteBufferPtr(-1)->scale((myRank)*100-5); // this vector is being sent to left
      p_pint->getRemoteBufferPtr( 3)->scale((myRank)*100+5); // this vector is being sent to right

      p_pint->getVectorPtr(-1)->scale((myRank)*100-3);
      p_pint->getVectorPtr( 0)->scale((myRank)*100-2); // this vector is recv (from left) 
      p_pint->getVectorPtr( 1)->scale((myRank)*100+1);
      p_pint->getVectorPtr( 2)->scale((myRank)*100+2); // this vector is recv (from right)
      p_pint->getVectorPtr( 3)->scale((myRank)*100+3);

      p_pint->boundaryExchangeSumInto();

      if(myRank!=2) {
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector();

        for(auto v : p_std) {
          double correct_value = ((myRank)*100+2) + ((myRank+1)*100-5); 
          bool correct = (v== correct_value); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << correct_value << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

      if(myRank!=0) {
        const std::vector<RealT> & p_std = *dynamic_cast<ROL::StdVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector();

        for(auto v : p_std) {
          double correct_value = ((myRank)*100-2) + ((myRank-1)*100+5); 
          bool correct = (v== correct_value); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << correct_value << " found " << v << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left/rigth

    *outStream << "Passed left/right sum boundary exchange" << std::endl;
#endif
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
