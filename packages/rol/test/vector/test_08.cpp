// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_PinTVector.hpp"

typedef double RealT;

int main(int argc, char* argv[]) 
{

  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
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
    
    {
      // allocate state vector
      std::vector<RealT> x_data(2); x_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); x_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
      std::vector<RealT> y_data(2); y_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); y_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
      std::vector<RealT> z_data(2); z_data[0] = ( (RealT)rand() / (RealT)RAND_MAX ); z_data[1] = ( (RealT)rand() / (RealT)RAND_MAX );
   
      PtrVector x_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(x_data));
      PtrVector y_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(y_data));
      PtrVector z_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(z_data));
  
      std::vector<int> stencil = {-1,0};
  
      ROL::Ptr<ROL::PinTVector<RealT>> x_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,x_vec,3*numRanks,stencil);
      ROL::Ptr<ROL::PinTVector<RealT>> y_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,y_vec,3*numRanks,stencil);
      ROL::Ptr<ROL::PinTVector<RealT>> z_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,z_vec,3*numRanks,stencil);
  
      /*
      *outStream << "X = " << std::endl;
      x_pint->print(*outStream);
      *outStream << "\n\nY = " << std::endl;
      y_pint->print(*outStream);
      *outStream << "\n\nZ = " << std::endl;
      z_pint->print(*outStream);
      */

      if(x_pint->numOwnedVectors()!=4) {
        std::stringstream ss;
        ss << procStr << "Number owned vectors is " << x_pint->numOwnedVectors() << " is not 4!";
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

    // test boundary exchange (insert into buffer from vector)
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    *outStream << "Testing boundary exchange (left stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {-1,0};

      std::vector<RealT> p_data(2); p_data[0] = 1.0; p_data[1] = 1.0;
      PtrVector p_vec = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtrFromRef(p_data));
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(-1)!=Teuchos::null);  // backwards time is owned for this stencil
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 2)!=Teuchos::null);

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
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr(1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr(2)!=Teuchos::null);

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
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr(1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr(2)!=Teuchos::null);

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
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(-1)!=Teuchos::null);  // backwards time is owned for this stencil
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 2)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 3)==Teuchos::null);

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
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(-1)==Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 2)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 3)!=Teuchos::null);

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
      ROL::Ptr<ROL::PinTVector<RealT>> p_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,p_vec,3*numRanks,stencil);

      TEUCHOS_ASSERT(p_pint->getVectorPtr(-1)!=Teuchos::null); // this comes from the stencil
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 0)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 1)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 2)!=Teuchos::null);
      TEUCHOS_ASSERT(p_pint->getVectorPtr( 3)!=Teuchos::null); // this comes from the stencil

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
  }
  catch (std::logic_error err) {
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
