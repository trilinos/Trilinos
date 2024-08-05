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
#include "Teuchos_DefaultMpiComm.hpp"

#include "ROL_PinTVector.hpp"
#include "ROL_PinTVectorCommunication_Tpetra.hpp"

#include "Tpetra_Core.hpp"

typedef double RealT;
typedef double ElementT;

typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<RealT, LO, GO, Node> MV;


int main(int argc, char* argv[]) 
{
  typedef ROL::Ptr<MV> MVP;
  typedef ROL::Ptr<ROL::Vector<RealT>> PtrVector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

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
    ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_Tpetra<RealT>>();

    ROL::Ptr<const Teuchos::Comm<int>> mpiSpaceComm = ROL::makePtr<Teuchos::MpiComm<int>>(pintComm->getSpaceCommunicator());

    int numEntries = 20;
    ROL::Ptr<Map> map = ROL::makePtr<Map>(numEntries,0,mpiSpaceComm);
 
    *outStream << "Testing checkVector" << std::endl; 

    {
      // allocate state vector
      MVP x_data = ROL::makePtr<MV>(map,1,true);
      MVP y_data = ROL::makePtr<MV>(map,1,true);
      MVP z_data = ROL::makePtr<MV>(map,1,true);
 
      x_data->randomize();
      y_data->randomize();
      z_data->randomize();
   
      PtrVector x_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(x_data);
      PtrVector y_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(y_data);
      PtrVector z_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_data);
  
      std::vector<int> stencil = {-1,0};
  
      ROL::Ptr<ROL::PinTVector<RealT>> x_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,x_vec,3*numRanks,stencil);
      ROL::Ptr<ROL::PinTVector<RealT>> y_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,y_vec,3*numRanks,stencil);
      ROL::Ptr<ROL::PinTVector<RealT>> z_pint = ROL::makePtr<ROL::PinTVector<RealT>>(pintComm,vectorComm,z_vec,3*numRanks,stencil);
  
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

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector()->getLocalView<Kokkos::HostSpace>();
          
        TEUCHOS_ASSERT(p_array.rank==2);
        TEUCHOS_ASSERT(p_array.dimension(0)==numEntries);
        TEUCHOS_ASSERT(p_array.dimension(1)==1);

        for(int i=0;i<p_array.dimension(0);i++) {
          bool correct = (p_array(i,0) == (myRank)*100+2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left boundary exchange failed: expected " << myRank*100+2 << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left

    *outStream << "Passed left boundary exchange" << std::endl;

    *outStream << "Testing boundary exchange (right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {1,0};

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector()->getLocalView<Kokkos::HostSpace>();

        TEUCHOS_ASSERT(p_array.rank==2);
        TEUCHOS_ASSERT(p_array.dimension(0)==numEntries);
        TEUCHOS_ASSERT(p_array.dimension(1)==1);

        for(int i=0;i<p_array.dimension(0);i++) {
          bool correct = (p_array(i,0)== (myRank+2)*100-2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of right boundary exchange failed: expected " << (myRank+2)*100-2 << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check right

    *outStream << "Passed right boundary exchange" << std::endl;

    *outStream << "Testing boundary exchange (left,right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {-1,0,1};

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(-1)).getVector()->getLocalView<Kokkos::HostSpace>();

        TEUCHOS_ASSERT(p_array.rank==2);
        TEUCHOS_ASSERT(p_array.dimension(0)==numEntries);
        TEUCHOS_ASSERT(p_array.dimension(1)==1);

        for(int i=0;i<p_array.dimension(0);i++) {
          bool correct = (p_array(i,0)== (myRank)*100+2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << myRank*100+2 << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

      if(myRank!=2) { // no right boundary exchange to check
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getRemoteBufferPtr(p_pint->numOwnedSteps())).getVector()->getLocalView<Kokkos::HostSpace>();

        for(int i=0;i<p_array.dimension(0);i++) {
          bool correct = (p_array(i,0)== (myRank+2)*100-2); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << (myRank+2)*100-2 << " found " << p_array(i,0) << std::endl;
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

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector()->getLocalView<Kokkos::HostSpace>();

        TEUCHOS_ASSERT(p_array.rank==2);
        TEUCHOS_ASSERT(p_array.dimension(0)==numEntries);
        TEUCHOS_ASSERT(p_array.dimension(1)==1);

        for(int i=0;i<p_array.dimension(0);i++) {
          double correct_value = ((myRank)*100+2) + ((myRank+1)*100-5); 
          bool correct = (p_array(i,0)== correct_value);
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left boundary exchange failed: expected " << correct_value << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left

    *outStream << "Passed left sum boundary exchange" << std::endl;

    *outStream << "Testing sum boundary exchange (right stencil)" << std::endl;
 
    {
      std::vector<int> stencil = {1,0};

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector()->getLocalView<Kokkos::HostSpace>();

        for(int i=0;i<p_array.dimension(0);i++) {
          double correct_value = ((myRank)*100-2) + ((myRank-1)*100+5); 
          bool correct = (p_array(i,0) == correct_value);
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of right boundary exchange failed: expected " << correct_value << " found " << p_array(i,0) << std::endl;
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

      MVP p_data = ROL::makePtr<MV>(map,1,true);
      p_data->putScalar(1.0);
      PtrVector p_vec = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_data);
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
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(2)).getVector()->getLocalView<Kokkos::HostSpace>();

        for(int i=0;i<p_array.dimension(0);i++) {
          double correct_value = ((myRank)*100+2) + ((myRank+1)*100-5); 
          bool correct = (p_array(i,0)== correct_value); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << correct_value << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

      if(myRank!=0) {
        dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector()->sync<Kokkos::HostSpace>();
        auto p_array = dynamic_cast<ROL::TpetraMultiVector<RealT>&>(*p_pint->getVectorPtr(0)).getVector()->getLocalView<Kokkos::HostSpace>();

        for(int i=0;i<p_array.dimension(0);i++) {
          double correct_value = ((myRank)*100-2) + ((myRank-1)*100+5); 
          bool correct = (p_array(i,0)== correct_value); 
          if(not correct) { 
            std::stringstream ss;
            ss << procStr << "Checking of left/right boundary exchange failed: expected " << correct_value << " found " << p_array(i,0) << std::endl;
            throw std::logic_error("Rank " + ss.str());
          }
        }
      } // end if myRank

    } // end check left/rigth

    *outStream << "Passed left/right sum boundary exchange" << std::endl;
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
