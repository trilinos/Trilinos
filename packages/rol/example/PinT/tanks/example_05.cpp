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

#include <iostream>
#include <iomanip>
#include <random>
#include <utility>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_PinTConstraint.hpp"
#include "ROL_PinTVectorCommunication_StdVector.hpp"

#include "Tanks_DynamicConstraint.hpp"
#include "Tanks_SerialConstraint.hpp"
#include "Tanks_ConstraintCheck.hpp"

using RealT = double;
using size_type = std::vector<RealT>::size_type;

void run_test_getTimeStampsByLevel(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void run_test_allocateVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void run_test_restrictVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void run_test_prolongVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void run_test_buildCommunicators(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);

int main( int argc, char* argv[] ) 
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  
  int myRank = Teuchos::GlobalMPISession::getRank();
//  int numProcs = Teuchos::GlobalMPISession::getNProc();

//  int iprint     = argc - 1;
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  try {

/*
    (*outStream) << "getTimeStampsByLevel" << std::endl;
    (*outStream) << "**************************************************" << std::endl;
    run_test_getTimeStampsByLevel(MPI_COMM_WORLD, outStream);

    (*outStream) << "allocateVectors" << std::endl;
    (*outStream) << "**************************************************" << std::endl;
    run_test_allocateVectors(MPI_COMM_WORLD, outStream);

    (*outStream) << "restrictVectors" << std::endl;
    (*outStream) << "**************************************************" << std::endl;
    run_test_restrictVectors(MPI_COMM_WORLD, outStream);

    (*outStream) << "prolongVectors" << std::endl;
    (*outStream) << "**************************************************" << std::endl;
    run_test_prolongVectors(MPI_COMM_WORLD, outStream);
*/

    (*outStream) << "buildCommunicators " << myRank << std::endl;
    (*outStream) << "**************************************************" << std::endl;
    run_test_buildCommunicators(MPI_COMM_WORLD, outStream);
  }
  catch (std::logic_error &err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }
  catch (...) {
    *outStream << "ERROR: Unknown exception on " << myRank << "\n";
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

struct ConstraintData {
  ROL::Ptr<ROL::PinTConstraint<RealT>> constraint;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;
  ROL::Ptr<const ROL::PinTCommunicators> communicators;
  ROL::Ptr<std::vector<ROL::TimeStamp<RealT>>> timeStamp;
  RealT totalTime;
};

ConstraintData
buildPinTConstraint(int local_Nt,MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
  using Bounds            = ROL::Bounds<RealT>;
  using PartitionedVector = ROL::PartitionedVector<RealT>;
  using State             = Tanks::StateVector<RealT>;
  using Control           = Tanks::ControlVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<RealT>>();
  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  ROL::Ptr< ROL::PinTVector<RealT>> state;
  ROL::Ptr< ROL::PinTVector<RealT>> control;

  auto  pl_ptr = ROL::getParametersFromXmlFile("tank-parameters.xml");
  auto& pl     = *pl_ptr;
  auto  con    = Tanks::DynamicConstraint<RealT>::create(pl);
  auto  height = pl.get("Height of Tank",              10.0  );
  auto  Qin00  = pl.get("Corner Inflow",               100.0 );
  auto  h_init = pl.get("Initial Fluid Level",         2.0   );
  auto  nrows  = static_cast<size_type>( pl.get("Number of Rows"   ,3) );
  auto  ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );
  auto  Nt     = numRanks*local_Nt; // use a power of two
  auto  T      = pl.get("Total Time", 20.0);

  RealT dt = T/Nt;

  ROL::Ptr<ROL::Vector<RealT>> initial_cond;
  {
    // control
    auto  z      = Control::create( pl, "Control (z)"     );    
    auto  vz     = z->clone( "Control direction (vz)"     );
    auto  z_lo   = z->clone( "Control Lower Bound (z_lo)" );
    auto  z_bnd  = makePtr<Bounds>( *z_lo );
    z_lo->zero();

    // State
    auto u_new     = State::create( pl, "New state (u_new)"   );
    auto u_old     = u_new->clone( "Old state (u_old)"        );
    auto u_initial = u_new->clone( "Initial conditions"       );
    auto u_new_lo  = u_new->clone( "State lower bound (u_lo)" );
    auto u_new_up  = u_new->clone( "State upper bound (u_up)" );
    auto u         = PartitionedVector::create( { u_old,    u_new    } );
    auto u_lo      = PartitionedVector::create( { u_new_lo, u_new_lo } );
    auto u_up      = PartitionedVector::create( { u_new_up, u_new_up } );
  
    u_lo->zero();
    u_up->setScalar( height );
    auto u_bnd = makePtr<Bounds>(u_new_lo,u_new_up);

    (*z)(0,0) = Qin00;

    for( size_type i=0; i<nrows; ++i ) {
      for( size_type j=0; j<ncols; ++j ) {
        u_old->h(i,j) = h_init;
        u_initial->h(i,j) = h_init;
      }
    }

    state        = ROL::buildStatePinTVector<RealT>(   communicators, vectorComm, Nt,     u_old);
    control      = ROL::buildControlPinTVector<RealT>( communicators, vectorComm, Nt,         z);

    initial_cond = u_initial;
    state->getVectorPtr(-1)->set(*u_initial);   // set the initial condition
  }

  assert(local_Nt==control->numOwnedSteps());

  auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(local_Nt);

  // do a scan to get the starting time
  double my_final_time = dt*local_Nt;
  double time_offset  = 0.0;
  MPI_Exscan(&my_final_time,&time_offset,1,MPI_DOUBLE,MPI_SUM,communicators->getTimeCommunicator());

  // (*outStream) << "Local_Nt = " << local_Nt << std::endl;
  // (*outStream) << "Dt = " << dt << std::endl;
  // (*outStream) << "TIME OFFSET = " << time_offset << std::endl;

  for( size_type k=0; k< static_cast<size_type>(local_Nt); ++k ) {
    timeStamp->at(k).t.resize(2);
    timeStamp->at(k).t.at(0) = k*dt+time_offset;
    timeStamp->at(k).t.at(1) = (k+1)*dt+time_offset;
  }

  ConstraintData cd;
  cd.constraint = ROL::makePtr<ROL::PinTConstraint<RealT>>(con,initial_cond,timeStamp);
  cd.state = state;
  cd.control = control;
  cd.communicators = communicators;
  cd.timeStamp = timeStamp;
  cd.totalTime = T;

  return cd;
}

void run_test_getTimeStampsByLevel(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
//  using Bounds            = ROL::Bounds<RealT>;
//  using PartitionedVector = ROL::PartitionedVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  int local_Nt = 16;
  Ptr<ROL::PinTConstraint<RealT>> pint_constraint = buildPinTConstraint(local_Nt,comm,outStream).constraint;

  auto level0_stamps = pint_constraint->getTimeStampsByLevel(0);
  double dt = level0_stamps->at(0).t[1]- level0_stamps->at(0).t[0];

  for(int l=0;l<3;l++) {
    auto level_stamps = pint_constraint->getTimeStampsByLevel(l);
  
    // make sure that the number of time steps is right (note the +1 is a result of initial condition weridness)
    int expected_size = 1+local_Nt/std::pow(2,l);
    if( level_stamps->size() != static_cast<size_type>(expected_size) ) {
      std::stringstream ss;
      ss << "Number of levels " << l << " stamps are incorrect, expected " << expected_size << " but got " << level_stamps->size() << std::endl;
      throw std::logic_error(ss.str());
    }

    // check the initial and final time node to ensure that they match on all levels
    if(level0_stamps->at(1).t[0]!=level_stamps->at(1).t[0]) {
      std::stringstream ss;
      ss << "Zero time node extrema does not match at level " << l << std::endl;
      throw std::logic_error(ss.str());
    }
    if(level0_stamps->at(level0_stamps->size()-1).t[1]!=level_stamps->at(level_stamps->size()-1).t[1]) {
      std::stringstream ss;
      ss << "Final time node extrema does not match at level " << l << std::endl;
      throw std::logic_error(ss.str());
    }

    // check that all the time step sizes are correct at this level
    for( size_type k=1; k<level_stamps->size(); k++ ) {
      double level_dt = level_stamps->at(k).t[1]- level_stamps->at(k).t[0];
      
      if(std::fabs(dt * std::pow(2,l) - level_dt)>1e-14) {
        std::stringstream ss;
        ss << "Time step size is not correct on level " << l << std::endl;
        throw std::logic_error(ss.str());
      }
    }

    // print out some stuff for grins
    std::stringstream ss;
    for( size_type k=1; k<level_stamps->size(); k++ ) {
      ss << "level " << l << " times = " << level_stamps->at(k).t[0] << " " << level_stamps->at(k).t[1] << std::endl;
    }
    (*outStream) << ss.str() << std::endl;

  }
}

void run_test_allocateVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
//  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
//  using Bounds            = ROL::Bounds<RealT>;
//  using PartitionedVector = ROL::PartitionedVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  int local_Nt = 16;
  ConstraintData cd = buildPinTConstraint(local_Nt,comm,outStream);
  Ptr<ROL::PinTConstraint<RealT>> pint_constraint = cd.constraint;

  auto level0_stamps = pint_constraint->getTimeStampsByLevel(0);
//  double dt = level0_stamps->at(0).t[1]- level0_stamps->at(0).t[0];

  // check allocation of the control and state vectors...note the need for a reference
  for(int level=0;level<3;level++) {
    auto state = pint_constraint->allocateSimVector(*cd.state,level);
    auto control = pint_constraint->allocateOptVector(*cd.control,level);
    ROL::PinTVector<RealT> & pint_state = dynamic_cast<ROL::PinTVector<RealT>&>(*state);
    ROL::PinTVector<RealT> & pint_control = dynamic_cast<ROL::PinTVector<RealT>&>(*control);

    std::stringstream ss;

    if(pint_state.stencil().size()!=2) {
      ss << "Stencil for the state on level " << level << " is incorrect" << std::endl;
      throw std::logic_error(ss.str());
    }

    if(pint_control.stencil().size()!=1) {
      ss << "Stencil for the control on level " << level << " is incorrect" << std::endl;
      throw std::logic_error(ss.str());
    }
    
    if(pint_state.numOwnedSteps()*std::pow(2,level)!=local_Nt) {
      ss << "Owned steps for the state on level " << level << " is incorrect" << std::endl;
      throw std::logic_error(ss.str());
    }

    if(pint_control.numOwnedSteps()*std::pow(2,level)!=local_Nt) {
      ss << "Owned steps for the control on level " << level << " is incorrect" << std::endl;
      throw std::logic_error(ss.str());
    }
      
  }
}

void run_test_restrictVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
//  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
//  using Bounds            = ROL::Bounds<RealT>;
//  using PartitionedVector = ROL::PartitionedVector<RealT>;

  auto & out = *outStream;
  std::stringstream ss;  // for errors

  RealT tolerance = 1e-14;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  int local_Nt = 16;
  ConstraintData cd = buildPinTConstraint(local_Nt,comm,outStream);
  Ptr<ROL::PinTConstraint<RealT>> pint_constraint = cd.constraint;

  for(int level=0;level<3;level++) {
    (*outStream) << "Checking level " << level << std::endl;

    auto state       = pint_constraint->allocateSimVector(*cd.state,level); // fine vectors
    auto control     = pint_constraint->allocateOptVector(*cd.control,level);
    auto state_crs   = pint_constraint->allocateSimVector(*cd.state,level+1); // coarse vectors (note the fine level can be use for reference)
    auto control_crs = pint_constraint->allocateOptVector(*cd.control,level+1);

    {
      ROL::PinTVector<RealT> & pint_state   = dynamic_cast<ROL::PinTVector<RealT>&>(*state); 
      ROL::PinTVector<RealT> & pint_control = dynamic_cast<ROL::PinTVector<RealT>&>(*control); 

      out << "Compare steps in control and state: " << pint_state.numOwnedSteps() << " " << pint_control.numOwnedSteps() << std::endl;
      if(pint_state.numOwnedSteps()!=pint_control.numOwnedSteps()) {
        ss << "restrictVectors: control and state time step count's don't match" << std::endl;
        throw std::logic_error(ss.str());
      }
  
      // fill up the fine level vectors
      pint_state.getVectorPtr(-1)->setScalar(-myRank-1.0); 
      for( size_type k=0;k < static_cast<size_type>(pint_state.numOwnedSteps()); k++ ) {
        pint_state.getVectorPtr(k)->setScalar(myRank*100.0+k); 
        pint_control.getVectorPtr(k)->setScalar(myRank*100.0+k); 
      }
    }

    out << "Calling state restriction" << std::endl;
    pint_constraint->restrictSimVector(*state,*state_crs);

    out << "Calling control restriction" << std::endl;
    pint_constraint->restrictOptVector(*control,*control_crs);

    // check the restriction operators
    out << "Checking restriction methods" << std::endl;
    { 
      ROL::PinTVector<RealT> & pint_state   = dynamic_cast<ROL::PinTVector<RealT>&>(*state); 
//      ROL::PinTVector<RealT> & pint_control = dynamic_cast<ROL::PinTVector<RealT>&>(*control); 

      ROL::PinTVector<RealT> & pint_state_crs   = dynamic_cast<ROL::PinTVector<RealT>&>(*state_crs); 
      ROL::PinTVector<RealT> & pint_control_crs = dynamic_cast<ROL::PinTVector<RealT>&>(*control_crs); 

      if(pint_state.dimension()!=pint_state_crs.dimension()*2) {
        ss << "restrictVectors: coarse state is the wrong size: found " << pint_state_crs.dimension() 
           << ", expected " << pint_state.dimension()/2 << std::endl;
        throw std::logic_error(ss.str());
      }

      // check the coarse state (read comment for restrict to understand what this is checking)
      //////////////////////////////////////////////
      
      int dim_state = pint_state_crs.getVectorPtr(-1)->dimension();

      // check that the virtual state is just copied from the fine vector
      if(std::fabs(pint_state.getVectorPtr(-1)->norm()-pint_state_crs.getVectorPtr(-1)->norm())>tolerance) {
        ss << "restrictVectors: norm of restricted state vector virtual variable is incorrect. " << std::endl;
        throw std::logic_error(ss.str());
      }

      // check the interior vectors
      for( size_type k=0; k<static_cast<size_type>(pint_state_crs.numOwnedSteps()-1); k++ ) {
        RealT exact_value = ((myRank*100.0+(2.0*k)) + (myRank*100.0+(2.0*k+1)) + (myRank*100.0+(2.0*k+2)))/3.0;
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_state);

        RealT coarse_norm = pint_state_crs.getVectorPtr(k)->norm();
        if(std::fabs(coarse_norm-exact_norm)>tolerance) {
          ss << "restrictVectors: norm of restricted state vector " << k << " is incorrect. "
             << "Expected " << exact_norm << ", found " << coarse_norm << std::endl;
          throw std::logic_error(ss.str());
        }
      }

      // check the final vector
      {
        int k = pint_state_crs.numOwnedSteps()-1; 
        
        RealT exact_value = 0.0;
        if(myRank==numRanks-1)
          exact_value = ((myRank*100.0+(2.0*k)) + (myRank*100.0+(2.0*k+1)))/2.0;
                                // on the last processor
        else
          exact_value = ((myRank*100.0+(2.0*k)) + (myRank*100.0+(2.0*k+1)) + ((myRank+1)*100.0))/3.0;
                                // get the neighobr
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_state);

        RealT coarse_norm = pint_state_crs.getVectorPtr(k)->norm();
        out << "Restriction of static vector on processor boundary: expected " << exact_norm << ", found " << coarse_norm << std::endl;
        if(std::fabs(coarse_norm-exact_norm)>tolerance) {
          ss << "restrictVectors: norm of restricted state vector " << k << " is incorrect (last on proc):"
             << "Expected " << exact_norm << ", found " << coarse_norm << std::endl;
        }
      }

      // check the coarse contro (read comment for restrict to understand what this is checking)
      //////////////////////////////////////////////
      
      int dim_control = pint_control_crs.getVectorPtr(0)->dimension();

      // check the interior vectors
      for(int k=0;k<pint_control_crs.numOwnedSteps();k++) {
        RealT exact_value = ((myRank*100.0+(2.0*k)) + (myRank*100.0+(2.0*k+1)))/2.0;
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_control);

        RealT coarse_norm = pint_control_crs.getVectorPtr(k)->norm();
        out << "restrictVectors: norm of restricted control vector " << k << ": "
            << "Expected " << exact_norm << ", found " << coarse_norm << std::endl;
        if(std::fabs(coarse_norm-exact_norm)>tolerance) {
          ss << "restrictVectors: norm of restricted control vector " << k << " is incorrect. "
             << "Expected " << exact_norm << ", found " << coarse_norm << std::endl;
          throw std::logic_error(ss.str());
        }
      }
    }
  }
}

void run_test_prolongVectors(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
//  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
//  using Bounds            = ROL::Bounds<RealT>;
//  using PartitionedVector = ROL::PartitionedVector<RealT>;

  auto & out = *outStream;
  std::stringstream ss;  // for errors

  RealT tolerance = 1e-14;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  int local_Nt = 16;
  ConstraintData cd = buildPinTConstraint(local_Nt,comm,outStream);
  Ptr<ROL::PinTConstraint<RealT>> pint_constraint = cd.constraint;

  for(int level=0;level<3;level++) {
    (*outStream) << "Checking level " << level << std::endl;

    auto state       = pint_constraint->allocateSimVector(*cd.state,level); // fine vectors
    auto control     = pint_constraint->allocateOptVector(*cd.control,level);
    auto state_crs   = pint_constraint->allocateSimVector(*cd.state,level+1); // coarse vectors (note the fine level can be use for reference)
    auto control_crs = pint_constraint->allocateOptVector(*cd.control,level+1);

    {
      ROL::PinTVector<RealT> & pint_state_crs   = dynamic_cast<ROL::PinTVector<RealT>&>(*state_crs); 
      ROL::PinTVector<RealT> & pint_control_crs = dynamic_cast<ROL::PinTVector<RealT>&>(*control_crs); 

      out << "Compare steps in control and state: " << pint_state_crs.numOwnedSteps() << " " << pint_control_crs.numOwnedSteps() << std::endl;
      if(pint_state_crs.numOwnedSteps()!=pint_control_crs.numOwnedSteps()) {
        ss << "prolongVectors: control and state time step count's don't match" << std::endl;
        throw std::logic_error(ss.str());
      }
  
      // fill up the fine level vectors
      pint_state_crs.getVectorPtr(-1)->setScalar(myRank*100.0-1); 
      for(int k=0;k<pint_state_crs.numOwnedSteps();k++) {
        pint_state_crs.getVectorPtr(k)->setScalar(myRank*100.0+k); 
        pint_control_crs.getVectorPtr(k)->setScalar(myRank*100.0+k); 
      }
    }

    out << "Calling state prolongation" << std::endl;
    pint_constraint->prolongSimVector(*state_crs,*state);

    out << "Calling control prolongation" << std::endl;
    pint_constraint->prolongOptVector(*control_crs,*control);

    // check the prolongion operators
    out << "Checking prolongation methods" << std::endl;
    { 
      ROL::PinTVector<RealT> & pint_state   = dynamic_cast<ROL::PinTVector<RealT>&>(*state); 
      ROL::PinTVector<RealT> & pint_control = dynamic_cast<ROL::PinTVector<RealT>&>(*control); 

      ROL::PinTVector<RealT> & pint_state_crs   = dynamic_cast<ROL::PinTVector<RealT>&>(*state_crs); 
      ROL::PinTVector<RealT> & pint_control_crs = dynamic_cast<ROL::PinTVector<RealT>&>(*control_crs); 

      if(pint_state.dimension()!=pint_state_crs.dimension()*2) {
        ss << "prolongVectors: coarse state is the wrong size: found " << pint_state_crs.dimension() 
           << ", expected " << pint_state.dimension()/2 << std::endl;
        throw std::logic_error(ss.str());
      }

      // check the coarse state (read comment for prolong to understand what this is checking)
      //////////////////////////////////////////////
      
      int dim_state = pint_state_crs.getVectorPtr(-1)->dimension();

      // check that the initial state is the sum of the previous processor
      {
        int num_previous_steps = pint_state_crs.numOwnedSteps()-1;
        RealT exact_value = 0.0;
        if(myRank==0)
          exact_value = (myRank*100.0+0)/2.0;
        else
          exact_value = (myRank*100.0+0+(myRank-1)*100.0+num_previous_steps)/2.0;

        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_state);
        RealT fine_norm = pint_state.getVectorPtr(0)->norm();
        if(std::fabs(fine_norm-exact_norm)>tolerance) {
          ss << "prolongVectors: norm of prolonged inital state vector is incorrect. " << std::endl;
          throw std::logic_error(ss.str());
        }
      }

      // check the injection vectors, note this includes the virtual variable
      for(int k=-1;k<pint_state_crs.numOwnedSteps();k++) {
        RealT exact_value = myRank*100.0+k;
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_state);

        RealT fine_norm = pint_state.getVectorPtr(2*k+1)->norm();
        if(std::fabs(fine_norm-exact_norm)>tolerance) {
          ss << "prolongVectors: norm of prolonged state vector " << 2*k+1 << " is incorrect. "
             << "Expected " << exact_norm << ", found " << fine_norm << std::endl;
          throw std::logic_error(ss.str());
        }
      }

      // check the injection vectors, note this includes the virtual variable
      for(int k=0;k<pint_state_crs.numOwnedSteps()-1;k++) {
        RealT exact_value = (myRank*100.0+k + myRank*100.0+k+1)/2.0;
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_state);

        RealT fine_norm = pint_state.getVectorPtr(2*k+2)->norm();
        if(std::fabs(fine_norm-exact_norm)>tolerance) {
          ss << "ERROR prolongVectors: norm of prolonged state vector " << 2*k+2 << " is incorrect. "
             << "Expected " << exact_norm << ", found " << fine_norm << std::endl;
          throw std::logic_error(ss.str());
        }
      }

      // check the coarse control (read comment for prolong to understand what this is checking)
      //////////////////////////////////////////////
      
      int dim_control = pint_control_crs.getVectorPtr(0)->dimension();

      // check the interior vectors
      for(int k=0;k<pint_control_crs.numOwnedSteps();k++) {
        RealT exact_value = myRank*100.0+k; 
        RealT exact_norm  = std::sqrt(std::pow(exact_value,2.0)*dim_control);

        RealT fine_norm_0 = pint_control.getVectorPtr(2*k+0)->norm();
        RealT fine_norm_1 = pint_control.getVectorPtr(2*k+1)->norm();
        out << "prolongVectors: norm of prolonged control vector " << 2*k << " and " << 2*k+1 << ": "
            << "Expected " << exact_norm << ", found " << fine_norm_0 << " and " << fine_norm_1 << std::endl;
        if(std::fabs(fine_norm_0-exact_norm)>tolerance) {
          ss << "prolongVectors: norm of prolonged control vector " << 2*k << " is incorrect. "
             << "Expected " << exact_norm << ", found " << fine_norm_0 << std::endl;
          throw std::logic_error(ss.str());
        }
        if(std::fabs(fine_norm_1-exact_norm)>tolerance) {
          ss << "prolongVectors: norm of prolonged control vector " << 2*k+1 << " is incorrect. "
             << "Expected " << exact_norm << ", found " << fine_norm_1 << std::endl;
          throw std::logic_error(ss.str());
        }
      }
    }
  }
}

void run_test_buildCommunicators(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;

  // auto & out = *outStream; // Unused
  std::stringstream ss;  // for errors

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  int local_Nt = 16;
  ConstraintData cd = buildPinTConstraint(local_Nt,comm,outStream);
  Ptr<ROL::PinTConstraint<RealT>> pint_constraint = cd.constraint;
  pint_constraint->applyMultigrid(3);
  pint_constraint->buildLevels(*cd.state);

  auto comm_0 = pint_constraint->getLevelCommunicators(0);
  auto comm_1 = pint_constraint->getLevelCommunicators(1);
  auto comm_2 = pint_constraint->getLevelCommunicators(2);

  auto ts_0 = pint_constraint->getTimeStampsByLevel_repart(0);
  auto ts_1 = pint_constraint->getTimeStampsByLevel_repart(1);
  auto ts_2 = pint_constraint->getTimeStampsByLevel_repart(2);

  if(myRank==0) {
    if(ROL::is_nullPtr(comm_0) || 
       ROL::is_nullPtr(comm_1) || 
       ROL::is_nullPtr(comm_2)) {
      ss << "Rank " << myRank << " has wrong communicators." << std::endl;
      throw std::logic_error(ss.str());
    }

    if(ROL::is_nullPtr(ts_0) || 
       ROL::is_nullPtr(ts_1) || 
       ROL::is_nullPtr(ts_2)) {
      ss << "Rank " << myRank << " has wrong time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }
  }
  else if(myRank==1) {
    if( ROL::is_nullPtr(comm_0) || 
        ROL::is_nullPtr(comm_1) || 
       !ROL::is_nullPtr(comm_2)) {
      ss << "Rank " << myRank << " has wrong communicators." << std::endl;
      throw std::logic_error(ss.str());
    }

    if( ROL::is_nullPtr(ts_0) || 
        ROL::is_nullPtr(ts_1) || 
       !ROL::is_nullPtr(ts_2)) {
      ss << "Rank " << myRank << " has wrong time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }
  }
  else if(myRank==2) {
    if( ROL::is_nullPtr(comm_0) || 
       !ROL::is_nullPtr(comm_1) || 
       !ROL::is_nullPtr(comm_2)) {
      ss << "Rank " << myRank << " has wrong communicators." << std::endl;
      throw std::logic_error(ss.str());
    }

    if( ROL::is_nullPtr(ts_0) || 
       !ROL::is_nullPtr(ts_1) || 
       !ROL::is_nullPtr(ts_2)) {
      ss << "Rank " << myRank << " has wrong time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }
  }
  else if(myRank==3) {
    if( ROL::is_nullPtr(comm_0) || 
       !ROL::is_nullPtr(comm_1) || 
       !ROL::is_nullPtr(comm_2)) {
      ss << "Rank " << myRank << " has wrong communicators." << std::endl;
      throw std::logic_error(ss.str());
    }

    if( ROL::is_nullPtr(ts_0) || 
       !ROL::is_nullPtr(ts_1) || 
       !ROL::is_nullPtr(ts_2)) {
      ss << "Rank " << myRank << " has wrong time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  // check the communicators are build correctly
  if(not ROL::is_nullPtr(comm_0)) {
    if(comm_0->getTimeSize()!=4 || comm_0->getTimeRank()!=myRank) {
      ss << "Rank " << myRank << " has wrong level 0 communicator." << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  if(not ROL::is_nullPtr(comm_1)) {
    if(comm_1->getTimeSize()!=2 || comm_1->getTimeRank()!=myRank) {
      ss << "Rank " << myRank << " has wrong level 1 communicator." << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  if(not ROL::is_nullPtr(comm_2)) {
    if(comm_2->getTimeSize()!=1 || comm_2->getTimeRank()!=myRank) {
      ss << "Rank " << myRank << " has wrong level 2 communicator." << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  // check the time stamps are built correctly
  if(not ROL::is_nullPtr(ts_0)) {
    // account for extra time stamp
    if(ts_0->size()!=cd.timeStamp->size()+1) { 
      ss << "Rank " << myRank << " has incorrect level 0 computed time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }

    for(size_t i = 1; i<ts_0->size();i++) {
      bool value = ts_0->at(i).t == cd.timeStamp->at(i-1).t;
      // *outStream << "checking time stamp " << i << " on rank " << myRank << ": " << value << std::endl;

      if(not value) {
        ss << "Rank " << myRank << " has incorrect level 0 computed time stamps." << std::endl;
        throw std::logic_error(ss.str());
      }
    }
  }

  if(not ROL::is_nullPtr(ts_1)) {
    // account for extra time stamp
    if(ts_1->size()!=cd.timeStamp->size()+1) { 
      ss << "Rank " << myRank << " has incorrect level 1 computed time stamps. " << ts_1->size() << " " << cd.timeStamp->size() << std::endl;
      throw std::logic_error(ss.str());
    }
  }

  if(not ROL::is_nullPtr(ts_2)) {
    if(myRank!=0)
      throw std::logic_error("Only expecting rank 0, found something else...something is seriously wrong!");

    // account for extra time stamp
    if(ts_2->size()!=cd.timeStamp->size()+1) { 
      ss << "Rank " << myRank << " has incorrect level 2 computed time stamps." << std::endl;
      throw std::logic_error(ss.str());
    }
 
    bool result = true;
    // note that each level contributes to this one, and this is the last one!    
    auto stamp_m1 = ts_2->at(0);
    result &= (stamp_m1.t.size() == 2);

    RealT dt = stamp_m1.t[1] - stamp_m1.t[0];
    result &= (std::fabs(dt - cd.totalTime/(ts_2->size()-1.0)) <= 1e-14); // make sure step size is correct
    *outStream << std::endl << dt << " or " << cd.totalTime/(ts_2->size()-1.0) << std::endl;

    for(size_type i=1;i<ts_2->size();i++) {
      auto stamp = ts_2->at(i);

      result &= (stamp.t.size() == 2);
      result &= (stamp.t[0] == (i-1) * dt);
      result &= (stamp.t[1] == (i+0) * dt);

      *outStream << "Stamp " << i << " " << stamp.t[0] << " " << stamp.t[1] << std::endl;
    }

    if(not result) {
      throw std::logic_error("BUZZZZ!");
    }
  }
}
