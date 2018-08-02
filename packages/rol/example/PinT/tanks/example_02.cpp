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

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_PinTConstraint.hpp"

#include "Tanks_DynamicConstraint.hpp"
#include "Tanks_ConstraintCheck.hpp"

typedef double Real;

/**
 * Run the test using a different communicator. This returns the last value of
 * the last step for the 'u' variable defined in the ODEConstraint_TimeSimOpt class.
 */
double run_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int numSteps);

int main(int argc, char* argv[]) 
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  MPI_Comm subComm;

  // split the comm so if three processors are used you can run a 3, 2 and 1 processor
  // case all at once
  if(numRanks!=3) {
    *outStream << "WARNING: This test is most effective with three processors" << std::endl;
  }
  else {
    MPI_Comm_split(MPI_COMM_WORLD,myRank==0,myRank,&subComm); 
  }
 
  int numSteps = 10;
  int errorFlag  = 0;

  try {
    double all_final = run_test(MPI_COMM_WORLD, outStream,numSteps);

    MPI_Barrier(MPI_COMM_WORLD);

    // run the serial and 2 processor cases as well
    if(numRanks==3) {

      // because of the splitting this actually runs two jobs at once
      double sub_final = run_test(subComm, outStream,numSteps);

      *outStream << "ALL processors result = " << all_final << std::endl;
      *outStream << "SUB processors result = " << sub_final << std::endl;

      if(sub_final!=all_final) {
        int subRanks = 0;
        MPI_Comm_size(subComm, &subRanks);

        std::stringstream ss;
        ss << "ALL processors result does not match SUB(p=" << subRanks << ") result: " << all_final << " != " << sub_final;
        throw std::logic_error(ss.str());
      }
    }
    else {
      *outStream << "ALL processors result = " << all_final << std::endl;
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  // cleanup split comm
  if(numRanks==3)
    MPI_Comm_free(&subComm); // cleanup

  int errors = std::abs(errorFlag);
  MPI_Allreduce(&errors,&errorFlag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}

double run_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream,int numSteps)
{
  typedef ROL::Ptr<ROL::Vector<Real>> PtrVector;

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

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  ROL::Ptr< ROL::PinTVector<Real>> state, state_unit;
  ROL::Ptr< ROL::PinTVector<Real>> control, control_unit;

  auto  pl_ptr = ROL::getParametersFromXmlFile("tank-parameters.xml");
  auto& pl     = *pl_ptr;
  auto  con    = Tanks::DynamicConstraint<RealT>::create(pl);
  auto  height = pl.get("Height of Tank",              10.0  );
  auto  Qin00  = pl.get("Corner Inflow",               100.0 );
  auto  h_init = pl.get("Initial Fluid Level",         2.0   );
  auto  nrows  = static_cast<size_type>( pl.get("Number of Rows"   ,3) );
  auto  ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );
  auto  Nt     = static_cast<size_type>( pl.get("Number of Time Steps",100) );
  auto  T       = pl.get("Total Time", 20.0);

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

    state        = ROL::buildStatePinTVector<Real>(   communicators, Nt,     u_old);
    control      = ROL::buildControlPinTVector<Real>( communicators, Nt,         z);

    initial_cond = u_initial;
  }

  auto timeStamp = makePtr<std::vector<ROL::TimeStamp<Real>>>(state->numOwnedSteps());
  for( size_type k=0; k<timeStamp->size(); ++k ) {
    timeStamp->at(k).t.resize(2);
    timeStamp->at(k).t.at(0) = k*dt;
    timeStamp->at(k).t.at(1) = (k+1)*dt;
  }

  ROL::PinTConstraint<RealT> pint_constraint(con,initial_cond,timeStamp);

  PtrVector u    = state;
  PtrVector z    = control->clone();
  PtrVector c    = state->clone();
  PtrVector jv   = c->clone();
  PtrVector v_u  = state->clone();
  PtrVector v_z  = control->clone();
  PtrVector w_u  = state->dual().clone();
  PtrVector w_z  = control->dual().clone();
  PtrVector w_c  = c->dual().clone();
  PtrVector hv_u = state->dual().clone();
  PtrVector hv_z = z->dual().clone();

  ROL::RandomizeVector<RealT>(*u); // this randomization doesn't really matter here as 'u' is complete determine by the solve
  ROL::RandomizeVector<RealT>(*z); 
  ROL::RandomizeVector<RealT>(*v_u); 
  ROL::RandomizeVector<RealT>(*v_z); 
  ROL::RandomizeVector<RealT>(*w_u); 
  ROL::RandomizeVector<RealT>(*w_z); 
  ROL::RandomizeVector<RealT>(*w_c); 

  // dynamic_cast<ROL::PinTVector<RealT>&>(*v_u).getVectorPtr(-1)->zero();
  // v_u->zero(); 
  
  double tol = 5e-14;

  // check the solve
  //////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking solve" << std::endl;

    double solveNorm = pint_constraint.checkSolve(*u,*z,*c,true,*outStream);

    if(solveNorm>tol) {
      std::stringstream ss;
      ss << "Solve failed on problem with " << numRanks << " processors and " << numSteps << ": residual = " << solveNorm
                                            << " > " << tol << std::endl;
      throw std::logic_error(ss.str());
    }
    else if(myRank==0) {
      *outStream << "Solve checked out for p=" << numRanks << " with residual = " << solveNorm << std::endl; 
    }
  }

  // check the Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Jacobian 1" << std::endl;

    auto errors = pint_constraint.checkApplyJacobian_1(*u,*z,*v_u,*jv,true,*outStream);
    if(errors[6][3]/errors[6][1] >= 1e-6)
      throw std::logic_error("Constraint apply jacobian 1 is incorrect");
  }

  // check the Jacobian_2
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Jacobian 2" << std::endl;

    auto errors = pint_constraint.checkApplyJacobian_2(*u,*z,*v_z,*jv,true,*outStream);
    if(errors[6][3]/errors[6][1] >= 1e-6)
      throw std::logic_error("Constraint apply jacobian 2 is incorrect");
  }

  // check the Adjoint Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Jacobian 1" << std::endl;

    auto error = pint_constraint.checkAdjointConsistencyJacobian_1(*w_u,*v_u,*u,*z,true,*outStream);
    if(error >= 1e-8)
      throw std::logic_error("Constraint apply adjoint jacobian 1 is incorrect");
  }

  // check the Adjoint Jacobian_2
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Jacobian 2" << std::endl;

    auto error = pint_constraint.checkAdjointConsistencyJacobian_2(*w_u,*v_z,*u,*z,true,*outStream);
    if(error >= 1e-8)
      throw std::logic_error("Constraint apply adjoint jacobian 2 is incorrect");
  }

#if 0
  // check the Adjoint Hessian 11
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Hessian_11" << std::endl;

    auto errors = pint_constraint->checkApplyAdjointHessian_11(*u,*z,*w_c,*v_u,*hv_u,true,*outStream);
    if(errors[6][3]/errors[6][1] >= 1e-4) // wow, this is really not very good accuracy
      throw std::logic_error("Constraint apply Adjoint Hessian 11 is incorrect: " + std::to_string( errors[6][3]/errors[6][1]));
  }

  // check the Adjoint Hessian_12
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Hessian_12" << std::endl;

    auto errors = pint_constraint->checkApplyAdjointHessian_12(*u,*z,*w_c,*v_u,*hv_z,true,*outStream);
    if(errors[6][3] >= 1e-5)
      throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
  }

  // check the Adjoint Hessian_21
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Hessian_21" << std::endl;

    auto errors = pint_constraint->checkApplyAdjointHessian_21(*u,*z,*w_c,*v_z,*hv_u,true,*outStream);
    if(errors[6][3] >= 1e-5)
      throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
  }

  // check the Adjoint Hessian_22
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Hessian_22" << std::endl;

    auto errors = pint_constraint->checkApplyAdjointHessian_22(*u,*z,*w_c,*v_z,*hv_z,true,*outStream);
    if(errors[6][3] >= 1e-5)
      throw std::logic_error("Constraint apply Adjoint Hessian 12 is incorrect");
  }

#endif

  // This computes and returns the last value of the 'u' component and shares
  // it with all processors.
  /////////////////////////////////////////////////////////////////////////////
  double final_result = 0.0;
  {
    PtrVector z_unit = control->clone();
    z->setScalar(3.0); 
 
    *outStream << "ZNORM = " << z->norm() << std::endl;
  
    pint_constraint.solve(*c,*u,*z,tol);
  
    // pull out the last value
    {
      PtrVector final_state = state->getVectorPtr(state->numOwnedSteps()-1);
      // *outStream << "Final state = " << final_state << std::endl;
      // final_state->print(*outStream);
  
      const std::vector<Real> & fs_stdvec = *dynamic_cast<const ROL::StdVector<Real>&>(*final_state).getVector();
      final_result = fs_stdvec[0];
    }
  
    MPI_Bcast(&final_result,1,MPI_DOUBLE,numRanks-1,comm);
  }

  return final_result;
}
