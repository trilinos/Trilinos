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

#include "ODEConstraint_TimeSimOpt.hpp"

#include "ROL_Constraint_PinTSimOpt.hpp"

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
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

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
 
  int numSteps = 9;
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
  typedef ROL::Ptr<const ROL::Vector<Real>> CPtrVector;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  double totaltime = 1.75; // we will do 1.75 periods, this way the final result is non-trivial
  double dt = totaltime/numSteps;
  double omega = 2.0*M_PI;
  double tol = 1e-15;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  ROL::Ptr< ROL::PinTVector<Real>> state;
  ROL::Ptr< ROL::PinTVector<Real>> control, control_unit;

  {
    // allocate local state vector
    ROL::Ptr<std::vector<Real>> u_data_ptr = ROL::makePtr<std::vector<Real>>(2);
    std::vector<Real> & u_data = *u_data_ptr;
    u_data[0] = 0.0;
    u_data[1] = omega;
    PtrVector u = ROL::makePtr<ROL::StdVector<Real>>(u_data_ptr);
      // this requirement to use a copy constructor tripped me up for a while

    // allocate control vector
    std::vector<Real> z_data(1);
    z_data[0] = 1.0;
    PtrVector z = ROL::makePtr<ROL::StdVector<Real>>(ROL::makePtr<std::vector<Real>>(z_data));

    state        = buildPinTVector(communicators,numSteps,{-1,0}   /* state stencil */,u);
    control      = buildPinTVector(communicators,numSteps,{0}    /* control stencil */,z);
    control_unit = buildPinTVector(communicators,numSteps,{0}    /* control stencil */,z);

    state->getVectorPtr(-1)->set(*u);   // set the initial condition
  }

  ROL::Ptr<ROL::Constraint_TimeSimOpt<Real>> step_constraint = ROL::makePtr<ODE_Constraint<Real>>(dt,omega);

  ROL::Ptr<ROL::Constraint_PinTSimOpt<Real>> pint_constraint = ROL::makePtr<ROL::Constraint_PinTSimOpt<Real>>(step_constraint);

  PtrVector u = state;
  PtrVector c = state->clone();
  PtrVector z = control->clone();

  ROL::RandomizeVector<RealT>(*u); // this randomization doesn't really matter here as 'u' is complete determine by the solve
  ROL::RandomizeVector<RealT>(*z); 
  z->scale(1.0);

  // check the solve
  //////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking solve" << std::endl;

    double solveNorm = pint_constraint->checkSolve(*u,*z,*c,true,*outStream);

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

  // This computes and returns the last value of the 'u' component and shares
  // it with all processors.
  /////////////////////////////////////////////////////////////////////////////
  double final_result = 0.0;
  {
    PtrVector z_unit = control_unit->clone();
    z->scale(0.0); // z can't be random as the different processor counts won't give the same values
    z->axpy(3.0,*z_unit); 
 
    *outStream << "ZNORM = " << z->norm() << std::endl;
  
    pint_constraint->solve(*c,*u,*z,tol);
  
    // pull out the last value
    if(myRank==numRanks-1) {
      PtrVector final_state = state->getVectorPtr(state->numOwnedSteps()-1);
      *outStream << "Final state = " << final_state << std::endl;
      final_state->print(*outStream);
  
      const std::vector<Real> & fs_stdvec = *dynamic_cast<const ROL::StdVector<Real>&>(*final_state).getVector();
      final_result = fs_stdvec[0];
    }
  
    MPI_Bcast(&final_result,1,MPI_DOUBLE,numRanks-1,comm);
  }

  return final_result;
}
