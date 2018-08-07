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

#include "test_05.hpp"
#include "ROL_PEBBL_Driver.hpp"
#include "ROL_StdBranchHelper_PEBBL.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int N = 10;
    ROL::Ptr<ROL::Vector<RealT>> x    = ROL::makePtr<ROL::StdVector<RealT>>(N,0.0);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    std::vector<RealT> alpha(N,0.0), beta(N,0.0);
    *outStream << std::endl << std::endl;
    for (int i = 0; i < N; ++i) {
      alpha[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      beta[i]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    }
    *outStream << std::endl << std::endl;
    ROL::Ptr<ROL::Objective<RealT>> obj = ROL::makePtr<Objective_SimpleBinary<RealT>>(alpha,beta);
    /**********************************************************************************************/
    /************************* CONSTRUCT CONSTRAINT ***********************************************/
    /**********************************************************************************************/
    RealT budget  = static_cast<RealT>(parlist->get("Budget", 3));
    bool  useIneq = parlist->get("Use Inequality", true);
    ROL::Ptr<ROL::Constraint<RealT>> icon;
    if (useIneq) {
      icon = ROL::makePtr<Constraint_SimpleBinary<RealT>>(0.0);
    }
    else {
      icon = ROL::makePtr<Constraint_SimpleBinary<RealT>>(static_cast<int>(budget));
    }
    ROL::Ptr<ROL::Vector<RealT>>          imul = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    ROL::Ptr<ROL::Vector<RealT>>          ilo  = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    ROL::Ptr<ROL::Vector<RealT>>          iup  = ROL::makePtr<ROL::StdVector<RealT>>(1,budget);
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilo,iup);
    /**********************************************************************************************/
    /************************* CONSTRUCT BOUND CONSTRAINT *****************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::Vector<RealT>>          xl  = ROL::makePtr<ROL::StdVector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>>          xu  = ROL::makePtr<ROL::StdVector<RealT>>(N,1.0);
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd = ROL::makePtr<ROL::Bounds<RealT>>(xl,xu);
    /**********************************************************************************************/
    /************************* SOLVE **************************************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    if (useIneq) {
      problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(obj,x,bnd,icon,imul,ibnd);
    }
    else {
      problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(obj,x,bnd,icon,imul);
    }
    problem->check(*outStream);
    ROL::Ptr<ROL::BranchHelper_PEBBL<RealT>> bHelper
      = ROL::makePtr<ROL::StdBranchHelper_PEBBL<RealT>>();
    ROL::ROL_PEBBL_Driver<RealT> pebbl(problem,parlist,bHelper,3,outStream);
    pebbl.solve(argc,argv,*outStream);
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
