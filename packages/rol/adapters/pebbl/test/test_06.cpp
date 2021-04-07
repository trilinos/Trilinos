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

#include "test_06.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"
#include "ROL_PEBBL_StdBranchHelper.hpp"

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
    /************************* GET PROBLEM PARAMETERS *********************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input_06.xml");
    /**********************************************************************************************/
    /************************* CONSTRUCT PROBLEM FACTORY ******************************************/
    /**********************************************************************************************/
    int M = parlist->sublist("Problem").get("Number of Facilities",5);
    int N = parlist->sublist("Problem").get("Number of Customers",10);
    std::vector<RealT> facX(M), facY(M), c(M);
    for (int i = 0; i < M; ++i) {
      facX[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      facY[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      c[i]    = static_cast<RealT>(100)*static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    }
    std::vector<RealT> cusX(N), cusY(N);
    for (int j = 0; j < N; ++j) {
      cusX[j] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      cusY[j] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    }
    std::vector<RealT> q(M*N,1.0);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        q[i + j*M] = static_cast<RealT>(50)*std::sqrt(std::pow(facX[i]-cusX[j],2) + std::pow(facY[i]-cusY[j],2));
      }
    }
    ROL::Ptr<FacilityLocationFactory<RealT>> factory
      = ROL::makePtr<FacilityLocationFactory<RealT>>(c,q,*parlist);
    /**********************************************************************************************/
    /************************* SOLVE **************************************************************/
    /**********************************************************************************************/
    RealT intTol = parlist->sublist("Problem").get("Integrality Tolerance",1e-6);
    int method    = parlist->sublist("Problem").get("Branching Method",0);
    int verbosity = parlist->sublist("Problem").get("BB Output Level",1);
    ROL::Ptr<ROL::PEBBL::StdBranchHelper<RealT>> bHelper
      = ROL::makePtr<ROL::PEBBL::StdBranchHelper<RealT>>(intTol,method);
    ROL::Ptr<ROL::PEBBL::Branching<RealT>> branching
      = ROL::makePtr<FacilityLocationBranching<RealT>>(factory,parlist,bHelper,verbosity,outStream);
    ROL::PEBBL::BranchAndBound<RealT> pebbl(branching);
    pebbl.solve(argc,argv,*outStream);
    *outStream << "Solution Vector:" << std::endl << "  ";
    pebbl.getSolution()->print(*outStream);
    *outStream << std::endl;
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
