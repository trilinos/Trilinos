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

/*! \file  example_01.cpp
    \brief Shows how to solve the binary advection-diffusion control problem.
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
//#include <fenv.h>

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_PEBBL_Driver.hpp"
#include "ROL_TeuchosBranchHelper_PEBBL.hpp"
#include "branchHelper.hpp"

#include "opfactory.hpp"

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile(filename);
    Teuchos::RCP<Teuchos::ParameterList> tparlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, tparlist.ptr() );

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    bool useQP = parlist->sublist("Problem").get("Use Quadratic Program",false);
    ROL::Ptr<ROL::OptimizationProblemFactory<RealT>> factory;
    ROL::Ptr<ROL::BranchHelper_PEBBL<RealT>> bHelper;
    if (!useQP) {
      factory = ROL::makePtr<BinaryAdvDiffFactory<RealT>>(*parlist,comm,outStream);
      bHelper = ROL::makePtr<PDEOPT_BranchHelper_PEBBL<RealT>>();
    }
    else {
      factory = ROL::makePtr<BinaryAdvDiffQPFactory<RealT>>(*parlist,comm,outStream);
      bHelper = ROL::makePtr<ROL::TeuchosBranchHelper_PEBBL<int,RealT>>();
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::ROL_PEBBL_Driver<RealT> pebbl(factory,parlist,bHelper,3,outStream);
    pebbl.solve(argc,argv,*outStream);

    //*outStream << "OPTIMAL CONTROLS" << std::endl;
    //for (int i=0; i<controlDim; ++i) {
    //  *outStream << (*z_ptr)[i] << "  ";
    //}
    //*outStream << std::endl;
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
