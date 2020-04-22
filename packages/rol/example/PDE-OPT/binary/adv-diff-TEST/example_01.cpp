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
#include "ROL_OptimizationSolver.hpp"
#include "ROL_PEBBL_Driver.hpp"
#include "ROL_TeuchosBranchHelper_PEBBL.hpp"
#include "opfactory.hpp"
#include "extractQP.hpp"
#include "branchHelper.hpp"
#include "branching.hpp"

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
    ROL::Ptr<BinaryAdvDiffFactory<RealT>>     factory;
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    ROL::Ptr<ROL::OptimizationSolver<RealT>>  solver;
    ROL::Ptr<ROL::Vector<RealT>>              z, u;
    int nx = parlist->sublist("Problem").get("Number of X-Cells",64);
    int ny = parlist->sublist("Problem").get("Number of Y-Cells",32);
    factory = ROL::makePtr<BinaryAdvDiffFactory<RealT>>(*parlist,comm,outStream);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false); 
    if (checkDeriv) factory->check(*outStream);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    problem = factory->build();
    if (checkDeriv) problem->check(*outStream);
    factory->getAssembler()->printMeshData(*outStream);
    factory->print(*outStream);
    bool solveQP = parlist->sublist("Problem").get("Solve as QP",false);
    bool useBB   = parlist->sublist("Problem").get("Solve using BB",true);
    solveQP = (useBB ? false : solveQP);
    if (useBB) {
      RealT intTol = parlist->sublist("Problem").get("Integrality Tolerance",1e-6);
      int method    = parlist->sublist("Problem").get("Branching Method",0);
      int incheur   = parlist->sublist("Problem").get("Incumbent Heuristic",0);
      int verbosity = parlist->sublist("Problem").get("BB Output Level",0);
      ROL::Ptr<ROL::BranchHelper_PEBBL<RealT>> bHelper
        = ROL::makePtr<PDEOPT_BranchHelper_PEBBL<RealT>>(intTol,method);
      ROL::Ptr<ADVDIFF_Branching<RealT>> branching
        = ROL::makePtr<ADVDIFF_Branching<RealT>>(factory,parlist,bHelper,verbosity,outStream,incheur);
      ROL::ROL_PEBBL_Driver<RealT> pebbl(branching);
      pebbl.solve(argc,argv,*outStream);
      z = pebbl.getSolution();
    }
    else {
      if (solveQP) {
        extractQP<RealT> qp(problem->getObjective(), problem->getSolutionVector(), problem->getBoundConstraint());
        problem = qp();
      }
      solver = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem, *parlist);
      solver->solve(*outStream);
      z = problem->getSolutionVector();
    }
    factory->getState(u,z);

    // Print
    bool usePC = parlist->sublist("Problem").get("Piecewise Constant Controls", true);
    std::stringstream uname, zname, xname, yname;
    if (!solveQP) {
      if (!usePC) {
        uname << "state.txt";
        zname << "control.txt";
        factory->getAssembler()->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<RealT>>(u)->getVector(),uname.str());
        factory->getAssembler()->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<RealT>>(z)->getVector(),zname.str());
      }
      else {
        uname << "state_"   << nx << "_" << ny << ".txt";
        zname << "control_" << nx << "_" << ny << ".txt";
        xname << "X_"       << nx << "_" << ny << ".txt";
        yname << "Y_"       << nx << "_" << ny << ".txt";
        factory->getAssembler()->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<RealT>>(u)->getVector(),uname.str());
        std::ofstream zfile, xfile, yfile;
        zfile.open(zname.str());
        xfile.open(xname.str());
        yfile.open(yname.str());
        std::vector<RealT> &zdata = *ROL::dynamicPtrCast<ROL::StdVector<RealT>>(z)->getVector();
        zfile << std::scientific << std::setprecision(18);
        for (int j = 0; j < nx; ++j) {
          for (int k = 0; k < ny; ++k) {
            zfile << zdata[j+k*nx] << std::endl;
            xfile << j << std::endl;
            yfile << k << std::endl;
          }
        }
        zfile.close();
        xfile.close();
        yfile.close();
      }
    }
    else {
      zname << "control_" << nx << "_" << ny << ".txt";
      xname << "X_"       << nx << "_" << ny << ".txt";
      yname << "Y_"       << nx << "_" << ny << ".txt";
      std::ofstream zfile, xfile, yfile;
      zfile.open(zname.str());
      xfile.open(xname.str());
      yfile.open(yname.str());
      Teuchos::SerialDenseVector<int,RealT> &zdata = *ROL::dynamicPtrCast<ROL::TeuchosVector<int,RealT>>(z)->getVector();
      for (int j = 0; j < nx; ++j) {
        for (int k = 0; k < ny; ++k) {
          zfile << zdata[j+k*nx] << std::endl;
          xfile << j << std::endl;
          yfile << k << std::endl;
        }
      }
      zfile.close();
      xfile.close();
      yfile.close();
    }
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
