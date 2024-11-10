// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "ROL_Solver.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"
#include "ROL_PEBBL_TeuchosBranchHelper.hpp"
#include "opfactory.hpp"
#include "extractQP.hpp"
#include "branchHelper.hpp"
#include "branching.hpp"

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
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
    ROL::Ptr<BinaryAdvDiffFactory<RealT>>       factory;
    ROL::Ptr<ROL::PEBBL::IntegerProblem<RealT>> problem;
    ROL::Ptr<ROL::Solver<RealT>> solver;
    ROL::Ptr<ROL::Vector<RealT>>                z, u;
    int nx = parlist->sublist("Problem").get("Number of X-Cells",64);
    int ny = parlist->sublist("Problem").get("Number of Y-Cells",32);
    factory = ROL::makePtr<BinaryAdvDiffFactory<RealT>>(*parlist,comm,outStream);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false); 
    if (checkDeriv) factory->check(*outStream);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    problem = factory->build();
    problem->finalize(false,true,*outStream);
    if (checkDeriv) problem->check(true,*outStream);
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
      ROL::Ptr<ROL::PEBBL::BranchHelper<RealT>> bHelper;
      if (factory->controlType())
        bHelper = ROL::makePtr<TpetraAdvDiffBranchHelper<RealT>>(intTol,method);
      else
        bHelper = ROL::makePtr<StdAdvDiffBranchHelper<RealT>>(intTol,method);
      ROL::Ptr<AdvDiffBranching<RealT>> branching
        = ROL::makePtr<AdvDiffBranching<RealT>>(factory,parlist,bHelper,verbosity,outStream,incheur);
      ROL::PEBBL::BranchAndBound<RealT> pebbl(branching);
      pebbl.solve(argc,argv,*outStream);
      z = pebbl.getSolution()->clone();
      z->set(*pebbl.getSolution());
    }
    else {
      if (solveQP) {
        extractQP<RealT> qp(problem->getObjective(), problem->getPrimalOptimizationVector(), problem->getBoundConstraint());
        problem = qp();
      }
      solver = ROL::makePtr<ROL::Solver<RealT>>(problem, *parlist);
      solver->solve(*outStream);
      z = problem->getPrimalOptimizationVector();
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
