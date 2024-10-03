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
#include "branchHelper.hpp"
#include "branching.hpp"

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
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
    bool print = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    print = (myRank==0 ? print : false);
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",print);

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Ptr<BinaryStefanBoltzmannFactory<RealT>> factory;
    ROL::Ptr<ROL::PEBBL::IntegerProblem<RealT>>   problem;
    ROL::Ptr<ROL::Solver<RealT>>   solver;
    ROL::Ptr<ROL::Vector<RealT>>                  z, u;
    //factory = ROL::makePtr<BinaryStefanBoltzmannFactory<RealT>>(*parlist,comm,outStream);
    factory = ROL::makePtr<BinaryStefanBoltzmannFactory<RealT>>(*parlist,outStream);
#ifdef HAVE_MPI
    ROL::Ptr<MPI_Comm> mpicomm = ROL::makePtr<MPI_Comm>(*ROL::staticPtrCast<const Teuchos::MpiComm<int>>(comm)->getRawMpiComm());
    factory->setCommunicator(mpicomm);
#endif
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
    bool useBB = parlist->sublist("Problem").get("Solve using BB",true);
    if (useBB) {
      RealT intTol = parlist->sublist("Problem").get("Integrality Tolerance",1e-6);
      int method    = parlist->sublist("Problem").get("Branching Method",0);
      int incheur   = parlist->sublist("Problem").get("Incumbent Heuristic",0);
      int verbosity = parlist->sublist("Problem").get("BB Output Level",0);
      bool useParam = parlist->sublist("Problem").get("Use Parametric Control", true);
      ROL::Ptr<ROL::PEBBL::BranchHelper<RealT>> bHelper;
      ROL::Ptr<ROL::PEBBL::Branching<RealT>> branching;
      if (useParam) {
        bHelper   = ROL::makePtr<ROL::PEBBL::StdBranchHelper<RealT>>(intTol,method);
        branching = ROL::makePtr<StdStefanBoltzmannBranching<RealT>>(factory,parlist,bHelper,verbosity,outStream,incheur);
      }
      else {
        bHelper   = ROL::makePtr<TpetraStefanBoltzmannBranchHelper<RealT>>(intTol,method);
        branching = ROL::makePtr<TpetraStefanBoltzmannBranching<RealT>>(factory,parlist,bHelper,verbosity,outStream,incheur);
      }
      ROL::PEBBL::BranchAndBound<RealT> pebbl(branching);
      pebbl.solve(argc,argv,*outStream);
      z = pebbl.getSolution()->clone();
      z->set(*pebbl.getSolution());
    }
    else {
      solver = ROL::makePtr<ROL::Solver<RealT>>(problem, *parlist);
      solver->solve(*outStream);
      z = problem->getPrimalOptimizationVector();
    }
    factory->getState(u,z);

    // Print
    std::stringstream uname, zname;
    uname << "state.txt";
    zname << "control.txt";
    factory->printTpetraVector(u,uname.str());
    factory->printControl(z,zname.str());
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
