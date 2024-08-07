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
#include "opfactory.hpp"
#include "hilbert.hpp"

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
    ROL::Ptr<BinaryAdvDiffFactory<RealT>> factory;
    ROL::Ptr<ROL::PEBBL::IntegerProblem<RealT>> problem;
    ROL::Ptr<ROL::Solver<RealT>> solver;
    ROL::Ptr<ROL::Vector<RealT>> z, omega, uomega, uz, du;
    RealT err(0);
    int order = parlist->sublist("Problem").get("Hilbert Curve Order",6);
    for (int i = 0; i < order; ++i) {
      parlist->sublist("Problem").set("Hilbert Curve Order",i+1);
      factory = ROL::makePtr<BinaryAdvDiffFactory<RealT>>(*parlist,comm,outStream);
      bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false); 
      if (checkDeriv) factory->check(*outStream);

      /*************************************************************************/
      /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
      /*************************************************************************/
      problem = factory->build();
      problem->finalize(false,true,*outStream);
      if (checkDeriv) problem->check(true,*outStream);
      solver = ROL::makePtr<ROL::Solver<RealT>>(problem, *parlist);
      solver->solve(*outStream);
      z = problem->getPrimalOptimizationVector();
      factory->getState(uz,z);

      // Sum Up Rounding
      omega = z->clone(); omega->zero();
      std::vector<RealT> &zdata = *ROL::dynamicPtrCast<PDE_OptVector<RealT>>(z)->getParameter()->getVector();
      std::vector<RealT> &odata = *ROL::dynamicPtrCast<PDE_OptVector<RealT>>(omega)->getParameter()->getVector();
      int n = std::pow(2,i+1), d(0);
      RealT g1(0), g2(0);
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          hilbert::xy2d(i+1,j,k,d);
          g1 += zdata[d];
          g2 += static_cast<RealT>(1) - zdata[d];
          if (g1 >= g2) {
            odata[d] = static_cast<RealT>(1);
            g1 -= static_cast<RealT>(1);
          }
          else {
            g2 -= static_cast<RealT>(1);
          }
        }
      }

      // Solve PDE
      factory->getState(uomega,omega);

      // Print
      std::stringstream uname, zname, oname, xname, yname;
      uname << "state_" << i+1 << ".txt";
      zname << "control_relaxed_" << i+1 << ".txt";
      oname << "control_binary_" << i+1 << ".txt";
      xname << "X_" << i+1 << ".txt";
      yname << "Y_" << i+1 << ".txt";
      factory->getAssembler()->outputTpetraVector(ROL::dynamicPtrCast<ROL::TpetraMultiVector<RealT>>(uomega)->getVector(),uname.str());
      std::ofstream zfile, ofile, xfile, yfile;
      zfile.open(zname.str());
      ofile.open(oname.str());
      xfile.open(xname.str());
      yfile.open(yname.str());
      int x(0), y(0);
      for (int j = 0; j < n*n; ++j) {
        zfile << zdata[j] << std::endl;
        ofile << odata[j] << std::endl;
        hilbert::d2xy(i+1, j, x, y);
        xfile << x << std::endl;
        yfile << y << std::endl;
      }
      zfile.close();
      ofile.close();
      xfile.close();
      yfile.close();

      du = uz->clone();
      du->set(*uz); du->axpy(static_cast<RealT>(-1),*uomega);
      err = du->norm();
      *outStream << "State Error: " << err << std::endl;
    }
    factory->print(*outStream);

    //*outStream << "OPTIMAL CONTROLS" << std::endl;
    //for (int i=0; i<controlDim; ++i) {
    //  *outStream << (*z_ptr)[i] << "  ";
    //}
    //*outStream << std::endl;
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
