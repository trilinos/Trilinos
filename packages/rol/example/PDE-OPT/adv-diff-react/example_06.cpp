// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve the stochastic advection-diffusion problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_PrimalDualRisk.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_stoch_adv_diff.hpp"
#include "obj_stoch_adv_diff.hpp"
#include "mesh_stoch_adv_diff.hpp"

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int>> serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();

  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input_ex06.xml");

    // Problem dimensions
    const int controlDim = 9, stochDim = 37;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshManager_stoch_adv_diff<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<PDE_stoch_adv_diff<RealT>> pde
      = ROL::makePtr<PDE_stoch_adv_diff<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT>> pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    pdeCon->getAssembler()->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>>  u_ptr = pdeCon->getAssembler()->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  p_ptr = pdeCon->getAssembler()->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  r_ptr = pdeCon->getAssembler()->createResidualVector();
    ROL::Ptr<std::vector<RealT>>     z_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    ROL::Ptr<ROL::Vector<RealT>> u, p, r, z;
    u = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,pdeCon->getAssembler());
    p = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,pdeCon->getAssembler());
    r = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,pdeCon->getAssembler());
    z = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(z_ptr));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_stoch_adv_diff<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_stoch_adv_diff<RealT>>();
    RealT stateCost   = parlist->sublist("Problem").get("State Cost",1.e5);
    RealT controlCost = parlist->sublist("Problem").get("Control Cost",1.e0);
    std::vector<RealT> wts = {stateCost, controlCost};
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,pdeCon->getAssembler());
    bool storage = parlist->sublist("Problem").get("Use State and Adjoint Storage",true);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, u, z, p, storage, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<ROL::Vector<RealT>> zlo, zup;
    zlo = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim, 0.0));
    zup = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim, 1.0));
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlo, zup);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT>> bounds(stochDim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** SOLVE STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    z->zero();
    ROL::Ptr<ROL::Problem<RealT>> problem
      = ROL::makePtr<ROL::Problem<RealT>>(objReduced, z);
    problem->addBoundConstraint(bnd);
    ROL::PrimalDualRisk<RealT> solver(problem, sampler, *parlist);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      problem->check(true,*outStream);
      solver.check(*outStream);
    }
    solver.run(*outStream);

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    //std::clock_t timer_print = std::clock();
    // Output control to file
    if ( myRank == 0 ) {
      std::ofstream zfile;
      zfile.open("control.txt");
      for (int i = 0; i < controlDim; i++) {
        zfile << std::scientific << std::setprecision(15);
        zfile << (*z_ptr)[i] << std::endl;
      }
      zfile.close();
    }

    RealT tol(1e-12);
    Teuchos::Array<RealT> res(1,0);
    pdeCon->solve(*r,*u,*z,tol);
    r_ptr->norm2(res.view(0,1));

    /*************************************************************************/
    /***************** CHECK RESIDUAL NORM ***********************************/
    /*************************************************************************/
    *outStream << "Residual Norm: " << res[0] << std::endl << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
