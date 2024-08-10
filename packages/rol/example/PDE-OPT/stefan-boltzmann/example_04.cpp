// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_04.cpp
    \brief Shows how to solve the Stefan-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_stoch_stefan_boltzmann.hpp"
#include "obj_stoch_stefan_boltzmann.hpp"
#include "mesh_stoch_stefan_boltzmann.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
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
    std::string filename = "input_ex04.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    const int controlDim = 1;
    RealT tol(1e-8);

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<StochasticStefanBoltzmannPDE<RealT> > pde
      = ROL::makePtr<StochasticStefanBoltzmannPDE<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  u_ptr, p_ptr, r_ptr, zbc_ptr;
    u_ptr   = assembler->createStateVector();    u_ptr->putScalar(0.0);
    p_ptr   = assembler->createStateVector();    p_ptr->putScalar(0.0);
    r_ptr   = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    zbc_ptr = assembler->createControlVector();  zbc_ptr->putScalar(280.0);
    ROL::Ptr<std::vector<RealT> >  zp_ptr;
    zp_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,0.0);
    ROL::Ptr<ROL::Vector<RealT> > up, pp, rp, zp, xp;
    up     = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    pp     = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    rp     = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zbc;
    zbc    = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zbc_ptr,pde,assembler);
    ROL::Ptr<ROL::StdVector<RealT>> zparam;
    zparam = ROL::makePtr<ROL::StdVector<RealT>>(zp_ptr);
    zp     = ROL::makePtr<PDE_OptVector<RealT>>(zbc, zparam);
    xp     = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(3,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateCost<RealT>>(pde->getVolFE(),*parlist);
    qoi_vec[1] = ROL::makePtr<QoI_ControlCost<RealT>>(
      pde->getVolFE(),pde->getBdryFE(0),pde->getBdryCellLocIds(0),*parlist);
    qoi_vec[2] = ROL::makePtr<QoI_AdvectionCost<RealT>>();
    ROL::Ptr<StochasticStefanBoltzmannStdObjective3<RealT> > std_obj
      = ROL::makePtr<StochasticStefanBoltzmannStdObjective3<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Bounds for boundary control
    RealT lower_bc = parlist->sublist("Problem").get("Lower Control Bound", 280.0);
    RealT upper_bc = parlist->sublist("Problem").get("Upper Control Bound", 370.0);
    ROL::Ptr<Tpetra::MultiVector<> >  zlo_bc_ptr
      = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> >  zhi_bc_ptr
      = assembler->createControlVector();
    zlo_bc_ptr->putScalar(static_cast<RealT>(lower_bc));
    zhi_bc_ptr->putScalar(static_cast<RealT>(upper_bc));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zlo_bc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_bc_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zhi_bc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_bc_ptr,pde,assembler);
    // Bounds for advection control
    RealT lower = parlist->sublist("Problem").get("Lower Advection Bound",-100.0);
    RealT upper = parlist->sublist("Problem").get("Upper Advection Bound", 100.0);
    ROL::Ptr<std::vector<RealT> > zlo_param_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim,lower);
    ROL::Ptr<std::vector<RealT> > zhi_param_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim,upper);
    ROL::Ptr<ROL::StdVector<RealT> > zlo_adv
      = ROL::makePtr<ROL::StdVector<RealT>>(zlo_param_ptr);
    ROL::Ptr<ROL::StdVector<RealT> > zhi_adv
      = ROL::makePtr<ROL::StdVector<RealT>>(zhi_param_ptr);
    // Combined bounds
    ROL::Ptr<ROL::Vector<RealT> > zlop
      = ROL::makePtr<PDE_OptVector<RealT>>(zlo_bc, zlo_adv);
    ROL::Ptr<ROL::Vector<RealT> > zhip
      = ROL::makePtr<PDE_OptVector<RealT>>(zhi_bc, zhi_adv);
    ROL::Ptr<ROL::BoundConstraint<RealT> > zbnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    ROL::Ptr<Tpetra::MultiVector<> > ulo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::Vector<RealT> > ulop
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ulo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > uhip
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uhi_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ubnd
      = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT> >(ubnd,zbnd);

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    bool useReducedSpace = parlist->sublist("Problem").get("Use Reduced Space",true);
    if (useReducedSpace) {
      problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(objReduced,zp,zbnd);
    }
    else {
      problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(obj,xp,bnd,con,pp);
    }
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      problem->check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::OptimizationSolver<RealT> solver(*problem,*parlist);
    (*zp_ptr)[0] = parlist->sublist("Problem").get("Advection Magnitude",0.0);
    u_ptr->putScalar(450.0);

    bool solveFS = parlist->sublist("Problem").get("Initial Solve for Full Space",true); 
    if (solveFS) {
      con->solve(*rp,*up,*zp,tol);
    }

    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    pdecon->outputTpetraVector(zbc_ptr,"control.txt");
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    *outStream << std::endl << "Advection value: " << (*zp_ptr)[0] << std::endl;
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
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
