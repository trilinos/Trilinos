// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the filtered Darcy porosity optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Solver.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_ConstraintFromObjective.hpp"

#include "../../../../TOOLS/meshreader.hpp"
#include "../../../../TOOLS/pdeconstraint.hpp"
#include "../../../../TOOLS/pdeobjective.hpp"
#include "../../../../TOOLS/pdevector.hpp"
#include "../../../../TOOLS/integralconstraint.hpp"

#include "pde_darcy.hpp"
#include "obj_darcy.hpp"
#include "pde_filter.hpp"
#include "filtered_obj.hpp"


typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
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
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    bool output = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    output = (iprint > 0) && (myRank==0) && output;
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",output);

    ROL::Ptr<MeshManager<RealT>>            meshMgr;
    ROL::Ptr<PDE<RealT>>                    pde, pdeFilter;
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con;
    ROL::Ptr<PDE_Constraint<RealT>>         pdecon;
    ROL::Ptr<Assembler<RealT>>              assembler;
    /*** Initialize main data structure. ***/
    meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing filtered Darcy equations.
    pde = ROL::makePtr<PDE_Darcy<RealT>>(*parlist);
    con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    pdeFilter = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    // Cast the constraint and get the assembler.
    pdecon = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<>>         u_ptr, p_ptr, f_ptr, r_ptr;
    ROL::Ptr<std::vector<RealT>>            z0_ptr;
    ROL::Ptr<ROL::StdVector<RealT>>         z0p;
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> f1p;
    ROL::Ptr<ROL::Vector<RealT>>            up, pp, fp, rp, xp;
    bool useParamVar = parlist->sublist("Problem").get("Use Optimal Constant Velocity",false);
    int dim = 2;
    u_ptr = assembler->createStateVector();    u_ptr->randomize();
    p_ptr = assembler->createStateVector();    p_ptr->randomize();
    f_ptr = assembler->createControlVector();  f_ptr->randomize();
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    if (useParamVar) {
      z0_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      z0p = ROL::makePtr<ROL::StdVector<RealT>>(z0_ptr);
      f1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(f_ptr,pde,assembler,*parlist);
    }
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    if (useParamVar)
      fp = ROL::makePtr<PDE_OptVector<RealT>>(f1p,z0p,myRank);
    else
      fp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(f_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,fp);

    // Initialize quadratic objective function.
    ROL::Ptr<QoI<RealT>>                           qoi;
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>         obj;
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> robj;
    ROL::Ptr<FilteredObjective<RealT>>             fobj;
    //qoi  = ROL::makePtr<QoI_Velocity_Darcy<RealT>>(*parlist,
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getPressureFE(),
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getControlFE(),
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getPressureBdryFE(4),
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getControlBdryFE(4),
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getBdryCellLocIds(4),
    //         ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getPermeability());
    qoi  = ROL::makePtr<QoI_VelocityTracking_Darcy<RealT>>(*parlist,
             ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getPressureFE(),
             ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getControlFE(),
             ROL::staticPtrCast<PDE_Darcy<RealT>>(pde)->getPermeability());
    obj  = ROL::makePtr<PDE_Objective<RealT>>(qoi,assembler);
    robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,fp,pp,true,false);
    fobj = ROL::makePtr<FilteredObjective<RealT>>(robj,pdeFilter,meshMgr,comm,*parlist,*outStream);

    // Build density vector
    ROL::Ptr<Tpetra::MultiVector<>> z_ptr;
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> z1p;
    ROL::Ptr<ROL::Vector<RealT>> zp;
    z_ptr = fobj->getAssembler()->createControlVector();  z_ptr->randomize();
    if (useParamVar) {
      z1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pdeFilter,fobj->getAssembler(),*parlist);
      zp = ROL::makePtr<PDE_OptVector<RealT>>(z1p,z0p,myRank);
    }
    else {
      zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pdeFilter,fobj->getAssembler(),*parlist);
    }

    // Build bound constraint
    ROL::Ptr<ROL::Vector<RealT>>          lp, hp;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd;
    if (useParamVar) {
      ROL::Ptr<ROL::StdVector<RealT>>         l0p, h0p;
      ROL::Ptr<Tpetra::MultiVector<>>         l1_ptr, h1_ptr;
      ROL::Ptr<ROL::TpetraMultiVector<RealT>> l1p, h1p;
      l0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_NINF<RealT>());
      h0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_INF<RealT>());
      l1_ptr = fobj->getAssembler()->createControlVector();
      h1_ptr = fobj->getAssembler()->createControlVector();
      l1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(l1_ptr,pdeFilter,fobj->getAssembler(),*parlist);
      h1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(h1_ptr,pdeFilter,fobj->getAssembler(),*parlist);
      l1p->setScalar(0.0);
      h1p->setScalar(1.0);
      lp = ROL::makePtr<PDE_OptVector<RealT>>(l1p,l0p,myRank);
      hp = ROL::makePtr<PDE_OptVector<RealT>>(h1p,h0p,myRank);
    }
    else {
      lp = zp->clone(); lp->setScalar(0.0);
      hp = zp->clone(); hp->setScalar(1.0);
    }
    bnd = ROL::makePtr<ROL::Bounds<RealT>>(lp, hp);
    // Build optimization problem
    ROL::Ptr<ROL::Problem<RealT>> optProb;
    optProb = ROL::makePtr<ROL::Problem<RealT>>(fobj, zp);
    optProb->addBoundConstraint(bnd);
    optProb->finalize(false,true,*outStream);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      //ROL::Ptr<ROL::Vector<RealT>> rup = up->clone(); rup->randomize(-1.0,1.0);
      //ROL::Ptr<ROL::Vector<RealT>> rpp = pp->clone(); rpp->randomize(-1.0,1.0);
      //ROL::Ptr<ROL::Vector<RealT>> dup = up->clone(); dup->randomize(-1.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> rzp = zp->clone(); rzp->randomize( 0.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> dzp = zp->clone(); dzp->randomize( 0.0,1.0);
      //con->checkApplyJacobian_1(*rup,*rzp,*dup,*rup,true,*outStream);
      //con->checkApplyJacobian_2(*rup,*rzp,*dzp,*rup,true,*outStream);
      //con->checkInverseJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      //con->checkInverseAdjointJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      //con->checkApplyAdjointHessian_11(*rup,*rzp,*rpp,*dup,*rup,true,*outStream);
      //con->checkApplyAdjointHessian_21(*rup,*rzp,*rpp,*dzp,*rup,true,*outStream);
      //con->checkApplyAdjointHessian_12(*rup,*rzp,*rpp,*dup,*rzp,true,*outStream);
      //con->checkApplyAdjointHessian_22(*rup,*rzp,*rpp,*dzp,*rzp,true,*outStream);
      //obj->checkGradient_1(*rup,*rzp,*dup,true,*outStream);
      //obj->checkGradient_2(*rup,*rzp,*dzp,true,*outStream);
      //obj->checkHessVec_11(*rup,*rzp,*dup,true,*outStream);
      //obj->checkHessVec_12(*rup,*rzp,*dzp,true,*outStream);
      //obj->checkHessVec_21(*rup,*rzp,*dup,true,*outStream);
      //obj->checkHessVec_22(*rup,*rzp,*dzp,true,*outStream);
      fobj->checkGradient(*rzp,*dzp,true,*outStream);
      fobj->checkHessVec(*rzp,*dzp,true,*outStream);
      //optProb->check(*outStream);
    }

    // Solve optimization problem
    zp->setScalar(0.5);
    up->zero(); pp->zero();
    bool opt = parlist->sublist("Problem").get("Solve Optimization Problem",true);
    if (opt) {
      std::ifstream infile("control.txt");
      if (infile.good())
        fobj->getAssembler()->inputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ifstream infile0; infile0.open("target.txt");
        if (infile0.good()) {
          for (int i = 0; i < dim; ++i) infile0 >> (*z0_ptr)[i];
          infile0.close();
        }
      }
      ROL::Solver<RealT> optSolver(optProb, *parlist);
      optSolver.solve(*outStream);
      pdecon->outputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ofstream outfile; outfile.open("target.txt");
        for (const auto x : *z0_ptr) {
          outfile << std::scientific << std::setprecision(16);
          outfile << x << std::endl;
        }
        outfile.close();
      }
    }

    // Output
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    fobj->applyFilter(*fp,*zp,false);
    con->solve(*rp,*up,*fp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(f_ptr,"filtered_control.txt");
    con->value(*rp,*up,*fp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    assembler->printDataPDE(pde,u_ptr,f_ptr);
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //pdecon->outputTpetraData();

    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> conFilter;
    ROL::Ptr<PDE_Constraint<RealT>>         pdeconFilter;
    ROL::Ptr<Assembler<RealT>>              assemblerFilter;
    conFilter = ROL::makePtr<PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
    pdeconFilter = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(conFilter);
    assemblerFilter = pdeconFilter->getAssembler();
    ROL::Ptr<Tpetra::MultiVector<>> filteru_ptr, filterz_ptr;
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> filteru, filterz;
    filteru_ptr = assemblerFilter->createStateVector();
    filterz_ptr = assemblerFilter->createControlVector();
    filteru = ROL::makePtr<PDE_PrimalSimVector<RealT>>(filteru_ptr,pdeFilter,assemblerFilter,*parlist);
    filterz = ROL::makePtr<PDE_PrimalOptVector<RealT>>(filterz_ptr,pdeFilter,assemblerFilter,*parlist);
    filteru->setScalar(0.0);
    filterz->setScalar(0.0);
    pdeconFilter->solve(*filteru, *filteru, *filterz, tol);
    pdeconFilter->outputTpetraData();

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
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
