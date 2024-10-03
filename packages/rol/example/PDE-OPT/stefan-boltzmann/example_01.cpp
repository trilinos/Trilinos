// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the stochastic Stefan-Boltzmann problem.
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

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Solver.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_stoch_stefan_boltzmann.hpp"
#include "obj_stoch_stefan_boltzmann.hpp"
#include "mesh_stoch_stefan_boltzmann.hpp"

typedef double RealT;

template<class Real>
Real random(const Teuchos::Comm<int> &comm,
            const Real a = -1, const Real b = 1) {
  Real val(0), u(0);
  if ( Teuchos::rank<int>(comm)==0 ) {
    u   = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    val = (b-a)*u + a;
  }
  Teuchos::broadcast<int,Real>(comm,0,1,&val);
  return val;
}

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
  ROL::Ptr<const Teuchos::Comm<int> > serial_comm
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
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    const int stochDim = 32;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT>>(*parlist);
    //  = ROL::makePtr<StochasticStefanBoltzmannMeshManager<RealT>>(*parlist);
    // Initialize Stefan-Boltzmann PDE.
    ROL::Ptr<StochasticStefanBoltzmannPDE<RealT> > pde
      = ROL::makePtr<StochasticStefanBoltzmannPDE<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    con->setSolveParameters(*parlist);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdeCon->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > pp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler);
    // Create residual vectors
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->randomize(); //r_ptr->putScalar(static_cast<RealT>(1));
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    // Create control vector and set to ones
    ROL::Ptr<Tpetra::MultiVector<> >  z_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > yz_ptr = assembler->createControlVector();
    z_ptr->randomize();  z_ptr->putScalar(static_cast<RealT>(280));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(0));
    yz_ptr->randomize(); //yz_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > dzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > yzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(yz_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtr<PDE_OptVector<RealT>>(zpde);
    ROL::Ptr<ROL::Vector<RealT> > dzp = ROL::makePtr<PDE_OptVector<RealT>>(dzpde);
    ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtr<PDE_OptVector<RealT>>(yzpde);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateCost<RealT>>(pde->getVolFE(),*parlist);
    qoi_vec[1] = ROL::makePtr<QoI_ControlCost<RealT>>(
      pde->getVolFE(),pde->getBdryFE(0),pde->getBdryCellLocIds(0),*parlist);
    ROL::Ptr<StochasticStefanBoltzmannStdObjective<RealT> > std_obj
      = ROL::makePtr<StochasticStefanBoltzmannStdObjective<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> >  zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(static_cast<RealT>(280));
    zhi_ptr->putScalar(static_cast<RealT>(370));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zlopde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zhipde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > zlop = ROL::makePtr<PDE_OptVector<RealT>>(zlopde);
    ROL::Ptr<ROL::Vector<RealT> > zhip = ROL::makePtr<PDE_OptVector<RealT>>(zhipde);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT> > bounds(stochDim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic", one);
    opt.setStochasticObjective(*parlist,sampler);
    ROL::Ptr<ROL::Problem<RealT>> newopt
      = ROL::makePtr<ROL::Problem<RealT>>(opt.getObjective(),opt.getSolutionVector());
    newopt->addBoundConstraint(opt.getBoundConstraint());

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      up->checkVector(*pp,*dup,true,*outStream);
      zp->checkVector(*yzp,*dzp,true,*outStream);
      std::vector<RealT> param(stochDim,0);
      objReduced->setParameter(param);
      obj->checkGradient(x,d,true,*outStream);
      obj->checkHessVec(x,d,true,*outStream);
      pdeCon->setParameter(param);
      pdeCon->checkApplyJacobian(x,d,*up,true,*outStream);
      pdeCon->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      pdeCon->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      pdeCon->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      pdeCon->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);
      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::Solver<RealT> solver(newopt,*parlist);
    //zp->zero();
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
    pdeCon->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    ROL::Ptr<ROL::BatchManager<RealT> > bman_Eu
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    std::vector<RealT> sample(stochDim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    file_samp << std::scientific << std::setprecision(15);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*rp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << std::setw(25) << std::left << sample[j];
      }
      file_samp << std::endl;
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    pdeCon->outputTpetraVector(p_ptr,"mean_state.txt");
    // Build objective function distribution
    RealT val1(0), val2(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > stateCost
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > redStateCost
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(stateCost, con, up, zp, pp, true, false);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > ctrlCost
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > redCtrlCost
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(ctrlCost, con, up, zp, pp, true, false);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      redStateCost->setParameter(sample);
      val1 = redStateCost->value(*zp,tol);
      redCtrlCost->setParameter(sample);
      val2 = redCtrlCost->value(*zp,tol);
      for (int j = 0; j < stochDim; ++j) {
        file << std::setw(25) << std::left << sample[j];
      }
      file << std::setw(25) << std::left << val1;
      file << std::setw(25) << std::left << val2 << std::endl;
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    pdeCon->solve(*rp,*up,*zp,tol);
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
