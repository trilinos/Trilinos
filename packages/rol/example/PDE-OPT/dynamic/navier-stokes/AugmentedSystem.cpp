// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_Bounds.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_PinTConstraint.hpp"
#include "ROL_PinTVectorCommunication.hpp"
//#include "ROL_PinTVectorCommunication_Tpetra.hpp"
#include "ROL_PinTVectorCommunication_StdTpetraComposite.hpp"

#include "../../TOOLS/dynconstraint.hpp"
#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "dynpde_navier-stokes.hpp"
#include "obj_navier-stokes.hpp"
#include "initial_condition.hpp"

using RealT = double;

class KKTOperator : public ROL::LinearOperator<RealT> {
public:
  ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;
  int myRank;
  ROL::Ptr<std::ostream> outStream;

  void apply( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    pint_con->applyAugmentedKKT(Hv,v,*state,*control,tol);
  }
};

class MGRITKKTOperator : public ROL::LinearOperator<RealT> {
public:
  ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;
  int myRank;
  ROL::Ptr<std::ostream> outStream;

  void apply( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    assert(false);
    pint_con->applyMultigridAugmentedKKT(Hv,v,*state,*control,tol);
  }
  void applyInverse( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    pint_con->applyMultigridAugmentedKKT(Hv,v,*state,*control,tol);
  }
};

class WathenKKTOperator : public ROL::LinearOperator<RealT> {
public:
  ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;
  int myRank;
  ROL::Ptr<std::ostream> outStream;

  void apply( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    assert(false);
  }
  void applyInverse( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    pint_con->applyWathenInverse(Hv,v,*state,*control,tol,false,0);
    // pint_con->applyLocalInverse(Hv,v,*state,*control,tol,0,false); // do this to exactly solve
  }
};

ROL::Ptr<ROL::Vector<RealT>> 
solveKKTSystem(const std::string & prefix,
               bool useWathenPrec,int myRank,
               double absTol,double relTol,
               const Teuchos::RCP<Teuchos::StackedTimer> & timer,
               const ROL::Ptr<std::ostream> & outStream,
               const ROL::Ptr<ROL::PinTConstraint<RealT>> & pint_con,
               const ROL::Ptr<ROL::Vector<RealT>> & kkt_b,
               const ROL::Ptr<ROL::PinTVector<RealT>> & state,
               const ROL::Ptr<ROL::PinTVector<RealT>> & control);


template<class Real>
void computeInitialCondition(const ROL::Ptr<ROL::Vector<Real>>       &u0,
                             const ROL::Ptr<ROL::Vector<Real>>       &ck,
                             const ROL::Ptr<ROL::Vector<Real>>       &uo,
                             const ROL::Ptr<ROL::Vector<Real>>       &un,
                             const ROL::Ptr<ROL::Vector<Real>>       &zk,
                             const ROL::Ptr<DynConstraint<Real>>     &con,
                             const Real                               dt,
                             std::ostream                            &outStream);

int main(int argc, char *argv[]) {

  using Teuchos::RCP;

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
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input.xml");
    int nt           = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
    RealT T          = parlist->sublist("Time Discretization").get("End Time",             1.0);
    RealT dt         = T/static_cast<RealT>(nt);
    bool useParametricControl = parlist->sublist("Problem").get("Use Parametric Control", false);
    int verbosity    = parlist->sublist("General").get("Print Verbosity", 0);
    verbosity        = (myRank==0 ? verbosity : 0);
    parlist->sublist("General").set("Print Verbosity", verbosity);
    bool solveOutput = parlist->sublist("Dynamic Constraint").sublist("Solve").get("Output Iteration History", false);
    solveOutput      = (myRank==0 ? solveOutput : false);
    parlist->sublist("Dynamic Constraint").sublist("Solve").set("Output Iteration History", solveOutput);
    solveOutput      = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History", false);
    solveOutput      = (myRank==0 ? solveOutput : false);
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History", solveOutput);

    /*** Add MGRIT parameter list options. ***/
    int sweeps         = parlist->get("MGRIT Sweeps", 1);
    RealT omega        = parlist->get("MGRIT Relaxation",2.0/3.0);
    int coarseSweeps   = parlist->get("MGRIT Coarse Sweeps", 1);
    RealT coarseOmega  = parlist->get("MGRIT Coarse Relaxation",1.0);
    int numLevels      = parlist->get("MGRIT Levels",3);
    double relTol      = parlist->get("MGRIT Krylov Relative Tolerance",1e-4);
    double absTol      = parlist->get("MGRIT Krylov Absolute Tolerance",1e-4);
    int spaceProc      = parlist->get("MGRIT Spatial Procs",1);
    double globalScale = parlist->get("MGRIT Global Scale",1.0);
    bool rebalance     = parlist->get("MGRIT Rebalance",false);
    int cgIterations   = parlist->get("MGRIT CG Iterations",1);
    double cntrlRegPar = parlist->get("MGRIT Control Regularization Parameter",1.0);

    ROL::Ptr<const ROL::PinTCommunicators> communicators           = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spaceProc);
    ROL::Ptr<const Teuchos::Comm<int>> mpiSpaceComm = ROL::makePtr<Teuchos::MpiComm<int>>(communicators->getSpaceCommunicator());

    // for "serial" in time cases, use the wathen preconditioner
    bool useWathenPrec = communicators->getTimeSize()==spaceProc;

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize mesh data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      //= ROL::makePtr<MeshReader<RealT>>(*parlist, numProcs);
      = ROL::makePtr<MeshReader<RealT>>(*parlist, 0);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<DynamicPDE_NavierStokes<RealT>> pde
      = ROL::makePtr<DynamicPDE_NavierStokes<RealT>>(*parlist);

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    ROL::Ptr<DynConstraint<RealT>> dyn_con
      = ROL::makePtr<DynConstraint<RealT>>(pde,meshMgr,mpiSpaceComm,*parlist,*outStream);
    const ROL::Ptr<Assembler<RealT>> assembler = dyn_con->getAssembler();
    dyn_con->setSolveParameters(*parlist);
    dyn_con->getAssembler()->printMeshData(*outStream);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u0_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> uo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> un_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> ck_ptr = assembler->createResidualVector();
    ROL::Ptr<ROL::Vector<RealT>> u0, uo, un, ck, zk;
    u0 = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u0_ptr,pde,*assembler,*parlist);
    uo = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uo_ptr,pde,*assembler,*parlist);
    un = ROL::makePtr<PDE_PrimalSimVector<RealT>>(un_ptr,pde,*assembler,*parlist);
    ck = ROL::makePtr<PDE_DualSimVector<RealT>>(ck_ptr,pde,*assembler,*parlist);
    if (!useParametricControl) {
      ROL::Ptr<Tpetra::MultiVector<>> zk_ptr = assembler->createControlVector();
      zk = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zk_ptr,pde,*assembler,*parlist);
    }
    else {
      zk = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(1));
    }
    ROL::Ptr<ROL::PartitionedVector<RealT>> z
      = ROL::PartitionedVector<RealT>::create(*zk, nt);

    ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm;
    if (!useParametricControl) {
      vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_Tpetra<RealT>>();
    }
    else {
      vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdTpetraComposite<RealT>>();
    }
    ROL::Ptr<ROL::PinTVector<RealT>> state   = ROL::buildStatePinTVector<RealT>(   communicators, vectorComm, nt, u0); // for Euler, Crank-Nicolson, stencil = [-1,0]
    ROL::Ptr<ROL::PinTVector<RealT>> control = ROL::buildControlPinTVector<RealT>( communicators, vectorComm, nt, zk); // time discontinous, stencil = [0]

    // make sure we are globally consistent
    state->boundaryExchange();
    control->boundaryExchange();

    ROL::Ptr<ROL::Vector<RealT>> kkt_vector = ROL::makePtr<ROL::PartitionedVector<RealT>>({state->clone(),control->clone(),state->clone()});

    /*************************************************************************/
    /************************* Set initial condition *************************/
    /*************************************************************************/
    // ??? yanked code, reinsert

    /*************************************************************************/
    /***************** Set initial guess for optimization ********************/
    /*************************************************************************/
    // ??? yanked code, reinsert
    state->zero();
    control->zero();

    /*************************************************************************/
    /***************** Build PINT Constraint *********************************/
    /*************************************************************************/
    auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(control->numOwnedSteps());
    for( uint k=0; k<timeStamp->size(); ++k ) {
      timeStamp->at(k).t.resize(2);
      timeStamp->at(k).t.at(0) = k*dt;
      timeStamp->at(k).t.at(1) = (k+1)*dt;
    }

    if(myRank==0) {
      (*outStream) << "Sweeps = " << sweeps       << std::endl;
      (*outStream) << "Omega = "  << omega        << std::endl;
      (*outStream) << "Levels = " << numLevels    << std::endl;
      (*outStream) << "Coarse Sweeps = " << coarseSweeps << std::endl;
      (*outStream) << "Coarse Omega = "  << coarseOmega  << std::endl;
      (*outStream) << "Global Scale  = " << globalScale  << std::endl;
      (*outStream) << "Rebalance     = " << rebalance    << std::endl;
      (*outStream) << "CG Iterations = " << cgIterations;
      if(cgIterations==0)
        *(outStream) << " (apply Wathen smoothing)" << std::endl;
      else
        *(outStream) << " (apply CG fix to Wathen smoothing)" << std::endl;
      (*outStream) << "Cntrl Reg Par = " << cntrlRegPar << std::endl;
    }

    // build the parallel in time constraint from the user constraint
    ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con = ROL::makePtr<ROL::PinTConstraint<RealT>>(dyn_con,u0,timeStamp);
    pint_con->applyMultigrid(numLevels,communicators,vectorComm);
    pint_con->setSweeps(sweeps);
    pint_con->setRelaxation(omega);
    pint_con->setCoarseSweeps(coarseSweeps);
    pint_con->setCoarseRelaxation(coarseOmega);
    pint_con->setGlobalScale(globalScale);
    pint_con->setCGIterations(cgIterations);
    pint_con->setControlRegParam(cntrlRegPar);
    pint_con->setRecordResidualReductions(true);

    /*************************************************************************/
    /***************** Run KKT Solver ***************************************/
    /*************************************************************************/

    double tol = 1e-12;
    auto kkt_x_in  = kkt_vector->clone();
    auto kkt_b     = kkt_vector->clone();

    ROL::RandomizeVector(*kkt_x_in);
    kkt_b->zero();

    if(myRank==0) {
      (*outStream) << "Applying augmented KKT system" << std::endl;
    }

    pint_con->applyAugmentedKKT(*kkt_b,*kkt_x_in,*state,*control,tol); // b = A * x

    if(myRank==0) {
      (*outStream) << "Applying augmented KKT system - complete" << std::endl;
      (*outStream) << std::endl;
    }

    auto timer = Teuchos::TimeMonitor::getStackedTimer();

  
    {
      RealT res0 = kkt_b->norm();
      kkt_b->scale(1.0/res0);
      res0 = kkt_b->norm();

      if(myRank==0) (*outStream) << std::endl;

      solveKKTSystem("SOLVE: ",
          useWathenPrec,myRank,absTol,relTol,
          timer,outStream,
          pint_con,kkt_b,state,control);

      if(myRank==0) (*outStream) << std::endl;
    }

    #if 0
    {
      RealT res0 = kkt_b->norm();
      kkt_b->scale(1.0/res0);
      res0 = kkt_b->norm();

      assert(std::fabs(res0-1.0) <= 1.0e-14);

      auto kkt_x_out = kkt_vector->clone();
      kkt_x_out->zero();

      KKTOperator kktOperator;
      kktOperator.pint_con = pint_con;
      kktOperator.state = state;
      kktOperator.control = control;
      kktOperator.outStream = outStream;
      kktOperator.myRank = myRank;

      RCP<ROL::LinearOperator<RealT>> precOperator;

      // build the preconditioner
      if(useWathenPrec) {
        RCP<WathenKKTOperator> wathenOperator(new WathenKKTOperator);
        wathenOperator->pint_con = pint_con;
        wathenOperator->state = state;
        wathenOperator->control = control;
        wathenOperator->outStream = outStream;
        wathenOperator->myRank = myRank;

        precOperator = wathenOperator;
      }
      else {
        RCP<MGRITKKTOperator> mgOperator(new MGRITKKTOperator);
        mgOperator->pint_con = pint_con;
        mgOperator->state = state;
        mgOperator->control = control;
        mgOperator->outStream = outStream;
        mgOperator->myRank = myRank;

        precOperator = mgOperator;
      }

      Teuchos::ParameterList parlist;
      ROL::ParameterList &krylovList = parlist.sublist("General").sublist("Krylov");
      krylovList.set("Absolute Tolerance",absTol);
      krylovList.set("Relative Tolerance",relTol);

      if(myRank==0)
        (*outStream) << "RELATIVE TOLERANCE = " << relTol << std::endl;

      ROL::GMRES<RealT> krylov(parlist); // TODO: Do Belos
      krylov.enableOutput(*outStream);

      int flag = 0, iter = 0;

      MPI_Barrier(MPI_COMM_WORLD);
      double t0 = MPI_Wtime();

      if(myRank==0) {
        (*outStream) << "Solving KKT system: ";
        if(useWathenPrec)
          (*outStream) << "using Wathen preconditioner" << std::endl;
        else
          (*outStream) << "using MGRIT preconditioner" << std::endl;
      }

      timer->start("krylov");
      RealT finalTol = krylov.run(*kkt_x_out,kktOperator,*kkt_b,*precOperator,iter,flag);
      timer->stop("krylov");

      if(myRank==0)
        (*outStream) << "Solving KKT system - complete" << std::endl;

      MPI_Barrier(MPI_COMM_WORLD);
      double tf = MPI_Wtime();

      std::stringstream ss;
      timer->report(ss);

      if(myRank==0) {
        (*outStream) << "Krylov Iteration = " << iter << " " << (finalTol / res0) << " " << tf-t0 << std::endl;
        (*outStream) << std::endl;
        (*outStream) << ss.str();
      }
    }
    #endif

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

template<class Real>
void computeInitialCondition(const ROL::Ptr<ROL::Vector<Real>>       &u0,
                             const ROL::Ptr<ROL::Vector<Real>>       &ck,
                             const ROL::Ptr<ROL::Vector<Real>>       &uo,
                             const ROL::Ptr<ROL::Vector<Real>>       &un,
                             const ROL::Ptr<ROL::Vector<Real>>       &zk,
                             const ROL::Ptr<DynConstraint<Real>>     &con,
                             const Real                               dt,
                             std::ostream                            &outStream) {
  Real T  = 80.0;
  int  nt = static_cast<int>(T/dt);
  std::vector<ROL::TimeStamp<Real>> ts(nt);
  for( int k=0; k<nt; ++k ) {
    ts.at(k).t.resize(2);
    ts.at(k).t.at(0) = k*dt;
    ts.at(k).t.at(1) = (k+1)*dt;
  }
  // Solve Navier-Stokes equation to determine initial condition
  zk->zero(); uo->set(*u0); un->zero();
  Real unorm = uo->norm();
  outStream << std::scientific << std::setprecision(6);
  outStream << std::right << std::setw(8)  << "ts"
            << std::right << std::setw(16) << "||u(ts)||"
            << std::right << std::setw(16) << "avg time (sec)"
            << std::endl;
  outStream << std::right << std::setw(8)  << 0
            << std::right << std::setw(16) << unorm
            << std::right << std::setw(16) << "---"
            << std::endl;
  std::vector<Real> time(10);
  std::clock_t timer_step;
  Real time_avg(0);
  for (int k = 1; k < nt; ++k) {
    // Advance time stepper
    timer_step = std::clock();
    con->solve(*ck, *uo, *un, *zk, ts[k]);
    time[k%10] = static_cast<Real>(std::clock()-timer_step)/static_cast<Real>(CLOCKS_PER_SEC);
    uo->set(*un);
    if ( k%10==0 ) {
      unorm = uo->norm();
      time_avg = 0.0;
      for (int i = 0; i < 10; ++i) {
        time_avg += time[i];
      }
      time_avg *= 0.1;
      outStream << std::right << std::setw(8)  << k
                << std::right << std::setw(16) << unorm
                << std::right << std::setw(16) << time_avg
                << std::endl;
    }
  }
  u0->set(*uo);
}

ROL::Ptr<ROL::Vector<RealT>> 
solveKKTSystem(const std::string & prefix,
               bool useWathenPrec,int myRank,
               double absTol,double relTol,
               const Teuchos::RCP<Teuchos::StackedTimer> & timer,
               const ROL::Ptr<std::ostream> & outStream,
               const ROL::Ptr<ROL::PinTConstraint<RealT>> & pint_con,
               const ROL::Ptr<ROL::Vector<RealT>> & kkt_b,
               const ROL::Ptr<ROL::PinTVector<RealT>> & state,
               const ROL::Ptr<ROL::PinTVector<RealT>> & control)
{
  using Teuchos::RCP;

  RealT res0 = kkt_b->norm();
  res0 = kkt_b->norm();

  auto kkt_x_out = kkt_b->clone();
  kkt_x_out->zero();

  KKTOperator kktOperator;
  kktOperator.pint_con = pint_con;
  kktOperator.state = state;
  kktOperator.control = control;
  kktOperator.outStream = outStream;
  kktOperator.myRank = myRank;

  RCP<ROL::LinearOperator<RealT>> precOperator;

  // build the preconditioner 
  if(useWathenPrec) {
    RCP<WathenKKTOperator> wathenOperator(new WathenKKTOperator);
    wathenOperator->pint_con = pint_con;
    wathenOperator->state = state;
    wathenOperator->control = control;
    wathenOperator->outStream = outStream;
    wathenOperator->myRank = myRank;

    precOperator = wathenOperator;
  }
  else {
    RCP<MGRITKKTOperator> mgOperator(new MGRITKKTOperator);
    mgOperator->pint_con = pint_con;
    mgOperator->state = state;
    mgOperator->control = control;
    mgOperator->outStream = outStream;
    mgOperator->myRank = myRank;

    precOperator = mgOperator;
  }

  Teuchos::ParameterList parlist;
  ROL::ParameterList &krylovList = parlist.sublist("General").sublist("Krylov");
  krylovList.set("Absolute Tolerance",absTol);
  krylovList.set("Relative Tolerance",relTol);

  if(myRank==0)
    (*outStream) << "Relative tolerance = " << relTol << std::endl;

  ROL::GMRES<RealT> krylov(parlist); // TODO: Do Belos
  // ROL::MINRES<RealT> krylov(1e0, 1e-6, 200); // TODO: Do Belos
 
  krylov.enableOutput(*outStream);

  int flag = 0;
  int iter = 0;


  if(myRank==0) {
    (*outStream) << "Solving KKT system: ";
    if(useWathenPrec)
      (*outStream) << "using Wathen preconditioner" << std::endl;
    else
      (*outStream) << "using MGRIT preconditioner" << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();

  timer->start("krylov");
  RealT finalTol = krylov.run(*kkt_x_out,kktOperator,*kkt_b,*precOperator,iter,flag);
  timer->stop("krylov");

  MPI_Barrier(MPI_COMM_WORLD);
  double tf = MPI_Wtime();

  // Don't trust residual computation, so repeat.
  auto myb = kkt_b->clone();
  myb->set(*kkt_b);
  auto norm_myb = myb->norm();
  auto tmp = kkt_b->clone();
  RealT dummytol = 1e-8;
  kktOperator.apply(*tmp, *kkt_x_out, dummytol);
  tmp->axpy(-1, *myb);
  auto norm_myAxmb = tmp->norm();
  auto norm_mysoln = kkt_x_out->norm();

  if(myRank==0) {
     (*outStream) << prefix << "Krylov Iteration = " << iter << " " << (finalTol / res0) << " " << tf-t0 << std::endl;
     (*outStream) << "||x||=" << norm_mysoln << "  ||b||=" << norm_myb << "  ||Ax-b||=" << norm_myAxmb << std::endl;

     const std::map<int,std::vector<double>> & preSmooth = pint_con->getPreSmoothResidualReductions();
     const std::map<int,std::vector<double>> & postSmooth = pint_con->getPostSmoothResidualReduction();
     const std::vector<double> & coarseSmooth = pint_con->getCoarseResidualReduction();

     if(not useWathenPrec) {
       // loop over presmoother residual reductions
       *outStream << "Presmooth Reductions: " << std::endl;
       for(auto iter_l : preSmooth) {
         double average = std::accumulate( iter_l.second.begin(), iter_l.second.end(), 0.0)/coarseSmooth.size(); 
         double min = *std::min_element(iter_l.second.begin(),iter_l.second.end());
         double max = *std::max_element(iter_l.second.begin(),iter_l.second.end());
         *outStream << "  level ";
         *outStream << "min,avg,max = " << min << ", " << average << ", " << max << std::endl;
       }

       // loop over postsmoother residual reductions
       *outStream << "Postsmooth Reductions: " << std::endl;
       for(auto iter_l : postSmooth) {
         double average = std::accumulate( iter_l.second.begin(), iter_l.second.end(), 0.0)/coarseSmooth.size(); 
         double min = *std::min_element(iter_l.second.begin(),iter_l.second.end());
         double max = *std::max_element(iter_l.second.begin(),iter_l.second.end());
         *outStream << "  level ";
         *outStream << "min,avg,max = " << min << ", " << average << ", " << max << std::endl;
       }

       // loop over coarse residual reductions
       {
         *outStream << "Coarse Reductions: ";
         double average = std::accumulate( coarseSmooth.begin(), coarseSmooth.end(), 0.0)/coarseSmooth.size(); 
         double min = *std::min_element(coarseSmooth.begin(),coarseSmooth.end());
         double max = *std::max_element(coarseSmooth.begin(),coarseSmooth.end());
         *outStream << "min,avg,max = " << min << ", " << average << ", " << max << std::endl;
         *outStream << std::endl;
       }
     }
  }

  pint_con->clearResidualReduction();

  return kkt_x_out;
}


