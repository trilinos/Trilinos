// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve a semilinear parabolic problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>

#include <sys/types.h>
#include <unistd.h>

#include "ROL_Bounds.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"
#include "ROL_PinTConstraint.hpp"
#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_PinTVectorCommunication_Tpetra.hpp"
#include "ROL_PinTVectorCommunication_StdTpetraComposite.hpp"

#include "../../TOOLS/dynconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "dynpde_semilinear.hpp"
#include "obj_semilinear.hpp"

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
  }
};

inline void pauseToAttach(MPI_Comm mpicomm)
{
  Teuchos::RCP<Teuchos::Comm<int> > comm =
    Teuchos::createMpiComm<int>(Teuchos::rcp(new Teuchos::OpaqueWrapper<MPI_Comm>(mpicomm)));
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setShowProcRank(true);
  out.setOutputToRootOnly(-1);

  comm->barrier();

  // try to get them to print out all at once
  out << "PID = " << getpid() << std::endl;

  if (comm->getRank() == 0)
    getchar();
  comm->barrier();
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
               const ROL::Ptr<ROL::PinTVector<RealT>> & control);

int main(int argc, char *argv[]) 
{

  using Teuchos::RCP;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // pauseToAttach(MPI_COMM_WORLD);

  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input_ex01.xml");
    int nt         = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
    RealT T        = parlist->sublist("Time Discretization").get("End Time",             1.0);
    RealT dt       = T/static_cast<RealT>(nt);

    const Teuchos::ParameterList & geomlist = parlist->sublist("Geometry");
    RealT dx       = geomlist.get<double>("Width") / geomlist.get<int>("NX");
    RealT dy       = geomlist.get<double>("Height") / geomlist.get<int>("NY");

    // Add MGRIT parameter list stuff
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
      = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    // Initialize PDE describe semilinear equation
    ROL::Ptr<DynamicPDE_Semilinear<RealT>> pde
      = ROL::makePtr<DynamicPDE_Semilinear<RealT>>(*parlist);

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    ROL::Ptr<DynConstraint<RealT> > dyn_con
      = ROL::makePtr<DynConstraint<RealT>>(pde,meshMgr,mpiSpaceComm,*parlist,*outStream);
    const ROL::Ptr<Assembler<RealT>> assembler = dyn_con->getAssembler();
    dyn_con->setSolveParameters(*parlist);
    dyn_con->getAssembler()->printMeshData(*outStream);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u0_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> zk_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT>> u0, zk;
    u0 = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u0_ptr,pde,*assembler,*parlist);
    zk = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zk_ptr,pde,*assembler,*parlist);

    // ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_Tpetra<RealT>>();
    ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdTpetraComposite<RealT>>();
    ROL::Ptr<ROL::PinTVector<RealT>> state   = ROL::buildStatePinTVector<RealT>(   communicators, vectorComm, nt,     u0); // for Euler, Crank-Nicolson, stencil = [-1,0]
    ROL::Ptr<ROL::PinTVector<RealT>> control = ROL::buildControlPinTVector<RealT>( communicators, vectorComm, nt,     zk); // time discontinous, stencil = [0]

    // make sure we are globally consistent
    state->boundaryExchange();
    control->boundaryExchange();
  
    ROL::Ptr<ROL::Vector<RealT>> kkt_vector = ROL::makePtr<ROL::PartitionedVector<RealT>>({state->clone(),control->clone(),state->clone()});

    /*************************************************************************/
    /***************** BUILD PINT CONSTRAINT *********************************/
    /*************************************************************************/

    /*************************************************************************/
    /************************* BUILD TIME STAMP ******************************/
    /*************************************************************************/
    auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(control->numOwnedSteps());

    double myFinalTime = dt*timeStamp->size();
    double timeOffset  = 0.0;
    MPI_Exscan(&myFinalTime,&timeOffset,1,MPI_DOUBLE,MPI_SUM,communicators->getTimeCommunicator());
    
    for( uint k=0; k<timeStamp->size(); ++k ) {
      timeStamp->at(k).t.resize(2);
      timeStamp->at(k).t.at(0) = k*dt+timeOffset;
      timeStamp->at(k).t.at(1) = (k+1)*dt+timeOffset;
    }

    if(myRank==0) {
      (*outStream) << "Sweeps        = " << sweeps       << std::endl;
      (*outStream) << "Omega         = " << omega        << std::endl;
      (*outStream) << "Levels        = " << numLevels    << std::endl;
      (*outStream) << "Coarse Sweeps = " << coarseSweeps << std::endl;
      (*outStream) << "Coarse Omega  = " << coarseOmega  << std::endl;
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
    pint_con->applyMultigrid(numLevels,communicators,vectorComm,rebalance);
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

    Teuchos::RCP<Teuchos::StackedTimer> timer = Teuchos::TimeMonitor::getStackedTimer();
    double tol = 1e-12;

    state->setScalar(1.0);
    control->setScalar(1.0);

    // setup and solve the lagrange multiplier RHS of the SQP algorithm
    ///////////////////////////////////////////////////////////////////////////
    {
      auto kkt_b     = kkt_vector->clone();
      kkt_b->setScalar(1.0);
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(0)->scale(dt*dx*dy); // u
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(1)->scale(dt*dx*dy); // z
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(2)->scale(0.0);      // lambda

      if(myRank==0) (*outStream) << std::endl;

      solveKKTSystem("LAGR MULTIPLIER STEP: ",
          useWathenPrec,myRank,absTol,relTol,
          timer,outStream,
          pint_con,kkt_b,state,control);

      if(myRank==0) (*outStream) << std::endl;
    }

    // compute quasi normal step
    ///////////////////////////////////////////////////////////////////////////
    
    ROL::Ptr<ROL::Vector<RealT>> g_k; // this will be use in tangential step
    {
      ROL::Ptr<ROL::Vector<RealT>> residual = state->clone();
      ROL::Ptr<ROL::Vector<RealT>> cau_pt_1;
      ROL::Ptr<ROL::Vector<RealT>> cau_pt_2;

      // residual = c(x_k)
      pint_con->value(*residual,*state,*control,tol );

      {
        ROL::Ptr<ROL::Vector<RealT>> jv_1     = state->clone();
        ROL::Ptr<ROL::Vector<RealT>> jv_2     = state->clone();
        ROL::Ptr<ROL::Vector<RealT>> ajv_1    = state->clone();
        ROL::Ptr<ROL::Vector<RealT>> ajv_2    = control->clone();

        // ajv = c'(x_k)^* c(x_k)
        pint_con->applyAdjointJacobian_1(*ajv_1,*residual,*state,*control,tol);
        pint_con->applyAdjointJacobian_2(*ajv_2,*residual,*state,*control,tol);

        // jv = c'(x_k) c'(x_k)^* c(x_k)
        pint_con->applyJacobian_1(*jv_1,*ajv_1,*state,*control,tol);
        pint_con->applyJacobian_2(*jv_2,*ajv_2,*state,*control,tol);

        RealT norm_ajv_sqr = std::pow(ajv_1->norm(),2.0) + std::pow(ajv_2->norm(),2.0);
        RealT norm_jv_sqr = std::pow(jv_1->norm(),2.0) + std::pow(jv_2->norm(),2.0);

        // now we we will modify the jacobian components
        jv_1->scale(norm_ajv_sqr/norm_jv_sqr);
        jv_2->scale(norm_ajv_sqr/norm_jv_sqr);
        residual->plus(*jv_1);
        residual->plus(*jv_2);

        jv_1 = ROL::nullPtr;
        jv_2 = ROL::nullPtr;

        // cauchy point
        cau_pt_1 = ajv_1; cau_pt_1->scale(norm_ajv_sqr/norm_jv_sqr);
        cau_pt_2 = ajv_2; cau_pt_2->scale(norm_ajv_sqr/norm_jv_sqr);
      }

      // build RHS vector
      auto kkt_b     = kkt_vector->clone();
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(0)->set(*cau_pt_1); // u
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(1)->set(*cau_pt_2); // z
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(2)->set(*residual); // lambda

      residual = ROL::nullPtr;

      if(myRank==0) (*outStream) << std::endl;

      auto delta_nk = solveKKTSystem("QUASI-NORMAL STEP: ",
                      useWathenPrec,myRank,absTol,relTol,
                      timer,outStream,
                      pint_con,kkt_b,state,control);

      if(myRank==0) (*outStream) << std::endl;

      kkt_b->setScalar(1.0);
      kkt_b->plus(*delta_nk);

      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(0)->scale(dt*dx*dy); // u
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(1)->scale(dt*dx*dy); // z
      ROL::dynamicPtrCast<ROL::PartitionedVector<double>>(kkt_b)->get(2)->setScalar(0.0); // lambda

      g_k = kkt_b;
    }

    // compute tangential step step
    ///////////////////////////////////////////////////////////////////////////
    {
      if(myRank==0) (*outStream) << std::endl;

      solveKKTSystem("TANG STEP: ",
          useWathenPrec,myRank,absTol,relTol,
          timer,outStream,
          pint_con,g_k,state,control);

      if(myRank==0) (*outStream) << std::endl;
    }

    // setup and solve the quasi normal step

    #if 0
    {
      // This was useful for debugging, maintaing it so that we can recreate it if required.
      // However, its not consistent right now so I don't want to mislead.
 
      ROL::Ptr<ROL::Vector<RealT>> scratch;
      {
        /*************************************************************************/
        /***************** BUILD SERIAL-IN-TIME CONSTRAINT ***********************/
        /*************************************************************************/
        auto serial_con = ROL::make_SerialConstraint(dyn_con, *u0, timeStamp);

        /*************************************************************************/
        /******************* BUILD SERIAL-IN-TIME VECTORS ************************/
        /*************************************************************************/
        auto U  = ROL::PartitionedVector<RealT>::create(*uo, nt);
        auto Z  = ROL::PartitionedVector<RealT>::create(*zk, nt);
        auto C  = ROL::PartitionedVector<RealT>::create(*ck, nt);

        scratch = U->clone(); // this will be used for testing in the PINT vector

        auto UZ    = ROL::make_Vector_SimOpt(U, Z);
        auto X     = UZ->clone();
        auto V_x   = X->clone();
        auto V_c   = C->clone();
        auto V_l   = C->dual().clone();
        auto ajv   = X->dual().clone();
        auto ajvpx = X->dual().clone();
        auto jv    = C->clone();

        /*************************************************************************/
        /******************* APPLY AUGMENTED SYSTEM TO 1's ***********************/
        /*************************************************************************/
        X->zero();
        ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(V_x)->get_1()->setScalar(1.0);
        ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(V_x)->get_2()->setScalar(2.0);
        V_l->setScalar(3.0);
        ajvpx->set(V_x->dual());
        serial_con->applyAdjointJacobian(*ajv, *V_l, *X, dummytol);
        ajvpx->plus(*ajv);

        auto ajv_u = ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(ajv)->get_1();
        auto ajv_z = ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(ajv)->get_2();

        serial_con->applyJacobian(*jv, *V_x, *X, dummytol);

        (*outStream) << "\nSerial components of augmented system:\n";
        (*outStream) << "  Norm of Aug*1 1st component = " << ajvpx->norm() << "\n";
        (*outStream) << "  Norm of Aug*1 2nd component = " << jv->norm() << "\n";

        auto Hv_u = ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(ajvpx)->get_1();
        auto Hv_z = ROL::dynamicPtrCast<ROL::Vector_SimOpt<RealT>>(ajvpx)->get_2();
        auto Hv_l = jv;

        (*outStream) << "  Norm of Aug*1_u = " << Hv_u->norm() << "\n";
        (*outStream) << "  Norm of Aug*1_z = " << Hv_z->norm() << "\n";
        (*outStream) << "  Norm of Aug*1_l = " << Hv_l->norm() << "\n";

      }
      {
        /*************************************************************************/
        /************************* TEST THE PINT SYSTEM  *************************/
        /*************************************************************************/

        /*************************************************************************/
        /******************* BUILD PARALLEL-IN-TIME VECTORS **********************/
        /*************************************************************************/
        auto U = ROL::buildStatePinTVector<RealT>(  communicators, vectorComm, nt, uo); // for Euler, Crank-Nicolson, stencil = [-1,0]
        auto Z = ROL::buildControlPinTVector<RealT>(communicators, vectorComm, nt, zk); // time discontinous, stencil = [0]
        auto C = ROL::buildStatePinTVector<RealT>(  communicators, vectorComm, nt, ck); // for Euler, Crank-Nicolson, stencil = [-1,0]

        auto triple = ROL::makePtr<ROL::PartitionedVector<RealT>>({U->clone(),Z->clone(),C->clone()});
        auto V      = triple->clone();
        auto Hv     = triple->dual().clone();

        /*************************************************************************/
        /******************* APPLY AUGMENTED SYSTEM TO 1's ***********************/
        /*************************************************************************/

        U->zero();
        Z->zero();
        ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(V)->get(0)->setScalar(1.0);
        ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(V)->get(1)->setScalar(2.0);
        ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(V)->get(2)->setScalar(3.0);
        Hv->setScalar(0.0);

        pint_con->applyAugmentedKKT(*Hv,*V,*U,*Z,tol);

        auto Hv_u = ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(Hv)->get(0);
        auto Hv_z = ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(Hv)->get(1);
        auto Hv_l = ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(Hv)->get(2);

        (*outStream) << "\nParallel compnents of augmented system:\n";
        (*outStream) << "  Norm of Aug*1_u = " << Hv_u->norm() << "\n";
        (*outStream) << "  Norm of Aug*1_z = " << Hv_z->norm() << "\n";
        (*outStream) << "  Norm of Aug*1_l = " << Hv_l->norm() << "\n";

        auto pint_Hv_u = ROL::dynamicPtrCast<ROL::PinTVector<RealT>>(Hv_u);
        auto part_scratch = ROL::dynamicPtrCast<ROL::PartitionedVector<RealT>>(scratch);
        for(int i=0;i<pint_Hv_u->numOwnedSteps();i++) {
          
          part_scratch->get(i)->setScalar(-1.0);  
          part_scratch->get(i)->axpy(1.0,*pint_Hv_u->getVectorPtr(2*i));
          part_scratch->get(i)->axpy(1.0,*pint_Hv_u->getVectorPtr(2*i+1));
        }

        (*outStream) << "  Norm of Aug*1_u modified = " << scratch->norm() << "\n";
      }
    }
    #endif
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if(myRank==0) {
    if (errorFlag != 0)
      std::cout << "End Result: TEST FAILED " << myRank << "\n";
    else
      std::cout << "End Result: TEST PASSED " << myRank << "\n";
  }

  return 0;
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
