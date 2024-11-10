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

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include "../../TOOLS/dynconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "dynpde_thermal-fluids.hpp"
#include "obj_thermal-fluids.hpp"
#include "initial_condition.hpp"

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
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int numProcs = (comm->getSize() > 1) ? comm->getSize() : 0;
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
    int verbosity    = parlist->sublist("General").get("Print Verbosity", 0);
    verbosity        = (myRank==0 ? verbosity : 0);
    bool solveOutput = parlist->sublist("Dynamic Constraint").sublist("Solve").get("Output Iteration History", false);
    solveOutput      = (myRank==0 ? solveOutput : false);
    parlist->sublist("General").set("Print Verbosity", verbosity);
    parlist->sublist("Dynamic Constraint").sublist("Solve").set("Output Iteration History", solveOutput);
    RealT Re         = parlist->sublist("Problem").get("Reynolds Number",200.0);
    RealT Gr         = parlist->sublist("Problem").get("Grashof Number", 40000.0);
    RealT Pr         = parlist->sublist("Problem").get("Prandtl Number", 0.72);
    RealT Tc         = parlist->sublist("Problem").get("Cylinder Temperature", 1.0);

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize mesh data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshReader<RealT>>(*parlist, numProcs);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<DynamicPDE_ThermalFluids<RealT>> pde
      = ROL::makePtr<DynamicPDE_ThermalFluids<RealT>>(*parlist);

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    ROL::Ptr<DynConstraint<RealT>> dyn_con
      = ROL::makePtr<DynConstraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
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
    zk = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(1));
    ROL::Ptr<ROL::PartitionedVector<RealT>> z
      = ROL::PartitionedVector<RealT>::create(*zk, nt);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(4,ROL::nullPtr), qoi_T(1,ROL::nullPtr);
    RealT w1 = parlist->sublist("Problem").get("State Cost",1.0);
    RealT w2 = parlist->sublist("Problem").get("Control Cost",0.0);
    std::vector<RealT> wts = {w1, w1, w1, w2}, wts_T = {w1};
    qoi_vec[0] = ROL::makePtr<QoI_State_ThermalFluids<RealT>>("Dissipation",
                                                              *parlist,
                                                              pde->getVelocityFE(),
                                                              pde->getPressureFE(),
                                                              pde->getThermalFE(),
                                                              pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_State_ThermalFluids<RealT>>("Bouyancy",
                                                              *parlist,
                                                              pde->getVelocityFE(),
                                                              pde->getPressureFE(),
                                                              pde->getThermalFE(),
                                                              pde->getFieldHelper());
    qoi_vec[2] = ROL::makePtr<QoI_DownStreamPower_ThermalFluids<RealT>>(pde->getVelocityFE(),
                                                                        pde->getPressureFE(),
                                                                        pde->getThermalFE(),
                                                                        pde->getVelocityBdryFE(1),
                                                                        pde->getBdryCellLocIds(1),
                                                                        pde->getFieldHelper());
    qoi_vec[3] = ROL::makePtr<QoI_RotationControl_ThermalFluids<RealT>>();
    qoi_T[0]   = ROL::makePtr<QoI_State_ThermalFluids<RealT>>("Tracking",
                                                              *parlist,
                                                              pde->getVelocityFE(),
                                                              pde->getPressureFE(),
                                                              pde->getThermalFE(),
                                                              pde->getFieldHelper());
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_k
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,assembler);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_T
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_T,wts_T,assembler);
    ROL::Ptr<LTI_Objective<RealT>> dyn_obj
      = ROL::makePtr<LTI_Objective<RealT>>(*parlist,obj_k,obj_T);

    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
    std::vector<ROL::TimeStamp<RealT>> timeStamp(nt);
    for( int k=0; k<nt; ++k ) {
      timeStamp.at(k).t.resize(2);
      timeStamp.at(k).t.at(0) = k*dt;
      timeStamp.at(k).t.at(1) = (k+1)*dt;
    }
    // Compute initial condition
    std::clock_t timer_init = std::clock();
    std::stringstream file;
    file << "initial_condition"
         << "_Re" << static_cast<int>(Re)
         << "_Gr" << static_cast<int>(Gr)
         << "_Pr" << static_cast<int>(100.0*Pr)
         << "_Tc" << static_cast<int>(Tc) << ".txt";
    std::ifstream infile(file.str());
    if (infile.good()) {
      dyn_con->inputTpetraVector(u0_ptr, file.str());
    }
    else {
      PotentialFlow<RealT> pf(pde->getVelocityFE(),
                              pde->getPressureFE(),
                              pde->getThermalFE(),
                              pde->getCellNodes(),
                              assembler->getDofManager()->getCellDofs(),
                              assembler->getCellIds(),
                              pde->getFieldHelper(),
                              *parlist);
      pf.build(u0_ptr);
      computeInitialCondition<RealT>(u0,ck,uo,un,zk,dyn_con,dt,*outStream);
      dyn_con->outputTpetraVector(u0_ptr, file.str());
    }
    *outStream << "Initial condition time: "
               << static_cast<RealT>(std::clock()-timer_init)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
    // Construct reduce dynamic objective function
    ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
    ROL::Ptr<ROL::ReducedDynamicObjective<RealT>> obj
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<ROL::PartitionedVector<RealT>> dz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::PartitionedVector<RealT>> hz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      zk->randomize(); z->randomize(); dz->randomize(); hz->randomize();
      uo->randomize(); un->randomize();
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicObjectiveCheck<RealT>::check(*dyn_obj,validate,*uo,*un,*zk);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*zk);
      obj->checkGradient(*z,*dz,true,*outStream);
      obj->checkHessVec(*z,*dz,true,*outStream);
      obj->checkHessSym(*z,*dz,*hz,true,*outStream);
    }

    /*************************************************************************/
    /***************** OUTPUT UNCONTROLLED STATE *****************************/
    /*************************************************************************/
    bool printU0 = parlist->sublist("Problem").get("Print Uncontrolled State", false);
    if (printU0) {
      std::clock_t timer_print0 = std::clock();
      // Output state and control to file
      zk->zero(); uo->set(*u0); un->zero();
      for (int k = 1; k < nt; ++k) {
        // Print previous state to file
        std::stringstream u0file;
        u0file << "uncontrolled_state." << k-1 << ".txt";
        dyn_con->outputTpetraVector(uo_ptr, u0file.str());
        // Advance time stepper
        dyn_con->solve(*ck, *uo, *un, *zk, timeStamp[k]);
        uo->set(*un);
      }
      // Print previous state to file
      std::stringstream u0file;
      u0file << "uncontrolled_state." << nt-1 << ".txt";
      dyn_con->outputTpetraVector(uo_ptr, u0file.str());
      *outStream << "Output uncontrolled state time: "
                 << static_cast<RealT>(std::clock()-timer_print0)/static_cast<RealT>(CLOCKS_PER_SEC)
                 << " seconds." << std::endl << std::endl;
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> problem(obj,z);
    ROL::OptimizationSolver<RealT> solver(problem,*parlist);
    z->zero();
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    // Output state and control to file
    uo->set(*u0); un->zero();
    for (int k = 1; k < nt; ++k) {
      // Print previous state to file
      std::stringstream ufile;
      ufile << "state." << k-1 << ".txt";
      dyn_con->outputTpetraVector(uo_ptr, ufile.str());
      // Print current control
      if (myRank == 0) {
        std::stringstream zname;
        zname << "control." << k-1 << ".txt";
        std::ofstream zfile;
        zfile.open(zname.str());
        zfile << std::scientific << std::setprecision(15);
        ROL::Ptr<std::vector<RealT>> zn
          = ROL::dynamicPtrCast<PDE_OptVector<RealT>>(z->get(k-1))->getParameter()->getVector();
        zfile << std::right << std::setw(25) << (*zn)[0];
        zfile.close();
      }
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z->get(k), timeStamp[k]);
      uo->set(*un);
    }
    // Print previous state to file
    std::stringstream ufile;
    ufile << "state." << nt-1 << ".txt";
    dyn_con->outputTpetraVector(uo_ptr, ufile.str());
    // Print current control
    if (myRank == 0) {
      std::stringstream zname;
      zname << "control." << nt-1 << ".txt";
      std::ofstream zfile;
      zfile.open(zname.str());
      zfile << std::scientific << std::setprecision(15);
      ROL::Ptr<std::vector<RealT>> zn
        = ROL::dynamicPtrCast<PDE_OptVector<RealT>>(z->get(nt-1))->getParameter()->getVector();
      zfile << std::right << std::setw(25) << (*zn)[0];
      zfile.close();
    }

    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
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
