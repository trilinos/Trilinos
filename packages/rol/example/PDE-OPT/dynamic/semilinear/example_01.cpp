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
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Bounds.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Solver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include "../../TOOLS/lindynconstraint.hpp"
#include "../../TOOLS/dynconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "dynpde_semilinear.hpp"
#include "obj_semilinear.hpp"


int main(int argc, char *argv[]) {
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
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input_ex01.xml");
    int nt         = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
    RealT T        = parlist->sublist("Time Discretization").get("End Time",             1.0);
    RealT dt       = T/static_cast<RealT>(nt);

    parlist->sublist("Reduced Dynamic Objective").set("State Domain Seed",12321*(myRank+1));
    parlist->sublist("Reduced Dynamic Objective").set("State Range Seed", 32123*(myRank+1));
    parlist->sublist("Reduced Dynamic Objective").set("Adjoint Domain Seed",23432*(myRank+1));
    parlist->sublist("Reduced Dynamic Objective").set("Adjoint Range Seed", 43234*(myRank+1));
    parlist->sublist("Reduced Dynamic Objective").set("State Sensitivity Domain Seed",34543*(myRank+1));
    parlist->sublist("Reduced Dynamic Objective").set("State Sensitivity Range Seed", 54345*(myRank+1));

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
    ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con;
    ROL::Ptr<Assembler<RealT>> assembler;
    int nltype = parlist->sublist("Semilinear").get("Nonlinearity Type",0);
    if (nltype == 0) {
      dyn_con = ROL::makePtr<LinDynConstraint<RealT>>(pde,meshMgr,comm,*parlist,true,*outStream);
      assembler = ROL::staticPtrCast<LinDynConstraint<RealT>>(dyn_con)->getAssembler();
    }
    else {
      dyn_con = ROL::makePtr<DynConstraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
      assembler = ROL::staticPtrCast<DynConstraint<RealT>>(dyn_con)->getAssembler();
    }
    dyn_con->setSolveParameters(*parlist);
    assembler->printMeshData(*outStream);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u0_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> uo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> un_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> ck_ptr = assembler->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<>> zk_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT>> u0, uo, un, ck, zk;
    u0 = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u0_ptr,pde,*assembler,*parlist);
    uo = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uo_ptr,pde,*assembler,*parlist);
    un = ROL::makePtr<PDE_PrimalSimVector<RealT>>(un_ptr,pde,*assembler,*parlist);
    ck = ROL::makePtr<PDE_DualSimVector<RealT>>(ck_ptr,pde,*assembler,*parlist);
    zk = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zk_ptr,pde,*assembler,*parlist);
    ROL::Ptr<ROL::PartitionedVector<RealT>> z
      = ROL::PartitionedVector<RealT>::create(*zk, nt);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_Semilinear<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_Semilinear<RealT>>(pde->getFE());
    RealT w1 = parlist->sublist("Problem").get("State Cost",1.0);
    RealT w2 = parlist->sublist("Problem").get("Control Cost",1e-2);
    std::vector<RealT> wts = {w1, w2};
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_k
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,assembler);
    ROL::Ptr<LTI_Objective<RealT>> dyn_obj
      = ROL::makePtr<LTI_Objective<RealT>>(*parlist,obj_k,false);

    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
    std::vector<ROL::TimeStamp<RealT>> timeStamp(nt);
    for( int k=0; k<nt; ++k ) {
      timeStamp.at(k).t.resize(2);
      timeStamp.at(k).t.at(0) = k*dt;
      timeStamp.at(k).t.at(1) = (k+1)*dt;
    }
    ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
    ROL::Ptr<ROL::ReducedDynamicObjective<RealT>> obj
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl,outStream);
    // print uncontrolled state
    bool print_unc = true;
    if (print_unc) {
      std::clock_t timer_print_unc = std::clock();
      // Output state and control to file
      uo->set(*u0); un->zero(); z->zero();
      for (int k = 1; k < nt; ++k) {
        // Print previous state to file
        std::stringstream ufile;
        ufile << "uncontrolled_state." << k-1 << ".txt";
        assembler->outputTpetraVector(uo_ptr, ufile.str());
        // Advance time stepper
        dyn_con->solve(*ck, *uo, *un, *z->get(k), timeStamp[k]);
        uo->set(*un);
      }
      // Print previous state to file
      std::stringstream ufile;
      ufile << "uncontrolled_state." << nt-1 << ".txt";
      assembler->outputTpetraVector(uo_ptr, ufile.str());
      *outStream << "Output time: "
                 << static_cast<RealT>(std::clock()-timer_print_unc)/static_cast<RealT>(CLOCKS_PER_SEC)
                 << " seconds." << std::endl << std::endl;
    }
 
    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    bool deactivate = parlist->sublist("Problem").get("Deactivate Bound Constraints",false);
    RealT lbnd = parlist->sublist("Problem").get("Lower Bound",0.0);
    RealT ubnd = parlist->sublist("Problem").get("Upper Bound",1.0);
    ROL::Ptr<ROL::PartitionedVector<RealT>> zlo = ROL::PartitionedVector<RealT>::create(*zk, nt);
    ROL::Ptr<ROL::PartitionedVector<RealT>> zhi = ROL::PartitionedVector<RealT>::create(*zk, nt);
    zlo->setScalar(lbnd);
    zhi->setScalar(ubnd);
    ROL::Ptr<ROL::BoundConstraint<RealT>>   bnd = ROL::makePtr<ROL::Bounds<RealT>>(zlo,zhi);

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
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Ptr<ROL::Problem<RealT>> problem = ROL::makePtr<ROL::Problem<RealT>>(obj,z);
    if (!deactivate) problem->addBoundConstraint(bnd);
    problem->finalize(false,true,*outStream);
    ROL::Solver<RealT> solver(problem,*parlist);
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
      assembler->outputTpetraVector(uo_ptr, ufile.str());
      // Print current control
      std::stringstream zfile;
      zfile << "control." << k-1 << ".txt";
      assembler->outputTpetraVector(ROL::staticPtrCast<PDE_PrimalOptVector<RealT>>(z->get(k-1))->getVector(),zfile.str());
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z->get(k), timeStamp[k]);
      uo->set(*un);
    }
    // Print previous state to file
    std::stringstream ufile;
    ufile << "state." << nt-1 << ".txt";
    assembler->outputTpetraVector(uo_ptr, ufile.str());
    // Print current control
    std::stringstream zfile;
    zfile << "control." << nt-1 << ".txt";
    assembler->outputTpetraVector(ROL::staticPtrCast<PDE_PrimalOptVector<RealT>>(z->get(nt-1))->getVector(),zfile.str());
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
