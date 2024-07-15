// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Solves a source inversion problem governed by the
           advection-diffusion equation.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Solver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"
#include "ROL_TypeBIndicatorObjective.hpp"
#include <iostream>
//#include <fenv.h>

#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/lindynconstraint.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "dynpde_adv_diff.hpp"
#include "obj_adv_diff.hpp"
#include "mesh_adv_diff.hpp"
#include "l1penaltydynamic.hpp"

#include "ROL_TypeP_TrustRegionAlgorithm.hpp"


int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input.xml");
    int nt         = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
    RealT T        = parlist->sublist("Time Discretization").get("End Time",             1.0);
    RealT dt       = T/static_cast<RealT>(nt);
    int controlDim = 9;

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshManager_adv_diff<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<DynamicPDE_adv_diff<RealT>> pde
      = ROL::makePtr<DynamicPDE_adv_diff<RealT>>(*parlist);

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    bool isLTI = !parlist->sublist("Problem").get("Time Varying Coefficients",false);
    ROL::Ptr<LinDynConstraint<RealT>> dyn_con
      = ROL::makePtr<LinDynConstraint<RealT>>(pde,meshMgr,comm,*parlist,isLTI,*outStream);
    dyn_con->getAssembler()->printMeshData(*outStream);

    /*************************************************************************/
    /***************** BUILD STATE VECTORS ***********************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u0_ptr, uo_ptr, un_ptr, ck_ptr;
    u0_ptr = dyn_con->getAssembler()->createStateVector();
    uo_ptr = dyn_con->getAssembler()->createStateVector();
    un_ptr = dyn_con->getAssembler()->createStateVector();
    ck_ptr = dyn_con->getAssembler()->createResidualVector();
    ROL::Ptr<ROL::Vector<RealT>> u0, uo, un, ck, zk;
    u0 = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u0_ptr,pde,*dyn_con->getAssembler());
    uo = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uo_ptr,pde,*dyn_con->getAssembler());
    un = ROL::makePtr<PDE_PrimalSimVector<RealT>>(un_ptr,pde,*dyn_con->getAssembler());
    ck = ROL::makePtr<PDE_DualSimVector<RealT>>(ck_ptr,pde,*dyn_con->getAssembler());
    zk = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim));
    ROL::Ptr<ROL::PartitionedVector<RealT>> z
      = ROL::PartitionedVector<RealT>::create(*zk, nt);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_State_MY_adv_diff<RealT>>(pde->getFE());
    RealT stateCost   = parlist->sublist("Problem").get("State Cost",1e5);
    RealT myCost      = 1e0;
    std::vector<RealT> wts = {stateCost, myCost};
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_k
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,dyn_con->getAssembler());
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
      = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl, outStream);
 
    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT AND L1 PENALTY *****************/
    /*************************************************************************/
    ROL::Ptr<ROL::PartitionedVector<RealT>> zlo = ROL::PartitionedVector<RealT>::create(*zk, nt);
    ROL::Ptr<ROL::PartitionedVector<RealT>> zhi = ROL::PartitionedVector<RealT>::create(*zk, nt);
    zlo->setScalar(-static_cast<RealT>(1)); zhi->setScalar(static_cast<RealT>(1));   
    ROL::Ptr<L1_Dyn_Objective<RealT>> nobj
      = ROL::makePtr<L1_Dyn_Objective<RealT>>(*parlist, timeStamp, zlo, zhi); 

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<ROL::PartitionedVector<RealT>> dz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      ROL::Ptr<ROL::PartitionedVector<RealT>> hz = ROL::PartitionedVector<RealT>::create(*zk, nt);
      z->randomize(); dz->randomize(); hz->randomize();
      zk->randomize(); uo->randomize(); un->randomize();
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicObjectiveCheck<RealT>::check(*dyn_obj,validate,*uo,*un,*zk,timeStamp[0]);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*zk,timeStamp[0]);
      obj->checkGradient(*z,*dz,true,*outStream);
      obj->checkHessVec(*z,*dz,true,*outStream);
      obj->checkHessSym(*z,*dz,*hz,true,*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    RealT gtol(1e-2), zerr(0), ztol0(1e-6), ztol(1);
    z->zero();
    ROL::Ptr<ROL::Vector<RealT>> zprev = z->clone(), zdiff = z->clone();
    ROL::Ptr<ROL::TypeP::TrustRegionAlgorithm<RealT>> algo;
    parlist->sublist("Status Test").set("Use Relative Tolerances",false);
    for (unsigned i = 0; i < 20; ++i) {
      *outStream << std::endl << "Moreau-Yosida Parameter: "
                 << myCost << std::endl
                 << "Gradient Tolerance: "
                 << gtol << std::endl << std::endl;

      ztol = ztol0 / std::max(static_cast<RealT>(1),z->norm());
      zprev->set(*z);
      parlist->sublist("Status Test").set("Gradient Tolerance",gtol);
      parlist->sublist("Status Test").set("Step Tolerance",static_cast<RealT>(1e-4)*gtol);
      algo = ROL::makePtr<ROL::TypeP::TrustRegionAlgorithm<RealT>>(*parlist); 
      std::clock_t timer = std::clock();
      algo->run(*z, *obj, *nobj, *outStream); 
      *outStream << "Optimization time: "
                 << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
                 << " seconds." << std::endl << std::endl;
      parlist->sublist("Reduced Dynamic Objective").set("State Rank", static_cast<int>(obj->getStateRank()));
      parlist->sublist("Reduced Dynamic Objective").set("Adjoint Rank", static_cast<int>(obj->getAdjointRank()));
      parlist->sublist("Reduced Dynamic Objective").set("State Sensitivity Rank", static_cast<int>(obj->getStateSensitivityRank()));

      zdiff->set(*z); zdiff->axpy(static_cast<RealT>(-1),*zprev);
      zerr = zdiff->norm();
      *outStream << std::endl << "Control Difference: "
                 << zerr << std::endl << std::endl;
      if (zerr < ztol && algo->getState()->gnorm < static_cast<RealT>(1e-8)) break;

      myCost *= static_cast<RealT>(2e0);
      gtol   *= static_cast<RealT>(1e-1); gtol = std::max(gtol,static_cast<RealT>(1e-10));
      wts     = {stateCost, myCost};
      obj_k   = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,dyn_con->getAssembler());
      dyn_obj = ROL::makePtr<LTI_Objective<RealT>>(*parlist,obj_k,false);
      obj     = ROL::makePtr<ROL::ReducedDynamicObjective<RealT>>(dyn_obj, dyn_con, u0, zk, ck, timeStamp, rpl, outStream);
    }

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    // Output control to file
    if ( myRank == 0 ) {
      ROL::Ptr<std::vector<RealT>> zn;
      std::ofstream zfile;
      zfile.open("control.txt");
      zfile << std::scientific << std::setprecision(15);
      for (int n = 0; n < nt; ++n) {
        zn = ROL::dynamicPtrCast<PDE_OptVector<RealT>>(z->get(n))->getParameter()->getVector();
        for (int i = 0; i < controlDim; ++i) {
          zfile << std::right << std::setw(25) << (*zn)[i];
        }
        zfile << std::endl;
      }
      zfile.close();
    }
    // Output state to file
    uo->set(*u0); un->zero();
    for (int k = 1; k < nt; ++k) {
      // Print previous state to file
      std::stringstream ufile;
      ufile << "state." << k-1 << ".txt";
      dyn_con->outputTpetraVector(uo_ptr, ufile.str());
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z->get(k), timeStamp[k]);
      uo->set(*un);
    }
    // Print previous state to file
    std::stringstream ufile;
    ufile << "state." << nt-1 << ".txt";
    dyn_con->outputTpetraVector(uo_ptr, ufile.str());
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
