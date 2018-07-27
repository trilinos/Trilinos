// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_DefaultPlatform.hpp"
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

#include "../../TOOLS/dynconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "dynpde_navier-stokes.hpp"
#include "obj_navier-stokes.hpp"


int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

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
    int verbosity  = parlist->sublist("General").get("Print Verbosity", 0);
    verbosity      = (myRank==0 ? verbosity : 0);
    parlist->sublist("General").set("Print Verbosity", verbosity);
    bool useParametricControl = parlist->sublist("Problem").get("Use Parametric Control", false);

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize mesh data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<DynamicPDE_NavierStokes<RealT>> pde
      = ROL::makePtr<DynamicPDE_NavierStokes<RealT>>(*parlist);

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
    if (!useParametricControl) {
      ROL::Ptr<Tpetra::MultiVector<>> zk_ptr = assembler->createControlVector();
      zk = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zk_ptr,pde,*assembler,*parlist);
    }
    else {
      zk = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(1));
    }
    ROL::Ptr<ROL::PartitionedVector<RealT>> z
      = ROL::PartitionedVector<RealT>::create(*zk, nt);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(1,ROL::nullPtr);
    RealT w1 = parlist->sublist("Problem").get("State Cost",1.0);
    std::vector<RealT> wts = {w1};
    qoi_vec[0] = ROL::makePtr<QoI_State_NavierStokes<RealT>>(*parlist,
                                                              pde->getVelocityFE(),
                                                              pde->getPressureFE(),
                                                              pde->getFieldHelper());
    if (!useParametricControl) {
      qoi_vec.push_back(ROL::makePtr<QoI_L2Penalty_NavierStokes<RealT>>(pde->getVelocityFE(),
                                                                        pde->getPressureFE(),
                                                                        pde->getVelocityBdryFE(),
                                                                        pde->getBdryCellLocIds(),
                                                                        pde->getFieldHelper()));
      RealT w2 = parlist->sublist("Problem").get("Control Cost",1e-2);
      wts.push_back(w2);
    }
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_k
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,assembler);
    ROL::Ptr<LTI_Objective<RealT>> dyn_obj
      = ROL::makePtr<LTI_Objective<RealT>>(obj_k,*zk,*parlist);

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
    std::clock_t timer_print0 = std::clock();
    // Output state and control to file
    z->zero();
    uo->set(*u0); un->zero();
    for (int k = 1; k < nt; ++k) {
      // Print previous state to file
      std::stringstream u0file;
      u0file << "uncontrolled_state." << k-1 << ".txt";
      dyn_con->outputTpetraVector(uo_ptr, u0file.str());
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z->get(k), timeStamp[k]);
      uo->set(*un);
    }
    // Print previous state to file
    std::stringstream u0file;
    u0file << "uncontrolled_state." << nt-1 << ".txt";
    dyn_con->outputTpetraVector(uo_ptr, u0file.str());
    *outStream << "Output uncontrolled state time: "
               << static_cast<RealT>(std::clock()-timer_print0)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

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
      if (!useParametricControl) {
        std::stringstream zfile;
        zfile << "control." << k-1 << ".txt";
        dyn_con->outputTpetraVector(ROL::dynamicPtrCast<PDE_PrimalOptVector<RealT>>(z->get(k-1))->getVector(),zfile.str());
      }
      else {
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
    if (!useParametricControl) {
      std::stringstream zfile;
      zfile << "control." << nt-1 << ".txt";
      dyn_con->outputTpetraVector(ROL::dynamicPtrCast<PDE_PrimalOptVector<RealT>>(z->get(nt-1))->getVector(),zfile.str());
    }
    else {
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
    }
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
