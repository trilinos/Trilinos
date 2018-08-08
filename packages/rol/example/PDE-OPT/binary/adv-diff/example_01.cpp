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
    \brief Shows how to solve the binary advection-diffusion control problem.
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
//#include <fenv.h>

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PEBBL_Driver.hpp"

#include "ROL_TeuchosBranchHelper_PEBBL.hpp"
#include "extractQP.hpp"

#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "pde_adv_diff.hpp"
#include "qoi_adv_diff.hpp"
#include "mesh_adv_diff.hpp"
#include "branchHelper.hpp"

#include "ROL_HelperFunctions.hpp"


int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile(filename);
    Teuchos::RCP<Teuchos::ParameterList> tparlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, tparlist.ptr() );

    // Problem dimensions
    const int controlDim = 9;

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshManager_adv_diff<RealT>>(*tparlist);
    ROL::Ptr<PDE_adv_diff<RealT>> pde
      = ROL::makePtr<PDE_adv_diff<RealT>>(*tparlist);
    ROL::Ptr<Linear_PDE_Constraint<RealT>> con
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde,meshMgr,comm,*tparlist,*outStream);
    const ROL::Ptr<Assembler<RealT>> assembler = con->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr;
    u_ptr = assembler->createStateVector();
    p_ptr = assembler->createStateVector();
    r_ptr = assembler->createResidualVector();
    ROL::Ptr<ROL::Vector<RealT>> up, pp, rp, zp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*tparlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*tparlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*tparlist);
    zp = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(1,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<RealT>>(pde->getFE());
    RealT stateCost = parlist->sublist("Problem").get("State Cost",4.0);
    std::vector<RealT> wts = {stateCost};
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,assembler);
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    /*************************************************************************/
    /***************** BUILD KNAPSACK CONSTRAINT *****************************/
    /*************************************************************************/
    RealT ctrlCost = parlist->sublist("Problem").get("Control Cost", 4.0);
    bool  useIneq  = parlist->sublist("Problem").get("Use Inequality", false);
    RealT budget   = (useIneq ? static_cast<RealT>(0) : ctrlCost);
    ROL::Ptr<QoI<RealT>> knapsack_qoi
      = ROL::makePtr<QoI_Control_Cost_adv_diff<RealT>>(budget);
    ROL::Ptr<ROL::Constraint<RealT>> knapsack_con
      = ROL::makePtr<IntegralOptConstraint<RealT>>(knapsack_qoi,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT>> knapsack_bnd = ROL::nullPtr;
    if (useIneq) {
      ROL::Ptr<ROL::Vector<RealT>> klop, khip;
      klop = ROL::makePtr<ROL::StdVector<RealT>>(1,static_cast<RealT>(0));
      khip = ROL::makePtr<ROL::StdVector<RealT>>(1,ctrlCost);
      knapsack_bnd = ROL::makePtr<ROL::Bounds<RealT>>(klop,khip);
    }
    ROL::Ptr<ROL::Vector<RealT>> knapsack_mul
      = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<ROL::Vector<RealT>> zlop
      = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim,0.0));
    ROL::Ptr<ROL::Vector<RealT>> zhip
      = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim,1.0));
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    bool useQP = parlist->sublist("Problem").get("Use Quadratic Program",false);
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    ROL::Ptr<ROL::BranchHelper_PEBBL<RealT>> bHelper;
    if (!useQP) {
      problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(robj,zp,bnd,knapsack_con,knapsack_mul,knapsack_bnd);
      bHelper = ROL::makePtr<PDEOPT_BranchHelper_PEBBL<RealT>>();
    }
    else {
      extractQP<RealT> extract(robj,zp,bnd,knapsack_con,knapsack_mul,knapsack_bnd);
      problem = extract();
      bHelper = ROL::makePtr<ROL::TeuchosBranchHelper_PEBBL<int,RealT>>();
    }
    bool derivCheck = parlist->sublist("Problem").get("Check Derivatives",true);
    if (derivCheck) {
      ROL::Ptr<Tpetra::MultiVector<>> vu_ptr, wu_ptr;
      vu_ptr = assembler->createStateVector();
      wu_ptr = assembler->createStateVector();
      ROL::Ptr<ROL::Vector<RealT>> vup, wup, vzp, wzp, vp, wp, xp;
      vup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(vu_ptr,pde,assembler,*tparlist);
      wup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(wu_ptr,pde,assembler,*tparlist);
      vzp = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim));
      wzp = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(controlDim));
      xp  = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);
      vp  = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(vup,vzp);
      wp  = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(wup,wzp);
      xp->randomize(); vp->randomize(); wp->randomize();
      pp->randomize(); rp->randomize();
      
      obj->checkGradient(*xp,*vp,true,*outStream);
      obj->checkGradient_1(*up,*zp,*vup,true,*outStream);
      obj->checkGradient_2(*up,*zp,*vzp,true,*outStream);
      obj->checkHessVec(*xp,*vp,true,*outStream);
      obj->checkHessVec_11(*up,*zp,*vup,true,*outStream);
      obj->checkHessVec_12(*up,*zp,*vzp,true,*outStream);
      obj->checkHessVec_21(*up,*zp,*vup,true,*outStream);
      obj->checkHessVec_22(*up,*zp,*vzp,true,*outStream);
      obj->checkHessSym(*xp,*vp,*wp,true,*outStream);

      con->checkSolve(*up,*zp,*rp,true,*outStream);
      con->checkAdjointConsistencyJacobian(*pp,*vp,*xp,true,*outStream);
      con->checkAdjointConsistencyJacobian_1(*pp,*vup,*up,*zp,true,*outStream);
      con->checkAdjointConsistencyJacobian_2(*pp,*vzp,*up,*zp,true,*outStream);
      con->checkInverseJacobian_1(*rp,*vup,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*rp,*vup,*up,*zp,true,*outStream);
      con->checkApplyJacobian(*xp,*vp,*rp,true,*outStream);
      con->checkApplyJacobian_1(*up,*zp,*vup,*rp,true,*outStream);
      con->checkApplyJacobian_2(*up,*zp,*vzp,*rp,true,*outStream);
      con->checkApplyAdjointHessian(*xp,*pp,*vp,*xp,true,*outStream);
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*vup,wup->dual(),true,*outStream);
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*vup,wzp->dual(),true,*outStream);
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*vzp,wup->dual(),true,*outStream);
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*vzp,wzp->dual(),true,*outStream);
      problem->check(*outStream);
      xp->zero(); pp->zero(); rp->zero();
    }

    ROL::ROL_PEBBL_Driver<RealT> pebbl(problem,parlist,bHelper,3,outStream);
    pebbl.solve(argc,argv,*outStream);

    //*outStream << "OPTIMAL CONTROLS" << std::endl;
    //for (int i=0; i<controlDim; ++i) {
    //  *outStream << (*z_ptr)[i] << "  ";
    //}
    //*outStream << std::endl;
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
