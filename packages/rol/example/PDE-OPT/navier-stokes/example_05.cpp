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
    \brief Shows how to solve the stochastic Stefan-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_Algorithm.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_navier-stokes.hpp"
#include "obj_navier-stokes.hpp"

#include "ROL_OptimizationSolver.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"

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
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    std::string filename = "input_ex05.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    //const int stochDim = 2;
    //const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<PDE_NavierStokes<RealT> > pde
      = ROL::makePtr<PDE_NavierStokes<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

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
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > pp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    // Create residual vectors
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->randomize(); //r_ptr->putScalar(static_cast<RealT>(1));
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create control vector and set to ones
    ROL::Ptr<Tpetra::MultiVector<> >  z_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > yz_ptr = assembler->createControlVector();
    z_ptr->randomize();  z_ptr->putScalar(static_cast<RealT>(0));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(0));
    yz_ptr->randomize(); //yz_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > dzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > yzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(yz_ptr,pde,assembler,*parlist);
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
    qoi_vec[0] = ROL::makePtr<QoI_State_NavierStokes<RealT>>(*parlist,
                                                                pde->getVelocityFE(),
                                                                pde->getPressureFE(),
                                                                pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_NavierStokes<RealT>>(pde->getVelocityFE(),
                                                                    pde->getPressureFE(),
                                                                    pde->getVelocityBdryFE(),
                                                                    pde->getBdryCellLocIds(),
                                                                    pde->getFieldHelper());
    ROL::Ptr<StdObjective_NavierStokes<RealT> > std_obj
      = ROL::makePtr<StdObjective_NavierStokes<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> >  zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(static_cast<RealT>(0));
    zhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zlopde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zhipde
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > zlop = ROL::makePtr<PDE_OptVector<RealT>>(zlopde);
    ROL::Ptr<ROL::Vector<RealT> > zhip = ROL::makePtr<PDE_OptVector<RealT>>(zhipde);
    ROL::Ptr<ROL::BoundConstraint<RealT> > zbnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    bool useBounds = parlist->sublist("Problem").get("Use bounds", false);
    if (!useBounds) {
      zbnd->deactivate();
    }
    // State bounds
    ROL::Ptr<Tpetra::MultiVector<> > ulo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::Vector<RealT> > ulop
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ulo_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > uhip
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uhi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ubnd
      = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    ubnd->deactivate();

    // SimOpt bounds
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT>>(ubnd,zbnd);

    /*************************************************************************/
    /***************** BUILD            PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      up->checkVector(*pp,*dup,true,*outStream);
      zp->checkVector(*yzp,*dzp,true,*outStream);
      obj->checkGradient(x,d,true,*outStream);
      obj->checkHessVec(x,d,true,*outStream);
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);
      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    RealT tol(1.e-8);
    up->setScalar(RealT(1));
    RandomizeVector(*zp);//->setScalar(RealT(1));
    con->solve(*rp,*up,*zp,tol);    
    ROL::OptimizationProblem<RealT> optProb(obj, makePtrFromRef(x), bnd, con, pp);
    ROL::OptimizationSolver<RealT> optSolver(optProb, *parlist);
    optSolver.solve(*outStream);
    std::clock_t timer = std::clock();
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
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
