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
    \brief Shows how to solve the Poisson topology optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "pde_poisson_topOpt.hpp"
#include "obj_poisson_topOpt.hpp"
#include "mesh_poisson_topOpt.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    RealT tol(1.e-8), one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT domainWidth  = parlist->sublist("Geometry").get("Width", 1.0);
    const RealT domainHeight = parlist->sublist("Geometry").get("Height", 1.0);
    const RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-2);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Poisson_TopOpt<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Poisson_TopOpt<RealT> > pde
      = ROL::makePtr<PDE_Poisson_TopOpt<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    con->setSolveParameters(*parlist);
    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT> > pdeFilter
      = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
    u_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    // Create control vector and set to ones
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = assembler->createControlVector();
    z_ptr->putScalar(volFraction);
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    // Create Lagrange multiplier vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > l_ptr = assembler->createStateVector();
    l_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > lp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(l_ptr,pde,assembler,*parlist);
    // Create residual vector and set to zeros
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create state direction vector and set to random
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    du_ptr->randomize(); //du_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    // Create control direction vector and set to random
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    dz_ptr->randomize(); //dz_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > dzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);
    // Create control test vector.
    ROL::Ptr<Tpetra::MultiVector<> > rz_ptr = assembler->createControlVector();
    rz_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > rzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(rz_ptr,pde,assembler,*parlist);

    ROL::Ptr<Tpetra::MultiVector<> > dualu_ptr = assembler->createStateVector();
    ROL::Ptr<ROL::Vector<RealT> > dualup
      = ROL::makePtr<PDE_DualSimVector<RealT>>(dualu_ptr,pde,assembler,*parlist);
    ROL::Ptr<Tpetra::MultiVector<> > dualz_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > dualzp
      = ROL::makePtr<PDE_DualOptVector<RealT>>(dualz_ptr,pde,assembler,*parlist);

    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      pdeWithFilter = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conFilter, *rp, *rp, *up, *zp, *zp);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Initialize compliance objective function.
    Teuchos::ParameterList list(*parlist);
    list.sublist("Vector").sublist("Sim").set("Use Riesz Map",true);
    list.sublist("Vector").sublist("Sim").set("Lump Riesz Map",false);
    // Has state Riesz map enabled for mesh-independent compliance scaling.
    ROL::Ptr<Tpetra::MultiVector<> > f_ptr = assembler->createResidualVector();
    f_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > fp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(f_ptr,pde,assembler,list);
    up->zero();
    con->value(*fp, *up, *zp, tol);
    RealT objScaling = objFactor, fnorm2 = fp->dot(*fp);
    if (fnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
      objScaling /= fnorm2;
    }
    u_ptr->randomize();
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(1,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_Energy_Poisson_TopOpt<RealT>>(pde->getFE(),
                                                                   pde->getForce(),
                                                                   objScaling);
    ROL::Ptr<StdObjective_Poisson_TopOpt<RealT> > std_obj
      = ROL::makePtr<StdObjective_Poisson_TopOpt<RealT>>();
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Objective<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,lp,true,false);

    // Initialize volume constraint
    ROL::Ptr<QoI<RealT> > qoi_vol
      = ROL::makePtr<QoI_Volume_Poisson_TopOpt<RealT>>(pde->getFE(),*parlist);
    ROL::Ptr<IntegralOptConstraint<RealT> > vcon
      = ROL::makePtr<IntegralOptConstraint<RealT>>(qoi_vol,assembler);
    // Create volume constraint vector and set to zero
    RealT vecScaling = one / std::pow(domainWidth*domainHeight*(one-volFraction), 2);
    ROL::Ptr<std::vector<RealT> > scalevec_ptr = ROL::makePtr<std::vector<RealT>>(1,vecScaling);
    ROL::Ptr<std::vector<RealT> > c1_ptr = ROL::makePtr<std::vector<RealT>>(1,0);
    ROL::Ptr<ROL::Vector<RealT> > c1p = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(c1_ptr, scalevec_ptr);
    ROL::Ptr<std::vector<RealT> > c2_ptr = ROL::makePtr<std::vector<RealT>>(1,1);
    ROL::Ptr<ROL::Vector<RealT> > c2p = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(c2_ptr, scalevec_ptr);

    // Build bound constraint
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = assembler->createControlVector();
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > hip
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Initialize Augmented Lagrangian functional.
    ROL::AugmentedLagrangian<RealT> augLag(robj,vcon,*c2p,1,*zp,*c1p,*parlist);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Opt Vector\n";
      zp->checkVector(*dzp,*rzp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(x,*lp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of PDE Constraint\n";
      con->checkApplyAdjointHessian_11(*up,*zp,*lp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of PDE Constraint\n";
      con->checkApplyAdjointHessian_21(*up,*zp,*lp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of PDE Constraint\n";
      con->checkApplyAdjointHessian_12(*up,*zp,*lp,*dup,*dualzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of PDE Constraint\n";
      con->checkApplyAdjointHessian_22(*up,*zp,*lp,*dzp,*dualzp,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of Filter\n";
      conFilter->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of Filter\n";
      conFilter->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of Filter\n";
      conFilter->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of Filter\n";
      conFilter->checkApplyAdjointHessian(x,*lp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Filter\n";
      conFilter->checkApplyAdjointHessian_11(*up,*zp,*lp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of Filter\n";
      conFilter->checkApplyAdjointHessian_21(*up,*zp,*lp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of Filter\n";
      conFilter->checkApplyAdjointHessian_12(*up,*zp,*lp,*dup,*dualzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Filter\n";
      conFilter->checkApplyAdjointHessian_22(*up,*zp,*lp,*dzp,*dualzp,true,*outStream);
      *outStream << "\n";
      conFilter->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      conFilter->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      conFilter->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyAdjointHessian(x,*lp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyAdjointHessian_11(*up,*zp,*lp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyAdjointHessian_21(*up,*zp,*lp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyAdjointHessian_12(*up,*zp,*lp,*dup,*dualzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Filtered PDE Constraint\n";
      pdeWithFilter->checkApplyAdjointHessian_22(*up,*zp,*lp,*dzp,*dualzp,true,*outStream);
      *outStream << "\n";
      pdeWithFilter->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      pdeWithFilter->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      pdeWithFilter->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      robj->checkHessVec(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Full Jacobian of Volume Constraint\n";
      vcon->checkApplyJacobian(x,d,*c1p,true,*outStream);
      *outStream << "\n";
      vcon->checkAdjointConsistencyJacobian(*c1p,d,x,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of Volume Constraint\n";
      vcon->checkApplyAdjointHessian(x,*c2p,d,x,true,*outStream);
      *outStream << "\n\nCheck Gradient of Augmented Lagrangian Function\n";
      augLag.checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Augmented Lagrangian Function\n";
      augLag.checkHessVec(*zp,*dzp,true,*outStream);
      *outStream << "\n";
    }

    ROL::Algorithm<RealT> algo("Augmented Lagrangian",*parlist,false);
    Teuchos::Time algoTimer("Algorithm Time", true);
    algo.run(*zp,*c2p,augLag,*vcon,*bnd,true,*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"density.txt");

    Teuchos::Array<RealT> res(1,0);
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
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
