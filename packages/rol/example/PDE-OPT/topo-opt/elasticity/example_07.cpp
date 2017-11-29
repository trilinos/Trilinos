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
    \brief Shows how to solve the stuctural topology optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Reduced_Constraint_SimOpt.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "obj_topo-opt.hpp"
#include "mesh_topo-opt.hpp"
#include "pde_elasticity.hpp"
#include "pde_filter.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makeSharedFromRef(std::cout);
  }
  else {
    outStream = ROL::makeSharedFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {
    RealT tol(1e-8), one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex07.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT domainWidth  = parlist->sublist("Geometry").get("Width", 1.0);
    const RealT domainHeight = parlist->sublist("Geometry").get("Height", 1.0);
    const RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    ROL::SharedPointer<MeshManager<RealT> > meshMgr
      = ROL::makeShared<MeshManager_TopoOpt<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::SharedPointer<PDE_Elasticity<RealT> > pde
      = ROL::makeShared<PDE_Elasticity<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > conPDE
      = ROL::makeShared<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::SharedPointer<PDE_Filter<RealT> > pdeFilter
      = ROL::makeShared<PDE_Filter<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makeShared<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::SharedPointer<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPointerCast<PDE_Constraint<RealT> >(conPDE);
    ROL::SharedPointer<Assembler<RealT> > assembler = pdecon->getAssembler();
    conPDE->setSolveParameters(*parlist);

    /*** Initialize vector storage. ***/
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr, du_ptr, p_ptr, z_ptr, dz_ptr,
                                         rz_ptr, r_ptr, dualu_ptr, dualz_ptr,
                                         dp_ptr, dr_ptr, yp_ptr, yr_ptr;
    u_ptr  = assembler->createStateVector();    u_ptr->randomize();
    du_ptr = assembler->createStateVector();    du_ptr->randomize();
    p_ptr  = assembler->createStateVector();    p_ptr->randomize();
    z_ptr  = assembler->createControlVector();  z_ptr->putScalar(volFraction);
    dz_ptr = assembler->createControlVector();  dz_ptr->randomize(); dz_ptr->scale(0.01);
    rz_ptr = assembler->createControlVector();  rz_ptr->randomize();
    r_ptr  = assembler->createResidualVector(); r_ptr->randomize();
    dp_ptr = assembler->createStateVector();    dp_ptr->randomize();
    dr_ptr = assembler->createResidualVector(); dr_ptr->randomize();
    yp_ptr = assembler->createStateVector();    yp_ptr->randomize();
    yr_ptr = assembler->createResidualVector(); yr_ptr->randomize();
    dualu_ptr = assembler->createStateVector();
    dualz_ptr = assembler->createControlVector();
    ROL::SharedPointer<ROL::Vector<RealT> > up, dup, pp, zp, dzp,
                                      rzp, rp, dualup, dualzp,
                                      dpp, drp, ypp, yrp;
    up     = ROL::makeShared<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    dup    = ROL::makeShared<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    pp     = ROL::makeShared<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp     = ROL::makeShared<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    dzp    = ROL::makeShared<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);
    rzp    = ROL::makeShared<PDE_PrimalOptVector<RealT>>(rz_ptr,pde,assembler,*parlist);
    rp     = ROL::makeShared<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    dpp    = ROL::makeShared<PDE_PrimalSimVector<RealT>>(dp_ptr,pde,assembler,*parlist);
    drp    = ROL::makeShared<PDE_DualSimVector<RealT>>(dr_ptr,pde,assembler,*parlist);
    ypp    = ROL::makeShared<PDE_PrimalSimVector<RealT>>(yp_ptr,pde,assembler,*parlist);
    yrp    = ROL::makeShared<PDE_DualSimVector<RealT>>(yr_ptr,pde,assembler,*parlist);
    dualup = ROL::makeShared<PDE_DualSimVector<RealT>>(dualu_ptr,pde,assembler,*parlist);
    dualzp = ROL::makeShared<PDE_DualOptVector<RealT>>(dualz_ptr,pde,assembler,*parlist);
    ROL::SharedPointer<ROL::Vector<RealT> > x, d;
    x = ROL::makeShared<ROL::Vector_SimOpt<RealT>>(up,zp);
    d = ROL::makeShared<ROL::Vector_SimOpt<RealT>>(dup,dzp);

    /*** Initialize "filtered" or "unfiltered" constraint. ***/
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter  = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = ROL::makeShared<ROL::CompositeConstraint_SimOpt<RealT>>(
            conPDE, conFilter, *rp, *rp, *up, *zp, *zp, useStorage);
    }
    else {
      pdeWithFilter = conPDE;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    /*** Initialize volume constraint. ***/
    ROL::SharedPointer<QoI<RealT> > qoi_vol
      = ROL::makeShared<QoI_Volume_TopoOpt<RealT>>(pde->getFE(),pde->getFieldHelper(),*parlist);
    ROL::SharedPointer<IntegralConstraint<RealT> > vcon
      = ROL::makeShared<IntegralConstraint<RealT>>(qoi_vol,assembler);

    /*** Initialize combined constraint. ***/
    std::vector<ROL::SharedPointer<ROL::Constraint<RealT> > > convec(2);
    convec[0] = pdeWithFilter;
    convec[1] = vcon;
    ROL::SharedPointer<ROL::Constraint<RealT> > con
      = ROL::makeShared<ROL::Constraint_Partitioned<RealT>>(convec);

    /*** Initialize constraint and multiplier vectors ***/
    RealT vecScaling = one / std::pow(domainWidth*domainHeight*(one-volFraction), 2);
    ROL::SharedPointer<std::vector<RealT> > scalevec_ptr, volres_ptr, volmult_ptr,
                                      volres1_ptr, volmult1_ptr, volres2_ptr,
                                      volmult2_ptr;
    scalevec_ptr = ROL::makeShared<std::vector<RealT>>(1,vecScaling);
    volres_ptr   = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult_ptr  = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volres1_ptr  = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult1_ptr = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volres2_ptr  = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult2_ptr = ROL::makeShared<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    ROL::SharedPointer<ROL::Vector<RealT> > volres, volmult, volres1, volmult1,
                                      volres2, volmult2;
    volres   = ROL::makeShared<ROL::PrimalScaledStdVector<RealT>>(volres_ptr, scalevec_ptr);
    volmult  = ROL::makeShared<ROL::DualScaledStdVector<RealT>>(volmult_ptr, scalevec_ptr);
    volres1  = ROL::makeShared<ROL::PrimalScaledStdVector<RealT>>(volres1_ptr, scalevec_ptr);
    volmult1 = ROL::makeShared<ROL::DualScaledStdVector<RealT>>(volmult1_ptr, scalevec_ptr);
    volres2  = ROL::makeShared<ROL::PrimalScaledStdVector<RealT>>(volres2_ptr, scalevec_ptr);
    volmult2 = ROL::makeShared<ROL::DualScaledStdVector<RealT>>(volmult2_ptr, scalevec_ptr);
    std::vector<ROL::SharedPointer<ROL::Vector<RealT> > > multvec(2), multvec1(2),
                                                    resvec(2), resvec1(2),
                                                    multvec2(2), resvec2(2);
    multvec[0]  = pp;  multvec[1]  = volmult;
    multvec1[0] = dpp; multvec1[1] = volmult1;
    multvec2[0] = ypp; multvec2[1] = volmult2;
    resvec[0]   = rp;  resvec[1]   = volres;
    resvec1[0]  = drp; resvec1[1]  = volres1;
    resvec2[0]  = yrp; resvec2[1]  = volres2;
    ROL::SharedPointer<ROL::Vector<RealT> > multv, multv1, resv, resv1, multv2, resv2;
    multv  = ROL::makeShared<ROL::PartitionedVector<RealT>>(multvec);
    multv1 = ROL::makeShared<ROL::PartitionedVector<RealT>>(multvec1);
    multv2 = ROL::makeShared<ROL::PartitionedVector<RealT>>(multvec2);
    resv   = ROL::makeShared<ROL::PartitionedVector<RealT>>(resvec);
    resv1  = ROL::makeShared<ROL::PartitionedVector<RealT>>(resvec1);
    resv2  = ROL::makeShared<ROL::PartitionedVector<RealT>>(resvec2);

    /*** Initialize compliance objective function. ***/
    std::vector<ROL::SharedPointer<QoI<RealT> > > qoi_vec(1,ROL::nullPointer);
    bool useEnergy = parlist->sublist("Problem").get("Use Energy Objective",false);
    if (useEnergy) {
      RealT objScaling = parlist->sublist("Problem").get("Objective Scaling",1.0);
      qoi_vec[0] = ROL::makeShared<QoI_Energy_TopoOpt<RealT>>(pde->getFE(),
                                                              pde->getMaterialTensor(),
                                                              pde->getFieldHelper(),
                                                              objScaling);
    }
    else {
      Teuchos::ParameterList list(*parlist);
      list.sublist("Vector").sublist("Sim").set("Use Riesz Map",true);
      list.sublist("Vector").sublist("Sim").set("Lump Riesz Map",false);
      // Has state Riesz map enabled for mesh-independent compliance scaling.
      ROL::SharedPointer<Tpetra::MultiVector<> > f_ptr = assembler->createResidualVector();
      ROL::SharedPointer<ROL::Vector<RealT> > fp
        = ROL::makeShared<PDE_DualSimVector<RealT>>(f_ptr,pde,assembler,list);
      fp->zero(); up->zero();
      pdeWithFilter->value(*fp, *up, *zp, tol);
      RealT objScaling = objFactor, fnorm2 = fp->dot(*fp);
      if (fnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
        objScaling /= fnorm2;
      }
      u_ptr->randomize();
      qoi_vec[0] = ROL::makeShared<QoI_TopoOpt<RealT>>(pde->getFE(),
                                                       pde->getLoad(),
                                                       pde->getFieldHelper(),
                                                       objScaling);
    }
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makeShared<PDE_Objective<RealT>>(qoi_vec,assembler);

    /*** Initialize bound constraints. ***/
    // Control bounds
    ROL::SharedPointer<Tpetra::MultiVector<> > zlo_ptr = assembler->createControlVector();
    ROL::SharedPointer<Tpetra::MultiVector<> > zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(0.0); zhi_ptr->putScalar(1.0);
    ROL::SharedPointer<ROL::Vector<RealT> > zlop
      = ROL::makeShared<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler);
    ROL::SharedPointer<ROL::Vector<RealT> > zhip
      = ROL::makeShared<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > zbnd
      = ROL::makeShared<ROL::Bounds<RealT>>(zlop,zhip);
    // State bounds
    ROL::SharedPointer<Tpetra::MultiVector<> > ulo_ptr = assembler->createStateVector();
    ROL::SharedPointer<Tpetra::MultiVector<> > uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::SharedPointer<ROL::Vector<RealT> > ulop
      = ROL::makeShared<PDE_PrimalSimVector<RealT>>(ulo_ptr,pde,assembler);
    ROL::SharedPointer<ROL::Vector<RealT> > uhip
      = ROL::makeShared<PDE_PrimalSimVector<RealT>>(uhi_ptr,pde,assembler);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > ubnd
      = ROL::makeShared<ROL::Bounds<RealT>>(ulop,uhip);
    ubnd->deactivate();
    // SimOpt bounds
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > bnd
      = ROL::makeShared<ROL::BoundConstraint_SimOpt<RealT>>(ubnd,zbnd);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Primal Constraint Vector\n";
      resv->checkVector(*resv1,*resv2,true,*outStream);

      *outStream << "\n\nCheck Dual Constraint Vector\n";
      multv->checkVector(*multv1,*multv2,true,*outStream);

      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(*x,*d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of Full Objective Function\n";
      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Full Objective Function\n";
      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(*x,*d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Full Objective Function\n";
      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of Full Objective Function\n";
      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of Full Objective Function\n";
      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Full Objective Function\n";
      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(*x,*d,*resv,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(*x,*multv,*d,*x,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*resv,*d,*x,true,*outStream);
      *outStream << "\n";
    }

    bool useAugLag = parlist->sublist("Problem").get("Use Augmented Lagrangian", false);
    pdeWithFilter->solve(*rp,*up,*zp,tol);
    multv->zero();
    if (useAugLag) {
      /*** Solve using Augmented Lagrangian. ***/
      ROL::AugmentedLagrangian<RealT> augLag(obj,con,*multv,1,
                                             *x,*resv,*parlist);
      ROL::Algorithm<RealT> algo("Augmented Lagrangian",*parlist,false);
      Teuchos::Time algoTimer("Algorithm Time", true);
      algo.run(*x,*multv,augLag,*con,*bnd,true,*outStream);
      algoTimer.stop();
      *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";
    }
    else {
      /*** Solve using Moreau-Yosida. ***/
      ROL::MoreauYosidaPenalty<RealT> myPen(obj,bnd,*x,*parlist);
      ROL::Algorithm<RealT> algo("Moreau-Yosida Penalty",*parlist,false);
      Teuchos::Time algoTimer("Algorithm Time", true);
      algo.run(*x,*multv,myPen,*con,*bnd,true,*outStream);
      algoTimer.stop();
      *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";
    }

    // Output.
    pdecon->printMeshData(*outStream);
    pdeWithFilter->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"density.txt");

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
