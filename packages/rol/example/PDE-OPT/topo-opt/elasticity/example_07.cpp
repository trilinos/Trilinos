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
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = Teuchos::rcp(&std::cout, false);
  }
  else {
    outStream = Teuchos::rcp(&bhs, false);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {
    RealT tol(1e-8), one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex07.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT domainWidth  = parlist->sublist("Geometry").get("Width", 1.0);
    const RealT domainHeight = parlist->sublist("Geometry").get("Height", 1.0);
    const RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_TopoOpt<RealT>(*parlist));
    // Initialize PDE describing elasticity equations.
    Teuchos::RCP<PDE_Elasticity<RealT> > pde
      = Teuchos::rcp(new PDE_Elasticity<RealT>(*parlist));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > conPDE
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,comm,*parlist,*outStream));
    // Initialize the filter PDE.
    Teuchos::RCP<PDE_Filter<RealT> > pdeFilter
      = Teuchos::rcp(new PDE_Filter<RealT>(*parlist));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > conFilter
      = Teuchos::rcp(new Linear_PDE_Constraint<RealT>(pdeFilter,meshMgr,comm,*parlist,*outStream));
    // Cast the constraint and get the assembler.
    Teuchos::RCP<PDE_Constraint<RealT> > pdecon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(conPDE);
    Teuchos::RCP<Assembler<RealT> > assembler = pdecon->getAssembler();
    conPDE->setSolveParameters(*parlist);

    /*** Initialize vector storage. ***/
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp, du_rcp, p_rcp, z_rcp, dz_rcp,
                                         rz_rcp, r_rcp, dualu_rcp, dualz_rcp,
                                         dp_rcp, dr_rcp, yp_rcp, yr_rcp;
    u_rcp  = assembler->createStateVector();    u_rcp->randomize();
    du_rcp = assembler->createStateVector();    du_rcp->randomize();
    p_rcp  = assembler->createStateVector();    p_rcp->randomize();
    z_rcp  = assembler->createControlVector();  z_rcp->putScalar(volFraction);
    dz_rcp = assembler->createControlVector();  dz_rcp->randomize(); dz_rcp->scale(0.01);
    rz_rcp = assembler->createControlVector();  rz_rcp->randomize();
    r_rcp  = assembler->createResidualVector(); r_rcp->randomize();
    dp_rcp = assembler->createStateVector();    dp_rcp->randomize();
    dr_rcp = assembler->createResidualVector(); dr_rcp->randomize();
    yp_rcp = assembler->createStateVector();    yp_rcp->randomize();
    yr_rcp = assembler->createResidualVector(); yr_rcp->randomize();
    dualu_rcp = assembler->createStateVector();
    dualz_rcp = assembler->createControlVector();
    Teuchos::RCP<ROL::Vector<RealT> > up, dup, pp, zp, dzp,
                                      rzp, rp, dualup, dualzp,
                                      dpp, drp, ypp, yrp;
    up     = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler,*parlist));
    dup    = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,assembler,*parlist));
    pp     = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler,*parlist));
    zp     = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(z_rcp,pde,assembler,*parlist));
    dzp    = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(dz_rcp,pde,assembler,*parlist));
    rzp    = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(rz_rcp,pde,assembler,*parlist));
    rp     = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler,*parlist));
    dpp    = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(dp_rcp,pde,assembler,*parlist));
    drp    = Teuchos::rcp(new PDE_DualSimVector<RealT>(dr_rcp,pde,assembler,*parlist));
    ypp    = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(yp_rcp,pde,assembler,*parlist));
    yrp    = Teuchos::rcp(new PDE_DualSimVector<RealT>(yr_rcp,pde,assembler,*parlist));
    dualup = Teuchos::rcp(new PDE_DualSimVector<RealT>(dualu_rcp,pde,assembler,*parlist));
    dualzp = Teuchos::rcp(new PDE_DualOptVector<RealT>(dualz_rcp,pde,assembler,*parlist));
    Teuchos::RCP<ROL::Vector<RealT> > x, d;
    x = Teuchos::rcp(new ROL::Vector_SimOpt<RealT>(up,zp));
    d = Teuchos::rcp(new ROL::Vector_SimOpt<RealT>(dup,dzp));

    /*** Initialize "filtered" or "unfiltered" constraint. ***/
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter  = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = Teuchos::rcp(new ROL::CompositeConstraint_SimOpt<RealT>(
            conPDE, conFilter, *rp, *rp, *up, *zp, *zp, useStorage));
    }
    else {
      pdeWithFilter = conPDE;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    /*** Initialize volume constraint. ***/
    Teuchos::RCP<QoI<RealT> > qoi_vol
      = Teuchos::rcp(new QoI_Volume_TopoOpt<RealT>(pde->getFE(),pde->getFieldHelper(),*parlist));
    Teuchos::RCP<IntegralConstraint<RealT> > vcon
      = Teuchos::rcp(new IntegralConstraint<RealT>(qoi_vol,assembler));

    /*** Initialize combined constraint. ***/
    std::vector<Teuchos::RCP<ROL::Constraint<RealT> > > convec(2);
    convec[0] = pdeWithFilter;
    convec[1] = vcon;
    Teuchos::RCP<ROL::Constraint<RealT> > con
      = Teuchos::rcp(new ROL::Constraint_Partitioned<RealT>(convec));

    /*** Initialize constraint and multiplier vectors ***/
    RealT vecScaling = one / std::pow(domainWidth*domainHeight*(one-volFraction), 2);
    Teuchos::RCP<std::vector<RealT> > scalevec_rcp, volres_rcp, volmult_rcp,
                                      volres1_rcp, volmult1_rcp, volres2_rcp,
                                      volmult2_rcp;
    scalevec_rcp = Teuchos::rcp(new std::vector<RealT>(1,vecScaling));
    volres_rcp   = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    volmult_rcp  = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    volres1_rcp  = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    volmult1_rcp = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    volres2_rcp  = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    volmult2_rcp = Teuchos::rcp(new std::vector<RealT>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX)));
    Teuchos::RCP<ROL::Vector<RealT> > volres, volmult, volres1, volmult1,
                                      volres2, volmult2;
    volres   = Teuchos::rcp(new ROL::PrimalScaledStdVector<RealT>(volres_rcp, scalevec_rcp));
    volmult  = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(volmult_rcp, scalevec_rcp));
    volres1  = Teuchos::rcp(new ROL::PrimalScaledStdVector<RealT>(volres1_rcp, scalevec_rcp));
    volmult1 = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(volmult1_rcp, scalevec_rcp));
    volres2  = Teuchos::rcp(new ROL::PrimalScaledStdVector<RealT>(volres2_rcp, scalevec_rcp));
    volmult2 = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(volmult2_rcp, scalevec_rcp));
    std::vector<Teuchos::RCP<ROL::Vector<RealT> > > multvec(2), multvec1(2),
                                                    resvec(2), resvec1(2),
                                                    multvec2(2), resvec2(2);
    multvec[0]  = pp;  multvec[1]  = volmult;
    multvec1[0] = dpp; multvec1[1] = volmult1;
    multvec2[0] = ypp; multvec2[1] = volmult2;
    resvec[0]   = rp;  resvec[1]   = volres;
    resvec1[0]  = drp; resvec1[1]  = volres1;
    resvec2[0]  = yrp; resvec2[1]  = volres2;
    Teuchos::RCP<ROL::Vector<RealT> > multv, multv1, resv, resv1, multv2, resv2;
    multv  = Teuchos::rcp(new ROL::PartitionedVector<RealT>(multvec));
    multv1 = Teuchos::rcp(new ROL::PartitionedVector<RealT>(multvec1));
    multv2 = Teuchos::rcp(new ROL::PartitionedVector<RealT>(multvec2));
    resv   = Teuchos::rcp(new ROL::PartitionedVector<RealT>(resvec));
    resv1  = Teuchos::rcp(new ROL::PartitionedVector<RealT>(resvec1));
    resv2  = Teuchos::rcp(new ROL::PartitionedVector<RealT>(resvec2));

    /*** Initialize compliance objective function. ***/
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(1,Teuchos::null);
    bool useEnergy = parlist->sublist("Problem").get("Use Energy Objective",false);
    if (useEnergy) {
      RealT objScaling = parlist->sublist("Problem").get("Objective Scaling",1.0);
      qoi_vec[0] = Teuchos::rcp(new QoI_Energy_TopoOpt<RealT>(pde->getFE(),
                                                              pde->getMaterialTensor(),
                                                              pde->getFieldHelper(),
                                                              objScaling));
    }
    else {
      Teuchos::ParameterList list(*parlist);
      list.sublist("Vector").sublist("Sim").set("Use Riesz Map",true);
      list.sublist("Vector").sublist("Sim").set("Lump Riesz Map",false);
      // Has state Riesz map enabled for mesh-independent compliance scaling.
      Teuchos::RCP<Tpetra::MultiVector<> > f_rcp = assembler->createResidualVector();
      Teuchos::RCP<ROL::Vector<RealT> > fp
        = Teuchos::rcp(new PDE_DualSimVector<RealT>(f_rcp,pde,assembler,list));
      fp->zero(); up->zero();
      pdeWithFilter->value(*fp, *up, *zp, tol);
      RealT objScaling = objFactor, fnorm2 = fp->dot(*fp);
      if (fnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
        objScaling /= fnorm2;
      }
      u_rcp->randomize();
      qoi_vec[0] = Teuchos::rcp(new QoI_TopoOpt<RealT>(pde->getFE(),
                                                       pde->getLoad(),
                                                       pde->getFieldHelper(),
                                                       objScaling));
    }
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,assembler));

    /*** Initialize bound constraints. ***/
    // Control bounds
    Teuchos::RCP<Tpetra::MultiVector<> > zlo_rcp = assembler->createControlVector();
    Teuchos::RCP<Tpetra::MultiVector<> > zhi_rcp = assembler->createControlVector();
    zlo_rcp->putScalar(0.0); zhi_rcp->putScalar(1.0);
    Teuchos::RCP<ROL::Vector<RealT> > zlop
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(zlo_rcp,pde,assembler));
    Teuchos::RCP<ROL::Vector<RealT> > zhip
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(zhi_rcp,pde,assembler));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > zbnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(zlop,zhip));
    // State bounds
    Teuchos::RCP<Tpetra::MultiVector<> > ulo_rcp = assembler->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> > uhi_rcp = assembler->createStateVector();
    ulo_rcp->putScalar(ROL::ROL_NINF<RealT>()); uhi_rcp->putScalar(ROL::ROL_INF<RealT>());
    Teuchos::RCP<ROL::Vector<RealT> > ulop
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(ulo_rcp,pde,assembler));
    Teuchos::RCP<ROL::Vector<RealT> > uhip
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(uhi_rcp,pde,assembler));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > ubnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(ulop,uhip));
    ubnd->deactivate();
    // SimOpt bounds
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp(new ROL::BoundConstraint_SimOpt<RealT>(ubnd,zbnd));

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
    pdecon->outputTpetraVector(u_rcp,"state.txt");
    pdecon->outputTpetraVector(z_rcp,"density.txt");

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
