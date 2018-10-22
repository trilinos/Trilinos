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
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Constraint_Partitioned.hpp"
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
#include "obj_compliance.hpp"
#include "obj_volume.hpp"
#include "mesh_topo-opt.hpp"
#include "pde_elasticity.hpp"
#include "pde_filter.hpp"

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
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_TopoOpt<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT> > pde
      = ROL::makePtr<PDE_Elasticity<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conPDE
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT> > pdeFilter
      = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(conPDE);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    conPDE->setSolveParameters(*parlist);

    /*** Initialize vector storage. ***/
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, z_ptr, r_ptr, dp_ptr, dr_ptr, yp_ptr, yr_ptr;
    u_ptr  = assembler->createStateVector();    u_ptr->randomize();
    p_ptr  = assembler->createStateVector();    p_ptr->randomize();
    z_ptr  = assembler->createControlVector();  z_ptr->putScalar(volFraction);
    r_ptr  = assembler->createResidualVector(); r_ptr->randomize();
    dp_ptr = assembler->createStateVector();    dp_ptr->randomize();
    dr_ptr = assembler->createResidualVector(); dr_ptr->randomize();
    yp_ptr = assembler->createStateVector();    yp_ptr->randomize();
    yr_ptr = assembler->createResidualVector(); yr_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, zp, rp, dpp, drp, ypp, yrp;
    up     = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp     = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp     = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp     = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    dpp    = ROL::makePtr<PDE_PrimalSimVector<RealT>>(dp_ptr,pde,assembler,*parlist);
    drp    = ROL::makePtr<PDE_DualSimVector<RealT>>(dr_ptr,pde,assembler,*parlist);
    ypp    = ROL::makePtr<PDE_PrimalSimVector<RealT>>(yp_ptr,pde,assembler,*parlist);
    yrp    = ROL::makePtr<PDE_DualSimVector<RealT>>(yr_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > x;
    x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    /*** Initialize "filtered" or "unfiltered" constraint. ***/
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter  = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(
            conPDE, conFilter, *rp, *rp, *up, *zp, *zp, useStorage);
    }
    else {
      pdeWithFilter = conPDE;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    /*** Initialize volume constraint. ***/
    ROL::Ptr<QoI<RealT> > qoi_vol
      = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getFE(),pde->getFieldHelper(),volFraction);
    ROL::Ptr<IntegralConstraint<RealT> > vcon
      = ROL::makePtr<IntegralConstraint<RealT>>(qoi_vol,assembler);

    /*** Initialize combined constraint. ***/
    std::vector<ROL::Ptr<ROL::Constraint<RealT> > > convec(2);
    convec[0] = pdeWithFilter;
    convec[1] = vcon;
    ROL::Ptr<ROL::Constraint<RealT> > con
      = ROL::makePtr<ROL::Constraint_Partitioned<RealT>>(convec);

    /*** Initialize constraint and multiplier vectors ***/
    RealT vecScaling = one / std::pow(domainWidth*domainHeight*(one-volFraction), 2);
    ROL::Ptr<std::vector<RealT> > scalevec_ptr, volres_ptr, volmult_ptr,
                                      volres1_ptr, volmult1_ptr, volres2_ptr,
                                      volmult2_ptr;
    scalevec_ptr = ROL::makePtr<std::vector<RealT>>(1,vecScaling);
    volres_ptr   = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult_ptr  = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volres1_ptr  = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult1_ptr = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volres2_ptr  = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    volmult2_ptr = ROL::makePtr<std::vector<RealT>>(1,static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX));
    ROL::Ptr<ROL::Vector<RealT> > volres, volmult, volres1, volmult1,
                                      volres2, volmult2;
    volres   = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(volres_ptr, scalevec_ptr);
    volmult  = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(volmult_ptr, scalevec_ptr);
    volres1  = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(volres1_ptr, scalevec_ptr);
    volmult1 = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(volmult1_ptr, scalevec_ptr);
    volres2  = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(volres2_ptr, scalevec_ptr);
    volmult2 = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(volmult2_ptr, scalevec_ptr);
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > multvec(2), multvec1(2),
                                                    resvec(2), resvec1(2),
                                                    multvec2(2), resvec2(2);
    multvec[0]  = pp;  multvec[1]  = volmult;
    multvec1[0] = dpp; multvec1[1] = volmult1;
    multvec2[0] = ypp; multvec2[1] = volmult2;
    resvec[0]   = rp;  resvec[1]   = volres;
    resvec1[0]  = drp; resvec1[1]  = volres1;
    resvec2[0]  = yrp; resvec2[1]  = volres2;
    ROL::Ptr<ROL::Vector<RealT> > multv, multv1, resv, resv1, multv2, resv2;
    multv  = ROL::makePtr<ROL::PartitionedVector<RealT>>(multvec);
    multv1 = ROL::makePtr<ROL::PartitionedVector<RealT>>(multvec1);
    multv2 = ROL::makePtr<ROL::PartitionedVector<RealT>>(multvec2);
    resv   = ROL::makePtr<ROL::PartitionedVector<RealT>>(resvec);
    resv1  = ROL::makePtr<ROL::PartitionedVector<RealT>>(resvec1);
    resv2  = ROL::makePtr<ROL::PartitionedVector<RealT>>(resvec2);

    /*** Initialize compliance objective function. ***/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(1,ROL::nullPtr);
    bool useEnergy = parlist->sublist("Problem").get("Use Energy Objective",false);
    if (useEnergy) {
      RealT objScaling = parlist->sublist("Problem").get("Objective Scaling",1.0);
      qoi_vec[0] = ROL::makePtr<QoI_Energy_TopoOpt<RealT>>(pde->getFE(),
                                                              pde->getMaterialTensor(),
                                                              pde->getFieldHelper(),
                                                              objScaling);
    }
    else {
      Teuchos::ParameterList list(*parlist);
      list.sublist("Vector").sublist("Sim").set("Use Riesz Map",true);
      list.sublist("Vector").sublist("Sim").set("Lump Riesz Map",false);
      // Has state Riesz map enabled for mesh-independent compliance scaling.
      ROL::Ptr<Tpetra::MultiVector<> > f_ptr = assembler->createResidualVector();
      ROL::Ptr<ROL::Vector<RealT> > fp
        = ROL::makePtr<PDE_DualSimVector<RealT>>(f_ptr,pde,assembler,list);
      fp->zero(); up->zero();
      pdeWithFilter->value(*fp, *up, *zp, tol);
      RealT objScaling = objFactor, fnorm2 = fp->dot(*fp);
      if (fnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
        objScaling /= fnorm2;
      }
      u_ptr->randomize();
      qoi_vec[0] = ROL::makePtr<QoI_Compliance_TopoOpt<RealT>>(pde->getFE(),
                                                               pde->getLoad(),
                                                               pde->getFieldHelper(),
                                                               objScaling);
    }
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,assembler);

    /*** Initialize bound constraints. ***/
    // Control bounds
    ROL::Ptr<Tpetra::MultiVector<> > zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(0.0); zhi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > zlop
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > zhip
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > zbnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
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

    // Create optimization problem and solve
    pdeWithFilter->solve(*rp,*up,*zp,tol);
    multv->zero();
    ROL::OptimizationProblem<RealT> optProb(obj,x,bnd,con,multv);
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if (checkDeriv) {
      optProb.check(*outStream);
    }
    ROL::OptimizationSolver<RealT>  optSolver(optProb,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    optSolver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

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
