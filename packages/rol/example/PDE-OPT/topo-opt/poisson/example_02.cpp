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

/*! \file  example_02.cpp
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
#include <fenv.h>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
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
  feenableexcept(FE_DIVBYZERO | FE_INVALID); // | FE_OVERFLOW);

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
    std::string filename = "input_ex02.xml";
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

    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, l_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    z_ptr = assembler->createControlVector();  z_ptr->putScalar(volFraction);
    l_ptr = assembler->createStateVector();    l_ptr->putScalar(0.0);
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > up, zp, lp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    lp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(l_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > x = ROL::makePtr<ROL::Vector_SimOpt<RealT> >(up,zp);

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

    // Initialize volume constraint
    ROL::Ptr<QoI<RealT> > qoi_vol
      = ROL::makePtr<QoI_Volume_Poisson_TopOpt<RealT>>(pde->getFE(),*parlist);
    ROL::Ptr<IntegralConstraint<RealT> > vcon
      = ROL::makePtr<IntegralConstraint<RealT>>(qoi_vol,assembler);
    // Create volume constraint vector and set to zero
    RealT vecScaling = one / std::pow(domainWidth*domainHeight*(one-volFraction), 2);
    ROL::Ptr<std::vector<RealT> > scalevec_ptr = ROL::makePtr<std::vector<RealT>>(1,vecScaling);
    ROL::Ptr<std::vector<RealT> > c1_ptr = ROL::makePtr<std::vector<RealT>>(1,0);
    ROL::Ptr<ROL::Vector<RealT> > c1p = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(c1_ptr, scalevec_ptr);
    ROL::Ptr<std::vector<RealT> > c2_ptr = ROL::makePtr<std::vector<RealT>>(1,1);
    ROL::Ptr<ROL::Vector<RealT> > c2p = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(c2_ptr, scalevec_ptr);

    // Build bound constraint
    ROL::Ptr<Tpetra::MultiVector<> > zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(0.0); zhi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > zlop
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > zhip
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::BoundConstraint<RealT> > zbnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    ROL::Ptr<Tpetra::MultiVector<> > ulo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::Vector<RealT> > ulop
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ulo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > uhip
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uhi_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ubnd
      = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT> >(ubnd,zbnd);

    // Set up optimization problem
    std::vector<ROL::Ptr<ROL::Constraint<RealT> > > econ(2);
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > emul(2);
    econ[0] = pdeWithFilter;
    econ[1] = vcon;
    emul[0] = lp;
    emul[1] = c2p;
    ROL::OptimizationProblem<RealT> problem(obj,x,bnd,econ,emul);

    // Solve optimization problem
    pdeWithFilter->solve(*rp,*up,*zp,tol);
    ROL::OptimizationSolver<RealT> solver(problem,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    solver.solve(*outStream);
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
