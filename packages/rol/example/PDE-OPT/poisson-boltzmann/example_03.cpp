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

/*! \file  example_03.cpp
    \brief Shows how to solve the Poisson-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_UnaryFunctions.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"

#include "../TOOLS/linearpdeconstraint.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_poisson_boltzmann_ex02.hpp"
#include "mesh_poisson_boltzmann_ex02.hpp"
#include "obj_poisson_boltzmann_ex02.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

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

    /*** Read in XML input ***/
    std::string filename = "input_ex03.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Example02<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Poisson_Boltzmann_ex02<RealT> > pde
      = ROL::makePtr<PDE_Poisson_Boltzmann_ex02<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT> > pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<PDE_Doping<RealT> > pdeDoping
      = ROL::makePtr<PDE_Doping<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conDoping
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeDoping,meshMgr,comm,*parlist,*outStream,true);
    const ROL::Ptr<Assembler<RealT> > assembler = pdeCon->getAssembler();
    assembler->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  r_ptr = assembler->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<> >  z_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, rp, zp;
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    r_ptr->randomize();  //r_ptr->putScalar(static_cast<RealT>(1));
    z_ptr->randomize();  //z_ptr->putScalar(static_cast<RealT>(1));
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::Vector<RealT> > xp
      = ROL::makePtr<ROL::Vector_SimOpt<RealT> >(up,zp);

    /*************************************************************************/
    /***************** BUILD REFERENCE DOPING AND POTENTIAL ******************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> > ru_ptr = assembler->createStateVector();
    ROL::Ptr<ROL::Vector<RealT> > rup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ru_ptr,pde,assembler,*parlist);
    ROL::Ptr<Tpetra::MultiVector<> > rz_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > rzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(rz_ptr,pde,assembler,*parlist);
    ROL::Ptr<Doping<RealT> > dope
      = ROL::makePtr<Doping<RealT>>(pde->getFE(), pde->getCellNodes(),
                                    assembler->getDofManager()->getCellDofs(),
                                    assembler->getCellIds(),*parlist);
    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithDoping
      = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conDoping,
        *rp, *rp, *up, *zp, *zp, true, true);
    pdeWithDoping->setSolveParameters(*parlist);
    dope->build(rz_ptr);
    RealT tol(1.e-8);
    pdeWithDoping->solve(*rp,*rup,*rzp,tol);
    pdeCon->outputTpetraVector(ru_ptr,"reference_state.txt");
    pdeCon->outputTpetraVector(rz_ptr,"reference_control.txt");

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(3,ROL::nullPtr);
    // Current flow over drain
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_1_Poisson_Boltzmann<RealT>>(pde->getFE(),
                                    pde->getBdryFE(),pde->getBdryCellLocIds(),*parlist);
    ROL::Ptr<IntegralObjective<RealT> > stateObj
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    // Deviation from reference doping
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_1_Poisson_Boltzmann<RealT>>(pde->getFE(),dope);
    ROL::Ptr<IntegralObjective<RealT> > ctrlObj1
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);
    // H1-Seminorm of doping
    qoi_vec[2] = ROL::makePtr<QoI_Control_Cost_2_Poisson_Boltzmann<RealT>>(pde->getFE());
    ROL::Ptr<IntegralObjective<RealT> > ctrlObj2
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[2],assembler);
    // Build standard vector objective function
    RealT currentWeight = parlist->sublist("Problem").get("Desired Current Scale",1.5);
    RealT J = stateObj->value(*rup,*rzp,tol); // Reference current flow over drain
    J *= currentWeight; // Increase current flow by 50%
    RealT w1 = parlist->sublist("Problem").get("State Cost Parameter",1e-3);
    RealT w2 = parlist->sublist("Problem").get("Control Misfit Parameter",1e-2);
    RealT w3 = parlist->sublist("Problem").get("Control Cost Parameter",1e-8);
    ROL::Ptr<ROL::StdObjective<RealT> > std_obj
      = ROL::makePtr<StdObjective_Poisson_Boltzmann<RealT>>(J,w1,w2,w3);
    // Build full-space objective
    ROL::Ptr<PDE_Objective<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
 
    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> > zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > zhi_ptr = assembler->createControlVector();
    ROL::Ptr<DopingBounds<RealT> > dopeBnd
      = ROL::makePtr<DopingBounds<RealT>>(pde->getFE(),pde->getCellNodes(),
                                             assembler->getDofManager()->getCellDofs(),
                                             assembler->getCellIds(),*parlist);
    dopeBnd->build(zlo_ptr,zhi_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zlop, zhip;
    zlop = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler);
    zhip = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler);
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
    bool deactivate = parlist->sublist("Problem").get("Deactivate Bound Constraints",false);
    if (deactivate) {
      bnd->deactivate();
    }

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(obj,xp,bnd,pdeWithDoping,pp);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::OptimizationSolver<RealT> solver(opt,*parlist);
    zp->set(*rzp);
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    // Output control to file
    pdeCon->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    pdeCon->outputTpetraVector(u_ptr,"state.txt");
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
