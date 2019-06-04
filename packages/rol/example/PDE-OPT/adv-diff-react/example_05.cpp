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
    \brief Shows how to solve the stochastic advection-diffusion problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_stoch_adv_diff.hpp"
#include "obj_stoch_adv_diff.hpp"
#include "mesh_stoch_adv_diff.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
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
    const int controlDim = 9;

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_stoch_adv_diff<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<PDE_stoch_adv_diff<RealT> > pde
      = ROL::makePtr<PDE_stoch_adv_diff<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT> > pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    const ROL::Ptr<Assembler<RealT> > assembler = pdeCon->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, r_ptr;
    ROL::Ptr<std::vector<RealT> > z_ptr;
    ROL::Ptr<ROL::Vector<RealT> > up, pp, rp, zp, xp;
    // Create state vectors
    u_ptr  = assembler->createStateVector(); u_ptr->putScalar(0.0);
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    p_ptr  = assembler->createStateVector(); p_ptr->putScalar(0.0);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    // Create residual vector
    r_ptr  = assembler->createResidualVector();
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create control vectors
    z_ptr  = ROL::makePtr<std::vector<RealT>>(controlDim,0);
    zp  = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(z_ptr));
    // Create SimOpt vector
    xp  = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_stoch_adv_diff<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_stoch_adv_diff<RealT>>();
    ROL::Ptr<StdObjective_stoch_adv_diff<RealT> > std_obj
      = ROL::makePtr<StdObjective_stoch_adv_diff<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT> > zlo_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,0);
    ROL::Ptr<std::vector<RealT> > zhi_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,1);
    ROL::Ptr<ROL::Vector<RealT> > zlop, zhip;
    zlop = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(zlo_ptr));
    zhip = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(zhi_ptr));
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

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(obj,xp,bnd,con,pp);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if (checkDeriv) {
      opt.check(*outStream);
    }
    ROL::OptimizationSolver<RealT>  solver(opt,*parlist);
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    // Output control to file
    if ( myRank == 0 ) {
      std::stringstream zname;
      zname << "control.txt";
      std::ofstream zfile;
      zfile.open(zname.str());
      for (int i = 0; i < controlDim; i++) {
        zfile << std::scientific << std::setprecision(15)
              << std::setw(25) << (*z_ptr)[i]
              << std::endl;
      }
      zfile.close();
    }
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
