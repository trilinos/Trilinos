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
    \brief Shows how to solve the Poisson control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "../TOOLS/meshreader.hpp"
#include "../TOOLS/linearpdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_poisson.hpp"
#include "obj_poisson.hpp"

typedef double RealT;

template<class Real>
class stateSolution : public Solution<Real> {
public:
  stateSolution(void) {}
  Real evaluate(const std::vector<Real> &x, const int fieldNumber) const {
    const Real pi(M_PI);
    Real u(1);
    int dim = x.size();
    for (int i=0; i<dim; ++i) {
      u *= std::sin(pi*x[i]);
    }
    return u;
  }
};

template<class Real>
class controlSolution : public Solution<Real> {
private:
  const Real alpha_;
public:
  controlSolution(const Real alpha) : alpha_(alpha) {}
  Real evaluate(const std::vector<Real> &x, const int fieldNumber) const {
    const Real eight(8), pi(M_PI);
    Real z(1);
    int dim = x.size();
    for (int i=0; i<dim; ++i) {
      z *= std::sin(eight*pi*x[i]);
    }
    Real coeff1(64), coeff2(dim);
    return -z/(alpha_*coeff1*coeff2*pi*pi);
  }
};

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

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    int probDim = parlist->sublist("Problem").get("Problem Dimension",2);
    ROL::SharedPointer<MeshManager<RealT> > meshMgr;
    if (probDim == 1) {
      meshMgr = ROL::makeShared<MeshManager_Interval<RealT>>(*parlist);
    } else if (probDim == 2) {
      meshMgr = ROL::makeShared<MeshManager_Rectangle<RealT>>(*parlist);
    } else if (probDim == 3) {
      meshMgr = ROL::makeShared<MeshReader<RealT>>(*parlist);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/example_01.cpp: Problem dim is not 1, 2 or 3!");
    }
    // Initialize PDE describe Poisson's equation
    ROL::SharedPointer<PDE_Poisson<RealT> > pde
      = ROL::makeShared<PDE_Poisson<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makeShared<Linear_PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::SharedPointer<Linear_PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPointerCast<Linear_PDE_Constraint<RealT> >(con);
    ROL::SharedPointer<Assembler<RealT> > assembler = pdecon->getAssembler();
    // Initialize quadratic objective function
    std::vector<ROL::SharedPointer<QoI<RealT> > > qoi_vec(2,ROL::nullPointer);
    qoi_vec[0] = ROL::makeShared<QoI_L2Tracking_Poisson<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makeShared<QoI_L2Penalty_Poisson<RealT>>(pde->getFE());
    RealT alpha = parlist->sublist("Problem").get("Control penalty parameter",1e-2);
    std::vector<RealT> wt(2); wt[0] = static_cast<RealT>(1); wt[1] = alpha;
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makeShared<PDE_Objective<RealT>>(qoi_vec,wt,assembler);

    // Create state vector and set to zeroes
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
    ROL::SharedPointer<ROL::Vector<RealT> > up, zp, pp, rp;
    u_ptr  = assembler->createStateVector();   u_ptr->putScalar(0.0);
    z_ptr  = assembler->createControlVector(); z_ptr->putScalar(0.0);
    p_ptr  = assembler->createStateVector();   p_ptr->putScalar(0.0);
    r_ptr  = assembler->createStateVector();   r_ptr->putScalar(0.0);
    up  = ROL::makeShared<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    zp  = ROL::makeShared<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);
    pp  = ROL::makeShared<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    rp  = ROL::makeShared<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);

    // Initialize reduced objective function
    ROL::SharedPointer<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);

    // Build optimization problem and check derivatives
    ROL::OptimizationProblem<RealT> optProb(robj,zp);
    optProb.check(*outStream);

    // Build optimization solver and solve
    zp->zero(); up->zero(); pp->zero();
    ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist);
    std::clock_t timerTR = std::clock();
    optSolver.solve(*outStream);
    *outStream << "Trust Region Time: "
               << static_cast<RealT>(std::clock()-timerTR)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    // Compute solution error
    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    ROL::SharedPointer<Solution<RealT> > usol
      = ROL::makeShared<stateSolution<RealT>>();
    RealT uerr = assembler->computeStateError(u_ptr,usol);
    *outStream << "State Error: " << uerr << std::endl;
    ROL::SharedPointer<Solution<RealT> > zsol
      = ROL::makeShared<controlSolution<RealT>>(alpha);
    RealT zerr = assembler->computeControlError(z_ptr,zsol);
    *outStream << "Control Error: " << zerr << std::endl;

    // Output.
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    pdecon->outputTpetraData();
    assembler->printMeshData(*outStream);

    errorFlag += (uerr > 10. ? 1 : 0);

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
