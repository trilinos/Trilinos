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
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_UnaryFunctions.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<StefanBoltzmannData<RealT> > data = ROL::makePtr<StefanBoltzmannData<RealT>>(comm, parlist, outStream);

    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::Ptr<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > c_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > Eu_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > Vu_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    // Set all values to 1 in u, z and c.
    u_ptr->putScalar(1.0);
    z_ptr->putScalar(1.0);
    c_ptr->putScalar(1.0);
    // Randomize d vectors.
    du_ptr->randomize();
    dz_ptr->randomize();
    // Create ROL::TpetraMultiVectors.
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > cp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dzp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    ROL::Ptr<ROL::Vector<RealT> > Eup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(Eu_ptr);
    ROL::Ptr<ROL::Vector<RealT> > Vup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(Vu_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj =
      ROL::makePtr<Objective_PDEOPT_StefanBoltzmann<RealT>>(data, parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con =
      ROL::makePtr<EqualityConstraint_PDEOPT_StefanBoltzmann<RealT>>(data, parlist);
    con->setSolveParameters(*parlist);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced =
      ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, up);

    /*** Build stochastic functionality. ***/
    int sdim = parlist->sublist("Problem").get("Stochastic Dimension",4);
    // Build sampler
    int nsamp = parlist->sublist("Problem").get("Number of Monte Carlo Samples",100);
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman,false,false,100);
    // Build stochastic problem
    ROL::OptimizationProblem<RealT> opt(objReduced,zp);
    opt.setStochasticObjective(*parlist,sampler);

    /*** Check functional interface. ***/
    std::vector<RealT> par(sdim,1.0);
    obj->setParameter(par);
    con->setParameter(par);
    objReduced->setParameter(par);
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

    /*** Solve optimization problem. ***/
    ROL::Algorithm<RealT> algo_tr("Trust Region",*parlist,false);
    zp->zero(); // set zero initial guess
    algo_tr.run(opt, true, *outStream);

    *outStream << " Solution Statistic: S(z) = " << opt.getSolutionStatistic() << "\n";

    data->outputTpetraVector(z_ptr, "control.txt");

    RealT w = 0.0, tol = 1.e-8;
    ROL::Elementwise::Power<RealT> sqr(2.0);
    Eup->zero(); Vup->zero();
    ROL::Ptr<ROL::Vector<RealT> > up2 = up->clone();
    for (int i = 0; i < sampler->numMySamples(); i++) {
      // Get samples and weights
      par = sampler->getMyPoint(i);
      w   = sampler->getMyWeight(i);
      // Solve state equation at current sample
      con->setParameter(par);
      con->solve(*cp,*up,*zp,tol);
      // Accumulate expected value
      Eup->axpy(w,*up);
      // Accumulate variance
      up2->set(*up); up2->applyUnary(sqr);
      Vup->axpy(w,*up2);
    }
    up2->set(*Eup); up2->applyUnary(sqr);
    Vup->axpy(-1.0,*up2);
    if (sampler->numMySamples() > 1) {
      Vup->scale((RealT)sampler->numMySamples()/(RealT)(sampler->numMySamples()-1));
    }
    data->outputTpetraVector(Eu_ptr, "expected_value_state.txt");
    data->outputTpetraVector(Vu_ptr, "variance_state.txt");
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
