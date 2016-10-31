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
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
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

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    Teuchos::RCP<const Teuchos::Comm<int> > FEMcomm
      = Teuchos::rcp(new Teuchos::SerialComm<int>());
    Teuchos::RCP<StefanBoltzmannData<RealT> > data
      = Teuchos::rcp(new StefanBoltzmannData<RealT>(FEMcomm, parlist, outStream));

    /*** Build vectors and dress them up as ROL vectors. ***/
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > c_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_c, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > pEu_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > pVu_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > Eu_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > Vu_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    // Set all values to 1 in u, z and c.
    u_rcp->putScalar(1.0);
    z_rcp->putScalar(1.0);
    c_rcp->putScalar(1.0);
    // Randomize d vectors.
    du_rcp->randomize();
    dz_rcp->randomize();
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > cp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(c_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > pEup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(pEu_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > pVup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(pVu_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > Eup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(Eu_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > Vup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(Vu_rcp));
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<RealT> > obj =
      Teuchos::rcp(new Objective_PDEOPT_StefanBoltzmann<RealT>(data, parlist));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > con =
      Teuchos::rcp(new EqualityConstraint_PDEOPT_StefanBoltzmann<RealT>(data, parlist));
    Teuchos::RCP<ROL::Reduced_ParametrizedObjective_SimOpt<RealT> > objReduced =
      Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(obj, con, up, up));

    /*** Build stochastic functionality. ***/
    int sdim = parlist->sublist("Problem").get("Stochastic Dimension",4);
    // Build sampler
    int nsamp = parlist->sublist("Problem").get("Number of Monte Carlo Samples",100);
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,bounds,bman,false,false,100));
    // Build stochastic problem
    ROL::StochasticProblem<RealT> opt(*parlist,objReduced,sampler,zp);

    /*** Check functional interface. ***/
    std::vector<RealT> par(sdim,1.0);
/*
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
    opt.checkObjectiveGradient(*dzp,true,*outStream);
    opt.checkObjectiveHessVec(*dzp,true,*outStream);
*/

    /*** Solve optimization problem. ***/
    ROL::Algorithm<RealT> algo_tr("Trust Region",*parlist,false);
    zp->zero(); // set zero initial guess
    algo_tr.run(opt, true, *outStream);

    *outStream << " Solution Statistic: S(z) = " << opt.getSolutionStatistic() << "\n";

    data->outputTpetraVector(z_rcp, "control.txt");

    RealT w = 0.0, tol = 1.e-8;
    ROL::Elementwise::Power<RealT> sqr(2.0);
    pEup->zero(); pVup->zero();
    Teuchos::RCP<ROL::Vector<RealT> > up2 = up->clone();
    for (int i = 0; i < sampler->numMySamples(); i++) {
      // Get samples and weights
      par = sampler->getMyPoint(i);
      w   = sampler->getMyWeight(i);
      // Solve state equation at current sample
      con->setParameter(par);
      con->solve(*cp,*up,*zp,tol);
      // Accumulate expected value
      pEup->axpy(w,*up);
      // Accumulate variance
      up2->set(*up); up2->applyUnary(sqr);
      pVup->axpy(w,*up2);
    }
    sampler->sumAll(*pEup,*Eup);
    sampler->sumAll(*pVup,*Vup);

    up2->set(*Eup); up2->applyUnary(sqr);
    Vup->axpy(-1.0,*up2);
    if (sampler->numMySamples() > 1) {
      Vup->scale((RealT)nsamp/(RealT)(nsamp-1));
    }
    data->outputTpetraVector(Eu_rcp, "expected_value_state.txt");
    data->outputTpetraVector(Vu_rcp, "variance_state.txt");

    // Build new sampler to evaluate random variable objective
    //int nsampCDF = parlist->sublist("Problem").get("Number of Samples for CDF Computation",10*nsamp);
    //Teuchos::RCP<ROL::SampleGenerator<RealT> > samplerCDF
    //  = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsampCDF,bounds,bman,false,false,100));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > samplerCDF = sampler;
    // Temporary vector storage to evaluate objective
    Teuchos::RCP<Tpetra::MultiVector<> > diff_rcp
      = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > Mdiff_rcp
      = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    // Open file
    std::stringstream name;
    name << "objective_samples." << comm->getRank() << ".txt";
    std::ofstream file;
    file.open(name.str().c_str());
    // Sample objective
    std::vector<RealT> fval(1,0.0);
    for (int i = 0; i < samplerCDF->numMySamples(); i++) {
      // Get sample and weight
      par = samplerCDF->getMyPoint(i);
      w   = samplerCDF->getMyWeight(i);
      // Solve PDE at sample
      con->setParameter(par);
      con->solve(*cp,*up,*zp,tol);
      // Compute tracking term
      diff_rcp->scale(1.0,*u_rcp);
      diff_rcp->update(-1.0,*(data->getVecUd()),1.0);
      data->getMatM()->apply(*diff_rcp,*Mdiff_rcp);
      diff_rcp->dot(*Mdiff_rcp,fval);
      fval[0] *= 0.5;
      // Print to file
      file << w << "  " << fval[0] << "\n"; 
    }
    file.close();
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
