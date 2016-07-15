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
    \brief Solves a source inversion problem governed by the
           advection-diffusion equation.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

template<class Real>
Real random(const Teuchos::RCP<const Teuchos::Comm<int> > &commptr) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*commptr)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*commptr,0,1,&val);
  return val;
}

template<class Real>
void randomize(std::vector<Real> &x,
               const Teuchos::RCP<const Teuchos::Comm<int> > &commptr) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = random<Real>(commptr);
  }
}

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm
    = Teuchos::rcp(new Teuchos::SerialComm<int>());
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
    Teuchos::RCP<PoissonData<RealT> > data
      = Teuchos::rcp(new PoissonData<RealT>(serial_comm, parlist, outStream));

    /*** Build vectors and dress them up as ROL vectors. ***/
    const RealT zero(0), one(1);
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
//    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
//    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > p_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > c_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_c, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
//    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp(new std::vector<RealT>(9,one));
    Teuchos::RCP<std::vector<RealT> > dz_rcp = Teuchos::rcp(new std::vector<RealT>(9,zero));
    // Set all values to 1 in u, z and c.
    u_rcp->putScalar(one);
//    z_rcp->putScalar(one);
    p_rcp->putScalar(one);
    c_rcp->putScalar(one);
    // Randomize d vectors.
    du_rcp->randomize();
    //dz_rcp->randomize();
    randomize<RealT>(*dz_rcp,comm);
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
//    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > pp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(p_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > cp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(c_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
//    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::StdVector<RealT>(dz_rcp));
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<RealT> > obj
      = Teuchos::rcp(new Objective_PDEOPT_Poisson<RealT>(data, parlist));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > con
      = Teuchos::rcp(new EqualityConstraint_PDEOPT_Poisson<RealT>(data, parlist));
    Teuchos::RCP<ROL::Reduced_ParametrizedObjective_SimOpt<RealT> > objReduced
      = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(obj, con, up, pp));
    Teuchos::RCP<std::vector<RealT> > zlo_rcp = Teuchos::rcp(new std::vector<RealT>(9,zero));
    Teuchos::RCP<std::vector<RealT> > zup_rcp = Teuchos::rcp(new std::vector<RealT>(9,one));
    Teuchos::RCP<ROL::Vector<RealT> > zlop = Teuchos::rcp(new ROL::StdVector<RealT>(zlo_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zupp = Teuchos::rcp(new ROL::StdVector<RealT>(zup_rcp));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp(new ROL::BoundConstraint<RealT>(zlop,zupp));

    /*** Build sampler ***/
    int sdim  = 37;
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one, one};
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::StdTeuchosBatchManager<RealT,int>(comm));
      //= Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
      //= Teuchos::rcp(new ROL::BatchManager<RealT>());
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,bounds,bman));
    // Build stochastic problem
    ROL::StochasticProblem<RealT> opt(*parlist,objReduced,sampler,zp,bnd);
    opt.setSolutionStatistic(zero);

    bool printMeanValueState = parlist->sublist("Problem").get("Print Mean Value State",false);
    if ( printMeanValueState ) {
      RealT tol = 1.e-8;
      std::vector<RealT> my_sample(sdim), mev_sample(sdim), gev_sample(sdim);
      for (int i = 0; i < sampler->numMySamples(); ++i) {
        my_sample = sampler->getMyPoint(i);
        for (int j = 0; j < sdim; ++j) {
          mev_sample[j] += sampler->getMyWeight(i)*my_sample[j];
        }
      }
      bman->sumAll(&mev_sample[0],&gev_sample[0],sdim);
      con->setParameter(gev_sample);
      zp->zero();
      con->solve(*cp,*up,*zp,tol);
      data->outputTpetraVector(u_rcp, "mean_value_state.txt");
    }

    /*** Check functional interface. ***/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      std::vector<RealT> param(sdim,1);
      objReduced->setParameter(param);
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
    }

    /*** Solve optimization problem. ***/
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    bool useBundle = parlist->sublist("Problem").get("Is problem nonsmooth?",false);
    if ( useBundle ) {
      algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Bundle",*parlist,false));
    }
    else {
      algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    }
    zp->zero(); // set zero initial guess
    std::clock_t timer = std::clock();
    algo->run(opt,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl;

    // Output control to file
    //data->outputTpetraVector(z_rcp, "control.txt");
    std::clock_t timer_print = std::clock();
    if ( myRank == 0 ) {
      std::ofstream zfile;
      zfile.open("control.txt");
      for (int i = 0; i < 9; i++) {
        zfile << (*z_rcp)[i] << "\n";
      }
      zfile.close();
    }

    // Output expected state to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman_Eu
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      con->setParameter(sampler->getMyPoint(i));
      con->solve(*cp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
    }
    bman_Eu->sumAll(*up,*pp);
    data->outputTpetraVector(p_rcp, "mean_state.txt");

    // Build objective function distribution
    RealT val(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
      //= Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,bounds,bman));
    std::vector<RealT> sample(sdim);
    std::stringstream name;
    name << "samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      objReduced->setParameter(sample);
      val = objReduced->value(*zp,tol);
      for (int j = 0; j < sdim; ++j) {
        file << sample[j] << "  ";
      }
      file << val << "\n";
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl;
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
