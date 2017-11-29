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
#include "ROL_Bounds.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

template<class Real>
Real random(const ROL::SharedPointer<const Teuchos::Comm<int> > &commptr) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*commptr)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*commptr,0,1,&val);
  return val;
}

template<class Real>
void randomize(std::vector<Real> &x,
               const ROL::SharedPointer<const Teuchos::Comm<int> > &commptr) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = random<Real>(commptr);
  }
}

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  ROL::SharedPointer<const Teuchos::Comm<int> > serial_comm
    = ROL::makeShared<Teuchos::SerialComm<int>>();
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
    ROL::SharedPointer<PoissonData<RealT> > data
      = ROL::makeShared<PoissonData<RealT>>(serial_comm, parlist, outStream);

    /*** Build vectors and dress them up as ROL vectors. ***/
    const RealT zero(0), one(1);
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
//    ROL::SharedPointer<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
//    ROL::SharedPointer<Tpetra::MultiVector<> > z_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > p_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > c_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > du_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
//    ROL::SharedPointer<Tpetra::MultiVector<> > dz_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::SharedPointer<std::vector<RealT> > z_ptr  = ROL::makeShared<std::vector<RealT>>(9,one);
    ROL::SharedPointer<std::vector<RealT> > dz_ptr = ROL::makeShared<std::vector<RealT>>(9,zero);
    // Set all values to 1 in u, z and c.
    u_ptr->putScalar(one);
//    z_ptr->putScalar(one);
    p_ptr->putScalar(one);
    c_ptr->putScalar(one);
    // Randomize d vectors.
    du_ptr->randomize();
    //dz_ptr->randomize();
    randomize<RealT>(*dz_ptr,comm);
    // Create ROL::TpetraMultiVectors.
    ROL::SharedPointer<ROL::Vector<RealT> > up = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(u_ptr);
//    ROL::SharedPointer<ROL::Vector<RealT> > zp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > pp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > cp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dup = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(du_ptr);
//    ROL::SharedPointer<ROL::Vector<RealT> > dzp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zp = ROL::makeShared<ROL::StdVector<RealT>>(z_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dzp = ROL::makeShared<ROL::StdVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makeShared<Objective_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makeShared<EqualityConstraint_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::SharedPointer<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);
    ROL::SharedPointer<std::vector<RealT> > zlo_ptr = ROL::makeShared<std::vector<RealT>>(9,zero);
    ROL::SharedPointer<std::vector<RealT> > zup_ptr = ROL::makeShared<std::vector<RealT>>(9,one);
    ROL::SharedPointer<ROL::Vector<RealT> > zlop = ROL::makeShared<ROL::StdVector<RealT>>(zlo_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zupp = ROL::makeShared<ROL::StdVector<RealT>>(zup_ptr);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > bnd
      = ROL::makeShared<ROL::Bounds<RealT>>(zlop,zupp);

    /*** Build sampler ***/
    int sdim  = 37;
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one, one};
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    ROL::SharedPointer<ROL::BatchManager<RealT> > bman
      = ROL::makeShared<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
      //= ROL::makeShared<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
      //= ROL::makeShared<ROL::BatchManager<RealT>>();
    ROL::SharedPointer<ROL::SampleGenerator<RealT> > sampler
      = ROL::makeShared<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);
    // Build stochastic problem
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic",zero);
    opt.setStochasticObjective(*parlist,sampler);

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
      data->outputTpetraVector(u_ptr, "mean_value_state.txt");
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
      opt.check(*outStream);
    }

    /*** Solve optimization problem. ***/
    ROL::SharedPointer<ROL::Algorithm<RealT> > algo;
    bool useBundle = parlist->sublist("Problem").get("Is problem nonsmooth?",false);
    if ( useBundle ) {
      algo = ROL::makeShared<ROL::Algorithm<RealT>>("Bundle",*parlist,false);
    }
    else {
      algo = ROL::makeShared<ROL::Algorithm<RealT>>("Trust Region",*parlist,false);
    }
    zp->zero(); // set zero initial guess
    std::clock_t timer = std::clock();
    algo->run(opt,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl;

    // Output control to file
    //data->outputTpetraVector(z_ptr, "control.txt");
    std::clock_t timer_print = std::clock();
    if ( myRank == 0 ) {
      std::ofstream zfile;
      zfile.open("control.txt");
      for (int i = 0; i < 9; i++) {
        zfile << (*z_ptr)[i] << "\n";
      }
      zfile.close();
    }

    // Output expected state to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    ROL::SharedPointer<ROL::BatchManager<RealT> > bman_Eu
      = ROL::makeShared<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    std::vector<RealT> sample(sdim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*cp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < sdim; ++j) {
        file_samp << sample[j] << "  ";
      }
      file_samp << "\n";
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    data->outputTpetraVector(p_ptr, "mean_state.txt");

    // Build objective function distribution
    RealT val(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
      //= ROL::makeShared<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::SharedPointer<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makeShared<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
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
