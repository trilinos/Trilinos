// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Solves a source inversion problem governed by the
           advection-diffusion equation.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

template<class Real>
Real random(const ROL::Ptr<const Teuchos::Comm<int>> &commptr) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*commptr)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*commptr,0,1,&val);
  return val;
}

template<class Real>
void randomize(std::vector<Real> &x,
               const ROL::Ptr<const Teuchos::Comm<int>> &commptr) {
  unsigned dim = x.size();
  for ( unsigned i = 0; i < dim; i++ ) {
    x[i] = random<Real>(commptr);
  }
}

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int>> serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();
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
    ROL::Ptr<PoissonData<RealT>> data
      = ROL::makePtr<PoissonData<RealT>>(serial_comm, parlist, outStream);

    /*** Build vectors and dress them up as ROL vectors. ***/
    const RealT zero(0), one(1);
    ROL::Ptr<const Tpetra::Map<>> vecmap_u = data->getMatA()->getDomainMap();
//    ROL::Ptr<const Tpetra::Map<>> vecmap_z = data->getMatB()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<>> vecmap_c = data->getMatA()->getRangeMap();
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
//    ROL::Ptr<Tpetra::MultiVector<>> z_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<>> p_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<>> c_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::Ptr<Tpetra::MultiVector<>> du_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
//    ROL::Ptr<Tpetra::MultiVector<>> dz_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<std::vector<RealT>> z_ptr  = ROL::makePtr<std::vector<RealT>>(9,one);
    ROL::Ptr<std::vector<RealT>> dz_ptr = ROL::makePtr<std::vector<RealT>>(9,zero);
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
    ROL::Ptr<ROL::Vector<RealT>> up = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
//    ROL::Ptr<ROL::Vector<RealT>> zp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT>> pp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::Ptr<ROL::Vector<RealT>> cp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::Ptr<ROL::Vector<RealT>> dup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
//    ROL::Ptr<ROL::Vector<RealT>> dzp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    ROL::Ptr<ROL::Vector<RealT>> zp = ROL::makePtr<ROL::StdVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT>> dzp = ROL::makePtr<ROL::StdVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj
      = ROL::makePtr<Objective_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con
      = ROL::makePtr<EqualityConstraint_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);
    ROL::Ptr<std::vector<RealT>> zlo_ptr = ROL::makePtr<std::vector<RealT>>(9,zero);
    ROL::Ptr<std::vector<RealT>> zup_ptr = ROL::makePtr<std::vector<RealT>>(9,one);
    ROL::Ptr<ROL::Vector<RealT>> zlop = ROL::makePtr<ROL::StdVector<RealT>>(zlo_ptr);
    ROL::Ptr<ROL::Vector<RealT>> zupp = ROL::makePtr<ROL::StdVector<RealT>>(zup_ptr);
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zupp);

    /*** Build sampler ***/
    int sdim  = 37;
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one, one};
    std::vector<std::vector<RealT>> bounds(sdim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
      //= ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
      //= ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);
    // Build stochastic problem
    ROL::Ptr<ROL::StochasticProblem<RealT>>
      opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objReduced,zp);
    opt->addBoundConstraint(bnd);
    parlist->sublist("SOL").sublist("Objective").set("Initial Statistic",zero);
    opt->makeObjectiveStochastic(*parlist,sampler);

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
      opt->check(true,*outStream);
    }

    /*** Solve optimization problem. ***/
    bool useBundle = parlist->sublist("Problem").get("Is problem nonsmooth?",false);
    if ( useBundle ) parlist.sublist("Step").set("Type","Bundle");
    else             parlist.sublist("Step").set("Type","Trust Region");{
    ROL::Solver<RealT> solver(opt,*parlist);
    zp->zero(); // set zero initial guess
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
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
    ROL::Ptr<ROL::BatchManager<RealT>> bman_Eu
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
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
      //= ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
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
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
