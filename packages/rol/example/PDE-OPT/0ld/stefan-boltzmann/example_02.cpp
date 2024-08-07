// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_SROMGenerator.hpp"
#include "ROL_DistributionFactory.hpp"
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
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
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
    filename = "SROMinput.xml";
    Teuchos::RCP<Teuchos::ParameterList> SROMlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, SROMlist.ptr() );

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
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced =
      ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, up);

    /*** Build stochastic functionality. ***/
    int sdim = parlist->sublist("Problem").get("Stochastic Dimension",4);
    // Build batch manager
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::BatchManager<RealT>>();
    // Build sampler
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(sdim);
    Teuchos::ParameterList distList;
    distList.sublist("SOL").sublist("Distribution").set("Name","Uniform");
    distList.sublist("SOL").sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    distList.sublist("SOL").sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    for (int d = 0; d < sdim; d++) {
      distVec[d] = ROL::DistributionFactory<RealT>(distList);
    }
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::SROMGenerator<RealT>>(*SROMlist,bman,distVec);
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
