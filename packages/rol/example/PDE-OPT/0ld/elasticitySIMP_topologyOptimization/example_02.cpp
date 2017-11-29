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

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationProblem.hpp"

#include "ROL_ScaledTpetraMultiVector.hpp"
#include "ROL_ScaledStdVector.hpp"

#include "ROL_Bounds.hpp"
#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Elementwise_Reduce.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "filter.hpp"
#include "parametrized_objective.hpp"
#include "parametrized_constraint.hpp"
#include "volume_constraint.hpp"
#include "build_sampler.hpp"

#include <fenv.h>

typedef double RealT;

ROL::SharedPointer<Tpetra::MultiVector<> > createTpetraVector(const ROL::SharedPointer<const Tpetra::Map<> > &map) {
  return ROL::makeShared<Tpetra::MultiVector<>>(map, 1, true);
}

int main(int argc, char *argv[]) {
//  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    std::string stoch_filename = "stochastic.xml";
    ROL::SharedPointer<Teuchos::ParameterList> stoch_parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( stoch_filename, stoch_parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::SharedPointer<const Teuchos::Comm<int> > serial_comm = ROL::makeShared<Teuchos::SerialComm<int>>();
    ROL::SharedPointer<ElasticitySIMPOperators<RealT> > data
      = ROL::makeShared<ElasticitySIMPOperators<RealT>>(serial_comm, parlist, outStream);
    /*** Initialize density filter. ***/
    ROL::SharedPointer<DensityFilter<RealT> > filter
      = ROL::makeShared<DensityFilter<RealT>>(serial_comm, parlist, outStream);
    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_u = data->getDomainMapA();
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_z = data->getCellMap();
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr      = createTpetraVector(vecmap_u);
    ROL::SharedPointer<Tpetra::MultiVector<> > z_ptr      = createTpetraVector(vecmap_z);
    ROL::SharedPointer<Tpetra::MultiVector<> > du_ptr     = createTpetraVector(vecmap_u);
    ROL::SharedPointer<Tpetra::MultiVector<> > dw_ptr     = createTpetraVector(vecmap_u);
    ROL::SharedPointer<Tpetra::MultiVector<> > dz_ptr     = createTpetraVector(vecmap_z);
    ROL::SharedPointer<Tpetra::MultiVector<> > dz2_ptr    = createTpetraVector(vecmap_z);
    ROL::SharedPointer<std::vector<RealT> >    vc_ptr     = ROL::makeShared<std::vector<RealT>>(1, 0);
    ROL::SharedPointer<std::vector<RealT> >    vc_lam_ptr = ROL::makeShared<std::vector<RealT>>(1, 0);
    ROL::SharedPointer<std::vector<RealT> >    vscale_ptr = ROL::makeShared<std::vector<RealT>>(1, 0);
    // Set all values to 1 in u, z.
    RealT one(1), two(2);
    u_ptr->putScalar(one);
    // Set z to gray solution.
    RealT volFrac = parlist->sublist("ElasticityTopoOpt").get<RealT>("Volume Fraction");
    z_ptr->putScalar(volFrac);
    // Set scaling vector for constraint
    RealT W = parlist->sublist("Geometry").get<RealT>("Width");
    RealT H = parlist->sublist("Geometry").get<RealT>("Height");
    (*vscale_ptr)[0] = one/std::pow(W*H*(one-volFrac),two);
    // Set Scaling vector for density
    bool  useZscale = parlist->sublist("Problem").get<bool>("Use Scaled Density Vectors");
    RealT densityScaling = parlist->sublist("Problem").get<RealT>("Density Scaling");
    ROL::SharedPointer<Tpetra::MultiVector<> > scaleVec = createTpetraVector(vecmap_z);
    scaleVec->putScalar(densityScaling);
    if ( !useZscale ) {
      scaleVec->putScalar(one);
    }
    ROL::SharedPointer<const Tpetra::Vector<> > zscale_ptr = scaleVec->getVector(0);

    // Randomize d vectors.
    du_ptr->randomize(); //du_ptr->scale(0);
    dw_ptr->randomize();
    dz_ptr->randomize(); //dz_ptr->scale(0);
    dz2_ptr->randomize();
    // Create ROL::TpetraMultiVectors.
    ROL::SharedPointer<ROL::Vector<RealT> > up
      = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dup
      = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dwp
      = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(dw_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zp 
      = ROL::makeShared<ROL::PrimalScaledTpetraMultiVector<RealT>>(z_ptr,zscale_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dzp
      = ROL::makeShared<ROL::PrimalScaledTpetraMultiVector<RealT>>(dz_ptr,zscale_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dz2p
      = ROL::makeShared<ROL::PrimalScaledTpetraMultiVector<RealT>>(dz2_ptr,zscale_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > vcp
      = ROL::makeShared<ROL::PrimalScaledStdVector<RealT>>(vc_ptr,vscale_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > vc_lamp
      = ROL::makeShared<ROL::DualScaledStdVector<RealT>>(vc_lam_ptr,vscale_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
    ROL::Vector_SimOpt<RealT> d2(dwp,dz2p);

    /*** Build sampler. ***/
    BuildSampler<RealT> buildSampler(comm,*stoch_parlist,*parlist);
    buildSampler.print("samples");

    /*** Compute compliance objective function scaling. ***/
    RealT min(ROL::ROL_INF<RealT>()), gmin(0), max(0), gmax(0), sum(0), gsum(0), tmp(0);
    Teuchos::Array<RealT> dotF(1, 0);
    RealT minDensity = parlist->sublist("ElasticitySIMP").get<RealT>("Minimum Density");
    for (int i = 0; i < buildSampler.get()->numMySamples(); ++i) {
      data->updateF(buildSampler.get()->getMyPoint(i));
      (data->getVecF())->dot(*(data->getVecF()),dotF.view(0,1));
      tmp = minDensity/dotF[0];
      min = ((min < tmp) ? min : tmp);
      max = ((max > tmp) ? max : tmp);
      sum += buildSampler.get()->getMyWeight(i) * tmp;
    }
    ROL::Elementwise::ReductionMin<RealT> ROLmin;
    buildSampler.getBatchManager()->reduceAll(&min,&gmin,1,ROLmin);
    ROL::Elementwise::ReductionMax<RealT> ROLmax;
    buildSampler.getBatchManager()->reduceAll(&max,&gmax,1,ROLmax);
    buildSampler.getBatchManager()->sumAll(&sum,&gsum,1);
    bool useExpValScale = stoch_parlist->sublist("Problem").get("Use Expected Value Scaling",false);
    RealT scale = (useExpValScale ? gsum : gmin);

    /*** Build objective function, constraint and reduced objective function. ***/
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj
       = ROL::makeShared<ParametrizedObjective_PDEOPT_ElasticitySIMP<RealT>>(data, filter, parlist,scale);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con
       = ROL::makeShared<ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP<RealT>>(data, filter, parlist);
    ROL::SharedPointer<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
       = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, dwp);
    ROL::SharedPointer<ROL::Constraint<RealT> > volcon
       = ROL::makeShared<EqualityConstraint_PDEOPT_ElasticitySIMP_Volume<RealT>>(data, parlist);

    /*** Build bound constraint ***/
    ROL::SharedPointer<Tpetra::MultiVector<> > lo_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > hi_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    ROL::SharedPointer<ROL::Vector<RealT> > lop
      = ROL::makeShared<ROL::PrimalScaledTpetraMultiVector<RealT>>(lo_ptr, zscale_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > hip
      = ROL::makeShared<ROL::PrimalScaledTpetraMultiVector<RealT>>(hi_ptr, zscale_ptr);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > bnd = ROL::makeShared<ROL::Bounds<RealT>>(lop,hip);

    /*** Build Stochastic Functionality. ***/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    opt.setStochasticObjective(*stoch_parlist,buildSampler.get());

    /*** Check functional interface. ***/
    bool checkDeriv = parlist->sublist("Problem").get("Derivative Check",false);
    if ( checkDeriv ) {
//      *outStream << "Checking Objective:" << "\n";
//      obj->checkGradient(x,d,true,*outStream);
//      obj->checkHessVec(x,d,true,*outStream);
//      obj->checkHessSym(x,d,d2,true,*outStream);
//      *outStream << "Checking Constraint:" << "\n";
//      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
//      con->checkAdjointConsistencyJacobian_1(*dwp, *dup, *up, *zp, true, *outStream);
//      con->checkAdjointConsistencyJacobian_2(*dwp, *dzp, *up, *zp, true, *outStream);
//      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
//      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
//      con->checkApplyJacobian(x,d,*up,true,*outStream);
//      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      *outStream << "Checking Reduced Objective:" << "\n";
      opt.check(*outStream);
      *outStream << "Checking Volume Constraint:" << "\n";
      volcon->checkAdjointConsistencyJacobian(*vcp,*dzp,*zp,true,*outStream);
      volcon->checkApplyJacobian(*zp,*dzp,*vcp,true,*outStream);
      volcon->checkApplyAdjointHessian(*zp,*vcp,*dzp,*zp,true,*outStream);
    }

    /*** Run optimization ***/
    ROL::AugmentedLagrangian<RealT> augLag(opt.getObjective(),volcon,*vc_lamp,1.0,
                                          *opt.getSolutionVector(),*vcp,*parlist);
    ROL::Algorithm<RealT> algo("Augmented Lagrangian",*parlist,false);
    std::clock_t timer = std::clock();
    algo.run(*opt.getSolutionVector(),*vc_lamp,augLag,*volcon,*opt.getBoundConstraint(),true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl;

    data->outputTpetraVector(z_ptr, "density.txt");
    data->outputTpetraVector(u_ptr, "state.txt");
    data->outputTpetraVector(zscale_ptr, "weights.txt");

    // Build objective function distribution
    RealT val(0);
    int nSamp = stoch_parlist->sublist("Problem").get("Number of Output Samples",10);
    stoch_parlist->sublist("Problem").set("Number of Samples",nSamp);
    BuildSampler<RealT> buildSampler_dist(comm,*stoch_parlist,*parlist);
    std::stringstream name;
    name << "samples_" << buildSampler_dist.getBatchManager()->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    std::vector<RealT> sample;
    RealT tol = 1.e-8;
    std::clock_t timer_print = std::clock();
    for (int i = 0; i < buildSampler_dist.get()->numMySamples(); ++i) {
      sample = buildSampler_dist.get()->getMyPoint(i);
      objReduced->setParameter(sample);
      val = objReduced->value(*zp,tol);
      for (int j = 0; j < static_cast<int>(sample.size()); ++j) {
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
