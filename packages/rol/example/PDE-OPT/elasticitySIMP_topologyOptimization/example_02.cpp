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

#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"
#include "ROL_StochasticProblem.hpp"

#include "ROL_ScaledTpetraMultiVector.hpp"
#include "ROL_ScaledStdVector.hpp"

#include "ROL_BoundConstraint.hpp"
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

Teuchos::RCP<Tpetra::MultiVector<> > createTpetraVector(const Teuchos::RCP<const Tpetra::Map<> > &map) {
  return Teuchos::rcp(new Tpetra::MultiVector<>(map, 1, true));
}

int main(int argc, char *argv[]) {
//  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    std::string stoch_filename = "stochastic.xml";
    Teuchos::RCP<Teuchos::ParameterList> stoch_parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( stoch_filename, stoch_parlist.ptr() );

    /*** Initialize main data structure. ***/
    Teuchos::RCP<const Teuchos::Comm<int> > serial_comm = Teuchos::rcp(new Teuchos::SerialComm<int>());
    Teuchos::RCP<ElasticitySIMPOperators<RealT> > data
      = Teuchos::rcp(new ElasticitySIMPOperators<RealT>(serial_comm, parlist, outStream));
    /*** Initialize density filter. ***/
    Teuchos::RCP<DensityFilter<RealT> > filter
      = Teuchos::rcp(new DensityFilter<RealT>(serial_comm, parlist, outStream));
    /*** Build vectors and dress them up as ROL vectors. ***/
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getDomainMapA();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getCellMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp      = createTpetraVector(vecmap_u);
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp      = createTpetraVector(vecmap_z);
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp     = createTpetraVector(vecmap_u);
    Teuchos::RCP<Tpetra::MultiVector<> > dw_rcp     = createTpetraVector(vecmap_u);
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp     = createTpetraVector(vecmap_z);
    Teuchos::RCP<Tpetra::MultiVector<> > dz2_rcp    = createTpetraVector(vecmap_z);
    Teuchos::RCP<std::vector<RealT> >    vc_rcp     = Teuchos::rcp(new std::vector<RealT>(1, 0));
    Teuchos::RCP<std::vector<RealT> >    vc_lam_rcp = Teuchos::rcp(new std::vector<RealT>(1, 0));
    Teuchos::RCP<std::vector<RealT> >    vscale_rcp = Teuchos::rcp(new std::vector<RealT>(1, 0));
    // Set all values to 1 in u, z.
    RealT one(1), two(2);
    u_rcp->putScalar(one);
    // Set z to gray solution.
    RealT volFrac = parlist->sublist("ElasticityTopoOpt").get<RealT>("Volume Fraction");
    z_rcp->putScalar(volFrac);
    // Set scaling vector for constraint
    RealT W = parlist->sublist("Geometry").get<RealT>("Width");
    RealT H = parlist->sublist("Geometry").get<RealT>("Height");
    (*vscale_rcp)[0] = one/std::pow(W*H*(one-volFrac),two);
    // Set Scaling vector for density
    bool  useZscale = parlist->sublist("Problem").get<bool>("Use Scaled Density Vectors");
    RealT densityScaling = parlist->sublist("Problem").get<RealT>("Density Scaling");
    Teuchos::RCP<Tpetra::MultiVector<> > scaleVec = createTpetraVector(vecmap_z);
    scaleVec->putScalar(densityScaling);
    if ( !useZscale ) {
      scaleVec->putScalar(one);
    }
    Teuchos::RCP<const Tpetra::Vector<> > zscale_rcp = scaleVec->getVector(0);

    // Randomize d vectors.
    du_rcp->randomize(); //du_rcp->scale(0);
    dw_rcp->randomize();
    dz_rcp->randomize(); //dz_rcp->scale(0);
    dz2_rcp->randomize();
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dwp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dw_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp 
      = Teuchos::rcp(new ROL::PrimalScaledTpetraMultiVector<RealT>(z_rcp,zscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new ROL::PrimalScaledTpetraMultiVector<RealT>(dz_rcp,zscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dz2p
      = Teuchos::rcp(new ROL::PrimalScaledTpetraMultiVector<RealT>(dz2_rcp,zscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > vcp
      = Teuchos::rcp(new ROL::PrimalScaledStdVector<RealT>(vc_rcp,vscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > vc_lamp
      = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(vc_lam_rcp,vscale_rcp));
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
    buildSampler.getBatchManager()->reduceAll(&min,&gmin,ROLmin);
    ROL::Elementwise::ReductionMax<RealT> ROLmax;
    buildSampler.getBatchManager()->reduceAll(&max,&gmax,ROLmax);
    buildSampler.getBatchManager()->sumAll(&sum,&gsum,1);
    bool useExpValScale = stoch_parlist->sublist("Problem").get("Use Expected Value Scaling",false);
    RealT scale = (useExpValScale ? gsum : gmin);

    /*** Build objective function, constraint and reduced objective function. ***/
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<RealT> > obj
       = Teuchos::rcp(new ParametrizedObjective_PDEOPT_ElasticitySIMP<RealT>(data, filter, parlist,scale));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > con
       = Teuchos::rcp(new ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP<RealT>(data, filter, parlist));
    Teuchos::RCP<ROL::Reduced_ParametrizedObjective_SimOpt<RealT> > objReduced
       = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(obj, con, up, dwp));
    Teuchos::RCP<ROL::EqualityConstraint<RealT> > volcon
       = Teuchos::rcp(new EqualityConstraint_PDEOPT_ElasticitySIMP_Volume<RealT>(data, parlist));

    /*** Build bound constraint ***/
    Teuchos::RCP<Tpetra::MultiVector<> > lo_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > hi_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    lo_rcp->putScalar(0.0); hi_rcp->putScalar(1.0);
    Teuchos::RCP<ROL::Vector<RealT> > lop
      = Teuchos::rcp(new ROL::PrimalScaledTpetraMultiVector<RealT>(lo_rcp, zscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > hip
      = Teuchos::rcp(new ROL::PrimalScaledTpetraMultiVector<RealT>(hi_rcp, zscale_rcp));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd = Teuchos::rcp(new ROL::BoundConstraint<RealT>(lop,hip));

    /*** Build Stochastic Functionality. ***/
    ROL::StochasticProblem<RealT> opt(*stoch_parlist,objReduced,buildSampler.get(),zp,bnd);

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
      opt.checkObjectiveGradient(*dzp,true,*outStream);
      opt.checkObjectiveHessVec(*dzp,true,*outStream);
      *outStream << "Checking Volume Constraint:" << "\n";
      volcon->checkAdjointConsistencyJacobian(*vcp,*dzp,*zp,true,*outStream);
      volcon->checkApplyJacobian(*zp,*dzp,*vcp,true,*outStream);
      volcon->checkApplyAdjointHessian(*zp,*vcp,*dzp,*zp,true,*outStream);
    }

    /*** Run optimization ***/
    ROL::AugmentedLagrangian<RealT> augLag(opt.getObjective(),volcon,*vc_lamp,1.0,
                                          *opt.getSolutionVector(),*vcp,*parlist);
    ROL::Algorithm<RealT> algo("Augmented Lagrangian",*parlist,false);
    algo.run(*opt.getSolutionVector(),*vc_lamp,augLag,*volcon,*opt.getBoundConstraint(),true,*outStream);    

    data->outputTpetraVector(z_rcp, "density.txt");
    data->outputTpetraVector(u_rcp, "state.txt");
    data->outputTpetraVector(zscale_rcp, "weights.txt");
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
