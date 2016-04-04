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
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

// ROL Vector Inputs
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_ScaledStdVector.hpp"

// Rol Algorithm
#include "ROL_Algorithm.hpp"

// SOL Inputs
#include "ROL_StochasticProblem.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "parametrized_objective.hpp"
#include "parametrized_constraint.hpp"
#include "volume_constraint.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

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
    Teuchos::RCP<ElasticitySIMPOperators<RealT> > data
      = Teuchos::rcp(new ElasticitySIMPOperators<RealT>(comm, parlist, outStream));
    /*** Build vectors and dress them up as ROL vectors. ***/
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getDomainMapA();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getCellMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dw_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz2_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<std::vector<RealT> >    vscale_rcp = Teuchos::rcp(new std::vector<RealT>(1, 0));
    Teuchos::RCP<std::vector<RealT> >    vc_rcp = Teuchos::rcp(new std::vector<RealT>(1, 0));
    Teuchos::RCP<std::vector<RealT> >    vc_lam_rcp = Teuchos::rcp(new std::vector<RealT>(1, 0));
    // Set all values to 1 in u, z.
    u_rcp->putScalar(1.0);
    // Set z to gray solution.
    RealT volFrac = parlist->sublist("ElasticityTopoOpt").get("Volume Fraction", 0.5), one(1), two(2);
    z_rcp->putScalar(volFrac);
    (*vscale_rcp)[0] = one/std::pow(static_cast<RealT>(z_rcp->getGlobalLength())*(one-volFrac),two);

    // Randomize d vectors.
    du_rcp->randomize();
    dw_rcp->randomize();
    dz_rcp->randomize();
    dz2_rcp->randomize();
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dwp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dw_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dz2p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz2_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > vcp
      = Teuchos::rcp(new ROL::PrimalScaledStdVector<RealT>(vc_rcp,vscale_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > vc_lamp
      = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(vc_lam_rcp,vscale_rcp));
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
    ROL::Vector_SimOpt<RealT> d2(dwp,dz2p);

    /*** Build objective function, constraint and reduced objective function. ***/
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<RealT> > obj
       = Teuchos::rcp(new ParametrizedObjective_PDEOPT_ElasticitySIMP<RealT>(data, parlist));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > con
       = Teuchos::rcp(new ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP<RealT>(data, parlist));
    Teuchos::RCP<ROL::Reduced_ParametrizedObjective_SimOpt<RealT> > objReduced
       = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(obj, con, up, dwp));
    Teuchos::RCP<ROL::EqualityConstraint<RealT> > volcon
       = Teuchos::rcp(new EqualityConstraint_PDEOPT_ElasticitySIMP_Volume<RealT>(data, parlist));

    /*** Build bound constraint ***/
    Teuchos::RCP<Tpetra::MultiVector<> > lo_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > hi_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    lo_rcp->putScalar(0.0); hi_rcp->putScalar(1.0);
    Teuchos::RCP<ROL::Vector<RealT> > lop = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(lo_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > hip = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(hi_rcp));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd = Teuchos::rcp(new ROL::BoundConstraint<RealT>(lop,hip));

    /*** Build stochastic functionality. ***/
    int sdim = 2;
    // Build sampler
    int nsamp = stoch_parlist->sublist("Problem").get("Number of Monte Carlo Samples",100);
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(sdim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::BatchManager<RealT>());
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,bounds,bman,false,false,100));
    // Build stochastic problem
    ROL::StochasticProblem<RealT> opt(*parlist,objReduced,sampler,zp,bnd);

    /*** Run optimization ***/
    ROL::AugmentedLagrangian<RealT> augLag(opt.getObjective(),volcon,*vc_lamp,1.0,*opt.getSolutionVector(),*vcp,*parlist);
    ROL::Algorithm<RealT> algo("Augmented Lagrangian",*parlist,false);
    algo.run(*opt.getSolutionVector(),*vc_lamp,augLag,*volcon,*bnd,true,*outStream);

    data->outputTpetraVector(z_rcp, "density.txt");
    data->outputTpetraVector(u_rcp, "state.txt");
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
