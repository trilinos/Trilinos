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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"

#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

typedef double RealT;

template<class Real>
class ObjectiveEx11 : public ROL::StdObjective<Real> {
private:
  const Real alpha_;
public:
  ObjectiveEx11(const Real alpha = 1e-4) : alpha_(alpha) {}

  Real value( const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    Real val(0);
    for ( unsigned i = 0; i < size; i++ ) {
      val += static_cast<Real>(0.5)*alpha_*x[i]*x[i] + x[i];
    }
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    for ( unsigned i = 0; i < size; i++ ) {
      g[i] = alpha_*x[i] + static_cast<Real>(1);
    }
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    for ( unsigned i = 0; i < size; i++ ) {
      hv[i] = alpha_*v[i];
    }
  }

  void precond( std::vector<Real> &pv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    for ( unsigned i = 0; i < size; i++ ) {
      pv[i] = v[i]/alpha_;
    }
  }
};

template<class Real>
class ConstraintEx11 : public ROL::StdConstraint<Real> {
private:
  const std::vector<Real> zeta_;

public:
  ConstraintEx11(const std::vector<Real> &zeta) : zeta_(zeta) {}

  void value( std::vector<Real> &c,
              const std::vector<Real> &x,
              Real &tol ) {
    unsigned size = x.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      c[i] = std::exp(param[7])/zeta_[i] - std::exp(param[i])*x[i];
    }
  }

  void applyJacobian( std::vector<Real> &jv,
                      const std::vector<Real> &v,
                      const std::vector<Real> &x,
                      Real &tol ) {
    unsigned size = x.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      jv[i] = -std::exp(param[i])*v[i];
    }
  }

   void applyAdjointJacobian( std::vector<Real> &ajv,
                              const std::vector<Real> &v,
                              const std::vector<Real> &x,
                              Real &tol ) {
    unsigned size = x.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      ajv[i] = -std::exp(param[i])*v[i];
    }
  }

  void applyAdjointHessian( std::vector<Real> &ahuv,
                            const std::vector<Real> &u,
                            const std::vector<Real> &v,
                            const std::vector<Real> &x,
                            Real &tol ) {
    unsigned size = x.size();
    for ( unsigned i = 0; i < size; ++i ) {
      ahuv[i] = static_cast<Real>(0);
    }
  }

  void applyPreconditioner( std::vector<Real> &pv,
                            const std::vector<Real> &v,
                            const std::vector<Real> &x,
                            const std::vector<Real> &g,
                            Real &tol ) {
    unsigned size = x.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      pv[i] = v[i]/std::pow(std::exp(param[i]),2);
    }
  }
};

template<class Real>
Real setUpAndSolve(Teuchos::ParameterList                    &list,
                   Teuchos::RCP<ROL::Objective<Real> >       &obj,
                   Teuchos::RCP<ROL::Vector<Real> >          &x,
                   Teuchos::RCP<ROL::BoundConstraint<Real> > &bnd,
                   Teuchos::RCP<ROL::Constraint<Real> >      &con,
                   Teuchos::RCP<ROL::Vector<Real> >          &mul,
                   Teuchos::RCP<ROL::BoundConstraint<Real> > &ibnd,
                   Teuchos::RCP<ROL::SampleGenerator<Real> > &sampler,
                   std::ostream & outStream) {
  ROL::OptimizationProblem<Real> optProblem(obj,x,bnd,con,mul,ibnd);
  optProblem.setMeanValueInequality(sampler);
  optProblem.check(outStream);
  // Run ROL algorithm
  ROL::OptimizationSolver<Real> optSolver(optProblem, list);
  optSolver.solve(outStream);
  Teuchos::RCP<ROL::Objective<Real> > robj = optProblem.getObjective();
  Real tol(1.e-8);
  return robj->value(*(optProblem.getSolutionVector()),tol);
}

template<class Real>
void printSolution(const std::vector<Real> &x,
                   std::ostream & outStream) {
  unsigned dim = x.size();
  outStream << "x = (";
  for ( unsigned i = 0; i < dim-1; i++ ) {
    outStream << x[i] << ", ";
  }
  outStream << x[dim-1] << ")\n";
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_11.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::ParameterList list = *parlist;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    const RealT zero(0), half(0.5), one(1), two(2);
    unsigned dim = 7;
    Teuchos::RCP<std::vector<RealT> > x_rcp;
    x_rcp = Teuchos::rcp( new std::vector<RealT>(dim,zero) );
    Teuchos::RCP<ROL::Vector<RealT> > x;
    x = Teuchos::rcp(new ROL::StdVector<RealT>(x_rcp));
    // Build samplers
    int nSamp = 50;
    unsigned sdim = dim + 1;
    std::vector<Teuchos::RCP<ROL::Distribution<RealT> > > dist(sdim);
    for ( unsigned i = 0; i < 2; ++i ) {
      //dist[i]
      //  = Teuchos::rcp(new ROL::Gaussian<RealT>(static_cast<RealT>(4.5856),  static_cast<RealT>(0.0392)));
      dist[i]
        = Teuchos::rcp(new ROL::Uniform<RealT>(std::log(80.),std::log(120.)));
    }
    for ( unsigned i = 2; i < sdim-1; ++i ) {
      //dist[i]
      //  = Teuchos::rcp(new ROL::Gaussian<RealT>(static_cast<RealT>(5.2787),  static_cast<RealT>(0.0392)));
      dist[i]
        = Teuchos::rcp(new ROL::Uniform<RealT>(std::log(160.),std::log(240.)));
    }
    //dist[sdim-1]
    //  = Teuchos::rcp(new ROL::Gaussian<RealT>(static_cast<RealT>(13.7413), static_cast<RealT>(0.1484)));
    dist[sdim-1]
      = Teuchos::rcp(new ROL::Uniform<RealT>(std::log(10.),std::log(40.)));
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::BatchManager<RealT>());
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nSamp,dist,bman));
    // Build objective function
    Teuchos::RCP<ROL::Objective<RealT> > obj
      = Teuchos::rcp(new ObjectiveEx11<RealT>);
    // Build bound constraints
    std::vector<RealT> lx(dim,half);
    std::vector<RealT> ux(dim,two);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp( new ROL::StdBoundConstraint<RealT>(lx,ux) );
    // Build inequality constraint
    std::vector<RealT> zeta(dim,one/static_cast<RealT>(std::sqrt(3.)));
    zeta[0] *= half;
    zeta[1] *= half;
    Teuchos::RCP<ROL::Constraint<RealT> > con
      = Teuchos::rcp(new ConstraintEx11<RealT>(zeta));
    // Build inequality bound constraint
    std::vector<RealT> lc(dim,ROL::ROL_NINF<RealT>());
    std::vector<RealT> uc(dim,zero);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > ibnd
      = Teuchos::rcp( new ROL::StdBoundConstraint<RealT>(lc,uc) );
    ibnd->deactivateLower();
    // Build multipliers
    std::vector<RealT> mean(sdim,zero), gmean(sdim,zero);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      std::vector<RealT> pt = sampler->getMyPoint(i);
      RealT              wt = sampler->getMyWeight(i);
      for (unsigned j = 0; j < sdim; ++j) {
        mean[j] += wt*std::exp(pt[j]);
      }
    }
    sampler->sumAll(&mean[0],&gmean[0],sdim);
    Teuchos::RCP<std::vector<RealT> > scaling_vec
      = Teuchos::rcp( new std::vector<RealT>(dim,zero) );
    for (unsigned i = 0; i < dim; ++i) {
      RealT cl = std::abs(gmean[dim]/zeta[i] - gmean[i]*lx[i]);
      RealT cu = std::abs(gmean[dim]/zeta[i] - gmean[i]*ux[i]);
      RealT scale = static_cast<RealT>(1e2)/std::pow(std::max(cl,cu),two);
      (*scaling_vec)[i] = (scale > std::sqrt(ROL::ROL_EPSILON<RealT>()))
                            ? scale : one;
    }
    Teuchos::RCP<std::vector<RealT> > l_rcp
      = Teuchos::rcp( new std::vector<RealT>(dim,one) );
    Teuchos::RCP<ROL::Vector<RealT> > l
      = Teuchos::rcp(new ROL::DualScaledStdVector<RealT>(l_rcp,scaling_vec));

    /**********************************************************************************************/
    /************************* SUPER QUANTILE QUADRANGLE ******************************************/
    /**********************************************************************************************/
    RealT val = setUpAndSolve(list,obj,x,bnd,con,l,ibnd,sampler,*outStream);
    *outStream << "Computed Solution" << std::endl;
    printSolution<RealT>(*x_rcp,*outStream);

    // Compute exact solution
    Teuchos::RCP<std::vector<RealT> > ximean, gximean, xsol;
    ximean  = Teuchos::rcp(new std::vector<RealT>(dim+1,0));
    gximean = Teuchos::rcp(new std::vector<RealT>(dim+1));
    xsol    = Teuchos::rcp(new std::vector<RealT>(dim));
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      std::vector<RealT> pt = sampler->getMyPoint(i);
      RealT              wt = sampler->getMyWeight(i);
      for (unsigned j = 0; j < dim+1; ++j) {
        (*ximean)[j] += wt*pt[j];
      }
    }
    bman->reduceAll(&(*ximean)[0],&(*gximean)[0],dim+1,ROL::Elementwise::ReductionMax<RealT>());
    for (unsigned i = 0; i < dim; ++i) {
      (*xsol)[i] = std::max(half,(std::exp((*gximean)[dim])/std::exp((*gximean)[i]))/zeta[i]);
    }
    Teuchos::RCP<ROL::Vector<RealT> > xtrue
      = Teuchos::rcp(new ROL::StdVector<RealT>(xsol));
    *outStream << "True Solution" << std::endl;
    printSolution<RealT>(*xsol,*outStream);
    xtrue->axpy(-one,*x);

    RealT error = xtrue->norm();
    *outStream << std::right
               << std::setw(20) << "obj val"
               << std::setw(20) << "error"
               << std::endl;
    *outStream << std::fixed << std::setprecision(0) << std::right
               << std::scientific << std::setprecision(11) << std::right
               << std::setw(20) << val
               << std::setw(20) << error
               << std::endl;

    errorFlag = (error < std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 0 : 1;
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
