// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
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
#include "ROL_UserInputGenerator.hpp"

typedef double RealT;

template<class Real>
class ObjectiveEx10 : public ROL::StdObjective<Real> {
private:
  const Real alpha_;
public:
  ObjectiveEx10(const Real alpha = 1e-4) : alpha_(alpha) {}

  Real value( const std::vector<Real> &x, Real &tol ) {
    unsigned size = x.size();
    Real val(0);
    for ( unsigned i = 0; i < size; i++ ) {
      val += static_cast<Real>(0.5)*alpha_*x[i]*x[i] + x[i];
    }
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    unsigned size = g.size();
    for ( unsigned i = 0; i < size; i++ ) {
      g[i] = alpha_*x[i] + static_cast<Real>(1);
    }
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    unsigned size = hv.size();
    for ( unsigned i = 0; i < size; i++ ) {
      hv[i] = alpha_*v[i];
    }
  }

  void precond( std::vector<Real> &pv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    unsigned size = pv.size();
    for ( unsigned i = 0; i < size; i++ ) {
      pv[i] = v[i]/alpha_;
    }
  }
};

template<class Real>
class ConstraintEx10 : public ROL::StdConstraint<Real> {
private:
  const std::vector<Real> zeta_;

public:
  ConstraintEx10(const std::vector<Real> &zeta) : zeta_(zeta) {}

  void value( std::vector<Real> &c,
              const std::vector<Real> &x,
              Real &tol ) {
    unsigned size = c.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      c[i] = std::exp(param[7])/zeta_[i] - std::exp(param[i])*x[i];
    }
  }

  void applyJacobian( std::vector<Real> &jv,
                      const std::vector<Real> &v,
                      const std::vector<Real> &x,
                      Real &tol ) {
    unsigned size = jv.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      jv[i] = -std::exp(param[i])*v[i];
    }
  }

   void applyAdjointJacobian( std::vector<Real> &ajv,
                              const std::vector<Real> &v,
                              const std::vector<Real> &x,
                              Real &tol ) {
    unsigned size = ajv.size();
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
    unsigned size = ahuv.size();
    for ( unsigned i = 0; i < size; ++i ) {
      ahuv[i] = static_cast<Real>(0);
    }
  }

  void applyPreconditioner( std::vector<Real> &pv,
                            const std::vector<Real> &v,
                            const std::vector<Real> &x,
                            const std::vector<Real> &g,
                            Real &tol ) {
    unsigned size = pv.size();
    const std::vector<Real> param = ROL::Constraint<Real>::getParameter();
    for ( unsigned i = 0; i < size; ++i ) {
      pv[i] = v[i]/std::pow(std::exp(param[i]),2);
    }
  }
};

template<class Real>
Real setUpAndSolve(ROL::ParameterList                    &list,
                   ROL::Ptr<ROL::Objective<Real> >       &obj,
                   ROL::Ptr<ROL::Vector<Real> >          &x,
                   ROL::Ptr<ROL::BoundConstraint<Real> > &bnd,
                   ROL::Ptr<ROL::Constraint<Real> >      &con,
                   ROL::Ptr<ROL::Vector<Real> >          &mul,
                   ROL::Ptr<ROL::BoundConstraint<Real> > &ibnd,
                   ROL::Ptr<ROL::SampleGenerator<Real> > &sampler,
                   std::ostream & outStream) {
  ROL::OptimizationProblem<Real> optProblem(obj,x,bnd,con,mul,ibnd);
  optProblem.setAlmostSureInequality(sampler);
  optProblem.check(outStream);
  // Run ROL algorithm
  ROL::OptimizationSolver<Real> optSolver(optProblem, list);
  optSolver.solve(outStream);
  ROL::Ptr<ROL::Objective<Real> > robj = optProblem.getObjective();
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
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input_10.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    const RealT zero(0), half(0.5), one(1), two(2);
    unsigned dim = 7;
    ROL::Ptr<std::vector<RealT> > x_ptr;
    x_ptr = ROL::makePtr<std::vector<RealT>>(dim,zero);
    ROL::Ptr<ROL::Vector<RealT> > x;
    x = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    // Build samplers
    int nSamp = 50;
    unsigned sdim = dim + 1;
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points.txt","weights.txt",nSamp,sdim,bman);
    // Build objective function
    ROL::Ptr<ROL::Objective<RealT> > obj
      = ROL::makePtr<ObjectiveEx10<RealT>>();
    // Build bound constraints
    std::vector<RealT> lx(dim,half);
    std::vector<RealT> ux(dim,two);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::StdBoundConstraint<RealT>>(lx,ux);
    // Build inequality constraint
    std::vector<RealT> zeta(dim,one/static_cast<RealT>(std::sqrt(3.)));
    zeta[0] *= half;
    zeta[1] *= half;
    ROL::Ptr<ROL::Constraint<RealT> > con
      = ROL::makePtr<ConstraintEx10<RealT>>(zeta);
    // Build inequality bound constraint
    std::vector<RealT> lc(dim,ROL::ROL_NINF<RealT>());
    std::vector<RealT> uc(dim,zero);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ibnd
      = ROL::makePtr<ROL::StdBoundConstraint<RealT>>(lc,uc);
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
    ROL::Ptr<std::vector<RealT> > scaling_vec
      = ROL::makePtr<std::vector<RealT>>(dim,zero);
    for (unsigned i = 0; i < dim; ++i) {
      RealT cl = std::abs(gmean[dim]/zeta[i] - gmean[i]*lx[i]);
      RealT cu = std::abs(gmean[dim]/zeta[i] - gmean[i]*ux[i]);
      RealT scale = static_cast<RealT>(1e2)/std::pow(std::max(cl,cu),two);
      (*scaling_vec)[i] = (scale > std::sqrt(ROL::ROL_EPSILON<RealT>()))
                            ? scale : one;
    }
    ROL::Ptr<std::vector<RealT> > l_ptr
      = ROL::makePtr<std::vector<RealT>>(dim,one);
    ROL::Ptr<ROL::Vector<RealT> > l
      = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(l_ptr,scaling_vec);

    /**********************************************************************************************/
    /************************* SUPER QUANTILE QUADRANGLE ******************************************/
    /**********************************************************************************************/
    RealT val = setUpAndSolve(list,obj,x,bnd,con,l,ibnd,sampler,*outStream);
    *outStream << "Computed Solution" << std::endl;
    printSolution<RealT>(*x_ptr,*outStream);

    // Compute exact solution
    ROL::Ptr<std::vector<RealT> > xmax, gmax;
    xmax = ROL::makePtr<std::vector<RealT>>(dim,half);
    gmax = ROL::makePtr<std::vector<RealT>>(dim);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      std::vector<RealT> pt = sampler->getMyPoint(i);
      for (unsigned j = 0; j < dim; ++j) {
        (*xmax)[j] = std::max( (*xmax)[j],
                               std::exp(pt[dim])/(zeta[j]*std::exp(pt[j])) );
      }
    }
    bman->reduceAll(&(*xmax)[0],&(*gmax)[0],dim,ROL::Elementwise::ReductionMax<RealT>());
    ROL::Ptr<ROL::Vector<RealT> > xtrue
      = ROL::makePtr<ROL::StdVector<RealT>>(gmax);
    *outStream << "True Solution" << std::endl;
    printSolution<RealT>(*gmax,*outStream);
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
