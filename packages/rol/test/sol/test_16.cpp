// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"
#include "ROL_Solver.hpp"

#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_Bounds.hpp"

#include "ROL_StochasticProblem.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"

#include <numeric>
#include <algorithm>

//#include <fenv.h>

typedef double RealT;

template<typename Real>
class LossEx16 : public ROL::StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    std::vector<Real> p = ROL::StdObjective<Real>::getParameter();
    return std::inner_product(x.begin(),x.end(),p.begin(),static_cast<Real>(0));
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    std::vector<Real> p = ROL::StdObjective<Real>::getParameter();
    g.assign(p.begin(),p.end());
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(hv.begin(),hv.end(),static_cast<Real>(0));
  }
};

template<typename Real>
class BudgetEx16 : public ROL::StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    c[0] = std::accumulate(x.begin(),x.end(),static_cast<Real>(0)) - static_cast<Real>(1);
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    jv[0] = std::accumulate(v.begin(),v.end(),static_cast<Real>(0));
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(ajv.begin(),ajv.end(),v[0]);
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    std::fill(ahuv.begin(),ahuv.end(),static_cast<Real>(0));
  }
};

void printSolution(const std::vector<RealT> &x,
                   std::ostream & outStream) {
  size_t dim = x.size();
  outStream << std::endl;
  outStream << "x = (";
  for ( size_t i = 0; i < dim-1; i++ ) {
    outStream << x[i] << ", ";
  }
  outStream << x[dim-1] << ")\n";
  outStream << std::endl;
}

int main(int argc, char* argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  bool print     = iprint > 0;
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
    std::string filename = "input_16.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    ROL::ParameterList list = *parlist;
    list.sublist("General").set("Output Level",print ? 1 : 0);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build vectors
    size_t dim = 10;
    ROL::Ptr<ROL::StdVector<RealT>>
      x = ROL::makePtr<ROL::StdVector<RealT>>(dim,0);
    // Build samplers
    int nsamp = 1000;
    RealT offset(0);
    std::vector<std::vector<RealT>> bounds;
    for (size_t i = 0; i < dim; ++i) {
      RealT li = -static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      RealT ui =  static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      bounds.push_back({li,ui});
      *outStream << "  Asset "        << i << ":"
                 << "  Lower Bound: " << li
                 << "  Upper Bound: " << ui << std::endl;
      offset += ui / static_cast<RealT>(dim);
    }
    offset *= static_cast<RealT>(0.35);
    ROL::Ptr<ROL::BatchManager<RealT>>
      bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT>>
      sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);
    // Build loss function
    ROL::Ptr<ROL::Objective<RealT>>
      loss = ROL::makePtr<LossEx16<RealT>>();
    loss->setParameter(sampler->getMyPoint(0));
    // Build loss constraint
    ROL::Ptr<ROL::Constraint<RealT>>
      losscon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(loss);
    ROL::Ptr<ROL::SingletonVector<RealT>>
      mul_losscon = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    ROL::Ptr<ROL::SingletonVector<RealT>>
      u_losscon = ROL::makePtr<ROL::SingletonVector<RealT>>(offset);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd_losscon = ROL::makePtr<ROL::Bounds<RealT>>(*u_losscon,false);
    // Build budget constraint
    ROL::Ptr<ROL::Constraint<RealT>>
      budget = ROL::makePtr<BudgetEx16<RealT>>();
    ROL::Ptr<ROL::StdVector<RealT>>
      mul_budget = ROL::makePtr<ROL::StdVector<RealT>>(1,0);
    // Build optimization problem
    ROL::Ptr<ROL::StochasticProblem<RealT>>
      problem = ROL::makePtr<ROL::StochasticProblem<RealT>>(loss,x);
    problem->addLinearConstraint("Budget",budget,mul_budget);
    problem->addLinearConstraint("Loss",losscon,mul_losscon,bnd_losscon);

    // Check deterministic problem
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);

    std::vector<std::pair<std::string,std::vector<RealT>>> sol;

    // Markowitz portfolio selection
    problem->edit();
    ROL::ParameterList objlist;
    objlist.sublist("SOL").sublist("Objective").set("Type","Deviation");
    objlist.sublist("SOL").sublist("Objective").sublist("Deviation Measure").set("Name","Variance");
    problem->makeObjectiveStochastic(objlist,sampler);
    ROL::ParameterList conlist;
    conlist.sublist("SOL").sublist("Loss").set("Type","Risk Neutral");
    problem->makeLinearConstraintStochastic("Loss",conlist,sampler,bman);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    ROL::Ptr<ROL::Solver<RealT>>
      solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);
    ROL::Ptr<ROL::Vector<RealT>> xm = x->clone(); xm->set(*x);
    sol.push_back({"Risk Neutral",*x->getVector()});

    // Markwotiz portfolio selection again
    problem->edit();
    problem->resetStochasticLinearConstraint("Loss");
    conlist.sublist("SOL").sublist("Loss").set("Type","Mean Value");
    problem->makeLinearConstraintStochastic("Loss",conlist,sampler);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);
    xm->axpy(-1.0,*x);
    errorFlag += (xm->norm() > 100.0*std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 1 : 0);
    sol.push_back({"Mean Value",*x->getVector()});

    // bPOE constrained portfolio selection
    problem->edit();
    conlist.sublist("SOL").sublist("Loss").set("Type","Risk Averse");
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").set("Name","CVaR");
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.95);
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1e-4);
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
    conlist.sublist("SOL").sublist("Loss").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);
    problem->removeLinearConstraint("Loss");
    problem->addConstraint("Loss",losscon,mul_losscon,bnd_losscon);
    problem->makeConstraintStochastic("Loss",conlist,sampler);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    solver->solve(*outStream);
    errorFlag += (solver->getAlgorithmState()->statusFlag == ROL::EXITSTATUS_CONVERGED ? 0 : 1);
    printSolution(*x->getVector(),*outStream);
    sol.push_back({"bPOE",*x->getVector()});

    *outStream << std::endl << std::scientific << std::setprecision(6);
    for (auto it = sol.begin(); it != sol.end(); ++it) {
      *outStream << "  ";
      *outStream << std::setw(20) << std::left << std::get<0>(*it);
      for (const auto& xi : std::get<1>(*it)) *outStream << std::setw(18) << std::left << xi;
      *outStream << std::endl;
    }
    *outStream << std::endl;

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
