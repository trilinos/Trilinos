// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize a quadratic functional with linear constraints
*/

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Solver.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

template<typename Real>
class ROL_TestObjective : public ROL::StdObjective<Real> {
public:
  ROL_TestObjective() : numVars(2) {
    x0.resize(2*numVars);
    for (int i = 0; i < 2*numVars; ++i)
      x0[i] = static_cast<Real>(-1) + static_cast<Real>(2)*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX); 
  }

  Real value(const std::vector<Real> &x, Real &tol) final override {
    Real val(0);
    for (int i = 0; i < 2*numVars; ++i)
      val += static_cast<Real>(0.5)*(x[i]-x0[i]) * (x[i]-x0[i]);
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) final override {
    for (int i = 0; i < 2*numVars; ++i)
      g[i] = (x[i] - x0[i]);
  }

private:
  const int numVars;
  std::vector<Real> x0;
};

template<typename Real>
class ROL_LinearConstraint : public ROL::StdConstraint<Real> {
public:
  /*
  ROL_LinearConstraint() : V({0.5376671395461000,  1.8338850145950865, -2.2588468610036481, 0.8621733203681206,
                              0.3187652398589808, -1.3076882963052734, -0.4335920223056836, 0.3426244665386499}),
  */
  ROL_LinearConstraint() : V({0.0567007, 0.0567007, 0.0567007, 0.0567007,
                              0.0403800, 0.0700000, 0.0403800, 0.0403800}),
                           numConstr(4), numVars(2) {}

  void value(std::vector<Real> &c, const std::vector<Real> &x, Real &tol) final override {
    for (int ic = 0; ic < numConstr; ++ic) {
      c[ic] = static_cast<Real>(0);
      c[ic+numConstr] = static_cast<Real>(0);
      for (int iv = 0; iv < numVars; ++iv) {
        c[ic] += V[ic+numConstr*iv]*x[iv];
        c[ic+numConstr] += V[ic+numConstr*iv]*x[iv+numVars];
      }
    }
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) final override {
    for (int ic = 0; ic < numConstr; ++ic) {
      jv[ic] = static_cast<Real>(0);
      jv[ic+numConstr] = static_cast<Real>(0);
      for (int iv = 0; iv < numVars; ++iv) {
        jv[ic] += V[ic+numConstr*iv]*v[iv];
        jv[ic+numConstr] += V[ic+numConstr*iv]*v[iv+numVars];
      }
    }
  }

  void applyAdjointJacobian(std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) final override {
    for (int iv = 0; iv < numVars; ++iv) {
      ajv[iv] = static_cast<Real>(0);
      ajv[iv+numVars] = static_cast<Real>(0);
      for (int ic = 0; ic < numConstr; ++ic) {
        ajv[iv] += V[ic+numConstr*iv]*v[ic];
        ajv[iv+numVars] += V[ic+numConstr*iv]*v[ic+numConstr];
      }
    }
  }

private:
  const std::vector<Real> V;
  const int numConstr, numVars;
};

typedef double RealT;

int main(int argc, char *argv[]) {
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

  // *** Example body.
 
  try {

    const int numVars = 2*2;
    ROL::Ptr<ROL::StdVector<RealT>> x = ROL::makePtr<ROL::StdVector<RealT>>(numVars, 0);

    ROL::Ptr<ROL::Objective<RealT>> pobj = ROL::makePtr<ROL_TestObjective<RealT>>();

    const int numConstraints = 2*4;
    ROL::Ptr<ROL::StdVector<RealT>>  imul = ROL::makePtr<ROL::StdVector<RealT>>(numConstraints, 0);
    ROL::Ptr<ROL::Constraint<RealT>> icon = ROL::makePtr<ROL_LinearConstraint<RealT>>();

    const RealT boundSize = 1e-3;
    ROL::Ptr<ROL::StdVector<RealT>> ll = ROL::makePtr<ROL::StdVector<RealT>>(numConstraints, -boundSize);
    ROL::Ptr<ROL::StdVector<RealT>> lu = ROL::makePtr<ROL::StdVector<RealT>>(numConstraints,  boundSize);
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ll,lu);

    // Solve using Augmented Lagrangian from ROL version 1.0
    ROL::ParameterList parlist;
    parlist.sublist("Step").set("Type","Augmented Lagrangian");
    ROL::OptimizationProblem<RealT> problem1(pobj, x, icon, imul, ibnd);
    ROL::OptimizationSolver<RealT> solver1(problem1, parlist);
    solver1.solve(*outStream);
    x->print(*outStream);

    // Solve using Augmented Lagrangian from ROL version 2.0
    x->zero();
    parlist.sublist("Step").set("Type","Augmented Lagrangian");
    parlist.sublist("General").set("Output Level", 1);
    ROL::Ptr<ROL::Problem<RealT>> problem2 = ROL::makePtr<ROL::Problem<RealT>>(pobj,x);
    problem2->addLinearConstraint("Linear Inequality Constraint", icon, imul, ibnd);
    problem2->finalize(true,true,*outStream);
    ROL::Solver<RealT> solver2(problem2, parlist);
    solver2.solve(*outStream);
    x->print(*outStream);

    // Solve using Augmented Lagrangian from ROL version 2.0
    x->zero();
    parlist.sublist("Step").set("Type","Trust Region");
    parlist.sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    parlist.sublist("General").set("Output Level", 1);
    parlist.sublist("General").sublist("Polyhedral Projection").set("Type","Semismooth Newton");
    ROL::Ptr<ROL::Problem<RealT>> problem3 = ROL::makePtr<ROL::Problem<RealT>>(pobj,x);
    problem3->addLinearConstraint("Linear Inequality Constraint", icon, imul, ibnd);
    problem3->setProjectionAlgorithm(parlist);
    problem3->finalize(false,true,*outStream);
    ROL::Solver<RealT> solver3(problem3, parlist);
    solver3.solve(*outStream);
    x->print(*outStream);

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
