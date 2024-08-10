// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*  Example of solving a problem with bound and inequality constraints
 *
 */

#define OPTIMIZATION_PROBLEM_REFACTOR 

#include "ROL_OptimizationSolver.hpp"

#include "ROL_RandomVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


/* OBJECTIVE FUNCTION */

template<class Real> 
class ObjectiveQL : public ROL::StdObjective<Real> {
private:
  std::vector<Real> coeff;

public:

  ObjectiveQL() : coeff({-21.98,-1.26,61.39,5.3,101.3}) { 
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    Real result = 0;
    for( int i=0; i<5; ++i ) {
      result +=x[i]*(0.5*x[i]+coeff[i]);
    }  
    return result;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    for( int i=0; i<5; ++i ) {
      g[i] = x[i]+coeff[i];
    }  
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    hv = v;
  }

}; // class ObjectiveQL

/* INEQUALITY CONSTRAINT */

template<class Real> 
class InequalityQL : public ROL::StdConstraint<Real> {
private:
  std::vector<Real> coeff;
  Real offset;

public:
  InequalityQL() : coeff({-7.56,0.0,0.0,0.0,0.5}), offset(39.1) {}

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol) {
    c[0] = offset;
    for( int i=0; i<5; ++i ) {
      c[0] += coeff[i]*x[i];
    }
  }

  void applyJacobian(  std::vector<Real> &jv, const std::vector<Real> &v ,const std::vector<Real> &x, Real &tol ) {
    jv[0] = 0;
    for( int i=0; i<5; ++i ) {
      jv[0] += coeff[i]*v[i];
    }
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v ,const std::vector<Real> &x, Real &tol ) {
    for( int i=0; i<5; ++i ) {
      ajv[i] = v[0]*coeff[i];
    }
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                       const std::vector<Real> &v ,const std::vector<Real> &x, Real &tol ) {
    ahuv.assign(5,0.0);
  }

}; // class InequalityQL


int main(int argc, char *argv[]) {

   

  typedef double RealT;
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);


  int errorFlag   = 0;

  try {
 
    ROL::ParameterList parlist;
    parlist.sublist("Step").set("Type","Augmented Lagrangian");
    

    ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(5,-100.0);
    ROL::Ptr<std::vector<RealT> > u_ptr = ROL::makePtr<std::vector<RealT>>(5, 100.0);

    ROL::Ptr<ROL::Vector<RealT> > lower = ROL::makePtr<ROL::StdVector<RealT>>( l_ptr );
    ROL::Ptr<ROL::Vector<RealT> > upper = ROL::makePtr<ROL::StdVector<RealT>>( u_ptr ); 

    ROL::Ptr<std::vector<RealT> > x_ptr  = ROL::makePtr<std::vector<RealT>>(5,1.0);
    ROL::Ptr<std::vector<RealT> > li_ptr = ROL::makePtr<std::vector<RealT>>(1,0.0);
    ROL::Ptr<std::vector<RealT> > ll_ptr = ROL::makePtr<std::vector<RealT>>(1,0.0);
    ROL::Ptr<std::vector<RealT> > lu_ptr = ROL::makePtr<std::vector<RealT>>(1,ROL::ROL_INF<RealT>());

    ROL::Ptr<ROL::Vector<RealT> > x  = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    ROL::Ptr<ROL::Vector<RealT> > li = ROL::makePtr<ROL::StdVector<RealT>>(li_ptr);
    ROL::Ptr<ROL::Vector<RealT> > ll = ROL::makePtr<ROL::StdVector<RealT>>(ll_ptr);
    ROL::Ptr<ROL::Vector<RealT> > lu = ROL::makePtr<ROL::StdVector<RealT>>(lu_ptr);

    ROL::Ptr<ROL::Objective<RealT> >             obj  = ROL::makePtr<ObjectiveQL<RealT>>();
    ROL::Ptr<ROL::BoundConstraint<RealT> >       bnd  = ROL::makePtr<ROL::Bounds<RealT>>(lower,upper);
    ROL::Ptr<ROL::Constraint<RealT> >            ineq = ROL::makePtr<InequalityQL<RealT>>();
    ROL::Ptr<ROL::BoundConstraint<RealT> >       ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ll,lu);

    ROL::OptimizationProblem<RealT> problem( obj, x, bnd, ineq, li, ibnd);    

    /* checkAdjointJacobianConsistency fails for the OptimizationProblem if we don't do this first... why? */
    ROL::Ptr<ROL::Vector<RealT> > u = x->clone(); 
    RandomizeVector(*u);
    ineq->checkAdjointConsistencyJacobian(*li,*x,*u,true,*outStream);
    /*******************************************************************************************************/

    problem.check(*outStream);

    ROL::OptimizationSolver<RealT> solver( problem, parlist );

    solver.solve(*outStream); 



     


    *outStream << "x_opt = [";
    for(int i=0;i<4;++i) {
      *outStream << (*x_ptr)[i] << ", " ;
    } 
    *outStream << (*x_ptr)[4] << "]" << std::endl;
    

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


  return 0;
}





