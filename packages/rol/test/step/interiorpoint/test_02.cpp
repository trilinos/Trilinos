// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#define OPTIMIZATION_PROBLEM_REFACTOR

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_NonlinearProgram.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_InteriorPointPenalty.hpp"
#include "ROL_PrimalDualInteriorPointResidual.hpp"

#include "HS_Problem_041.hpp"

#include <iomanip>

/*! \file test_02.cpp 
    \brief Perform a finite difference check on for the Hessian of the Lagrangian
           resulting from a Type-EB problem. Note that this test has a required
           dependence on Sacado. 
           NOTE: The finite difference check can only be expected to pass if
           the system is not symmetrized and all bounds are finite. 

*/


template<class Real> 
void printVector( const ROL::Vector<Real> &x, std::ostream &outStream ) {

  try {
    ROL::Ptr<const std::vector<Real> > xp = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();

    outStream << "Standard Vector" << std::endl;
    for( size_t i=0; i<xp->size(); ++i ) {
      outStream << (*xp)[i] << std::endl;
    }
  }
  catch( const std::bad_cast& e ) {
    outStream << "Partitioned Vector" << std::endl;
    
    typedef ROL::PartitionedVector<Real>    PV;
    typedef typename PV::size_type          size_type;

    const PV &xpv = dynamic_cast<const PV&>(x);

    for( size_type i=0; i<xpv.numVectors(); ++i ) {
      outStream << "--------------------" << std::endl;
      printVector( *(xpv.get(i)), outStream );
    }
    outStream << "--------------------" << std::endl;
  }

}


// Exact residual for H&S Problem 41
template<class Real>
void value( ROL::Vector<Real> &c, const ROL::Vector<Real> &sol, const Real &mu ) {

  typedef std::vector<Real>            vector;
  typedef ROL::StdVector<Real>         SV;
  typedef ROL::PartitionedVector<Real> PV;

  typedef typename PV::size_type size_type;

  
  
  using ROL::dynamicPtrCast;

  const size_type OPT   = 0;
  const size_type EQUAL = 1;
  const size_type LOWER = 2;
  const size_type UPPER = 3;

  const PV &sol_pv = dynamic_cast<const PV&>(sol);
  const vector &x  = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(OPT))->getVector());
  const vector &l  = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(EQUAL))->getVector());
  const vector &zl = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(LOWER))->getVector());
  const vector &zu = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(UPPER))->getVector());

  PV &c_pv = dynamic_cast<PV&>(c);
  vector &cx  = *(ROL::dynamicPtrCast<SV>(c_pv.get(OPT))->getVector());
  vector &cl  = *(ROL::dynamicPtrCast<SV>(c_pv.get(EQUAL))->getVector());
  vector &czl = *(ROL::dynamicPtrCast<SV>(c_pv.get(LOWER))->getVector());
  vector &czu = *(ROL::dynamicPtrCast<SV>(c_pv.get(UPPER))->getVector());
 
  cx[0] = -x[1]*x[2] +   l[0] - zl[0] + zu[0];
  cx[1] = -x[0]*x[1] + 2*l[0] - zl[1] + zu[1];
  cx[2] = -x[0]*x[1] + 2*l[0] - zl[2] + zu[2];
  cx[3] =            -   l[0] - zl[3] + zu[3];

  cl[0] = x[0] + 2*x[1] + 2*x[2] - x[3];

  czl[0] = x[0]*zl[0] - mu;
  czl[1] = x[1]*zl[1] - mu;
  czl[2] = x[2]*zl[2] - mu;
  czl[3] = x[3]*zl[3] - mu;
 
  czu[0] = (1.0-x[0])*zu[0] - mu;
  czu[1] = (1.0-x[1])*zu[1] - mu;
  czu[2] = (1.0-x[2])*zl[2] - mu;
  czu[3] = (2.0-x[3])*zl[3] - mu;
}
 


// Exact residual for H&S Problem 41
template<class Real>
void applyJacobian( ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &sol ) { 

  typedef std::vector<Real>            vector;
  typedef ROL::StdVector<Real>         SV;
  typedef ROL::PartitionedVector<Real> PV;

  typedef typename PV::size_type size_type;

  
  
  using ROL::dynamicPtrCast;

  const size_type OPT   = 0;
  const size_type EQUAL = 1;
  const size_type LOWER = 2;
  const size_type UPPER = 3;

  const PV &sol_pv = dynamic_cast<const PV&>(sol);
  const vector &x  = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(OPT))->getVector());
//const vector &l  = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(EQUAL))->getVector());
  const vector &zl = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(LOWER))->getVector());
  const vector &zu = *(ROL::dynamicPtrCast<const SV>(sol_pv.get(UPPER))->getVector());

  const PV &v_pv = dynamic_cast<const PV&>(v);
  const vector &vx  = *(ROL::dynamicPtrCast<const SV>(v_pv.get(OPT))->getVector());
  const vector &vl  = *(ROL::dynamicPtrCast<const SV>(v_pv.get(EQUAL))->getVector());
  const vector &vzl = *(ROL::dynamicPtrCast<const SV>(v_pv.get(LOWER))->getVector());
  const vector &vzu = *(ROL::dynamicPtrCast<const SV>(v_pv.get(UPPER))->getVector());

  PV &jv_pv = dynamic_cast<PV&>(jv);
  vector &jvx  = *(ROL::dynamicPtrCast<SV>(jv_pv.get(OPT))->getVector());
  vector &jvl  = *(ROL::dynamicPtrCast<SV>(jv_pv.get(EQUAL))->getVector());
  vector &jvzl = *(ROL::dynamicPtrCast<SV>(jv_pv.get(LOWER))->getVector());
  vector &jvzu = *(ROL::dynamicPtrCast<SV>(jv_pv.get(UPPER))->getVector());
 
  jvx[0] = -x[1]*vx[2] - x[2]*vx[1] +   vl[0] - vzl[0] + vzu[0];
  jvx[1] = -x[0]*vx[2] - x[2]*vx[0] + 2*vl[0] - vzl[1] + vzu[1];
  jvx[2] = -x[0]*vx[1] - x[1]*vx[0] + 2*vl[0] - vzl[2] + vzu[2];
  jvx[3] =                          -   vl[0] - vzl[3] + vzu[3];

  jvl[0] = vx[0] + 2*vx[1] + 2*vx[2] - vx[3];

  jvzl[0] = zl[0]*vx[0] + vzl[0]*x[0];
  jvzl[1] = zl[1]*vx[1] + vzl[1]*x[1];
  jvzl[2] = zl[2]*vx[2] + vzl[2]*x[2];
  jvzl[3] = zl[3]*vx[3] + vzl[3]*x[3];

  jvzu[0] = -zu[0]*vx[0] + vzu[0]*(1.0-x[0]);
  jvzu[1] = -zu[1]*vx[1] + vzu[1]*(1.0-x[1]);
  jvzu[2] = -zu[2]*vx[2] + vzu[2]*(1.0-x[2]);
  jvzu[3] = -zu[3]*vx[3] + vzu[3]*(2.0-x[3]);

}
 

typedef double RealT;

int main(int argc, char *argv[]) {
 
//  typedef std::vector<RealT>                          vector;

  typedef ROL::ParameterList                      PL;

  typedef ROL::Vector<RealT>                          V;
  typedef ROL::PartitionedVector<RealT>               PV;
  typedef ROL::Objective<RealT>                       OBJ;
  typedef ROL::Constraint<RealT>                      CON;
  typedef ROL::BoundConstraint<RealT>                 BND;
  typedef ROL::OptimizationProblem<RealT>             OPT;
  typedef ROL::NonlinearProgram<RealT>                NLP;

  typedef ROL::InteriorPointPenalty<RealT>            PENALTY;
  typedef ROL::PrimalDualInteriorPointResidual<RealT> RESIDUAL;

   

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs;
  if( iprint > 0 ) 
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;
   
  try {

    RealT mu = 0.1;

    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

    PL parlist;

    PL &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");
    PL &lblist = iplist.sublist("Barrier Objective");

    iplist.set("Symmetrize Primal Dual System",false);

    lblist.set("Use Linear Damping", true);
    lblist.set("Linear Damping Coefficient",1.e-4);
    lblist.set("Initial Barrier Parameter", mu);

    // Need an example problem that satisfies the following criteria:
    // 1) Has an equality constraint 
    // 2) Has a bound constraint where all variables have finite upper and lower bounds

    ROL::Ptr<NLP> nlp = ROL::makePtr<HS::Problem_041<RealT>>(); 
    ROL::Ptr<OPT> opt = nlp->getOptimizationProblem();
 
    ROL::Ptr<V>   x   = opt->getSolutionVector();
    ROL::Ptr<V>   l   = opt->getMultiplierVector();
    ROL::Ptr<V>   zl  = x->clone();
    ROL::Ptr<V>   zu  = x->clone();

    ROL::Ptr<V>   scratch = x->clone();

    ROL::Ptr<PV>  x_pv = ROL::dynamicPtrCast<PV>(x);

    ROL::Ptr<V>   sol = CreatePartitionedVector(x,l,zl,zu);

    ROL::Ptr<V>   c = sol->clone();
    ROL::Ptr<V>   v = sol->clone();
    ROL::Ptr<V>  jv = sol->clone();

    ROL::Ptr<V>   c_exact = c->clone();
    ROL::Ptr<V>  jv_exact = jv->clone();

    ROL::RandomizeVector(*l, -1.0, 1.0);
    ROL::RandomizeVector(*v,  0.0, 1.0);


    ROL::Ptr<OBJ> obj = opt->getObjective();
    ROL::Ptr<CON> con = opt->getConstraint();
    ROL::Ptr<BND> bnd = opt->getBoundConstraint();

    PENALTY penalty(obj,bnd,parlist);
 
    ROL::Ptr<const V> maskL = penalty.getLowerMask();
    ROL::Ptr<const V> maskU = penalty.getUpperMask();

    zl->set(*maskL);
    zu->set(*maskU);

    ROL::Ptr<CON> res = ROL::makePtr<RESIDUAL>(obj,con,bnd,*sol,maskL,maskU,scratch,mu,false);


    *outStream << "\n[x|lambda|zl|zu]" << std::endl;
    printVector(*sol,*outStream);

    res->value(*c,*sol,tol);
    value(*c_exact,*sol,mu);

    *outStream << "\n[cx|clambda|czl|czu]" << std::endl;

    printVector(*c,*outStream);

    c->axpy(-1.0,*c_exact);

    RealT cerror = c->norm();

    if( cerror>tol ) {
      ++errorFlag;
    }

    *outStream << "\n\n||c-c_exact|| = " << cerror << std::endl;

    res->applyJacobian(*jv,*v,*sol,tol);
    applyJacobian(*jv_exact,*v,*sol);

    *outStream << "\n[vx|vlambda|vzl|vzu]" << std::endl;

    printVector(*v,*outStream);

    *outStream << "\n[jvx|jvlambda|jvzl|jvzu]" << std::endl;

    printVector(*jv,*outStream);

    jv->axpy(-1.0,*jv_exact);

    RealT jverror = jv->norm();
    
    if( jverror > tol ) {
      ++errorFlag;  
    }

    *outStream << "\n\n||jv-jv_exact|| = " << jverror << std::endl;

    *outStream << "Residual Jacobian Finite Difference Check" << std::endl;
    res->checkApplyJacobian(*sol,*v,*v);



  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;
}

