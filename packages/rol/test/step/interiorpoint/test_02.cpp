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

#define OPTIMIZATION_PROBLEM_REFACTOR

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
    Teuchos::RCP<const std::vector<Real> > xp = 
      Teuchos::dyn_cast<const ROL::StdVector<Real> >(x).getVector();

    outStream << "Standard Vector" << std::endl;
    for( size_t i=0; i<xp->size(); ++i ) {
      outStream << (*xp)[i] << std::endl;
    }
  }
  catch( const std::bad_cast& e ) {
    outStream << "Partitioned Vector" << std::endl;
    
    typedef ROL::PartitionedVector<Real>    PV;
    typedef typename PV::size_type          size_type;

    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

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

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

  const size_type OPT   = 0;
  const size_type EQUAL = 1;
  const size_type LOWER = 2;
  const size_type UPPER = 3;

  const PV &sol_pv = dyn_cast<const PV>(sol);
  const vector &x  = *(rcp_dynamic_cast<const SV>(sol_pv.get(OPT))->getVector());
  const vector &l  = *(rcp_dynamic_cast<const SV>(sol_pv.get(EQUAL))->getVector());
  const vector &zl = *(rcp_dynamic_cast<const SV>(sol_pv.get(LOWER))->getVector());
  const vector &zu = *(rcp_dynamic_cast<const SV>(sol_pv.get(UPPER))->getVector());

  PV &c_pv = dyn_cast<PV>(c);
  vector &cx  = *(rcp_dynamic_cast<SV>(c_pv.get(OPT))->getVector());
  vector &cl  = *(rcp_dynamic_cast<SV>(c_pv.get(EQUAL))->getVector());
  vector &czl = *(rcp_dynamic_cast<SV>(c_pv.get(LOWER))->getVector());
  vector &czu = *(rcp_dynamic_cast<SV>(c_pv.get(UPPER))->getVector());
 
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

  using Teuchos::RCP;
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

  const size_type OPT   = 0;
  const size_type EQUAL = 1;
  const size_type LOWER = 2;
  const size_type UPPER = 3;

  const PV &sol_pv = dyn_cast<const PV>(sol);
  const vector &x  = *(rcp_dynamic_cast<const SV>(sol_pv.get(OPT))->getVector());
//const vector &l  = *(rcp_dynamic_cast<const SV>(sol_pv.get(EQUAL))->getVector());
  const vector &zl = *(rcp_dynamic_cast<const SV>(sol_pv.get(LOWER))->getVector());
  const vector &zu = *(rcp_dynamic_cast<const SV>(sol_pv.get(UPPER))->getVector());

  const PV &v_pv = dyn_cast<const PV>(v);
  const vector &vx  = *(rcp_dynamic_cast<const SV>(v_pv.get(OPT))->getVector());
  const vector &vl  = *(rcp_dynamic_cast<const SV>(v_pv.get(EQUAL))->getVector());
  const vector &vzl = *(rcp_dynamic_cast<const SV>(v_pv.get(LOWER))->getVector());
  const vector &vzu = *(rcp_dynamic_cast<const SV>(v_pv.get(UPPER))->getVector());

  PV &jv_pv = dyn_cast<PV>(jv);
  vector &jvx  = *(rcp_dynamic_cast<SV>(jv_pv.get(OPT))->getVector());
  vector &jvl  = *(rcp_dynamic_cast<SV>(jv_pv.get(EQUAL))->getVector());
  vector &jvzl = *(rcp_dynamic_cast<SV>(jv_pv.get(LOWER))->getVector());
  vector &jvzu = *(rcp_dynamic_cast<SV>(jv_pv.get(UPPER))->getVector());
 
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

  typedef Teuchos::ParameterList                      PL;

  typedef ROL::Vector<RealT>                          V;
  typedef ROL::PartitionedVector<RealT>               PV;
  typedef ROL::Objective<RealT>                       OBJ;
  typedef ROL::Constraint<RealT>                      CON;
  typedef ROL::BoundConstraint<RealT>                 BND;
  typedef ROL::OptimizationProblem<RealT>             OPT;
  typedef ROL::NonlinearProgram<RealT>                NLP;

  typedef ROL::InteriorPointPenalty<RealT>            PENALTY;
  typedef ROL::PrimalDualInteriorPointResidual<RealT> RESIDUAL;

  using Teuchos::RCP; using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;
  RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs;
  if( iprint > 0 ) 
    outStream = rcp(&std::cout,false);
  else
    outStream = rcp(&bhs,false);

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

    RCP<NLP> nlp = rcp( new HS::Problem_041<RealT>() ); 
    RCP<OPT> opt = nlp->getOptimizationProblem();
 
    RCP<V>   x   = opt->getSolutionVector();
    RCP<V>   l   = opt->getMultiplierVector();
    RCP<V>   zl  = x->clone();
    RCP<V>   zu  = x->clone();

    RCP<V>   scratch = x->clone();

    RCP<PV>  x_pv = Teuchos::rcp_dynamic_cast<PV>(x);

    RCP<V>   sol = CreatePartitionedVector(x,l,zl,zu);

    RCP<V>   c = sol->clone();
    RCP<V>   v = sol->clone();
    RCP<V>  jv = sol->clone();

    RCP<V>   c_exact = c->clone();
    RCP<V>  jv_exact = jv->clone();

    ROL::RandomizeVector(*l, -1.0, 1.0);
    ROL::RandomizeVector(*v,  0.0, 1.0);


    RCP<OBJ> obj = opt->getObjective();
    RCP<CON> con = opt->getConstraint();
    RCP<BND> bnd = opt->getBoundConstraint();

    PENALTY penalty(obj,bnd,parlist);
 
    RCP<const V> maskL = penalty.getLowerMask();
    RCP<const V> maskU = penalty.getUpperMask();

    zl->set(*maskL);
    zu->set(*maskU);

    RCP<CON> res = rcp( new RESIDUAL(obj,con,bnd,*sol,maskL,maskU,scratch,mu,false)  );


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
  catch (std::logic_error err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;
}

