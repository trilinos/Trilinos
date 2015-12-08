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
//         
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_07.cpp
    \brief Interior Point test using Hock & Schittkowski problem 32.
*/

#include "ROL_HS32.hpp"
#include "ROL_LogBarrierObjective.hpp"
#include "ROL_InteriorPoint.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_OptimizationProblem.hpp"


template<class Real> 
void print_vector(const ROL::Vector<Real> &x) {
  typedef ROL::StdVector<Real> SV;
  using Teuchos::dyn_cast;
  const SV &xs = dyn_cast<const SV>(x);
  Teuchos::RCP<const std::vector<Real> > x_rcp = xs.getVector();
   
  for(int i=0;i<xs.dimension();++i) {
    std::cout << (*x_rcp)[i] << std::endl; 
  }
}

template<class Real> 
void print_subvector(const ROL::Vector<Real> &x, const int i){
  typedef ROL::PartitionedVector<Real> PV;
  typedef typename PV::size_type size_type;
  size_type n = static_cast<size_type>(i);
  using Teuchos::dyn_cast;
  const PV &xp = dyn_cast<const PV>(x);
  print_vector(*(xp.get(n)));
}

typedef double RealT;

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp; 

  typedef std::vector<RealT>            vec;
  typedef ROL::StdVector<RealT>         SV;
  typedef RCP<ROL::Vector<RealT> >      RCPV;

//  typedef ROL::PartitionedVector<RealT> PV;


  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = rcp(&std::cout, false);
  else
    outStream = rcp(&bhs, false);

  int errorFlag = 0;

  try {

    int xopt_dim  = 3; // Dimension of optimization vectors
    int ce_dim = 1;    // Dimension of equality constraint
    int ci_dim = 4;    // Dimension of inequality constraint

    RCP<vec> xopt_rcp = rcp( new vec(xopt_dim,0.0) );
    RCP<vec> dopt_rcp = rcp( new vec(xopt_dim,0.0) );
    RCP<vec> vopt_rcp = rcp( new vec(xopt_dim,0.0) );
   
    RCP<vec> vec_rcp  = rcp( new vec(ce_dim,1.0) );
    RCP<vec> vel_rcp  = rcp( new vec(ce_dim,1.0) );
   
    RCP<vec> vic_rcp  = rcp( new vec(ci_dim,0.0) );
    RCP<vec> vil_rcp  = rcp( new vec(ci_dim,0.0) );
     
    // Slack variables
    RCP<vec> xs_rcp = rcp( new vec(ci_dim,1.0) );    
    RCP<vec> vs_rcp = rcp( new vec(ci_dim,0.0) );
    RCP<vec> ds_rcp = rcp( new vec(ci_dim,0.0) );

    // Feasible initial guess
    (*xopt_rcp)[0] = 0.1;
    (*xopt_rcp)[1] = 0.7;
    (*xopt_rcp)[2] = 0.2;

    RealT left = -1e0, right = 1e0;
    for (int i=0; i<xopt_dim; i++) {
      (*dopt_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vopt_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }    

    for (int i=0; i<ci_dim; i++) {
      (*vic_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vil_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vs_rcp)[i]  = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*ds_rcp)[i]  = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }    

    RCPV xopt = rcp( new SV(xopt_rcp) );
    RCPV dopt = rcp( new SV(dopt_rcp) );
    RCPV vopt = rcp( new SV(vopt_rcp) );
    RCPV vec  = rcp( new SV(vec_rcp) );
    RCPV vel  = rcp( new SV(vel_rcp) );
    RCPV vic  = rcp( new SV(vic_rcp) );
    RCPV vil  = rcp( new SV(vil_rcp) );
    RCPV xs   = rcp( new SV(xs_rcp) );
    RCPV vs   = rcp( new SV(vs_rcp) );
    RCPV ds   = rcp( new SV(ds_rcp) );

    RCPV gopt = xopt->dual().clone();
    RCPV gs   = xs->dual().clone();

   

    // Partitioned vectors of optimization and slack variables
    RCPV x = CreatePartitionedVector(xopt,xs);
    RCPV v = CreatePartitionedVector(vopt,vs);
    RCPV d = CreatePartitionedVector(dopt,ds);
    RCPV vc = CreatePartitionedVector(vic,vec);
    RCPV vl = CreatePartitionedVector(vil,vel);

    RCPV g = CreatePartitionedVector(gopt,gs);

    // Original obective
    using ROL::ZOO::Objective_HS32;
    using ROL::ZOO::EqualityConstraint_HS32;
    using ROL::ZOO::InequalityConstraint_HS32;    

    RCP<ROL::Objective<RealT> > obj_hs32 = rcp( new Objective_HS32<RealT> ); 
    RCP<ROL::EqualityConstraint<RealT> > eqcon_hs32 = rcp( new EqualityConstraint_HS32<RealT> );
    RCP<ROL::InequalityConstraint<RealT> > incon_hs32 = rcp( new  InequalityConstraint_HS32<RealT> );
    

    RCP<Teuchos::ParameterList> parlist = rcp(new Teuchos::ParameterList);
    std::string stepname = "Interior Point";

    RealT mu = 0.1;            // Initial penalty parameter
    RealT factor = 0.2;        // Penalty reduction factor

    parlist->sublist("Step").sublist("Interior Point").set("Initial Barrier Penalty",mu);
    parlist->sublist("Step").sublist("Interior Point").set("Minimium Barrier Penalty",1e-8);
    parlist->sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor",factor);
    parlist->sublist("Step").sublist("Interior Point").set("Subproblem Iteration Limit",30);

    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist->sublist("Step").sublist("Composite Step").set("Output Level",0);

    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Iteration Limit",100);

    ROL::OptimizationProblem<RealT> problem( obj_hs32, xopt, eqcon_hs32, vel, incon_hs32, vil, parlist);  

/*
    *outStream << "\nChecking individual objectives and constraints separately\n" << std::endl;

    *outStream << "\nObjective\n" << std::endl;
    obj_hs32.checkGradient(*xopt,*dopt,true,*outStream);
    obj_hs32.checkHessVec(*xopt,*vopt,true,*outStream);

    
    *outStream << "\nEquality Constraint\n" << std::endl;
    eqcon_hs32.checkApplyJacobian(*xopt,*vopt,*vec,true,*outStream); 
    eqcon_hs32.checkApplyAdjointJacobian(*xopt,*vel,*vec,*xopt,true,*outStream); 

    *outStream << "\nInequality Constraint\n" << std::endl;
    incon_hs32.checkApplyJacobian(*xopt,*vopt,*vic,true,*outStream); 
    incon_hs32.checkApplyAdjointJacobian(*xopt,*vil,*vic,*xopt,true,*outStream); 
    incon_hs32.checkApplyAdjointHessian(*xopt,*vil,*dopt,*xopt,true,*outStream);


    *outStream << "\nCheck Interior Point objective\n" << std::endl;
    ipobj->checkGradient(*x,*d,true,*outStream);
    ipobj->checkHessVec(*x,*v,true,*outStream);

    *outStream << "\nCheck Interior Point constraints\n" << std::endl;
    ipcon->checkApplyJacobian(*x,*v,*vc,true,*outStream);
    ipcon->checkApplyAdjointJacobian(*x,*vl,*vc,*x,true,*outStream);
    ipcon->checkApplyAdjointHessian(*x,*vl,*d,*x,true,*outStream);    
*/
    // Define algorithm.

    RCP<ROL::Algorithm<RealT> > algo;    
    algo = rcp( new ROL::Algorithm<RealT>(stepname,*parlist) );

    algo->run(problem,true,*outStream);   

/*
    for(int iter = 0; iter<8; ++iter) {
      mu *= factor;
      parlist.sublist("Step").sublist("Interior Point").set("Initial Barrier Penalty",mu);
      ipobj->updatePenalty(mu);     

      algo = rcp( new ROL::Algorithm<RealT>(stepname,parlist) );

      algo->run(*x,*g,*vl,*vc,*ipobj,*ipcon,true,*outStream);   


    }

    RealT tol = 0.0; 

    ipobj->gradient(*g,*x,tol);
    ipcon->value(*vc,*x,tol);

    *outStream << "------------" << std::endl;
    *outStream << "c_I(x)" << std::endl;
    print_subvector(*vc,0);
    *outStream << "------------" << std::endl;

    *outStream << "c_E(x)" << std::endl;
    print_subvector(*vc,1);
    *outStream << "------------" << std::endl;

    *outStream << "x_opt" << std::endl;
    print_subvector(*x,0);
    *outStream << "-----------" << std::endl;

    *outStream << "x_slack" << std::endl;
    print_subvector(*x,1);
    *outStream << "------------" << std::endl;

    *outStream << "g_opt" << std::endl;
    print_subvector(*g,0);
    *outStream << "-----------" << std::endl;

    *outStream << "g_slack" << std::endl;
    print_subvector(*g,1);
    *outStream << "------------" << std::endl;

    *outStream << g->norm() << std::endl;
*/

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

