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
#include "ROL_LinearOperatorFromConstraint.hpp"
#include "ROL_KrylovFactory.hpp"

#include "HS_ProblemFactory.hpp"

#include <iomanip>

/*! \file test_03.cpp 
    \brief Verify that the symmetrized version of a primal dual
           system is indeed symmetric and that the solution to 
           the unsymmetrized version satisfies the symmetrized version.
 
           Note: CG will almost certainly fail with exit flag 2 (negative
           eigenvalues)
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

template<class Real> 
void printMatrix( const std::vector<ROL::Ptr<ROL::Vector<Real> > > &A,
                  const std::vector<ROL::Ptr<ROL::Vector<Real> > > &I,
                  std::ostream &outStream ) {
  typedef typename std::vector<Real>::size_type uint;
  uint dim = A.size();
   
  for( uint i=0; i<dim; ++i ) {
    for( uint j=0; j<dim; ++j ) {
      outStream << std::setw(6) << A[j]->dot(*(I[i])); 
    }
    outStream << std::endl;
  }
}


template<class Real> 
class IdentityOperator : public ROL::LinearOperator<Real> {
public:
  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    Hv.set(v);
  }
}; // IdentityOperator


typedef double RealT;

int main(int argc, char *argv[]) {
 
//  typedef std::vector<RealT>                             vector;

  typedef ROL::ParameterList                           PL;

  typedef ROL::Vector<RealT>                               V;
  typedef ROL::PartitionedVector<RealT>                    PV;
  typedef ROL::Objective<RealT>                            OBJ;
  typedef ROL::Constraint<RealT>                           CON;
  typedef ROL::BoundConstraint<RealT>                      BND;
  typedef ROL::OptimizationProblem<RealT>                  OPT;
  typedef ROL::NonlinearProgram<RealT>                     NLP;
  typedef ROL::LinearOperator<RealT>                       LOP;  
  typedef ROL::LinearOperatorFromConstraint<RealT>         LOPEC;
  typedef ROL::Krylov<RealT>                               KRYLOV;


  typedef ROL::InteriorPointPenalty<RealT>                 PENALTY;
  typedef ROL::PrimalDualInteriorPointResidual<RealT>      RESIDUAL;

   

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

    lblist.set("Use Linear Damping", true);         // Not used in this test
    lblist.set("Linear Damping Coefficient",1.e-4); // Not used in this test 
    lblist.set("Initial Barrier Parameter", mu);

    PL &krylist = parlist.sublist("General").sublist("Krylov");
   
    krylist.set("Absolute Tolerance", 1.e-6);
    krylist.set("Relative Tolerance", 1.e-6);
    krylist.set("Iteration Limit", 50);

    // Create a Conjugate Gradients solver 
    krylist.set("Type","Conjugate Gradients"); 
    ROL::Ptr<KRYLOV> cg = ROL::KrylovFactory<RealT>(parlist);
    HS::ProblemFactory<RealT> problemFactory;

    // Choose an example problem with inequality constraints and
    // a mixture of finite and infinite bounds
    ROL::Ptr<NLP> nlp = problemFactory.getProblem(16);
    ROL::Ptr<OPT> opt = nlp->getOptimizationProblem();
 
    ROL::Ptr<V>   x   = opt->getSolutionVector();
    ROL::Ptr<V>   l   = opt->getMultiplierVector();
    ROL::Ptr<V>   zl  = x->clone(); zl->zero();
    ROL::Ptr<V>   zu  = x->clone(); zu->zero();

    ROL::Ptr<V>   scratch = x->clone();

    ROL::Ptr<PV>  x_pv = ROL::dynamicPtrCast<PV>(x);
    // New slack variable initialization does not guarantee strict feasibility.
    // This ensures that the slack variables are the same as the previous
    // implementation.
    (*ROL::dynamicPtrCast<ROL::StdVector<RealT> >(x_pv->get(1))->getVector())[0] = 1.0;

    ROL::Ptr<V>   sol = CreatePartitionedVector(x,l,zl,zu);   

    std::vector< ROL::Ptr<V> > I;
    std::vector< ROL::Ptr<V> > J;

    for( int k=0; k<sol->dimension(); ++k ) {
      I.push_back(sol->basis(k));
      J.push_back(sol->clone());
    }

    ROL::Ptr<V>   u = sol->clone();
    ROL::Ptr<V>   v = sol->clone();

    ROL::Ptr<V>   rhs = sol->clone();
    ROL::Ptr<V>   symrhs = sol->clone();

    ROL::Ptr<V>   gmres_sol = sol->clone();   gmres_sol->set(*sol);
    ROL::Ptr<V>   cg_sol = sol->clone();      cg_sol->set(*sol);
 
    IdentityOperator<RealT> identity;

    RandomizeVector(*u,-1.0,1.0);
    RandomizeVector(*v,-1.0,1.0);

    ROL::Ptr<OBJ> obj = opt->getObjective();
    ROL::Ptr<CON> con = opt->getConstraint();
    ROL::Ptr<BND> bnd = opt->getBoundConstraint();

    PENALTY penalty(obj,bnd,parlist);
 
    ROL::Ptr<const V> maskL = penalty.getLowerMask();
    ROL::Ptr<const V> maskU = penalty.getUpperMask();

    zl->set(*maskL);
    zu->set(*maskU);

    /********************************************************************************/
    /* Nonsymmetric representation test                                             */
    /********************************************************************************/

    int gmres_iter = 0;
    int gmres_flag = 0;

    // Form the residual's Jacobian operator
    ROL::Ptr<CON> res = ROL::makePtr<RESIDUAL>(obj,con,bnd,*sol,maskL,maskU,scratch,mu,false);
    ROL::Ptr<LOP> lop = ROL::makePtr<LOPEC>( sol, res );

    // Evaluate the right-hand-side
    res->value(*rhs,*sol,tol);

    // Create a GMRES solver
    krylist.set("Type","GMRES");
    ROL::Ptr<KRYLOV> gmres = ROL::KrylovFactory<RealT>(parlist);

     for( int k=0; k<sol->dimension(); ++k ) {
      res->applyJacobian(*(J[k]),*(I[k]),*sol,tol);
    }

    *outStream << "Nonsymmetric Jacobian" << std::endl;
    printMatrix(J,I,*outStream);

   // Solve the system 
    gmres->run( *gmres_sol, *lop, *rhs, identity, gmres_iter, gmres_flag );

    errorFlag += gmres_flag;

    *outStream << "GMRES terminated after " << gmres_iter << " iterations "
               << "with exit flag " << gmres_flag << std::endl;


    /********************************************************************************/
    /* Symmetric representation test                                                */
    /********************************************************************************/

    int cg_iter = 0;
    int cg_flag = 0;

    ROL::Ptr<V> jv = v->clone();
    ROL::Ptr<V> ju = u->clone();

    iplist.set("Symmetrize Primal Dual System",true);
    ROL::Ptr<CON> symres = ROL::makePtr<RESIDUAL>(obj,con,bnd,*sol,maskL,maskU,scratch,mu,true);
    ROL::Ptr<LOP> symlop = ROL::makePtr<LOPEC>( sol, res );
    symres->value(*symrhs,*sol,tol);

    symres->applyJacobian(*jv,*v,*sol,tol);
    symres->applyJacobian(*ju,*u,*sol,tol);
    *outStream << "Symmetry check |u.dot(jv)-v.dot(ju)| = "
               << std::abs(u->dot(*jv)-v->dot(*ju)) << std::endl;
 
    cg->run( *cg_sol, *symlop, *symrhs, identity, cg_iter, cg_flag );

    *outStream << "CG terminated after " << cg_iter << " iterations "
               << "with exit flag " << cg_flag << std::endl;

    *outStream << "Check that GMRES solution also is a solution to the symmetrized system"
               << std::endl;

    symres->applyJacobian(*ju,*gmres_sol,*sol,tol);
    ju->axpy(-1.0,*symrhs);
    RealT mismatch = ju->norm();
    if(mismatch>tol) {
      errorFlag++;
    }
    *outStream << "||J_sym*sol_nonsym-rhs_sym|| = " << mismatch << std::endl;

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

