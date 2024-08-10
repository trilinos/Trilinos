// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the equality constrained NLP
           from Nocedal/Wright, 2nd edition, page 574, example 18.2.
*/

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
//#include <fenv.h>

typedef double RealT;

/*** Declare four vector spaces. ***/

// Forward declarations:

template <class Real, class Element=Real>
class OptStdVector;  // Optimization space.

template <class Real, class Element=Real>
class OptDualStdVector;  // Dual optimization space.

template <class Real, class Element=Real>
class ConStdVector;  // Constraint space.

template <class Real, class Element=Real>
class ConDualStdVector;  // Dual constraint space.

// Vector space definitions:

// Optimization space.
template <class Real, class Element>
class OptStdVector : public ROL::Vector<Real> {

typedef std::vector<Element>       vector;
typedef ROL::Vector<Real>          V;
typedef typename vector::size_type uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<OptDualStdVector<Real> >  dual_vec_;

public:

OptStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
    
  
  ROL::Ptr<const vector> xp = dynamic_cast<const OptStdVector&>(x).getVector();
 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xp)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {

    

  ROL::Ptr<const vector> xp = dynamic_cast<const OptStdVector&>(x).getVector();
  Real val = 0;
 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xp)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

ROL::Ptr<ROL::Vector<Real> > clone() const {
  return ROL::makePtr<OptStdVector>( ROL::makePtr<std::vector<Element>>(std_vec_->size()) );
}

ROL::Ptr<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::Ptr<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(),0.0);
  ROL::Ptr<V> e = ROL::makePtr<OptStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makePtr<OptDualStdVector<Real>>( ROL::makePtr<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

Real apply( const ROL::Vector<Real> &x ) const {
  Real val = 0;
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptDualStdVector<Real,Element>&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

typedef std::vector<Element>       vector;
typedef ROL::Vector<Real>          V;
typedef typename vector::size_type uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<OptStdVector<Real> >  dual_vec_;

public:

OptDualStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const OptDualStdVector&>(x).getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xp)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const OptDualStdVector&>(x).getVector();
  Real val = 0;
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xp)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

ROL::Ptr<ROL::Vector<Real> > clone() const {
  return ROL::makePtr<OptDualStdVector>( ROL::makePtr<std::vector<Element>>(std_vec_->size()) );
}

ROL::Ptr<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::Ptr<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>( std_vec_->size(), 0.0 );
  ROL::Ptr<V> e = ROL::makePtr<OptDualStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makePtr<OptStdVector<Real>>( ROL::makePtr<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

Real apply( const ROL::Vector<Real> &x ) const {
  Real val = 0;
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptStdVector<Real,Element>&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

}; // class OptDualStdVector


// Constraint space.
template <class Real, class Element>
class ConStdVector : public ROL::Vector<Real> {

typedef std::vector<Element>       vector;
typedef ROL::Vector<Real>          V;
typedef typename vector::size_type uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<ConDualStdVector<Real> >  dual_vec_;

public:

ConStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const ConStdVector&>(x).getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xp)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const ConStdVector&>(x).getVector();
  Real val = 0;
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xp)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

ROL::Ptr<ROL::Vector<Real> > clone() const {
  return ROL::makePtr<ConStdVector>( ROL::makePtr<std::vector<Element>>(std_vec_->size()));
}

ROL::Ptr<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::Ptr<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(), 0.0);
  ROL::Ptr<V> e = ROL::makePtr<ConStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makePtr<ConDualStdVector<Real>>( ROL::makePtr<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

Real apply( const ROL::Vector<Real> &x ) const {
  Real val = 0;
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const ConDualStdVector<Real,Element>&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>       vector;
  typedef ROL::Vector<Real>          V;
  typedef typename vector::size_type uint;

private:

ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<ConStdVector<Real> >  dual_vec_;

public:

ConDualStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const ConDualStdVector&>(x).getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xp)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {
    
  ROL::Ptr<const vector> xp = dynamic_cast<const ConDualStdVector&>(x).getVector();
  Real val = 0;
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xp)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

ROL::Ptr<ROL::Vector<Real> > clone() const {
  return ROL::makePtr<ConDualStdVector>( ROL::makePtr<std::vector<Element>>(std_vec_->size()));
}

ROL::Ptr<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::Ptr<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(),0.0);
  ROL::Ptr<V> e = ROL::makePtr<ConDualStdVector>(e_ptr);
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makePtr<ConStdVector<Real>>( ROL::makePtr<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

Real apply( const ROL::Vector<Real> &x ) const {
  Real val = 0;
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const ConStdVector<Real,Element>&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

}; // class ConDualStdVector

/*** End of declaration of four vector spaces. ***/


int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  typedef std::vector<RealT> vector;
  typedef vector::size_type  uint;

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

    uint dim = 5;
    uint nc = 3;
    ROL::Ptr<ROL::Objective<RealT> > obj;
    ROL::Ptr<ROL::Constraint<RealT> > constr;
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> sol_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<OptStdVector<RealT>> x   = ROL::makePtr<OptStdVector<RealT>>(x_ptr);   // Iteration vector.
    ROL::Ptr<OptStdVector<RealT>> sol = ROL::makePtr<OptStdVector<RealT>>(sol_ptr); // Reference solution vector.

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained <RealT, OptStdVector<RealT>, OptDualStdVector<RealT>, ConStdVector<RealT>, ConDualStdVector<RealT> > SEC;
    obj = SEC.getObjective();
    constr = SEC.getEqualityConstraint();
    x->set(*SEC.getInitialGuess());
    sol->set(*SEC.getSolution());

    // Run derivative checks, etc.
    RealT left = -1e0, right = 1e0;
    ROL::Ptr<vector> xtest_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> g_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> d_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> gd_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> vc_ptr = ROL::makePtr<vector>(nc, 0.0);
    ROL::Ptr<vector> vl_ptr = ROL::makePtr<vector>(nc, 0.0);
    ROL::Ptr<OptStdVector<RealT>>  xtest = ROL::makePtr<OptStdVector<RealT>>(xtest_ptr);
    ROL::Ptr<OptDualStdVector<RealT>>  g = ROL::makePtr<OptDualStdVector<RealT>>(g_ptr);
    ROL::Ptr<OptStdVector<RealT>>      d = ROL::makePtr<OptStdVector<RealT>>(d_ptr);
    ROL::Ptr<OptDualStdVector<RealT>> gd = ROL::makePtr<OptDualStdVector<RealT>>(gd_ptr);
    ROL::Ptr<OptStdVector<RealT>>      v = ROL::makePtr<OptStdVector<RealT>>(v_ptr);
    ROL::Ptr<ConStdVector<RealT>>     vc = ROL::makePtr<ConStdVector<RealT>>(vc_ptr);
    ROL::Ptr<ConDualStdVector<RealT>> vl = ROL::makePtr<ConDualStdVector<RealT>>(vl_ptr);
    // set xtest, d, v
    for (uint i=0; i<dim; i++) {
      (*xtest_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*gd_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    // set vc, vl
    for (uint i=0; i<nc; i++) {
      (*vc_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vl_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    obj->checkGradient(*xtest, *g, *d, true, *outStream);                      *outStream << std::endl;
    obj->checkHessVec(*xtest, *g, *v, true, *outStream);                       *outStream << std::endl;
    obj->checkHessSym(*xtest, *g, *d, *v, true, *outStream);                   *outStream << std::endl;
    constr->checkApplyJacobian(*xtest, *v, *vc, true, *outStream);             *outStream << std::endl;
    constr->checkApplyAdjointJacobian(*xtest, *vl, *vc, *g, true, *outStream); *outStream << std::endl;
    constr->checkApplyAdjointHessian(*xtest, *vl, *d, *g, true, *outStream);   *outStream << std::endl;

    ROL::Ptr<vector> v1_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v2_ptr = ROL::makePtr<vector>(nc, 0.0);
    ROL::Ptr<OptStdVector<RealT>>     v1 = ROL::makePtr<OptStdVector<RealT>>(v1_ptr);
    ROL::Ptr<ConDualStdVector<RealT>> v2 = ROL::makePtr<ConDualStdVector<RealT>>(v2_ptr);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(*v1, *v2, *gd, *vc, *xtest, augtol);

    // Define problem.
    vl->zero(); vc->zero(); g->zero();
    ROL::Ptr<ROL::Problem<RealT>>
      problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x,g);
    problem->addConstraint("Equality",constr,ROL::staticPtrCast<ROL::Vector<RealT>>(vl),ROL::staticPtrCast<ROL::Vector<RealT>>(vc));
    problem->finalize(false,true,*outStream);

    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Augmented Lagrangian";
    parlist.sublist("General").set("Output Level",1);
    parlist.sublist("Step").set("Type",stepname);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Ptr<ROL::Solver<RealT>> solver
      = ROL::makePtr<ROL::Solver<RealT>>(problem,parlist);

    // Run Algorithm
    //(*x_ptr)[0] = 3.0; (*x_ptr)[1] = 2.0; (*x_ptr)[2] = 2.0; (*x_ptr)[3] = 1.0; (*x_ptr)[4] = 1.0;
    //(*x_ptr)[0] = -5.0; (*x_ptr)[1] = -5.0; (*x_ptr)[2] = -5.0; (*x_ptr)[3] = -6.0; (*x_ptr)[4] = -6.0;
    solver->solve(*outStream);

    // Compute Error
    *outStream << "\nReference solution x_r =\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[4] << "\n";
    *outStream << "\nOptimal solution x =\n";
    *outStream << std::scientific << "  " << (*x_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[4] << "\n";
    x->axpy(-1.0, *sol);
    RealT abserr = x->norm();
    RealT relerr = abserr/sol->norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
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

