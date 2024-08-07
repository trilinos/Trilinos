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
#include "ROL_Algorithm.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

/*** Declare four vector spaces. ***/

// Forward declarations:

template <class Real, class Element=Real>
class OptStdVector;      // Optimization space.

template <class Real, class Element=Real>
class OptDualStdVector;  // Dual optimization space.

template <class Real, class Element=Real>
class ConStdVector;      // Constraint space.

template <class Real, class Element=Real>
class ConDualStdVector;  // Dual constraint space.

// Vector space definitions:

// Optimization space.
template <class Real, class Element>
class OptStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>           vector;
  typedef ROL::Vector<Real>              V;
  typedef OptStdVector<Real,Element>     PV;
  typedef OptDualStdVector<Real,Element> DV;  

  typedef typename vector::size_type     uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;

public:

OptStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::Ptr<const vector> xp;
  try {
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) {
    xp = dynamic_cast<const DV&>(x).getVector();
 }
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
    
  Real val = 0;
  ROL::Ptr<const vector> xp;
  try {  // duality pairing
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) { // inner product
    xp = dynamic_cast<const DV&>(x).getVector();
  }
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
  
  return ROL::makePtr<OptStdVector>( ROL::makePtr<vector>(std_vec_->size()) );
}

ROL::Ptr<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::Ptr<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(),0.0);
  ROL::Ptr<V> e = ROL::makePtr<PV>(e_ptr);
  (*e_ptr)[i] = 1.0; 
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>           vector;
  typedef ROL::Vector<Real>              V;
  typedef OptStdVector<Real,Element>     DV;
  typedef OptDualStdVector<Real,Element> PV;  

  typedef typename vector::size_type     uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;

public:

OptDualStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
    

  ROL::Ptr<const vector> xp;
  try {
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) {
    xp = dynamic_cast<const DV&>(x).getVector();
  }
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
    
  Real val = 0;

  ROL::Ptr<const vector> xp;
 
  try {  // duality pairing
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) { // inner product
    xp = dynamic_cast<const DV&>(x).getVector(); 
  }

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
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(),0.0);
  ROL::Ptr<V> e = ROL::makePtr<PV>(e_ptr);
  (*e_ptr)[i] = 1.0; 
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

}; // class OptDualStdVector


// Constraint space.
template <class Real, class Element>
class ConStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>           vector;
  typedef ROL::Vector<Real>              V;
  typedef ConStdVector<Real,Element>     PV;
  typedef ConDualStdVector<Real,Element> DV;  

  typedef typename vector::size_type     uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;

public:

ConStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::Ptr<const vector> xp;
  try {
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) {
    xp = dynamic_cast<const DV&>(x).getVector();
  }

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
    
  ROL::Ptr<const vector> xp; 
  Real val = 0;
  try {  // duality pairing
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) { // inner product
    xp = dynamic_cast<const DV&>(x).getVector();
  }
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
    
  ROL::Ptr<vector> e_ptr = ROL::makePtr<vector>(std_vec_->size(),0.0);
  ROL::Ptr<V> e = ROL::makePtr<PV>(e_ptr);
  (*e_ptr)[i] = 1.0; 
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>           vector;
  typedef ROL::Vector<Real>              V;
  typedef ConStdVector<Real,Element>     DV;
  typedef ConDualStdVector<Real,Element> PV;  

  typedef typename vector::size_type     uint;

private:

ROL::Ptr<std::vector<Element> >  std_vec_;

public:

ConDualStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
    

  ROL::Ptr<const vector> xp;

  try {
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) {
    xp = dynamic_cast<const DV&>(x).getVector();
  }
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
    
  Real val = 0;
  ROL::Ptr<const vector> xp;
  try {  // duality pairing
    xp = dynamic_cast<const PV&>(x).getVector();
  }
  catch (const std::bad_cast &e) { // inner product
    xp = dynamic_cast<const DV&>(x).getVector();
  }

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
  ROL::Ptr<V> e = ROL::makePtr<PV>(e_ptr);
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

}; // class ConDualStdVector

/*** End of declaration of four vector spaces. ***/


int main(int argc, char *argv[]) {

  typedef std::vector<RealT>         vector;
  typedef typename vector::size_type uint;

    

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
    OptStdVector<RealT> x(x_ptr);      // Iteration vector.
    OptStdVector<RealT> sol(sol_ptr);  // Reference solution vector.

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained <RealT, OptStdVector<RealT>, OptDualStdVector<RealT>, ConStdVector<RealT>, ConDualStdVector<RealT> > SEC;
    obj = SEC.getObjective();
    constr = SEC.getEqualityConstraint();
    x.set(*SEC.getInitialGuess());
    sol.set(*SEC.getSolution());

    // Run derivative checks, etc.
    RealT left = -1e0, right = 1e0;
    ROL::Ptr<vector> xtest_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> g_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> d_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> gd_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v_ptr  = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> vc_ptr = ROL::makePtr<vector>(nc, 0.0);
    ROL::Ptr<vector> vl_ptr = ROL::makePtr<vector>(nc, 0.0);
    OptStdVector<RealT> xtest(xtest_ptr);
    OptDualStdVector<RealT> g(g_ptr);
    OptStdVector<RealT> d(d_ptr);
    OptDualStdVector<RealT> gd(gd_ptr);
    OptStdVector<RealT> v(v_ptr);
    ConStdVector<RealT> vc(vc_ptr);
    ConDualStdVector<RealT> vl(vl_ptr);
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
    obj->checkGradient(xtest, g, d, true, *outStream);                      *outStream << "\n"; 
    obj->checkHessVec(xtest, g, v, true, *outStream);                       *outStream << "\n";
    obj->checkHessSym(xtest, g, d, v, true, *outStream);                    *outStream << "\n";
    constr->checkApplyJacobian(xtest, v, vc, true, *outStream);             *outStream << "\n";
    constr->checkApplyAdjointJacobian(xtest, vl, vc, g, true, *outStream);  *outStream << "\n";
    constr->checkApplyAdjointHessian(xtest, vl, d, g, true, *outStream);    *outStream << "\n";

    ROL::Ptr<vector> v1_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> v2_ptr = ROL::makePtr<vector>(nc, 0.0);
    OptStdVector<RealT> v1(v1_ptr);
    ConDualStdVector<RealT> v2(v2_ptr);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(v1, v2, gd, vc, xtest, augtol);
    

    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::CompositeStep<RealT>>(parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Run Algorithm
    vl.zero();
    //(*x_ptr)[0] = 3.0; (*x_ptr)[1] = 2.0; (*x_ptr)[2] = 2.0; (*x_ptr)[3] = 1.0; (*x_ptr)[4] = 1.0;
    //(*x_ptr)[0] = -5.0; (*x_ptr)[1] = -5.0; (*x_ptr)[2] = -5.0; (*x_ptr)[3] = -6.0; (*x_ptr)[4] = -6.0;

    std::vector<std::string> output = algo.run(x, g, vl, vc, *obj, *constr, false);
    for ( uint i = 0; i < output.size(); i++ ) {
      *outStream << output[i];
    }

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
    x.axpy(-1.0, sol);
    RealT abserr = x.norm();
    RealT relerr = abserr/sol.norm();
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

