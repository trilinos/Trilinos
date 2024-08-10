// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to minimize Rosenbrock's function using Newton-Krylov.
*/

#define USE_HESSVEC 1

#include "ROL_Rosenbrock.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <iostream>

typedef double RealT;


/*** Declare two vector spaces. ***/

// Forward declarations:

template <class Real, class Element=Real>
class OptStdVector;  // Optimization space.

template <class Real, class Element=Real>
class OptDualStdVector;  // Dual optimization space.


// Vector space definitions:

// Optimization space.
template <class Real, class Element>
class OptStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element> vector;
  typedef ROL::Vector<Real>    V;

  typedef typename vector::size_type uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<OptDualStdVector<Real> >  dual_vec_;

public:

OptStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptStdVector&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
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
    
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptStdVector&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
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
  (*e_ptr)[i] = 1.0;
  
  ROL::Ptr<V> e = ROL::makePtr<OptStdVector>( e_ptr );
  
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

  typedef std::vector<Element> vector;
  typedef ROL::Vector<Real>    V;

  typedef typename vector::size_type uint;

private:
ROL::Ptr<std::vector<Element> >  std_vec_;
mutable ROL::Ptr<OptStdVector<Real> >  dual_vec_;

public:

OptDualStdVector(const ROL::Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPtr) {}

void plus( const ROL::Vector<Real> &x ) {
      
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptDualStdVector&>(x).getVector(); 

  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
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
     
  ROL::Ptr<const vector> xvalptr = dynamic_cast<const OptDualStdVector&>(x).getVector(); 
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
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
  (*e_ptr)[i] = 1.0;
  
  ROL::Ptr<V> e = ROL::makePtr<OptDualStdVector>( e_ptr );
  return e;
}

int dimension() const {return std_vec_->size();}

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


/*** End of declaration of two vector spaces. ***/






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

    // Define objective function.
    ROL::Ptr<ROL::ZOO::Objective_Rosenbrock<RealT, OptStdVector<RealT>, OptDualStdVector<RealT>>> obj
      = ROL::makePtr<ROL::ZOO::Objective_Rosenbrock<RealT, OptStdVector<RealT>, OptDualStdVector<RealT>>>();
    int dim = 100; // Set problem dimension. Must be even.

    // Iteration Vector
    ROL::Ptr<std::vector<RealT>> x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> g_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    // Set Initial Guess
    for (int i=0; i<dim/2; i++) {
      (*x_ptr)[2*i]   = -1.2;
      (*x_ptr)[2*i+1] =  1.0;
      (*g_ptr)[2*i]   = 0;
      (*g_ptr)[2*i+1] = 0;
    }
    ROL::Ptr<OptStdVector<RealT>>     x = ROL::makePtr<OptStdVector<RealT>>(x_ptr);     // Iteration Vector
    ROL::Ptr<OptDualStdVector<RealT>> g = ROL::makePtr<OptDualStdVector<RealT>>(g_ptr); // zeroed gradient vector in dual space

    // Check vector.
    ROL::Ptr<std::vector<RealT> > aa_ptr = ROL::makePtr<std::vector<RealT>>(1, 1.0);
    OptDualStdVector<RealT> av(aa_ptr);
    ROL::Ptr<std::vector<RealT> > bb_ptr = ROL::makePtr<std::vector<RealT>>(1, 2.0);
    OptDualStdVector<RealT> bv(bb_ptr);
    ROL::Ptr<std::vector<RealT> > cc_ptr = ROL::makePtr<std::vector<RealT>>(1, 3.0);
    OptDualStdVector<RealT> cv(cc_ptr);
    std::vector<RealT> std_vec_err = av.checkVector(bv,cv,true,*outStream);

    // Build optimization problem.
    ROL::Ptr<ROL::Problem<RealT>> problem
      = ROL::makePtr<ROL::Problem<RealT>>(obj,x,g);
    problem->finalize(false,true,*outStream);

    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Trust Region";
    parlist.sublist("Step").set("Type",stepname);
    //parlist.sublist("Step").sublist(stepname).sublist("Descent Method").set("Type","Quasi-Newton Method");
    parlist.sublist("Step").sublist(stepname).set("Subproblem Solver", "Truncated CG");
    parlist.sublist("General").set("Output Level",1);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance",1e-2);
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit",10);
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance",1e-4);
    parlist.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");
    parlist.sublist("General").sublist("Secant").set("Use as Hessian",true);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Ptr<ROL::Solver<RealT>> solver
      = ROL::makePtr<ROL::Solver<RealT>>(problem,parlist);

    // Run Algorithm
    solver->solve(*outStream);

    // Get True Solution
    ROL::Ptr<std::vector<RealT> > xtrue_ptr = ROL::makePtr<std::vector<RealT>>(dim, 1.0);
    OptStdVector<RealT> xtrue(xtrue_ptr); 
   
    // Compute Errors
    x->axpy(-1.0, xtrue);
    RealT abserr = x->norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute solution error: " << abserr;
    *outStream << std::scientific << "\n   Relative solution error: " << relerr;
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
    ROL::Ptr<std::vector<RealT> > vec_err_ptr = ROL::makePtr<std::vector<RealT>>(std_vec_err);
    ROL::StdVector<RealT> vec_err(vec_err_ptr);
    *outStream << std::scientific << "\n   Linear algebra error: " << vec_err.norm() << std::endl;
    if ( vec_err.norm() > 1e2*ROL::ROL_EPSILON<RealT>() ) {
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

