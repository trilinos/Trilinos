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

/*! \file  example_01.cpp
    \brief Shows how to solve the equality constrained NLP
           from Nocedal/Wright, 2nd edition, page 574, example 18.2.
*/

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_CompositeStep.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

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
ROL::SharedPointer<std::vector<Element> >  std_vec_;
mutable ROL::SharedPointer<OptDualStdVector<Real> >  dual_vec_;

public:

OptStdVector(const ROL::SharedPointer<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPointer) {}

void plus( const ROL::Vector<Real> &x ) {
    
  
  ROL::SharedPointer<const vector> xp = dynamic_cast<const OptStdVector&>(x).getVector();
 
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

    

  ROL::SharedPointer<const vector> xp = dynamic_cast<const OptStdVector&>(x).getVector();
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

ROL::SharedPointer<ROL::Vector<Real> > clone() const {
  return ROL::makeShared<OptStdVector>( ROL::makeShared<std::vector<Element>>(std_vec_->size()) );
}

ROL::SharedPointer<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::SharedPointer<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::SharedPointer<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::SharedPointer<vector> e_ptr = ROL::makeShared<vector>(std_vec_->size(),0.0);
  ROL::SharedPointer<V> e = ROL::makeShared<OptStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makeShared<OptDualStdVector<Real>>( ROL::makeShared<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

typedef std::vector<Element>       vector;
typedef ROL::Vector<Real>          V;
typedef typename vector::size_type uint;

private:
ROL::SharedPointer<std::vector<Element> >  std_vec_;
mutable ROL::SharedPointer<OptStdVector<Real> >  dual_vec_;

public:

OptDualStdVector(const ROL::SharedPointer<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPointer) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const OptDualStdVector&>(x).getVector();
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
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const OptDualStdVector&>(x).getVector();
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

ROL::SharedPointer<ROL::Vector<Real> > clone() const {
  return ROL::makeShared<OptDualStdVector>( ROL::makeShared<std::vector<Element>>(std_vec_->size()) );
}

ROL::SharedPointer<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::SharedPointer<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::SharedPointer<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::SharedPointer<vector> e_ptr = ROL::makeShared<vector>( std_vec_->size(), 0.0 );
  ROL::SharedPointer<V> e = ROL::makeShared<OptDualStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makeShared<OptStdVector<Real>>( ROL::makeShared<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

}; // class OptDualStdVector


// Constraint space.
template <class Real, class Element>
class ConStdVector : public ROL::Vector<Real> {

typedef std::vector<Element>       vector;
typedef ROL::Vector<Real>          V;
typedef typename vector::size_type uint;

private:
ROL::SharedPointer<std::vector<Element> >  std_vec_;
mutable ROL::SharedPointer<ConDualStdVector<Real> >  dual_vec_;

public:

ConStdVector(const ROL::SharedPointer<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPointer) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const ConStdVector&>(x).getVector();
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
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const ConStdVector&>(x).getVector();
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

ROL::SharedPointer<ROL::Vector<Real> > clone() const {
  return ROL::makeShared<ConStdVector>( ROL::makeShared<std::vector<Element>>(std_vec_->size()));
}

ROL::SharedPointer<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::SharedPointer<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::SharedPointer<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::SharedPointer<vector> e_ptr = ROL::makeShared<vector>(std_vec_->size(), 0.0);
  ROL::SharedPointer<V> e = ROL::makeShared<ConStdVector>( e_ptr );
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makeShared<ConDualStdVector<Real>>( ROL::makeShared<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element>       vector;
  typedef ROL::Vector<Real>          V;
  typedef typename vector::size_type uint;

private:

ROL::SharedPointer<std::vector<Element> >  std_vec_;
mutable ROL::SharedPointer<ConStdVector<Real> >  dual_vec_;

public:

ConDualStdVector(const ROL::SharedPointer<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(ROL::nullPointer) {}

void plus( const ROL::Vector<Real> &x ) {
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const ConDualStdVector&>(x).getVector();
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
    
  ROL::SharedPointer<const vector> xp = dynamic_cast<const ConDualStdVector&>(x).getVector();
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

ROL::SharedPointer<ROL::Vector<Real> > clone() const {
  return ROL::makeShared<ConDualStdVector>( ROL::makeShared<std::vector<Element>>(std_vec_->size()));
}

ROL::SharedPointer<const std::vector<Element> > getVector() const {
  return std_vec_;
}

ROL::SharedPointer<std::vector<Element> > getVector() {
  return std_vec_;
}

ROL::SharedPointer<ROL::Vector<Real> > basis( const int i ) const {
    
  ROL::SharedPointer<vector> e_ptr = ROL::makeShared<vector>(std_vec_->size(),0.0);
  ROL::SharedPointer<V> e = ROL::makeShared<ConDualStdVector>(e_ptr);
  (*e_ptr)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = ROL::makeShared<ConStdVector<Real>>( ROL::makeShared<std::vector<Element>>(*std_vec_) );
  return *dual_vec_;
}

}; // class ConDualStdVector

/*** End of declaration of four vector spaces. ***/


int main(int argc, char *argv[]) {

  typedef std::vector<RealT> vector;
  typedef vector::size_type  uint;

    

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makeSharedFromRef(std::cout);
  else
    outStream = ROL::makeSharedFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {

    ROL::SharedPointer<ROL::Objective<RealT> > obj;
    ROL::SharedPointer<ROL::Constraint<RealT> > constr;
    ROL::SharedPointer<vector> x_ptr = ROL::makeShared<vector>(0, 0.0);
    ROL::SharedPointer<vector> sol_ptr = ROL::makeShared<vector>(0, 0.0);
    OptStdVector<RealT> x(x_ptr);      // Iteration vector.
    OptStdVector<RealT> sol(sol_ptr);  // Reference solution vector.

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained <RealT, OptStdVector<RealT>, OptDualStdVector<RealT>, ConStdVector<RealT>, ConDualStdVector<RealT> > (obj, constr, x, sol);

    // Run derivative checks, etc.
    uint dim = 5;
    uint nc = 3;
    RealT left = -1e0, right = 1e0;
    ROL::SharedPointer<vector> xtest_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> g_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> d_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> gd_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> v_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> vc_ptr = ROL::makeShared<vector>(nc, 0.0);
    ROL::SharedPointer<vector> vl_ptr = ROL::makeShared<vector>(nc, 0.0);
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

    ROL::SharedPointer<vector> v1_ptr = ROL::makeShared<vector>(dim, 0.0);
    ROL::SharedPointer<vector> v2_ptr = ROL::makeShared<vector>(nc, 0.0);
    OptStdVector<RealT> v1(v1_ptr);
    ConDualStdVector<RealT> v2(v2_ptr);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(v1, v2, gd, vc, xtest, augtol);
    

    // Define algorithm.
    Teuchos::ParameterList parlist;
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Algorithm<RealT> algo(stepname, parlist);

    // Run Algorithm
    vl.zero();
    //(*x_ptr)[0] = 3.0; (*x_ptr)[1] = 2.0; (*x_ptr)[2] = 2.0; (*x_ptr)[3] = 1.0; (*x_ptr)[4] = 1.0;
    //(*x_ptr)[0] = -5.0; (*x_ptr)[1] = -5.0; (*x_ptr)[2] = -5.0; (*x_ptr)[3] = -6.0; (*x_ptr)[4] = -6.0;
    algo.run(x, g, vl, vc, *obj, *constr, true, *outStream);

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

