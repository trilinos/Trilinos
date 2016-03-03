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
    \brief Shows how to minimize Rosenbrock's function using Newton-Krylov.
*/

#define USE_HESSVEC 1

#include "ROL_Rosenbrock.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "Teuchos_oblackholestream.hpp"
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
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptDualStdVector<Real> >  dual_vec_;

public:

OptStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const ROL::Vector<Real> &x ) {
  using Teuchos::RCP;  using Teuchos::dyn_cast; 
  RCP<const vector> xvalptr = dyn_cast<const OptStdVector>(x).getVector(); 
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
  using Teuchos::RCP;  using Teuchos::dyn_cast;
  RCP<const vector> xvalptr = dyn_cast<const OptStdVector>(x).getVector(); 
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

Teuchos::RCP<ROL::Vector<Real> > clone() const {
  return Teuchos::rcp( new OptStdVector( Teuchos::rcp( new std::vector<Element>(std_vec_->size()) ) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<std::vector<Element> > getVector() {
  return std_vec_;
}

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  using Teuchos::RCP;  using Teuchos::rcp;
  RCP<vector> e_rcp = rcp( new vector(std_vec_->size(),0.0) );
  (*e_rcp)[i] = 1.0;
  
  RCP<V> e = rcp( new OptStdVector( e_rcp ) );
  
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new OptDualStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

  typedef std::vector<Element> vector;
  typedef ROL::Vector<Real>    V;

  typedef typename vector::size_type uint;

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptStdVector<Real> >  dual_vec_;

public:

OptDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const ROL::Vector<Real> &x ) {
  using Teuchos::RCP;  using Teuchos::dyn_cast;  
  RCP<const vector> xvalptr = dyn_cast<const OptDualStdVector>(x).getVector(); 

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
  using Teuchos::RCP;  using Teuchos::dyn_cast; 
  RCP<const vector> xvalptr = dyn_cast<const OptDualStdVector>(x).getVector(); 
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

Teuchos::RCP<ROL::Vector<Real> > clone() const {
  return Teuchos::rcp( new OptDualStdVector( Teuchos::rcp( new std::vector<Element>(std_vec_->size()) ) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<std::vector<Element> > getVector() {
  return std_vec_;
}

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  using Teuchos::RCP;  using Teuchos::rcp;
  RCP<vector> e_rcp = rcp( new vector(std_vec_->size(),0.0) );
  (*e_rcp)[i] = 1.0;
  
  RCP<V> e = rcp( new OptDualStdVector( e_rcp ) );
  return e;
}

int dimension() const {return std_vec_->size();}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new OptStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class OptDualStdVector


/*** End of declaration of two vector spaces. ***/






int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {

    ROL::ZOO::Objective_Rosenbrock<RealT, OptStdVector<RealT>, OptDualStdVector<RealT> > obj;
    int dim = 100; // Set problem dimension. Must be even.

    // Define algorithm.
    Teuchos::ParameterList parlist;
    std::string stepname = "Trust Region";
    parlist.sublist("Step").sublist(stepname).set("Subproblem Solver", "Truncated CG");
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit",10);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance",1e-2);
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance",1e-4);
    parlist.sublist("General").sublist("Secant").set("Use as Hessian",true);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Algorithm<RealT> algo(stepname,parlist);

    // Iteration Vector
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    // Set Initial Guess
    for (int i=0; i<dim/2; i++) {
      (*x_rcp)[2*i]   = -1.2;
      (*x_rcp)[2*i+1] =  1.0;
      (*g_rcp)[2*i]   = 0;
      (*g_rcp)[2*i+1] = 0;
    }

    OptStdVector<RealT> x(x_rcp); // Iteration Vector
    OptDualStdVector<RealT> g(g_rcp); // zeroed gradient vector in dual space

    Teuchos::RCP<std::vector<RealT> > aa_rcp = Teuchos::rcp( new std::vector<RealT> (1, 1.0) );
    OptDualStdVector<RealT> av(aa_rcp);
    Teuchos::RCP<std::vector<RealT> > bb_rcp = Teuchos::rcp( new std::vector<RealT> (1, 2.0) );
    OptDualStdVector<RealT> bv(bb_rcp);
    Teuchos::RCP<std::vector<RealT> > cc_rcp = Teuchos::rcp( new std::vector<RealT> (1, 3.0) );
    OptDualStdVector<RealT> cv(cc_rcp);
    std::vector<RealT> std_vec_err = av.checkVector(bv,cv,true,*outStream);

    // Run Algorithm
    std::vector<std::string> output = algo.run(x,g, obj, true, *outStream);

    // Get True Solution
    Teuchos::RCP<std::vector<RealT> > xtrue_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 1.0) );
    OptStdVector<RealT> xtrue(xtrue_rcp); 
   
    // Compute Errors
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute solution error: " << abserr;
    *outStream << std::scientific << "\n   Relative solution error: " << relerr;
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
    Teuchos::RCP<std::vector<RealT> > vec_err_rcp = Teuchos::rcp( new std::vector<RealT> (std_vec_err) );
    ROL::StdVector<RealT> vec_err(vec_err_rcp);
    *outStream << std::scientific << "\n   Linear algebra error: " << vec_err.norm() << std::endl;
    if ( vec_err.norm() > 1e2*ROL::ROL_EPSILON<RealT>() ) {
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

