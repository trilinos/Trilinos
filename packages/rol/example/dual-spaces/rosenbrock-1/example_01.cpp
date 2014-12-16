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
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
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

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptDualStdVector<Real> >  dual_vec_;

public:

OptStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const ROL::Vector<Real> &x ) {
  OptStdVector &ex = Teuchos::dyn_cast<OptStdVector>(const_cast <ROL::Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {
  Real val = 0;
  OptStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
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

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<OptStdVector> e = Teuchos::rcp( new OptStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

const ROL::Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new OptDualStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptStdVector<Real> >  dual_vec_;

public:

OptDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const ROL::Vector<Real> &x ) {
  OptDualStdVector &ex = Teuchos::dyn_cast<OptDualStdVector>(const_cast <ROL::Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const ROL::Vector<Real> &x ) const {
  Real val = 0;
  OptDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
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

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<OptDualStdVector> e = Teuchos::rcp( new OptDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
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

    Teuchos::ParameterList parlist;
    // Enumerations
    parlist.set("Descent Type",                           "Quasi-Newton Method");
    parlist.set("Secant Type",                            "Limited-memory BFGS");
    parlist.set("Linesearch Type",                        "Cubic Interpolation");
    parlist.set("Linesearch Curvature Condition",         "Wolfe");
    // Linesearch Parameters
    parlist.set("Maximum Number of Function Evaluations", 20);
    parlist.set("Sufficient Decrease Parameter",          1.e-4);
    parlist.set("Curvature Conditions Parameter",         0.9);
    parlist.set("Backtracking Rate",                      0.5);
    parlist.set("Initial Linesearch Parameter",           1.0);
    parlist.set("User Defined Linesearch Parameter",      false);
    // Krylov Parameters
    parlist.set("Absolute Krylov Tolerance",              1.e-4);
    parlist.set("Relative Krylov Tolerance",              1.e-2);
    parlist.set("Maximum Number of Krylov Iterations",    10);
    // Trust Region Parameters
    parlist.set("Trust-Region Subproblem Solver Type","Truncated CG");
    parlist.set("Use Secant Hessian-Times-A-Vector",true);
    // Define Step
    //ROL::LineSearchStep<RealT> step(parlist);
    ROL::TrustRegionStep<RealT> step(parlist);


    // Define Status Test
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    

    // Define Algorithm
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);

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
    av.checkVector(bv,cv);


    // Run Algorithm
    std::vector<std::string> output = algo.run(x,g, obj, false);
    for ( unsigned i = 0; i < output.size(); i++ ) {
      std::cout << output[i];
    }

    // Get True Solution
    Teuchos::RCP<std::vector<RealT> > xtrue_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 1.0) );
    OptStdVector<RealT> xtrue(xtrue_rcp); 
   
    // Compute Error
    x.axpy(-1.0, xtrue);
    RealT abserr = x.norm();
    RealT relerr = abserr/xtrue.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON) ) {
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

