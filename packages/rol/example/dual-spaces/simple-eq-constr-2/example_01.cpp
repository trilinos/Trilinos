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
#include "ROL_CompositeStepSQP.hpp"
#include "ROL_Algorithm.hpp"
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

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

OptStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
  try {
    OptStdVector &ex = Teuchos::dyn_cast<OptStdVector>(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) {
    OptDualStdVector<Real, Element> &ex = Teuchos::dyn_cast<OptDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
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
  try {  // duality pairing
    OptDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) { // inner product
    OptStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
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

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public ROL::Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

OptDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
  try {
    OptDualStdVector &ex = Teuchos::dyn_cast<OptDualStdVector>(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) {
    OptStdVector<Real, Element> &ex = Teuchos::dyn_cast<OptStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
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
  try {  // duality pairing
    OptStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) { // inner product
    OptDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
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

}; // class OptDualStdVector


// Constraint space.
template <class Real, class Element>
class ConStdVector : public ROL::Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

ConStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
  try {
    ConStdVector &ex = Teuchos::dyn_cast<ConStdVector>(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) {
    ConDualStdVector<Real, Element> &ex = Teuchos::dyn_cast<ConDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
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
  try {  // duality pairing
    ConDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) { // inner product
    ConStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<ROL::Vector<Real> > clone() const {
  return Teuchos::rcp( new ConStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<ConStdVector> e = Teuchos::rcp( new ConStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public ROL::Vector<Real> {
private:

Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

ConDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

void plus( const ROL::Vector<Real> &x ) {
  try {
    ConDualStdVector &ex = Teuchos::dyn_cast<ConDualStdVector>(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) {
    ConStdVector<Real, Element> &ex = Teuchos::dyn_cast<ConStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*std_vec_)[i] += (*xvalptr)[i];
    }
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
  try {  // duality pairing
    ConStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  catch (const std::bad_cast &e) { // inner product
    ConDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConDualStdVector<Real, Element> >(const_cast <ROL::Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    unsigned dimension  = std_vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*(*xvalptr)[i];
    }
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<ROL::Vector<Real> > clone() const {
  return Teuchos::rcp( new ConDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<ConDualStdVector> e = Teuchos::rcp( new ConDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

}; // class ConDualStdVector

/*** End of declaration of four vector spaces. ***/


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

    Teuchos::RCP<ROL::Objective<RealT> > obj;
    Teuchos::RCP<ROL::EqualityConstraint<RealT> > constr;
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (0, 0.0) );
    Teuchos::RCP<std::vector<RealT> > sol_rcp = Teuchos::rcp( new std::vector<RealT> (0, 0.0) );
    OptStdVector<RealT> x(x_rcp);      // Iteration vector.
    OptStdVector<RealT> sol(sol_rcp);  // Reference solution vector.

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained <RealT, OptStdVector<RealT>, OptDualStdVector<RealT>, ConStdVector<RealT>, ConDualStdVector<RealT> > (obj, constr, x, sol);

    Teuchos::ParameterList parlist;
    // Define Step
    parlist.set("Nominal SQP Optimality Solver Tolerance", 1.e-2);
    ROL::CompositeStepSQP<RealT> step(parlist);

    // Run derivative checks, etc.
    int dim = 5;
    int nc = 3;
    RealT left = -1e0, right = 1e0;
    Teuchos::RCP<std::vector<RealT> > xtest_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > d_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > gd_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > v_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > vc_rcp = Teuchos::rcp( new std::vector<RealT> (nc, 0.0) );
    Teuchos::RCP<std::vector<RealT> > vl_rcp = Teuchos::rcp( new std::vector<RealT> (nc, 0.0) );
    OptStdVector<RealT> xtest(xtest_rcp);
    OptDualStdVector<RealT> g(g_rcp);
    OptStdVector<RealT> d(d_rcp);
    OptDualStdVector<RealT> gd(gd_rcp);
    OptStdVector<RealT> v(v_rcp);
    ConStdVector<RealT> vc(vc_rcp);
    ConDualStdVector<RealT> vl(vl_rcp);
    // set xtest, d, v
    for (int i=0; i<dim; i++) {
      (*xtest_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*gd_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    // set vc, vl
    for (int i=0; i<nc; i++) {
      (*vc_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vl_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    obj->checkGradient(xtest, g, d, true);  *outStream << "\n"; 
    obj->checkHessVec(xtest, g, v, true);  *outStream << "\n";
    obj->checkHessSym(xtest, g, d, v, true);  *outStream << "\n";
    constr->checkApplyJacobian(xtest, v, vc, true);  *outStream << "\n";
    constr->checkApplyAdjointJacobian(xtest, vl, vc, g, true);  *outStream << "\n";
    constr->checkApplyAdjointHessian(xtest, vl, d, g, true);  *outStream << "\n";

    Teuchos::RCP<std::vector<RealT> > v1_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > v2_rcp = Teuchos::rcp( new std::vector<RealT> (nc, 0.0) );
    OptStdVector<RealT> v1(v1_rcp);
    ConDualStdVector<RealT> v2(v2_rcp);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(v1, v2, gd, vc, xtest, augtol);
    
    // Define Status Test
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT ctol  = 1e-12;  // norm of constraint tolerance
    RealT stol  = 1e-18;  // norm of step tolerance
    int   maxit = 1000;    // maximum number of iterations
    ROL::StatusTestSQP<RealT> status(gtol, ctol, stol, maxit);    

    // Define Algorithm
    ROL::DefaultAlgorithm<RealT> algo(step, status, false);

    // Run Algorithm
    vl.zero();
    //(*x_rcp)[0] = 3.0; (*x_rcp)[1] = 2.0; (*x_rcp)[2] = 2.0; (*x_rcp)[3] = 1.0; (*x_rcp)[4] = 1.0;
    //(*x_rcp)[0] = -5.0; (*x_rcp)[1] = -5.0; (*x_rcp)[2] = -5.0; (*x_rcp)[3] = -6.0; (*x_rcp)[4] = -6.0;

    std::vector<std::string> output = algo.run(x, g, vl, vc, *obj, *constr, false);
    for ( unsigned i = 0; i < output.size(); i++ ) {
      std::cout << output[i];
    }

    // Compute Error
    *outStream << "\nReference solution x_r =\n";
    *outStream << std::scientific << "  " << (*sol_rcp)[0] << "\n";
    *outStream << std::scientific << "  " << (*sol_rcp)[1] << "\n";
    *outStream << std::scientific << "  " << (*sol_rcp)[2] << "\n";
    *outStream << std::scientific << "  " << (*sol_rcp)[3] << "\n";
    *outStream << std::scientific << "  " << (*sol_rcp)[4] << "\n";
    *outStream << "\nOptimal solution x =\n";
    *outStream << std::scientific << "  " << (*x_rcp)[0] << "\n";
    *outStream << std::scientific << "  " << (*x_rcp)[1] << "\n";
    *outStream << std::scientific << "  " << (*x_rcp)[2] << "\n";
    *outStream << std::scientific << "  " << (*x_rcp)[3] << "\n";
    *outStream << std::scientific << "  " << (*x_rcp)[4] << "\n";
    x.axpy(-1.0, sol);
    RealT abserr = x.norm();
    RealT relerr = abserr/sol.norm();
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

