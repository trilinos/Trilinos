// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for the equality constrained NLP:

            minimize  x^2 + y^2
          subject to  (x-2)^2 + y^2 = 1

    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_PARABOLOIDCIRCLE_HPP
#define ROL_PARABOLOIDCIRCLE_HPP

#include "ROL_TestProblem.hpp"
#include "ROL_StdVector.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Objective function:
             f(x,y) = x^2 + y^2
   */
  template< class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real> >
  class Objective_ParaboloidCircle : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;

  typedef typename vector::size_type uint;
   

  private:

    template<class VectorType>
    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const VectorType&>(x).getVector();
    }

    template<class VectorType>
    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<VectorType&>(x).getVector();
    }

  public:
    Objective_ParaboloidCircle() {}

    Real value( const Vector<Real> &x, Real &tol ) {
 
      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x); 

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective value): "
                                                                   "Primal vector x must be of length 2.");

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];

      Real val = x1*x1 + x2*x2;

      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<vector> gp = getVector<XDual>(g); 

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective gradient): "
                                                                   " Primal vector x must be of length 2."); 

      n = gp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective gradient): "
                                                                   "Gradient vector g must be of length 2."); 

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];

      Real two(2);

      (*gp)[0] = two*x1;
      (*gp)[1] = two*x2;
    }

    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<XPrim>(v);
      ROL::Ptr<vector> hvp = getVector<XDual>(hv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective hessVec): "
                                                                   "Primal vector x must be of length 2."); 

      n = vp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective hessVec): "
                                                                   "Input vector v must be of length 2."); 

      n = hvp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, objective hessVec): "
                                                                   "Output vector hv must be of length 2."); 

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];

      Real two(2);

      (*hvp)[0] = two*v1;
      (*hvp)[1] = two*v2;
    }

  };


  /** \brief constraint c(x,y) = (x-2)^2 + y^2 - 1.
   */
  template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
  class Constraint_ParaboloidCircle : public Constraint<Real> {

    typedef std::vector<Real> vector;
    typedef Vector<Real>      V;

    typedef typename vector::size_type uint;

  private:
    template<class VectorType>
    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const VectorType&>(x).getVector();
    }

    template<class VectorType> 
    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<VectorType&>(x).getVector(); 
    }

  public:
    Constraint_ParaboloidCircle() {}

    void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<vector> cp = getVector<CPrim>(c);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint value): "
                                                                   "Primal vector x must be of length 2.");

      uint m = cp->size();
      ROL_TEST_FOR_EXCEPTION( (m != 1), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint value): "
                                                                   "Constraint vector c must be of length 1.");

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
 
      Real one(1), two(2);

      (*cp)[0] = (x1-two)*(x1-two) + x2*x2 - one;
    }
  
    void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<XPrim>(v);
      ROL::Ptr<vector> jvp = getVector<CPrim>(jv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyJacobian): "
                                                                   "Primal vector x must be of length 2.");

      uint d = vp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyJacobian): "
                                                                   "Input vector v must be of length 2.");
      d = jvp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 1), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyJacobian): "
                                                                   "Output vector jv must be of length 1.");
      
      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];

      Real two(2);

      (*jvp)[0] = two*(x1-two)*v1 + two*x2*v2;
    } //applyJacobian

    void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<CDual>(v);
      ROL::Ptr<vector> ajvp = getVector<XDual>(ajv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointJacobian): "
                                                                   "Primal vector x must be of length 2.");

      uint d = vp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 1), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointJacobian): "
                                                                   "Input vector v must be of length 1.");

      d = ajvp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointJacobian): "
                                                                   "Output vector ajv must be of length 2.");
      
      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];

      Real v1 = (*vp)[0];

      Real two(2);

      (*ajvp)[0] = two*(x1-two)*v1;
      (*ajvp)[1] = two*x2*v1;

    } //applyAdjointJacobian

    void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      
      bool useFD = true;

      if (useFD) {
        Constraint<Real>::applyAdjointHessian( ahuv, u, v, x, tol );
      }
      else {
        
        ROL::Ptr<const vector> xp = getVector<XPrim>(x);
        ROL::Ptr<const vector> up = getVector<CDual>(u);
        ROL::Ptr<const vector> vp = getVector<XPrim>(v);
        ROL::Ptr<vector> ahuvp = getVector<XDual>(ahuv);

        uint n = xp->size();
        ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointHessian): "
                                                                     "Primal vector x must be of length 2.");

        n = vp->size();
        ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointHessian): "
                                                                     "Direction vector v must be of length 2.");

        n = ahuvp->size();
        ROL_TEST_FOR_EXCEPTION( (n != 2), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointHessian): "
                                                                     "Output vector ahuv must be of length 2.");
        uint d = up->size();
        ROL_TEST_FOR_EXCEPTION( (d != 1), std::invalid_argument, ">>> ERROR (ROL_ParaboloidCircle, constraint applyAdjointHessian): "
                                                                     "Dual constraint vector u must be of length 1.");
        
        Real v1 = (*vp)[0];
        Real v2 = (*vp)[1];

        Real u1 = (*up)[0];

        Real two(2);

        (*ahuvp)[0] = two*u1*v1;
        (*ahuvp)[1] = two*u1*v2;
      }
    } //applyAdjointHessian

  };


  template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
  class getParaboloidCircle : public TestProblem<Real> {
    typedef std::vector<Real> vector;
    typedef typename vector::size_type uint;
  public:
    getParaboloidCircle(void) {}

    Ptr<Objective<Real>> getObjective(void) const {
      // Instantiate objective function.
      return ROL::makePtr<Objective_ParaboloidCircle<Real,XPrim,XDual>>();
    }
     
    Ptr<Vector<Real>> getInitialGuess(void) const {
      uint n = 2;
      // Get initial guess.
      ROL::Ptr<vector> x0p  = makePtr<vector>(n,0.0);
      (*x0p)[0] = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      (*x0p)[1] = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      return makePtr<XPrim>(x0p);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      uint n = 2;
      // Get solution.
      Real zero(0), one(1);
      ROL::Ptr<vector> solp = makePtr<vector>(n,0.0);
      (*solp)[0] = one;
      (*solp)[1] = zero;
      return makePtr<XPrim>(solp);
    }

    Ptr<Constraint<Real>> getEqualityConstraint(void) const {
      // Instantiate constraints.
      return ROL::makePtr<Constraint_ParaboloidCircle<Real,XPrim,XDual,CPrim,CDual>>();
    }

    Ptr<Vector<Real>> getEqualityMultiplier(void) const {
      ROL::Ptr<vector> lp = makePtr<vector>(1,0.0);
      return makePtr<CDual>(lp);
    }
  };

} // End ZOO Namespace
} // End ROL Namespace

#endif
