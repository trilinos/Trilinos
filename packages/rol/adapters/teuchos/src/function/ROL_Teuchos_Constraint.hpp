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

#ifndef ROL_TEUCHOS_EQUALITYCONSTRAINT_H
#define ROL_TEUCHOS_EQUALITYCONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_TeuchosVector.hpp"

/** @ingroup func_group
    \class ROL::TeuchosConstraint
    \brief Defines the equality constraint operator interface for 
           Teuchos::SerialDenseVectors

*/

namespace ROL {

template <class Ordinal, class Real>
class TeuchosConstraint : public Constraint<Real> {

  template<class T> using RCP = Teuchos::RCP<T>;
  using SerialDenseVector = Teuchos::SerialDenseVector<Ordinal,Real>;

  RCP<const SerialDenseVector> getVector( const Vector<Real>& x ) {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(x).getVector();
  }
    
  RCP<SerialDenseVector> getVector( Vector<Real>& x ) {
    return dynamic_cast<TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

public:

  virtual ~TeuchosConstraint() {}

  virtual void update( const SerialDenseVector &x, bool flag = true, int iter = -1 ) {}  

  using Constraint<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    auto xp = getVector(x);
    update(*xp,flag,true);
  }


  virtual void value( SerialDenseVector &c, const SerialDenseVector &x, Real &tol ) = 0;

  using Constraint<Real>::value;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    auto cp = getVector(c);
    auto xp = getVector(x);
    value(*cp,*xp,tol);
  }


  virtual void applyJacobian( SerialDenseVector &jv, const SerialDenseVector &v, 
                              const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosConstraint): applyJacobian not implemented!");
  }


  using Constraint<Real>::applyJacobian;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                     const Vector<Real> &x, Real &tol) {

    auto jvp = getVector(jv);
    auto vp  = getVector(v);
    auto xp  = getVector(x);
    try {
      applyJacobian(*jvp,*vp,*xp,tol);      
    } 
    catch (std::exception &e ){
      Constraint<Real>::applyJacobian(jv,v,x,tol);
    }
  }

  virtual void applyAdjointJacobian( SerialDenseVector &ajv, const SerialDenseVector &v, 
                                      const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosConstraint): applyAdjointJacobian not implemented!");
  }


  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                            const Vector<Real> &x, Real &tol) {

    auto ajvp = getVector(ajv);
    auto vp   = getVector(v);
    auto xp   = getVector(x);

    try {
      applyAdjointJacobian(*ajvp,*vp,*xp,tol);      
    } 
    catch (std::exception &e ){
      Constraint<Real>::applyAdjointJacobian(ajv,v,x,tol);
    }
  }

  virtual void applyAdjointHessian( SerialDenseVector &ahuv, const SerialDenseVector &u,
                                    const SerialDenseVector &v, const SerialDenseVector &x,
                                    Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::TeuchosConstraint) : applyAdjointHessian not implemented!");
  }


  using Constraint<Real>::applyAdjointHessian;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) {

    auto ahuvp = getVector(ahuv);
    auto up    = getVector(u);
    auto vp    = getVector(v);
    auto xp    = getVector(x);

    try {
      applyAdjointHessian( *ahuvp, *up, *vp, *xp, tol );
    }
    catch (std::exception &e) {
      Constraint<Real>::applyAdjointHessian(ahuv,u,v,x,tol);
    }   

  }

  virtual std::vector<Real> solveAugmentedSystem( SerialDenseVector &v1, SerialDenseVector &v2,
                                                  const SerialDenseVector &b1, const SerialDenseVector &b2,
                                                  const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::TeuchosConstraint) : solveAugmentedSystem not implemented!");
  }


  using Constraint<Real>::solveAugmentedSystem;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                         const Vector<Real> &b1, const Vector<Real> &b2,
                                         const Vector<Real> &x, Real &tol) {
    auto v1p = getVector(v1);
    auto v2p = getVector(v2);
    auto b1p = getVector(b1);
    auto b2p = getVector(b2);
    auto xp  = getVector(x);

    try {
      return solveAugmentedSystem( *v1p, *v2p, *b1p, *b2p, *xp, tol );
    }
    catch (std::exception &e) {
      return Constraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,x,tol); 
    } 
  }

  virtual void applyPreconditioner( SerialDenseVector &pv, const SerialDenseVector &v,
                                    const SerialDenseVector &x, const SerialDenseVector &g, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::TeuchosConstraint) : applyPreconditioner not implemented!");
  }


  using Constraint<Real>::applyPreconditioner;
  void applyPreconditioner(Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x,
                           const Vector<Real> &g, Real &tol) {
    auto Pvp = getVector(Pv);
    auto vp  = getVector(v);
    auto xp  = getVector(x);
    auto gp  = getVector(g);

    try {
      applyPreconditioner( *Pvp, *vp, *xp, *gp, tol );
    }
    catch (std::exception &e) {
      Constraint<Real>::applyPreconditioner(Pv,v,x,g,tol);
    }
  }

}; // class TeuchosConstraint

} // namespace ROL


#endif // ROL_TEUCHOS_EQUALITYCONSTRAINT_H
