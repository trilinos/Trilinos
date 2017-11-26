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

#ifndef ROL_STDEQUALITY_CONSTRAINT_H
#define ROL_STDEQUALITY_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::StdConstraint
    \brief Defines the equality constraint operator interface for StdVectors

*/

namespace ROL {

template <class Real>
class StdConstraint : public virtual Constraint<Real> {
public:

  virtual ~StdConstraint() {}


  using Constraint<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    update(*(xs.getVector()),flag,true);
  }

  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}  


  using Constraint<Real>::value;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    StdVector<Real> cs = dynamic_cast<StdVector<Real>&>(c);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    value(*(cs.getVector()),*(xs.getVector()),tol);
  }

  virtual void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) = 0;




  using Constraint<Real>::applyJacobian;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                             const Vector<Real> &x, Real &tol) {
    StdVector<Real> jvs = dynamic_cast<StdVector<Real>&>(jv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    try {
      applyJacobian(*(jvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);      
    } 
    catch (std::exception &e ){
      Constraint<Real>::applyJacobian(jv,v,x,tol);
    }
  }

  virtual void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                              const std::vector<Real> &x, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdConstraint): applyJacobian not implemented!");
  }



  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv,     const Vector<Real> &v,
                                    const Vector<Real> &x, Real &tol) {
    StdVector<Real> ajvs = dynamic_cast<StdVector<Real>&>(ajv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    try {
       applyAdjointJacobian(*(ajvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);      
    } 
    catch (std::exception &e ){
      Constraint<Real>::applyAdjointJacobian(ajv,v,x,tol);
    }
  }

   virtual void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                                      const std::vector<Real> &x, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdConstraint): applyAdjointJacobian not implemented!");
  }



  using Constraint<Real>::applyAdjointHessian;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) {
    StdVector<Real> ahuvs = dynamic_cast<StdVector<Real>&>(ahuv);
    const StdVector<Real> us = dynamic_cast<const StdVector<Real>&>(u);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    try {
      applyAdjointHessian( *(ahuvs.getVector()), *(us.getVector()), *(vs.getVector()), 
                           *(xs.getVector()), tol );
    }
    catch (std::exception &e) {
      Constraint<Real>::applyAdjointHessian(ahuv,u,v,x,tol);
    }   

  }

  virtual void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                                    const std::vector<Real> &v, const std::vector<Real> &x,
                                    Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::StdConstraint) : applyAdjointHessian not implemented!");
  }




  using Constraint<Real>::solveAugmentedSystem;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                         const Vector<Real> &b1, const Vector<Real> &b2,
                                         const Vector<Real> &x, Real &tol) {
    StdVector<Real> v1s = dynamic_cast<StdVector<Real>&>(v1);
    StdVector<Real> v2s = dynamic_cast<StdVector<Real>&>(v2);
    const StdVector<Real> b1s = dynamic_cast<const StdVector<Real>&>(b1);
    const StdVector<Real> b2s = dynamic_cast<const StdVector<Real>&>(b2);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    try {
      return solveAugmentedSystem( *(v1s.getVector()), *(v2s.getVector()), *(b1s.getVector()),
                            *(b2s.getVector()), *(xs.getVector()), tol );
    }
    catch (std::exception &e) {
      return Constraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,x,tol); 
    }
  }

  virtual std::vector<Real> solveAugmentedSystem( std::vector<Real> &v1, std::vector<Real> &v2,
                                                  const std::vector<Real> &b1, const std::vector<Real> &b2,
                                                  const std::vector<Real> &x, Real tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::StdConstraint) : solveAugmentedSystem not implemented!");
    return std::vector<Real>();
  }


  using Constraint<Real>::applyPreconditioner;
  void applyPreconditioner(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x,
                           const Vector<Real> &g, Real &tol) {
    StdVector<Real> pvs = dynamic_cast<StdVector<Real>&>(pv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    const StdVector<Real> gs = dynamic_cast<const StdVector<Real>&>(g);
    try {
      applyPreconditioner( *(pvs.getVector()), *(vs.getVector()), *(xs.getVector()),
                           *(gs.getVector()), tol );
    }
    catch (std::exception &e) {
      Constraint<Real>::applyPreconditioner(pv,v,x,g,tol);
    }
  }

  virtual void applyPreconditioner( std::vector<Real> &pv, const std::vector<Real> &v,
                                    const std::vector<Real> &x, const std::vector<Real> &g, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
      ">>> ERROR (ROL::StdConstraint) : applyPreconditioner not implemented!");
  }


}; // class StdConstraint

} // namespace ROL


#endif
