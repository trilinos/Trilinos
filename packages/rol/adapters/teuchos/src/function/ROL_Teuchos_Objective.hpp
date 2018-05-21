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

#ifndef ROL_TEUCHOS_OBJECTIVE_H
#define ROL_TEUCHOS_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_TeuchosVector.hpp"

/** @ingroup func_group
    \class TeuchosObjective
    \brief Specializes the Objective interface for objective functions
    that operate on TeuchosVector
*/


namespace ROL {

template <class Ordinal,  class Real>
class TeuchosObjective : public Objective<Real> {

  template<class T> using RCP = Teuchos::RCP<T>;

  using SerialDenseVector = Teuchos::SerialDenseVector<Ordinal,Real>;

  RCP<const SerialDenseVector> getVector( const Vector<Real>& x ) {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(x).getVector();
  }
    
  RCP<SerialDenseVector> getVector( Vector<Real>& x ) {
    return dynamic_cast<TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

public:

  virtual ~TeuchosObjective() {}

  virtual void update( const SerialDenseVector &x, bool flag = true, int iter = -1 ) {}

  using Objective<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    auto xp = getVector(x);
    update(*xp,flag,true);
  }

  virtual Real value( const SerialDenseVector &x, Real &tol ) = 0;

  using Objective<Real>::value;
  Real value( const Vector<Real> &x, Real &tol ) {
    auto xp = getVector(x);
    return value(*xp,tol);
  }

  virtual void gradient( SerialDenseVector &g, const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (TeuchosObjective): gradient not implemented!");
  }

  using Objective<Real>::gradient;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    auto gp = getVector(g);
    auto xp = getVector(x);

    try {
      gradient(*gp,*xp,tol);
    }
    catch (std::exception &e) {
      Objective<Real>::gradient(g,x,tol);
    }
  }

  virtual Real dirDeriv( const SerialDenseVector &x, const SerialDenseVector &d, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (TeuchosObjective): dirDeriv not implemented!");
  }

  using Objective<Real>::dirDeriv;
  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    auto xp = getVector(x);
    auto dp = getVector(d);
    try {
      return dirDeriv(*xp,*dp,tol);
    }
    catch (std::exception &e) {
      return Objective<Real>::dirDeriv(x,d,tol);
    }
  }

  virtual void hessVec( SerialDenseVector &hv, const SerialDenseVector &v, 
                        const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (TeuchosObjective): hessVec not implemented!");
  }

  using Objective<Real>::hessVec;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    auto hvp = getVector(hv);
    auto vp  = getVector(v);
    auto xp  = getVector(x);
    try {
      hessVec(*hvp,*vp,*xp,tol);
    }
    catch (std::exception &e) {
      Objective<Real>::hessVec(hv,v,x,tol);
    }
  }

  virtual void invHessVec( SerialDenseVector &hv, const SerialDenseVector &v, 
                           const SerialDenseVector &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (TeuchosObjective): invHessVec not implemented!");
  }

  using Objective<Real>::invHessVec;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, 
                   const Vector<Real> &x, Real &tol ) {
    auto hvp = getVector(hv);
    auto vp  = getVector(v);
    auto xp  = getVector(x);
    invHessVec(*hvp,*vp,*xp,tol);
  }

  virtual void precond( SerialDenseVector &Pv, const SerialDenseVector &v, 
                        const SerialDenseVector &x, Real &tol ) {
    Pv.assign(v);
  }

  using Objective<Real>::precond;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, 
                const Vector<Real> &x, Real &tol ) {
    auto Pvp = getVector(Pv);
    auto vp  = getVector(v);
    auto xp  = getVector(x);
   
    precond(*Pvp,*vp,*xp,tol);
  }
}; // class TeuchosObjective

} // namespace ROL

#endif // ROL_TEUCHOS_OBJECTIVE_H
