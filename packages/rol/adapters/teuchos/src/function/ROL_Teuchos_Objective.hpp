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
    \class ROL::TeuchosObjective
    \brief Specializes the ROL::Objective interface for objective functions
    that operate on ROL::TeuchosVector
*/


namespace ROL {

template <class Ordinal,  class Real>
class TeuchosObjective : public Objective<Real> {

  template <typename T> using ROL::SharedPointer = ROL::SharedPointer<T>;

  typedef Teuchos::SerialDenseVector<Ordinal,Real> SDV;
  typedef TeuchosVector<Ordinal,Real>              TV;

public:

  virtual void update( const SDV &x, bool flag = true, int iter = -1 ) {}

  using Objective<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();
    update(*xp,flag,true);
  }

  virtual Real value( const SDV &x, Real &tol ) = 0;

  using Objective<Real>::value;
  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();
    return value(*xp,tol);
  }

  virtual void gradient( SDV &g, const SDV &x, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosObjective): gradient not implemented!");
  }

  using Objective<Real>::gradient;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::SharedPointer<SDV> gp = dynamic_cast<TV&>(g).getVector();
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();

    try {
      gradient(*gp,*xp,tol);
    }
    catch (std::exception &e) {
      Objective<Real>::gradient(g,x,tol);
    }
  }

  virtual Real dirDeriv( const SDV &x, const SDV &d, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosObjective): dirDeriv not implemented!");
  }

  using Objective<Real>::dirDeriv;
  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();
    ROL::SharedPointer<const SDV> dp = dynamic_cast<const TV&>(d).getVector();
   try {
      return dirDeriv(*xp,*dp,tol);
    }
    catch (std::exception &e) {
      return Objective<Real>::dirDeriv(x,d,tol);
    }
  }

  virtual void hessVec( SDV &hv, const SDV &v, const SDV &x, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosObjective): hessVec not implemented!");
  }

  using Objective<Real>::hessVec;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::SharedPointer<SDV> hvp = dynamic_cast<TV&>(hv).getVector();
    ROL::SharedPointer<const SDV> vp = dynamic_cast<const TV&>(v).getVector();
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();
    try {
      hessVec(*hvp,*vp,*xp,tol);
    }
    catch (std::exception &e) {
      Objective<Real>::hessVec(hv,v,x,tol);
    }
  }

  virtual void invHessVec( SDV &hv, const SDV &v, const SDV &x, Real &tol ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosObjective): invHessVec not implemented!");
  }

  using Objective<Real>::invHessVec;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::SharedPointer<SDV> hvp = dynamic_cast<TV&>(hv).getVector();
    ROL::SharedPointer<const SDV> vp = dynamic_cast<const TV&>(v).getVector();
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();
    invHessVec(*hvp,*vp,*xp,tol);
  }

  virtual void precond( SDV &Pv, const SDV &v, const SDV &x, Real &tol ) {
    Pv.assign(v.begin(),v.end());
  }

  using Objective<Real>::precond;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::SharedPointer<SDV> Pvp = dynamic_cast<TV&>(Pv).getVector();
    ROL::SharedPointer<const SDV> vp = dynamic_cast<const TV&>(v).getVector();
    ROL::SharedPointer<const SDV> xp = dynamic_cast<const TV&>(x).getVector();

    precond(*Pvp,*vp,*xp,tol);
  }
}; // class TeuchosObjective

} // namespace ROL

#endif // ROL_TEUCHOS_OBJECTIVE_H
