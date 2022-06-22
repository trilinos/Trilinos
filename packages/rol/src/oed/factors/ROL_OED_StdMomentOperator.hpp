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

#ifndef ROL_OED_STD_COVARIANCE_OPERATOR_HPP
#define ROL_OED_STD_COVARIANCE_OPERATOR_HPP

#include "ROL_OED_MomentOperator.hpp"
#include "ROL_LinearAlgebra.hpp"
#include "ROL_LAPACK.hpp"
#include "ROL_BLAS.hpp"

#include <algorithm>
#include <functional>
#include <sstream>

namespace ROL {
namespace OED {

template<typename Real>
class StdMomentOperator : public MomentOperator<Real> {
private:
  const Ptr<LAPACK<int,Real>>                                lapack_;
  const Ptr<BLAS<int,Real>>                                    blas_;
  LA::Matrix<Real>                     M_, Minv_, U_, V_, Xpred_, P_;
  std::vector<Real>                                    sval_, sfull_;
  bool isBuilt_, isFactorized_, isSet_, isFullSet_, useSVD_, isPset_;
  std::vector<LA::Matrix<Real>>                       Xdata_, Xfull_;
  int                                                          nobs_;

  void initialize(int nfactors);
  void build(const Vector<Real> &p);
  void factorize(const Vector<Real> &p);

protected:
  /***************************************************************************/
  /* Begin Accessor Functions                                                */
  /***************************************************************************/
  virtual std::vector<Real>& getData(Vector<Real> &x) const;
  virtual const std::vector<Real>& getConstData(const Vector<Real> &x) const;
  /***************************************************************************/
  /* End Accessor Functions                                                  */
  /***************************************************************************/

  using MomentOperator<Real>::get;
  using MomentOperator<Real>::set;
  using MomentOperator<Real>::getNoise;
  using MomentOperator<Real>::getLocalDesign;
  using MomentOperator<Real>::getConstLocalDesign;
  using MomentOperator<Real>::sumAll;
  using MomentOperator<Real>::applyPerturbation;
  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  StdMomentOperator(RegressionType regType = LEASTSQUARES,
                    bool homNoise = true,
                    const Ptr<Noise<Real>> &noise = nullPtr);

  Ptr<MomentOperator<Real>> clone() const;

  void update(const Vector<Real> &p, UpdateType type, int iter = -1);

  void applyInverse(Vector<Real> &Mx,
              const Vector<Real> &x,
              const Vector<Real> &p);

  void apply(Vector<Real> &Mc,
             const Vector<Real> &c,
             const Vector<Real> &p);

  void applyDeriv(Vector<Real> &Mc,
                  const Vector<Real> &c,
                  const Vector<Real> &p);

  void applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v);

  void setFactors(const Ptr<Factors<Real>> &factors);

  void setPerturbation(const Ptr<LinearOperator<Real>> &pOp);

  Real logDeterminant(const Vector<Real> &z);

}; // class StdMomentOperator

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_StdMomentOperator_Def.hpp"

#endif
