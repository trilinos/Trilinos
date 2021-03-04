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

#ifndef ROL2_TYPEU_NONLINEARCG_DEF_H
#define ROL2_TYPEU_NONLINEARCG_DEF_H

namespace ROL2 {
namespace TypeU {

template<typename Real>
NonlinearCG<Real>::NonlinearCG(       ParameterList&          parlist,
                                const Ptr<NonlinearCG<Real>>& nlcg = nullPtr)
  : nlcg_(nlcg), enlcg_(Type::UserDefined) {

  // Initialize secant object
  auto& llist = parlist.sublist("Step").sublist("Line Search");
  auto& dmlist = llist.sublist("Descent Method");
  if ( nlcg == nullPtr ) {
    ncgName_ = dmlist.get("Nonlinear CG Type","Oren-Luenberger");
    enlcg_   = StringToENonlinearCG(ncgName_);
    nlcg_    = makePtr<NonlinearCG<Real>>(enlcg_);
  }
  else 
    ncgName_ = dmlist.get("User Defined Nonlinear CG Name",
                          "Unspecified User Define Nonlinear CG Method");
}

template<typename Real>
void NonlinearCG<Real>::compute(       Vector<Real>&    s, 
                                       Real&            snorm, 
                                       Real&            sdotg, 
                                       int&             iter, 
                                       int&             flag,
                                 const Vector<Real>&    x, 
                                 const Vector<Real>&    g, 
                                       Objective<Real>& obj) {
  nlcg_->run(s,g,x,obj);
  sdotg = -s.apply(g);
  if (sdotg >= static_cast<Real>(0)) {
    s.set(g.dual());
    sdotg = -s.apply(g);
  }
  s.scale(static_cast<Real>(-1));
  snorm = s.norm();
  iter  = 0;
  flag  = 0;
}

template<typename Real>
void NonlinearCG<Real>::writeName( std::ostream& os ) const override {
  os << ncgName_ << " Nonlinear CG";
}

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_NONLINEARCG_DEF_H






