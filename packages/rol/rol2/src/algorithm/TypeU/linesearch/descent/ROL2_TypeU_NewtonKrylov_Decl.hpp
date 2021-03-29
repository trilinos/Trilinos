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

#pragma once
#ifndef ROL2_TYPEU_NEWTONKRYLOV_DECL_HPP
#define ROL2_TYPEU_NEWTONKRYLOV_DECL_HPP

/** @ingroup step_group
    \class ROL2::TypeU::NewtonKrylov
    \brief Provides the interface to compute optimization steps 
           with Krylov approximations to Newton's method using line search.

*/

namespace ROL2 {
namespace TypeU {

template <class Real>
class NewtonKrylov : public DescentDirection<Real>  {
public:

  using KrylovType = typename Krylov<Real>::Type;
  using SecantType = typename Secant<Real>::Type;

  //------------------------------------------------------------
  class HessianNK : public LinearOperator<Real> {
  public:
    HessianNK( const Ptr<Objective<Real>>&    obj,
               const Ptr<const Vector<Real>>& x );

    void apply(       Vector<Real>& Hv,
                const Vector<Real>& v,
                      Real&         tol ) const override;
  private:
    const Ptr<Objective<Real>>    obj_;
    const Ptr<const Vector<Real>> x_;
  }; // class HessianNK

  //------------------------------------------------------------
  class PrecondNK : public LinearOperator<Real> {
  public:
    PrecondNK( const Ptr<Objective<Real>>&    obj,
               const Ptr<const Vector<Real>>& x );

    void apply(       Vector<Real>& Hv,
                const Vector<Real>& v,
                      Real&         tol ) const override;

    void applyInverse(       Vector<Real>& Hv,
                       const Vector<Real>& v,
                             Real&         tol ) const override;
  private:
    const Ptr<Objective<Real>>    obj_;
    const Ptr<const Vector<Real>> x_;
  }; // class HessianNK
  //------------------------------------------------------------


  /** \brief Constructor.

      Standard constructor to build a NewtonKrylovStep object.  Algorithmic
      specifications are passed in through a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  NewtonKrylov( ParameterList& parlist );

  NewtonKrylov(       ParameterList&     parlist, 
                const Ptr<Krylov<Real>>& krylov,
                      Ptr<Secant<Real>>& secant,
                      bool               computeObj = true );
                
  virtual ~NewtonKrylov() = default;

  virtual void compute(       Vector<Real>&    s, 
                              Real&            snorm, 
                              Real&            sdotg, 
                              int&             iter, 
                              int&             flag,
                        const Vector<Real>&    x, 
                        const Vector<Real>&    g, 
                              Objective<Real>& obj) override;

  virtual void update( const Vector<Real>& x, 
                       const Vector<Real>& s,
                       const Vector<Real>& gold, 
                       const Vector<Real>& gnew,
                             Real          snorm, 
                             int           iter ) override;

  virtual void writeName( std::ostream& os ) const { 
    os << "Newton-Krylov Method using " << krylovName_; 
    if( useSecantPrecond_ ) {
      os << " with " << secantName_ << " preconditioning";
    }
  }

private:

  Ptr<Secant<Real>>         secant_;
  Ptr<Krylov<Real>>         krylov_;
  Ptr<LinearOperator<Real>> precond_;

  KrylovType krylovType_;
  SecantType secantType_;

  std::string krylovName_;
  std::string secantName_;

  bool useSecantPrecond_ = false;

}; // class NewtonKrylov

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_NEWTONKRYLOV_DECL_HPP
