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
#ifndef ROL_KOKKOS_UNARYFUNCTIONS_DECL_HPP
#define ROL_KOKKOS_UNARYFUNCTIONS_DECL_HPP

namespace ROL {

namespace Elementwise {

template<typename Real, typename Device>
class KokkosUnaryFunction {
public:
  virtual ~KokkosUnaryFunction() = default;
  virtual void apply( ::ROL::KokkosVector<Real,Device>& x ) const {}

  static std::unique_ptr<KokkosUnaryFunction> 
  create( const ::ROL::Elementwise::UnaryFunction<Real>& uf );

private:

  struct Factory : public UnaryFunction<Real>::Visitor {
    void visit( const ::ROL::Elementwise::AbsoluteValue<Real>&    ) override;
    void visit( const ::ROL::Elementwise::Fill<Real>&             ) override;
    void visit( const ::ROL::Elementwise::Heaviside<Real>&        ) override;
    void visit( const ::ROL::Elementwise::Logarithm<Real>&        ) override;
    void visit( const ::ROL::Elementwise::Power<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Reciprocal<Real>&       ) override;
    void visit( const ::ROL::Elementwise::Round<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Scale<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Shift<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Sign<Real>&             ) override;
    void visit( const ::ROL::Elementwise::SquareRoot<Real>&       ) override;
    void visit( const ::ROL::Elementwise::ThresholdLower<Real>&   ) override;
    void visit( const ::ROL::Elementwise::ThresholdUpper<Real>&   ) override;
    std::unique_ptr<KokkosUnaryFunction> ufun_;

  }; // Factory

}; // KokkosUnaryFunction


template<typename Real, typename Device>
class KokkosAbsoluteValue : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosAbsoluteValue() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosFill : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosFill( Real value ) : value_(value) {}
  virtual ~KokkosFill() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
  Real get_value() const { return value_; }
private:
  Real value_;
};

template<typename Real, typename Device>
class KokkosHeaviside : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosHeaviside() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosLogarithm : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosLogarithm() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosPower : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosPower( Real exponent ) : exponent_(exponent) {}
  virtual ~KokkosPower() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
private:
  Real exponent_;
};

template<typename Real, typename Device>
class KokkosReciprocal : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosReciprocal() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosRound : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosRound() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosScale : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosScale( Real value ) : value_(value) {}
  virtual ~KokkosScale() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
private:
  Real value_;
};

template<typename Real, typename Device>
class KokkosShift : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosShift( Real value ) : value_(value) {}
  virtual ~KokkosShift() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
private:
  Real value_;
};

template<typename Real, typename Device>
class KokkosSign : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosSign() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosSquareRoot : public KokkosUnaryFunction<Real,Device> {
public:
  virtual ~KokkosSquareRoot() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
};

template<typename Real, typename Device>
class KokkosThresholdLower : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosThresholdLower( Real threshold ) : threshold_(threshold) {}
  virtual ~KokkosThresholdLower() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
private:
  Real threshold_;
};

template<typename Real, typename Device>
class KokkosThresholdUpper : public KokkosUnaryFunction<Real,Device> {
public:
  KokkosThresholdUpper( Real threshold ) : threshold_(threshold) {}
  virtual ~KokkosThresholdUpper() = default;
  void apply( ::ROL::KokkosVector<Real,Device>& ) const override;
private:
  Real threshold_;
};


} // namespace Elementwise
} // namespace ROL

#endif  // ROL_KOKKOS_UNARYFUNCTIONS_DECL_HPP

