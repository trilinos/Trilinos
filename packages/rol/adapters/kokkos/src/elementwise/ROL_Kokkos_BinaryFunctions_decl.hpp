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
#ifndef ROL_KOKKOS_BINARYFUNCTIONS_DECL_HPP
#define ROL_KOKKOS_BINARYFUNCTIONS_DECL_HPP

namespace ROL {

namespace Elementwise {

template<typename Real, typename Device>
class KokkosBinaryFunction {
public:
  virtual ~KokkosBinaryFunction() = default;
  virtual void apply(       ::ROL::KokkosVector<Real,Device>&, 
                      const ::ROL::KokkosVector<Real,Device>&  ) const {}

  static std::unique_ptr<KokkosBinaryFunction> 
  create( const ::ROL::Elementwise::BinaryFunction<Real>& uf );

private:

  struct Factory : public BinaryFunction<Real>::Visitor {
    void visit( const ::ROL::Elementwise::Axpy<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Aypx<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Divide<Real>&          ) override;
    void visit( const ::ROL::Elementwise::DivideAndInvert<Real>& ) override;
    void visit( const ::ROL::Elementwise::Greater<Real>&         ) override;
    void visit( const ::ROL::Elementwise::Lesser<Real>&          ) override;
    void visit( const ::ROL::Elementwise::Max<Real>&             ) override;
    void visit( const ::ROL::Elementwise::Min<Real>&             ) override;
    void visit( const ::ROL::Elementwise::Multiply<Real>&        ) override;
    void visit( const ::ROL::Elementwise::Plus<Real>&            ) override;
    void visit( const ::ROL::Elementwise::Set<Real>&             ) override;

    std::unique_ptr<KokkosBinaryFunction> bfun_;
  }; // Factory
}; // KokkosBinaryFunction


template<typename Real, typename Device>
class KokkosAxpy : public KokkosBinaryFunction<Real,Device> {
public:
  KokkosAxpy( Real alpha ) : alpha_(alpha) {}
  virtual ~KokkosAxpy() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
private:
  Real alpha_;
};


template<typename Real, typename Device>
class KokkosAypx : public KokkosBinaryFunction<Real,Device> {
public:
  KokkosAypx( Real alpha ) : alpha_(alpha) {}
  virtual ~KokkosAypx() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
private:
  Real alpha_;
};

template<typename Real, typename Device>
class KokkosDivide : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosDivide() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosDivideAndInvert : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosDivideAndInvert() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosGreater : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosGreater() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosLesser : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosLesser() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
}; 

template<typename Real, typename Device>
class KokkosMax : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosMax() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosMin : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosMin() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosMultiply : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosMultiply() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosPlus : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosPlus() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

template<typename Real, typename Device>
class KokkosSet : public KokkosBinaryFunction<Real,Device> {
public:
  virtual ~KokkosSet() = default;
  void apply(       ::ROL::KokkosVector<Real,Device>&,
              const ::ROL::KokkosVector<Real,Device>&  ) const override;
};

} // namespace Elementwise
} // namespace ROL

#endif  // ROL_KOKKOS_BINARYFUNCTIONS_DECL_HPP

