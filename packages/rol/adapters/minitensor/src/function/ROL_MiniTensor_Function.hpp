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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_Function_hpp)
#define ROL_MiniTensor_Function_hpp

#include "Intrepid2_MiniTensor.h"

namespace ROL {

///
/// Function base class that defines the interface to Mini Solvers.
///
template<typename FunctionDerived, typename S>
struct Function_Base
{
public:
  Function_Base()
  {
  }

  ///
  /// By default use merit function 0.5 dot(residual,residual)
  /// as the target to optimize if only the residual is provided.
  ///
  template<typename T, Intrepid2::Index N>
  T
  value(FunctionDerived & f, Intrepid2::Vector<T, N> const & x);

  ///
  /// By default compute gradient with AD from value().
  ///
  template<typename T, Intrepid2::Index N>
  Intrepid2::Vector<T, N>
  gradient(FunctionDerived & f, Intrepid2::Vector<T, N> const & x);

  ///
  /// Defined explicitly.
  ///
  template<typename T, Intrepid2::Index N>
  Intrepid2::Vector<T, N>
  residual(FunctionDerived & f, Intrepid2::Vector<T, N> const & x);

  ///
  /// By default compute Hessian with AD from gradient().
  ///
  template<typename T, Intrepid2::Index N>
  Intrepid2::Tensor<T, N>
  hessian(FunctionDerived & f, Intrepid2::Vector<T, N> const & x);

  ///
  /// Signal that something has gone horribly wrong.
  ///
  bool
  failed{false};

};
} // namespace ROL

#endif // ROL_MiniTensor_Function_hpp
