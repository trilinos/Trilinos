//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOSSOLVERFACTORY_BELOS_HPP
#define BELOSSOLVERFACTORY_BELOS_HPP

#include "Belos_Details_registerSolverFactory.hpp"

// Note that this file is currently included by BelosSolverFactory.hpp
// to maintain backwards compatibility. We don't include it here because
// gcc won't resolve the circular includes for namespacing
// #include "BelosSolverFactory.hpp"

#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"

namespace Belos {

class BelosSolverFactory : public Impl::SolverFactoryParent<double,MultiVec<double>,Operator<double>>
{
  public:
    BelosSolverFactory() {
      Details::registerSolverFactory();
    };
};

namespace Impl {

template<>
class SolverFactorySelector<double,MultiVec<double>,Operator<double>> {
  public:
    typedef BelosSolverFactory type;
};

#ifdef HAVE_TEUCHOS_COMPLEX
template<>
class SolverFactorySelector<std::complex<double>,MultiVec<std::complex<double>>,Operator<std::complex<double>>> {
  public:
    typedef BelosSolverFactory type;
};
#endif

} // namespace Impl
} // namespace Belos

#endif // BELOSSOLVERFACTORY_BELOS_HPP
