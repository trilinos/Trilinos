// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_RANDOM_UQ_PCE_HPP
#define KOKKOS_RANDOM_UQ_PCE_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_KOKKOSALGORITHMS)

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE_Contiguous.hpp"
#include "Kokkos_Random.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos random functions for Sacado::UQ::PCE scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

  template<class Generator, class Storage>
  struct rand<Generator,Sacado::UQ::PCE<Storage> > {
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef typename Scalar::value_type BaseScalar;
    typedef rand<Generator,BaseScalar> BaseRand;

    KOKKOS_INLINE_FUNCTION
    static Scalar max() { return BaseRand::max(); }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen) {
      return BaseRand::draw(gen);
    }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen, const Scalar& range) {
      return BaseRand::draw(gen, range.coeff(0));
    }

    KOKKOS_INLINE_FUNCTION
    static Scalar draw(Generator& gen, const Scalar& start, const Scalar& end) {
      return BaseRand::draw(gen, start.coeff(0), end.coeff(0));
    }
  };

  template<class S, class L, class D, class M, class RandomPool>
  void fill_random( const View<Sacado::UQ::PCE<S>**,L,D,M>& a,
                    RandomPool g,
                    const Sacado::UQ::PCE<S>& begin,
                    const Sacado::UQ::PCE<S>& end ) {
    typedef View<Sacado::UQ::PCE<S>**,L,D,M> Vector;
    typename Vector::flat_array_type a_flat = a;
    fill_random( a_flat, g, begin.fastAccessCoeff(0), end.fastAccessCoeff(0) );
  }

}

#endif

#endif /* #ifndef KOKKOS_RANDOM_UQ_PCE_HPP */
