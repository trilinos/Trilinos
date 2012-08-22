//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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

#ifndef __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp

/// \file TsqrFactory_TbbTsqr.hpp
///
/// \warning Trilinos users should _not_ include this file directly.

#include "Kokkos_ConfigDefs.hpp" // HAVE_KOKKOSCLASSIC_TBB

#ifdef HAVE_KOKKOSCLASSIC_TBB
#  include "TbbTsqr.hpp"
#endif // HAVE_KOKKOSCLASSIC_TBB


namespace TSQR {
  namespace Trilinos {

#ifdef HAVE_KOKKOSCLASSIC_TBB
    /// \class TbbTsqrFactory
    /// \brief Subclass of TsqrFactory that uses \c TSQR::TBB::TbbTsqr.
    /// \author Mark Hoemmen
    ///
    /// \tparam LO "LocalOrdinal": the type of indices into the
    ///   node-local part of the matrix.
    ///
    /// \tparam S "Scalar": the type of entries in the node-local part
    ///   of the matrix.
    ///
    /// All of this class' public methods, other than the constructor
    /// and destructor, are implemented in the parent class.
    template<class LO, class S>
    class TbbTsqrFactory :
      public TsqrFactory<LO, S, TSQR::TBB::TbbTsqr<LO, S>, DistTsqr<LO, S> > {
    public:
      // Help C++ pull in the typedefs from the base class.  C++ needs
      // help when both the base and the derived classes are
      // templated.
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::dist_tsqr_type dist_tsqr_type;
      typedef typename base_type::tsqr_type tsqr_type;
      typedef typename base_type::scalar_messenger_type scalar_messenger_type;

      TbbTsqrFactory () {}
      virtual ~TbbTsqrFactory () {}
    };
#endif // HAVE_KOKKOSCLASSIC_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
