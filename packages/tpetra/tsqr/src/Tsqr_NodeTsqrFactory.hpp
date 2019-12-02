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
// ************************************************************************
//@HEADER

#ifndef TSQR_NODETSQRFACTORY_HPP
#define TSQR_NODETSQRFACTORY_HPP

#include "Tsqr_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

#ifdef HAVE_KOKKOSTSQR_TBB
#  include "TbbTsqr.hpp"
#endif // HAVE_KOKKOSTSQR_TBB

#include "Tsqr_KokkosNodeTsqr.hpp"
#include "Tsqr_SequentialTsqr.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include <stdexcept>


namespace TSQR {

  /// \class NodeTsqrFactory
  /// \brief Factory for creating an instance of the right NodeTsqr
  ///   subclass.
  /// \author Mark Hoemmen
  ///
  /// \tparam Scalar The type of entries in the matrices to factor.
  /// \tparam LocalOrdinal The type of local indices in the matrices
  ///   to factor.
  /// \tparam Device Kokkos::Device specialization used by the
  ///   matrices to factor.
  ///
  /// This class maps from (Scalar, LocalOrdinal, Device), to the
  /// corresponding NodeTsqr subclass.  It lets you construct a
  /// default ParameterList for that NodeTsqr subclass, as well as an
  /// instance of the NodeTsqr subclass.  It also provides type
  /// aliases for template metaprogramming.
  ///
  /// The "right" NodeTsqr subclass is a function of Device, and
  /// possibly also of the other template parameters.
  ///
  /// \note If this class does <i>not</i> have a partial
  ///   specialization for your Device type, it defaults to use
  ///   SequentialTsqr.  That class does <i>not</i> use threads, and
  ///   only knows how to deal with host data; it cannot handle GPU
  ///   device-resident data.  Thus, it may perform poorly.
  template<class Scalar, class LocalOrdinal, class Device>
  class NodeTsqrFactory {
  public:
    //! The NodeTsqr subclass corresponding to the Kokkos Node type.
    using node_tsqr_type = SequentialTsqr<LocalOrdinal, Scalar>;
  };
} // namespace TSQR

#endif // TSQR_NODETSQRFACTORY_HPP
