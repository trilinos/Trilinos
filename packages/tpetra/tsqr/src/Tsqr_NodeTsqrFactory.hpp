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
#ifdef HAVE_KOKKOSTSQR_TBB
#  include "TbbTsqr.hpp"
#endif // HAVE_KOKKOSTSQR_TBB
#include "Tsqr_KokkosNodeTsqr.hpp"
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr_CombineNodeTsqr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include <vector>
#include <string>

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
  private:
    using host_serial_node_tsqr_type =
      SequentialTsqr<LocalOrdinal, Scalar>;
    using host_parallel_node_tsqr_type =
      KokkosNodeTsqr<LocalOrdinal, Scalar>;
    using combine_node_tsqr_type =
      CombineNodeTsqr<LocalOrdinal, Scalar>;

  public:
    using node_tsqr_type = NodeTsqr<LocalOrdinal, Scalar>;

    static Teuchos::RCP<node_tsqr_type>
    getNodeTsqr ()
    {
      using execution_space = typename Device::execution_space;
#ifdef KOKKOS_ENABLE_CUDA
      constexpr bool is_cuda =
        std::is_same<execution_space, Kokkos::Cuda>::value;
#else
      constexpr bool is_cuda = false;
#endif // KOKKOS_ENABLE_CUDA
      if (is_cuda) {
        // FIXME (mfh 02 Dec 2019): We don't yet have a CUDA option.
        // Just run SequentialTsqr (on host) for now.  This need not
        // necessarily rely on UVM, since the adapter can access the
        // host version of the data.
        return Teuchos::rcp (new host_serial_node_tsqr_type);
      }

#ifdef HAVE_KOKKOSTSQR_COMPLEX
      constexpr bool is_complex =
        std::is_same<Scalar, std::complex<double>>::value ||
        std::is_same<Scalar, std::complex<float>>::value;
#else
      constexpr bool is_complex = false;
#endif // HAVE_KOKKOSTSQR_COMPLEX
      if (is_complex) {
        return Teuchos::rcp (new combine_node_tsqr_type);
      }

      execution_space execSpace;
      if (execSpace.concurrency () == 1) {
        return Teuchos::rcp (new host_serial_node_tsqr_type);
      }
      else {
        return Teuchos::rcp (new host_parallel_node_tsqr_type);
      }
    }

    static Teuchos::RCP<node_tsqr_type>
    getNodeTsqr (const std::string& name)
    {
      using Teuchos::rcp;
      if (name == "SequentialTsqr" || name == "Sequential") {
        return rcp (new SequentialTsqr<LocalOrdinal, Scalar>);
      }
      else if (name == "KokkosNodeTsqr" || name == "Kokkos") {
        return rcp (new KokkosNodeTsqr<LocalOrdinal, Scalar>);
      }
      else if (name == "CombineNodeTsqr" || name == "Combine") {
        return rcp (new CombineNodeTsqr<LocalOrdinal, Scalar>);
      }
      else if (name == "Default") {
        return getNodeTsqr ();
      }
      else {
        const char prefix[] = "TSQR::NodeTsqrFactory::getNodeTsqr: ";
        const std::vector<std::string> validNames
          {{"SequentialTsqr",
            "KokkosNodeTsqr",
            "CombineNodeTsqr",
            "Default"}};
        std::ostringstream os;
        os << prefix << "Invalid NodeTsqr subclass name \"" << name
           << "\".  Valid names are: {";
        for (size_t k = 0; k < validNames.size (); ++k) {
          os << "\"" << validNames[k] << "\"";
          if (k + size_t (1) < validNames.size ()) {
            os << ", ";
          }
        }
        os << "}.";
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, os.str ());
      }
    }
  };
} // namespace TSQR

#endif // TSQR_NODETSQRFACTORY_HPP
