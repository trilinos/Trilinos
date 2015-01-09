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

#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

/// \file Kokkos_SerialNode.hpp
/// \brief Declaration and definition of the (now DEPRECATED)
///   KokkosClassic::SerialNode Node type.
/// \warning KokkosClassic::SerialNode has been DEPRECATED.  For a
///   Node with comparable intent, please use
///   Kokkos::Compat::KokkosSerialWrapperNode instead.

#include "Kokkos_ConfigDefs.hpp"

// mfh 08 Jan 2015: Don't enable the contents of this file unless the
// appropriate CMake option is enabled.  This avoids deprecation
// warnings once we deprecate this Node type.
#ifdef HAVE_TPETRACLASSIC_SERIAL
#include <Kokkos_StandardNodeMemoryModel.hpp>
#include "Kokkos_NodeHelpers.hpp"

#ifdef HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT
#  include "KokkosCore_config.h"
#  ifdef KOKKOS_HAVE_SERIAL
#    include "Kokkos_Serial.hpp"
#  endif // KOKKOS_HAVE_SERIAL
#endif // HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT

namespace KokkosClassic {

  /// \brief Node API implementation that uses sequential execution
  ///   for "thread-level parallelism."
  /// \ingroup kokkos_node_api
  /// \warning This class has been DEPRECATED.  For a Node with
  ///   comparable intent, please use
  ///   Kokkos::Compat::KokkosSerialWrapperNode instead.
  class TPETRA_DEPRECATED SerialNode : public StandardNodeMemoryModel {
  public:
    /// \brief This is a "classic" Node type.
    ///
    /// That means we plan to deprecate it with the 11.14 release of
    /// Trilinos, and remove it entirely with the 12.0 release.
    static const bool classic = true;

    //! Constructor; sets default parameters.
    SerialNode ();

    /// \brief Constructor that takes a list of parameters.
    ///
    /// This class doesn't currently use any parameters.
    SerialNode (Teuchos::ParameterList& pl);

    //! Get default parameters for this class.
    static Teuchos::ParameterList getDefaultParameters ();

    /// \brief Skeleton for a parallel for, with a serial implementation.
    ///
    /// See \ref kokkos_node_api "Kokkos Node API"
    template <class WDP>
    static void parallel_for (const int beg, const int end, WDP wd) {
      for (int i = beg; i != end; ++i) {
        wd.execute (i);
      }
    }

    /// \brief Skeleton for a parallel reduction, with a serial implementation.
    ///
    /// See \ref kokkos_node_api "Kokkos Node API"
    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce (const int begin, const int end, WDP wd) {
      typename WDP::ReductionType result = wd.identity ();
      for (int i=begin; i != end; ++i) {
        result = wd.reduce (result, wd.generate (i));
      }
      return result;
    }

    //! No-op for SerialNode.
    inline void sync () const {}

    /// \brief Return the human-readable name of this Node.
    ///
    /// See \ref kokkos_node_api "Kokkos Node API"
    static std::string name ();
  };

#ifdef _MSC_VER
#  pragma warning(push)
// destructor could not be generated because a base class destructor is inaccessible
#  pragma warning(disable : 4624)
#endif // _MSC_VER

  template <> class ArrayOfViewsHelper<SerialNode> :
    public ArrayOfViewsHelperTrivialImpl<SerialNode>
  {};

#ifdef _MSC_VER
#  pragma warning(pop)
#endif // _MSC_VER

} // namespace KokkosClassic

#if defined(HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT) && defined(KOKKOS_HAVE_SERIAL)
namespace Kokkos {
  namespace Compat {
    template <>
    struct NodeDevice<KokkosClassic::SerialNode> {
      typedef Kokkos::Serial type;
    };
  } // namespace Compat
} // namespace Kokkos
#endif // HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT && KOKKOS_HAVE_SERIAL

#endif // HAVE_TPETRACLASSIC_SERIAL
#endif // KOKKOS_SERIALNODE_HPP_
