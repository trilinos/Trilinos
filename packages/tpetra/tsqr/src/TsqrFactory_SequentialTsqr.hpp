// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

/// \file TsqrFactory_SequentialTsqr.hpp
/// \brief Declaration and definition of SequentialTsqrFactory.

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

namespace TSQR {
  namespace Trilinos {

    /// \class SequentialTsqrFactory
    ///
    /// Subclass of \c TsqrFactory that knows how to instantiate
    /// \c SequentialTsqr as the intranode TSQR implementation.
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
    class SequentialTsqrFactory :
      public TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> > {
    public:
      // This class' parent class.
      typedef TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> > base_type;

      // Pull in the typedefs from the base class.  C++ doesn't do
      // this when both the base and the derived classes are
      // templated.
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::dist_tsqr_type dist_tsqr_type;
      typedef typename base_type::tsqr_type tsqr_type;
      typedef typename base_type::scalar_messenger_type scalar_messenger_type;

      SequentialTsqrFactory () {}
      virtual ~SequentialTsqrFactory () {}
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
