// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrCommFactory_hpp
#define __TSQR_Trilinos_TsqrCommFactory_hpp

/// \file TsqrCommFactory.hpp
/// \brief Factory for TSQR messenger objects.
///
/// \warning TSQR users should <i>not</i> include this file directly.

#include "Teuchos_RCP.hpp"
#include "Tsqr_MessengerBase.hpp"

namespace TSQR {
  namespace Trilinos {

    // Forward declaration. TsqrFactory.hpp includes
    // TsqrCommFactory.hpp, and TsqrTypeAdaptor.hpp includes
    // TsqrFactory.hpp.  We resolve the dependency with a forward
    // declaration of TsqrTypeAdaptor here.
    template< class S, class LO, class GO, class MV >
    class TsqrTypeAdaptor;

    /// \class CommFactory
    /// \brief Factory for TSQR messenger objects.
    ///
    template< class S, class LO, class GO, class MV >
    class CommFactory {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef MV multivector_type;

      typedef TsqrTypeAdaptor< S, LO, GO, MV >    type_adaptor;
      typedef typename type_adaptor::comm_type    comm_type;
      typedef typename type_adaptor::comm_ptr     comm_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >  scalar_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase< LO > > ordinal_messenger_ptr;

      /// Given a raw pointer to a communicator object (implementing
      /// the underlying distributed-memory communication protocol),
      /// return (through the last two reference arguments) wrappers
      /// suitable for TSQR.
      virtual void
      makeMessengers (const comm_ptr& comm,
                      scalar_messenger_ptr& scalarMessenger,
                      ordinal_messenger_ptr& ordinalMessenger) = 0;

      //! Virtual destructor for memory safety.
      virtual ~CommFactory () {}
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrCommFactory_hpp
