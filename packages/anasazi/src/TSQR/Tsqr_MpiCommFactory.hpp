// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __TSQR_MpiCommFactory_hpp
#define __TSQR_MpiCommFactory_hpp

#include <mpi.h>
#include <Tsqr_Config.hpp>
#include <Tsqr_MpiMessenger.hpp>
#include <Teuchos_RCP.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace MPI {

    namespace details {

      template< class Scalar >
      Teuchos::RCP< MessengerBase< Scalar > >
      makeMpiComm (MPI_Comm comm)
      {
	return Teuchos::rcp_implicit_cast< MessengerBase< Scalar > >(new MpiMessenger< Scalar > (comm));
      }
    } // namespace details

#ifdef HAVE_MPI_COMM_NETWORK
    /// \fn makeMpiCommNetwork
    /// \brief Return wrapper of MPI_COMM_NETWORK
    template< class Scalar >
    Teuchos::RCP< MessengerBase< Scalar > >
    makeMpiCommNetwork () 
    {
      makeMpiComm (MPI_COMM_NETWORK);
    }
#endif // HAVE_MPI_COMM_NETWORK
    
#ifdef HAVE_MPI_COMM_NODE
    /// \fn makeMpiCommNode
    /// \brief Return wrapper of MPI_COMM_NODE
    template< class Scalar >
    Teuchos::RCP< MessengerBase< Scalar > >
    makeMpiCommNode () 
    {
      makeMpiComm (MPI_COMM_NODE);
    }
#endif // HAVE_MPI_COMM_NODE

    /// \fn makeMpiCommWorld
    /// \brief Return wrapper of MPI_COMM_WORLD
    ///
    /// \warning NOT REENTRANT
    template< class Scalar >
    Teuchos::RCP< MessengerBase< Scalar > >
    makeMpiCommWorld () 
    {
      makeMpiComm (MPI_COMM_WORLD);
    }

  } // namespace MPI
} // namespace TSQR

#endif // __TSQR_MpiCommFactory_hpp

