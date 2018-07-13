// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
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

