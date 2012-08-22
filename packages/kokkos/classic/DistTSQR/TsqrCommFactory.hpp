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

#ifndef __TSQR_Trilinos_TsqrCommFactory_hpp
#define __TSQR_Trilinos_TsqrCommFactory_hpp

/// \file TsqrCommFactory.hpp
/// \brief Factory for TSQR messenger objects.
///
/// \warning TSQR users should <i>not</i> include this file directly.

#include <Teuchos_RCP.hpp>
#include <Tsqr_MessengerBase.hpp>

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
