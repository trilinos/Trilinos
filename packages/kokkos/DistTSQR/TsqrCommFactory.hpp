//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
