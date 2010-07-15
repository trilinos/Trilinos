// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
#define __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp

/// \file TsqrCommFactory_Tpetra.hpp
///
/// \warning Anasazi users should _not_ include this file directly.
///   Include "AnasaziTpetraAdaptor.hpp" instead.

#include "Tpetra_MultiVector.hpp"
#include "Tsqr_TpetraMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class Node >
    class TpetraCommFactory :
      public CommFactory< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > >
    {
    public:
      // MV def looks circular but is not, because C++ does not pass
      // typedefs from a base class to a derived class when both the
      // base and derived classes are templated.  MV here is just a
      // shorthand so that the base_type typedef fits in one line.
      typedef Tpetra::MultiVector< S, LO, GO, Node > MV;
      typedef CommFactory< S, LO, GO, MV > base_type;

      // C++ doesn't pass typedefs from a templated base class to a
      // templated derived class.  Here, we avoid repeating code from
      // the base class by pulling in its typedefs into the derived
      // class.
      typedef typename base_type::comm_ptr comm_ptr;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      TpetraCommFactory () {}
      virtual ~TpetraCommFactory () {}

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger)
      {
	using TSQR::Trilinos::TpetraMessenger;

	scalarMessenger = Teuchos::rcp (new TpetraMessenger< S > (comm));
	ordinalMessenger = Teuchos::rcp (new TpetraMessenger< LO > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrCommFactory_Tpetra_hpp
