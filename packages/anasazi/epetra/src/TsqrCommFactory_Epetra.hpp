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

#ifndef __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp
#define __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp

/// \file TsqrCommFactory_Epetra.hpp
///
/// \warning Anasazi users should _not_ include this file directly.
///   Include "AnasaziEpetraAdaptor.hpp" instead, if it exists.  If it
///   doesn't exist, you're out of luck; use Tpetra instead.

#include "Epetra_MultiVector.hpp"
#include "Tsqr_EpetraMessenger.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    class EpetraCommFactory :
      public CommFactory< double, int, int, Epetra_MultiVector >
    {
    public:
      typedef CommFactory< double, int, int, Epetra_MultiVector > base_type;
      typedef base_type::comm_ptr comm_ptr;
      typedef base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      EpetraCommFactory () {}
      virtual ~EpetraCommFactory () {}

      virtual void
      makeMessengers (const comm_ptr& comm,
		      scalar_messenger_ptr& scalarMessenger,
		      ordinal_messenger_ptr& ordinalMessenger)
      {
	using TSQR::Trilinos::EpetraMessenger;

	scalarMessenger = Teuchos::rcp (new EpetraMessenger< S > (comm));
	ordinalMessenger = Teuchos::rcp (new EpetraMessenger< LO > (comm));
      }
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrCommFactory_Epetra_hpp
