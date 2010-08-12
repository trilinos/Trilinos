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

#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp

/// \file TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp
///
/// \warning Anasazi users should _not_ include this file directly.
///   Include "AnasaziTpetraAdaptor.hpp" instead.

#include "Tpetra_MultiVector.hpp"
#include "Kokkos_SerialNode.hpp"
#include "TsqrFactory_SequentialTsqr.hpp"
#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO >
    class TsqrTypeAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Kokkos::SerialNode > >
    {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef Tpetra::MultiVector< S, LO, GO, Kokkos::SerialNode > multivector_type;

      typedef TSQR::SequentialTsqr< LO, S > node_tsqr_type;
      typedef TSQR::Tsqr< LO, S, node_tsqr_type > tsqr_type;
      typedef SequentialTsqrFactory< local_ordinal_type, scalar_type > factory_type;
      typedef Teuchos::Comm<int> comm_type;
      typedef Teuchos::RCP< const comm_type > comm_ptr;
    };

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode_hpp
