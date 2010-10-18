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

#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_hpp

/// \file TsqrTypeAdaptor.hpp
///
/// \warning Anasazi users should _not_ include this file directly.

#include <Teuchos_RCP.hpp>
#include <Tsqr_Config.hpp>
#include <TsqrFactory.hpp>
#include <Tsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    class UndefinedComm {};

    /// \class TsqrTypeAdaptor
    ///
    /// \brief Mapping between multivector class MV and appropriate
    /// {intra,inter}-node TSQR classes.
    ///
    /// TsqrAdaptor has to map a multivector type to two different
    /// classes:
    ///
    /// \li node_tsqr_type, responsible for the intranode part of the
    ///   TSQR computations
    /// \li tsqr_type, responsible for the internode part of the TSQR
    ///   computations
    ///
    /// TsqrTypeAdaptor maps from the multivector type MV, to these
    /// two classes.  It also gives the appropriate TsqrFactory type
    /// to use for constructing a TsqrAdaptor. 
    ///
    /// \note Implementers who want to port TSQR to a new MV class (by
    ///   mapping to an existing TSQR implementation) should first
    ///   specialize a new TsqrTypeAdaptor class for that MV.  They
    ///   should then implement the corresponding TsqrAdaptor class.
    template< class S, class LO, class GO, class MV >
    class TsqrTypeAdaptor {
    public:
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef MV multivector_type;

      ///
      /// Type representing the intranode part of TSQR.
      /// Defaults to sequential, cache-blocked TSQR.
      ///
      typedef TSQR::SequentialTsqr< LO, S >  node_tsqr_type;
      typedef Teuchos::RCP< node_tsqr_type > node_tsqr_ptr;

      ///
      /// Type representing the internode part of TSQR.
      ///
      typedef TSQR::DistTsqr< LO, S >        dist_tsqr_type;
      typedef Teuchos::RCP< dist_tsqr_type > dist_tsqr_ptr;

      ///
      /// Type representing the whole TSQR method.
      /// Depends on node_tsqr_type and dist_tsqr_type.
      ///
      typedef TSQR::Tsqr< LO, S, node_tsqr_type, dist_tsqr_type > tsqr_type;
      typedef Teuchos::RCP< tsqr_type >                           tsqr_ptr;

      ///
      /// Type of the TsqrFactory object that knows how to construct
      /// node_tsqr_type and dist_tsqr_type objects.
      ///
      typedef TsqrFactory< LO, S, node_tsqr_type, dist_tsqr_type > factory_type;

      ///
      /// Type of the (raw) communicator object used by the given
      /// multivector type.  Communicator objects are always handled
      /// via Teuchos::RCP.
      ///
      typedef UndefinedComm comm_type;
      typedef Teuchos::RCP< const comm_type > comm_ptr;
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrTypeAdaptor_hpp
