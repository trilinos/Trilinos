// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file TsqrTypeAdaptor.hpp
/// \brief Traits class mapping between multivector type and TSQR
///   implementation types.
///
/// \warning Trilinos users should <i>not</i> include this file directly.

#ifndef __TSQR_Trilinos_TsqrTypeAdaptor_hpp
#define __TSQR_Trilinos_TsqrTypeAdaptor_hpp

#include "Teuchos_RCP.hpp"
#include "TsqrFactory.hpp"
#include "Tsqr.hpp"

namespace TSQR {
namespace Trilinos {
/// \class UndefinedComm
/// \brief Class used to catch undefined specializations of \c TsqrTypeAdaptor.
class UndefinedComm {};

/// \class TsqrTypeAdaptor
/// \brief Traits class mapping between multivector type and TSQR implementation types.
///
/// This traits class maps between the specific multivector type
/// MV, and the corresponding appropriate intranode and internode
/// TSQR implementation classes.
///
/// \tparam S The Scalar type (the type of elements in the matrix
///   for which to compute a QR factorization)
/// \tparam LO The "local ordinal" type, as one would find in
///   Tpetra distributed linear algebra objects.  (In Epetra, the
///   local and global ordinal types are both the same, namely
///   "int".).
/// \tparam GO The "global ordinal" type, as one would find in
///   Tpetra distributed linear algebra objects.  (In Epetra, the
///   local and global ordinal types are both the same, namely
///   "int".).
/// \tparam MV The multivector type.
///
/// This class maps a multivector type to three different classes:
/// - node_tsqr_type, the type of the intranode part of TSQR
/// - dist_tsqr_type, the type of the internode part of TSQR
/// - tsqr_type, which is constrained by the above two
///   definitions; it is the type of the full TSQR implementation,
///   including both intra- and internode components.
///
/// It also gives the appropriate \c TsqrFactory type to use for
/// constructing a \c TsqrAdaptor.
///
/// \note Implementers who want to port TSQR to a new MV class (by
///   mapping to an existing TSQR implementation) should first
///   specialize a new TsqrTypeAdaptor class for that MV.  They
///   should then implement the corresponding \c TsqrAdaptor
///   class.
template <class S, class LO, class GO, class MV>
class TsqrTypeAdaptor {
 public:
  typedef S scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef MV multivector_type;

  /// \typedef node_tsqr_type
  /// \brief Type representing the intranode part of TSQR.
  ///
  /// Defaults to sequential, cache-blocked TSQR.
  typedef TSQR::SequentialTsqr<LO, S> node_tsqr_type;
  typedef Teuchos::RCP<node_tsqr_type> node_tsqr_ptr;

  /// \typedef dist_tsqr_type
  /// \brief Type representing the internode part of TSQR.
  typedef TSQR::DistTsqr<LO, S> dist_tsqr_type;
  typedef Teuchos::RCP<dist_tsqr_type> dist_tsqr_ptr;

  /// \typedef tsqr_type
  /// \brief Type representing the whole TSQR method.
  ///
  /// Depends on \c node_tsqr_type and \c dist_tsqr_type.
  using tsqr_type = TSQR::Tsqr<LO, S>;
  typedef Teuchos::RCP<tsqr_type> tsqr_ptr;

  /// \typedef factory_type
  ///
  /// Type of the \c TsqrFactory object that knows how to construct
  /// \c node_tsqr_type and \c dist_tsqr_type objects.
  typedef TsqrFactory<LO, S, node_tsqr_type, dist_tsqr_type> factory_type;

  /// \typedef comm_type
  ///
  /// Type of the (raw) communicator object used by the given
  /// multivector type.  Communicator objects are always handled
  /// via Teuchos::RCP.  The default is \c UndefinedComm, which
  /// catches missing or partially defined specializations of
  /// TsqrTypeAdaptor as syntax errors.
  typedef UndefinedComm comm_type;
  typedef Teuchos::RCP<const comm_type> comm_ptr;
};

}  // namespace Trilinos
}  // namespace TSQR

#endif  // __TSQR_Trilinos_TsqrTypeAdaptor_hpp
