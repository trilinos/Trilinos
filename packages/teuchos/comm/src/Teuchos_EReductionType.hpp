// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_EREDUCTIONTYPE_HPP
#define TEUCHOS_EREDUCTIONTYPE_HPP

/// \file Teuchos_EReductionType.hpp
/// \brief Declaration of Teuchos::EReductionType enum, and related functions

#include "Teuchos_config.h"
#include "Teuchos_DLLExportMacro.h"
#ifdef HAVE_TEUCHOS_MPI
#  include <mpi.h> // need this for MPI_Op (see below)
#endif // HAVE_TEUCHOS_MPI

namespace Teuchos {

/// \brief Predefined reduction operations that Teuchos::Comm understands
/// \relates Comm
///
/// Teuchos' Comm class wraps MPI (the Message Passing Interface for
/// distributed-memory parallel programming).  If you do not have MPI,
/// it imitates MPI's functionality, as if you were running with a
/// single "parallel" process.  This means that Teuchos must wrap a
/// subset of MPI's functionality, so that it can build without MPI.
///
/// Comm provides wrappers for \c MPI_Reduce, \c MPI_Allreduce, and
/// other collectives that take a reduction operator \c MPI_Op.
/// Teuchos wraps \c MPI_Op in two different ways.  The first way is
/// this enum, which lets users pick from a set of common predefined
/// \c MPI_Op.  The second way is through Teuchos' wrappers for custom
/// \c MPI_Op, namely ValueTypeReductionOp and ValueTypeReductionOp.
/// Most users should find the reduction operators below sufficient.
enum EReductionType {
  REDUCE_SUM, ///< Sum
  REDUCE_MIN, ///< Min
  REDUCE_MAX, ///< Max
  REDUCE_AND, ///< Logical AND
  REDUCE_BOR  ///< Bitwise OR
};

/// \brief Convert EReductionType to string representation.
/// \relates EReductionType
const char* toString (const EReductionType reductType);

#ifdef HAVE_TEUCHOS_MPI
namespace Details {

/// \brief Get the raw MPI_Op corresponding to the given reduction
///   type enum value.
///
/// \warning This is an implementation detail and not for public use.
///   It only exists when Trilinos was built with MPI.
TEUCHOSCOMM_LIB_DLL_EXPORT MPI_Op getMpiOpForEReductionType (const enum EReductionType reductionType);

} // namespace Details
#endif // HAVE_TEUCHOS_MPI

} // namespace Teuchos

#endif // TEUCHOS_EREDUCTIONTYPE_HPP
