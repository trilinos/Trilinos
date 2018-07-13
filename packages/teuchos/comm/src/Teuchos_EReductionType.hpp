// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef TEUCHOS_EREDUCTIONTYPE_HPP
#define TEUCHOS_EREDUCTIONTYPE_HPP

/// \file Teuchos_EReductionType.hpp
/// \brief Declaration of Teuchos::EReductionType enum, and related functions

#include "Teuchos_config.h"
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
MPI_Op getMpiOpForEReductionType (const enum EReductionType reductionType);

} // namespace Details
#endif // HAVE_TEUCHOS_MPI

} // namespace Teuchos

#endif // TEUCHOS_EREDUCTIONTYPE_HPP
