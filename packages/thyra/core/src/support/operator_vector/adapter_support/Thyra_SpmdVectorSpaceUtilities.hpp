// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SPMD_VECTOR_SPACE_UTILITIES_HPP
#define THYRA_SPMD_VECTOR_SPACE_UTILITIES_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Comm.hpp"

namespace Thyra {

namespace SpmdVectorSpaceUtilities {

/** \brief . */
Ordinal computeMapCode( const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim );

/** \brief . */
Ordinal computeLocalOffset( const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim );

/** \brief . */
Ordinal computeGlobalDim( const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim );

} // namespace SpmdVectorSpaceUtiltities

} // namespace Thyra

#endif // THYRA_SPMD_VECTOR_SPACE_UTILITIES_HPP
