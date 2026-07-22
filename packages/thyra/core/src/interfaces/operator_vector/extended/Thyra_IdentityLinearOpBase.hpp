// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_IDENTITY_LINEAR_OP_BASE_HPP
#define THYRA_IDENTITY_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

/** \brief Interface class for identity linear operators.
 *
 * This interface represents a identity linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = I
 \endverbatim
 *
 * In other words, subclasses define <tt>apply()</tt> as:
 *
 \verbatim

 y = alpha*M*x + beta*y
   = alpha * x + beta * y
 \endverbatim
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class IdentityLinearOpBase : virtual public LinearOpBase<Scalar>
{};

} // namespace Thyra

#endif	// THYRA_IDENTITY_LINEAR_OP_BASE_HPP
