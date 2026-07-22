// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ZERO_LINEAR_OP_BASE_HPP
#define THYRA_ZERO_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

/** \brief Interface class for zero linear operators.
 *
 * This interface represents a zero linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = 0
 \endverbatim
 *
 * In other words, subclasses define <tt>apply()</tt> as:
 *
 \verbatim

 y = alpha*M*x + beta*y
   = beta * y
 \endverbatim
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class ZeroLinearOpBase : virtual public LinearOpBase<Scalar>
{};

} // namespace Thyra

#endif	// THYRA_ZERO_LINEAR_OP_BASE_HPP
