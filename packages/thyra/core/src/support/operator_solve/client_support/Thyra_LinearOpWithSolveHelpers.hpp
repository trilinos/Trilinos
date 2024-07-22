// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_SUBCLASS_HELPERS_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_SUBCLASS_HELPERS_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_toString.hpp"


namespace Thyra {


/** \brief Assert that a LOWSB object supports a particular solve type.
 *
 * This function will throw an excetion with a very good error message if the
 * requested solve type is not supported.
 *
 * \relates LinearOpWithSolveBase
 */
template<class Scalar>
void assertSolveSupports(
  const LinearOpWithSolveBase<Scalar> &lows,
  const EOpTransp M_trans,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria = Teuchos::null
  );

/** \brief Assert that a LOWSB object supports a particular solve type.
 *
 * This function will throw an excetion with a very good error message if the
 * requested solve type is not supported.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
void assertSupportsSolveMeasureType(
  const LinearOpWithSolveBase<Scalar> &lows,
  const EOpTransp M_trans,
  const SolveMeasureType &solveMeasureType
  );
// 2010/08/22: rabartl: ToDo: Deprecate this bug 4915 is finished!

} // namespace Thyra


//
// Implementations
//


template<class Scalar>
void Thyra::assertSolveSupports(
  const LinearOpWithSolveBase<Scalar> &lows,
  const EOpTransp M_trans,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  )
{
  using Teuchos::toString;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !lows.solveSupports(M_trans, solveCriteria),
    std::logic_error,
    "Error, the LinearOpWithSolve object \"" << lows.description() << "\"\n"
    "for M_trans = " << toString(M_trans) << " does not support the solve"
    " criteria = "
    << ( nonnull(solveCriteria) ? toString(*solveCriteria) : std::string("null") )
    << "!"
    );
}
// 2010/08/22: rabartl: Bug 4915 ToDo: Move the above into the NIV function
// solve(...).

template<class Scalar>
void Thyra::assertSupportsSolveMeasureType(
  const LinearOpWithSolveBase<Scalar> &lows,
  const EOpTransp M_trans,
  const SolveMeasureType &solveMeasureType
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !solveSupportsSolveMeasureType(lows,M_trans,solveMeasureType),
    std::logic_error,
    "Error, the LinearOpWithSolve object \"" << lows.description() << "\"\n"
    "for M_trans = " << toString(M_trans) << " does not support the solve"
    " measure = "
    << toString(solveMeasureType.numerator)
    << "/"
    << toString(solveMeasureType.denominator)
    << "!"
    );
}

#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_SUBCLASS_HELPERS_HPP
