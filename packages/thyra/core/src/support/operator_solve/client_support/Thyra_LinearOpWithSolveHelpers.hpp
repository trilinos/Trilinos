// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
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
