// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_TESTER_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_TESTER_DECL_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Thyra {

/** \brief Testing class for <tt>LinearOpWithSolveBase</tt>.
 *
 * This testing class can be used in several different roles:
 *
 * <ul>
 *
 * <li> <b>Comprehensive unit testing class:</b> In this mode, all of the
 * capabilities of the solve functions are tested and the post conditions are
 * examined according to a set of tolerances.  This is the default mode.
 *
 * <li> <b>Check of accuracy for specific types of linear solves:</b> In this
 * mode, the client can check only one or more types of linear solves for the
 * sake of efficiency.
 *
 * </ul>
 *
 * This class can check single linear solves or multi-RHS (i.e. multi-vector)
 * linear solves.
 *
 * As a general rule, the client can specify tolerances to determine if a test
 * should succeed or fail.
 *
 * The client can pass in randomizer objects that create valid random vectors
 * and multi-vectors that are used by the testing class.  This can be very
 * important in some cases and this gives the client full control over what
 * data is used in the tests.
 *
 * This testing class is not designed to test the <tt>LinearOpBase</tt>
 * interface.  For that purpose, use the <tt>LinearOpTester</tt> class in
 * conjunction with this testing class to fully validate a
 * <tt>LinearOpWithSolveBase</tt> objects.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class LinearOpTester {
public:

  /** \brief Local typedef for scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

}; // class LinearOpTester

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_DECL_HPP
