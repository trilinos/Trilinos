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
#include "Teuchos_FancyOStream.hpp"

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
 * As a general rule, the client can specify tolerances that are used to
 * determine if a test is successful or not.
 *
 * The client can pass in randomizer objects that create "random" vectors and
 * multi-vectors that are used by the testing class.  This can be very
 * important in some cases and this gives the client full control over what
 * data is used in the tests.
 *
 * This testing class is not designed to test the <tt>LinearOpBase</tt>
 * interface.  For that purpose, use the <tt>LinearOpTester</tt> class in
 * conjunction with this testing class to fully validate a
 * <tt>LinearOpWithSolveBase</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveTester {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<RangeScalar>::magnitudeType              RangeScalarMag;
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<DomainScalar>::magnitudeType             DomainScalarMag;
  /** \brief . */
  typedef typename Teuchos::PromotionTraits<RangeScalar,DomainScalar>::promote    Scalar;
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType                   ScalarMag;

  /** \brief Default constructor which sets default parameter values.
   *
   * Note: It is not recommended that the client pass in values in this
   * constructor since the argument list may change in the near future but
   * instead use the below set functions to change an option after
   * construction.
   */
  LinearOpWithSolveTester(
    const bool                 check_forward_default                    = true
    ,const RangeScalarMag      forward_default_residual_warning_tol     = 1e-6
    ,const RangeScalarMag      forward_default_residual_error_tol       = 1e-5
    ,const DomainScalarMag     forward_default_solution_error_warning_tol= 1e-6
    ,const DomainScalarMag     forward_default_solution_error_error_tol = 1e-5
    ,const bool                check_forward_residual                   = true
    ,const RangeScalarMag      forward_residual_solve_tol               = 1e-5
    ,const RangeScalarMag      forward_residual_slack_warning_tol       = 1e-6
    ,const RangeScalarMag      forward_residual_slack_error_tol         = 1e-5
    ,const bool                check_forward_solution_error             = true
    ,const RangeScalarMag      forward_solution_error_solve_tol         = 1e-5
    ,const RangeScalarMag      forward_solution_error_slack_warning_tol = 1e-6
    ,const RangeScalarMag      forward_solution_error_slack_error_tol   = 1e-5
    ,const bool                check_adjoint_default                    = true
    ,const DomainScalarMag     adjoint_default_residual_warning_tol     = 1e-6
    ,const DomainScalarMag     adjoint_default_residual_error_tol       = 1e-5
    ,const RangeScalarMag      adjoint_default_solution_error_warning_tol= 1e-6
    ,const RangeScalarMag      adjoint_default_solution_error_error_tol = 1e-5
    ,const bool                check_adjoint_residual                   = true
    ,const DomainScalarMag     adjoint_residual_solve_tol               = 1e-5
    ,const DomainScalarMag     adjoint_residual_slack_warning_tol       = 1e-6
    ,const DomainScalarMag     adjoint_residual_slack_error_tol         = 1e-5
    ,const bool                check_adjoint_solution_error             = true
    ,const DomainScalarMag     adjoint_solution_error_solve_tol         = 1e-5
    ,const DomainScalarMag     adjoint_solution_error_slack_warning_tol = 1e-6
    ,const DomainScalarMag     adjoint_solution_error_slack_error_tol   = 1e-5
    ,const int                 num_random_vectors                       = 1
    ,const bool                show_all_tests                           = false
    ,const bool                dump_all                                 = false
    ,const int                 num_rhs                                  = 1
    );
  
  /** \brief Set if a default forward solve will be performed on not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_forward_default  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_default_residual_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_default_residual_error_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, forward_default_solution_error_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, forward_default_solution_error_error_tol  );

  /** \brief Set if a tolerance on the residual of the forward solve should checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_forward_residual  );

  /** \brief Set the relative tolerance that will be requested in the residual
   * for the forward solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_residual_solve_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_residual_slack_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_residual_slack_error_tol  );

  /** \brief Set if a tolerance on the solution error of the forward solve should checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_forward_solution_error  );

  /** \brief Set the relative tolerance that will be requested in the solution_error
   * for the forward solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_solution_error_solve_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_solution_error_slack_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, forward_solution_error_slack_error_tol  );

  /** \brief Set if a default forward solve will be performed on not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint_default  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, adjoint_default_residual_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, adjoint_default_residual_error_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, adjoint_default_solution_error_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, adjoint_default_solution_error_error_tol  );

  /** \brief Set if a tolerance on the residual of the adjoint solve should
   * checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint_residual  );

  /** \brief Set the relative tolerance that will be requested in the residual
   * in the adjoint solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, adjoint_residual_solve_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, adjoint_residual_slack_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( DomainScalarMag, adjoint_residual_slack_error_tol  );

  /** \brief Set if a tolerance on the solution error of the adjoint solve should checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint_solution_error  );

  /** \brief Set the relative tolerance that will be requested in the solution_error
   * for the adjoint solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, adjoint_solution_error_solve_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, adjoint_solution_error_slack_warning_tol  );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( RangeScalarMag, adjoint_solution_error_slack_error_tol  );

  /** \brief Set the number random vectors that is generated during each test.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_random_vectors );

  /** \brief Set if all tests are shown or just summaries.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, show_all_tests );

  /** \brief Set if all of the vectors are dumped or not (only relevant if
   * <tt>show_all_tests()==true</tt>).
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );

  /** \brief Set the number of right-hand-sides in the multivectors
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_rhs );

  /** \brief Turn off all tests so that individual tests can be set.
   *
   * <b>Postconditions:</b><ul>
   * <li>
   * </ul>
   */
  void turn_off_all_tests();

  /** \brief Set all the solve tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>
   * </ul>
   */
  void set_all_solve_tol( const ScalarMag solve_tol );

  /** \brief Set all the warning tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>
   * </ul>
   */
  void set_all_slack_warning_tol( const ScalarMag slack_warning_tol );

  /** \brief Set all the error tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>
   * </ul>
   */
  void set_all_slack_error_tol( const ScalarMag slack_error_tol );
  
  /** \brief Check a <tt>LinearOpWithSolveBase</tt> object.
   *
   * ToDo: Finish documentation!
   */
  bool check(
    const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &op
    ,Teuchos::FancyOStream                                  *out
    ) const;


}; // class LinearOpWithSolveTester

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_DECL_HPP
