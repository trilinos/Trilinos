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
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
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
template<class Scalar>
class LinearOpWithSolveTester
  : virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \name Public types . */
  //@{

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  //@}

  /** \name Constructors/initializers */
  //@{

  /** \brief Default constructor. */
  LinearOpWithSolveTester();
  
  /** \brief Set if a default forward solve will be performed on not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_forward_default );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_default_residual_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_default_residual_error_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_default_solution_error_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_default_solution_error_error_tol );

  /** \brief Set if a tolerance on the residual of the forward solve should checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_forward_residual );
  /** \brief Set the relative tolerance that will be requested in the residual
   * for the forward solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_residual_solve_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_residual_slack_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    forward_residual_slack_error_tol );

  /** \brief Set if a default forward solve will be performed on not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint_default );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_default_residual_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_default_residual_error_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_default_solution_error_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_default_solution_error_error_tol );

  /** \brief Set if a tolerance on the residual of the adjoint solve should
   * checked or not. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint_residual );
  /** \brief Set the relative tolerance that will be requested in the residual
   * in the adjoint solve . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_residual_solve_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_residual_slack_warning_tol );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag,
    adjoint_residual_slack_error_tol );

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
   * <li>???
   * </ul>
   */
  void turn_off_all_tests();

  /** \brief Set all the solve tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>???
   * </ul>
   */
  void set_all_solve_tol( const ScalarMag solve_tol );

  /** \brief Set all the warning tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>???
   * </ul>
   */
  void set_all_slack_warning_tol( const ScalarMag slack_warning_tol );

  /** \brief Set all the error tolerances to the same value.
   *
   * <b>Postconditions:</b><ul>
   * <li>???
   * </ul>
   */
  void set_all_slack_error_tol( const ScalarMag slack_error_tol );

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<ParameterList>& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name LOWS testing */
  //@{
  
  /** \brief Check a <tt>LinearOpWithSolveBase</tt> object.
   *
   * ToDo: Finish documentation!
   */
  bool check(
    const LinearOpWithSolveBase<Scalar> &op,
    Teuchos::FancyOStream *out
    ) const;

  //@}

private:

  static const bool check_forward_default_default_;
  static const bool check_forward_residual_default_;
  static const bool check_adjoint_default_default_;
  static const bool check_adjoint_residual_default_;

  static const ScalarMag warning_tol_default_;
  static const ScalarMag error_tol_default_;
  static const ScalarMag solve_tol_default_;
  static const ScalarMag slack_warning_tol_default_;
  static const ScalarMag slack_error_tol_default_;

  static const int num_random_vectors_default_;
  static const bool show_all_tests_default_;
  static const bool dump_all_default_;
  static const int num_rhs_default_;

  static const std::string AllSolveTol_name_;
  static const std::string AllSlackWarningTol_name_;
  static const std::string AllSlackErrorTol_name_;
  static const std::string ShowAllTests_name_;
  static const std::string DumpAll_name_;

}; // class LinearOpWithSolveTester


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_DECL_HPP
