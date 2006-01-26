/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief <tt>LinearOpWithSolveFactoryBase</tt> subclass implemented in terms
 * of <tt>AztecOO</tt>.
 *
 * This class creates objects of type <tt>AztecOOLinearOpWithSolve</tt>
 * (through the <tt>LinearOpWithSolveBase</tt> interface) using
 * <tt>AztecOO</tt> objects.
 *
 * The class can support both externally defined preconditioners and built-in
 * aztec preconditioners.  Then built-in aztec preconditioners are used (as
 * specified by the input parameter list), <tt>*this</tt> only supports very
 * limited functionality and does not support adjoint solves.  However, when
 * no preconditioning or externally defined preconditioners are used,
 * <tt>*this</tt> supports a wide range of features which include:
 *
 * <ul>
 * <li>Handling of implicitly scaled and transposed <tt>LinearOpBase</tt>
 * objects through the <tt>ScaledAdjointLinearOpBase</tt> interface.
 * <li>Supports forward and adjoint solves.
 * </ul>
 *
 * <b>Warning:</b> One must be very careful what options are set using the
 * parameter lists passed in using <tt>setFwdAztecSolveParameters()</tt> and
 * <tt>setAdjAztecSolveParameters()</tt> as some of these options will cause
 * great problems and may even result in <tt>exit()</tt> being called to
 * terminate your program!  In the future, a new parameter sublist will be
 * defined that will define a safer way to control the underlying aztec
 * solvers.
 *
 * Click on the above "examples" link at the top to see how this class is
 * used.
 *
 * \ingroup AztecOO_Thyra_adapters_grp
 */
class AztecOOLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

	/** Construct uninitialized but with default option values.
	 *
	 * Note, these defaults where taken from
	 * NOX::EpetraNew::LinearSystemAztecOO::applyJacobianInverse(...) on
	 * 2005/08/15.
	 */
 	AztecOOLinearOpWithSolveFactory(
	 	const int       fwdDefaultMaxIterations       = 400
		,const double   fwdDefaultTol                 = 1e-6
	 	,const int      adjDefaultMaxIterations       = 400
		,const double   adjDefaultTol                 = 1e-6
		);

	/** \brief The default maximum number of iterations for forward solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, fwdDefaultMaxIterations )
  /** \brief The default solution tolerance on the residual for forward solves. */
 	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, fwdDefaultTol )
	/** \brief The default maximum number of iterations for adjoint solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, adjDefaultMaxIterations )
  /** \brief The default solution tolerance on the residual for adjoint solves. */
 	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, adjDefaultTol )

  /** \brief Set the parameters that will be used for the forward aztec sovler.
   *
   * See <tt>AztecOO::SetParameters()</tt> for what options these are!
   */
  void setFwdAztecSolveParameters(
    const Teuchos::RefCountPtr<Teuchos::ParameterList>     &fwdSolveParamlist
    ,bool                                                  fwd_cerr_warning_if_unused = false
    );

  /** \brief Set the parameters that will be used for the adjoint aztec sovler.
   *
   * See <tt>AztecOO::SetParameters()</tt> for what options these are!
   */
  void setAdjAztecSolveParameters(
    const Teuchos::RefCountPtr<Teuchos::ParameterList>     &adjSolveParamlist
    ,bool                                                  adj_cerr_warning_if_unused = false
    );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
    ,LinearOpWithSolveBase<double>                             *Op
    ,const ESupportSolveUse                                    supportSolveUse
    ) const;

  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
    ,LinearOpWithSolveBase<double>                             *Op
    ) const;

  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<double> >    &precOp
    ,const EPreconditionerInputType                             precOpType
    ,LinearOpWithSolveBase<double>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<double>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *fwdOp
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *precOp
    ,EPreconditionerInputType                           *precOpType
    ,ESupportSolveUse                                   *supportSolveUse
    ) const;

  //@}


private:

  // /////////////////////////
  // Private data members

  Teuchos::RefCountPtr<Teuchos::ParameterList>     fwdSolveParamlist_;
  bool                                             fwd_cerr_warning_if_unused_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>     adjSolveParamlist_;
  bool                                             adj_cerr_warning_if_unused_;

  bool useAztecPrec_;

  // /////////////////////////
  // Private member functions

  void initializeOp_impl(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<double> >    &precOp
    ,const EPreconditionerInputType                             precOpType
    ,const bool                                                 reusePrec
    ,LinearOpWithSolveBase<double>                              *Op
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
