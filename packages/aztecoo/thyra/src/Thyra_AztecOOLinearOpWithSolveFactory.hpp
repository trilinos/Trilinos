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

#include "Thyra_PreconditionedLinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief <tt>PreconditionedLinearOpWithSolveFactoryBase</tt> subclass
 * implemented in terms of <tt>AztecOO</tt>.
 *
 * ToDo: Finish documentation!
 */
class AztecOOLinearOpWithSolveFactory : public PreconditionedLinearOpWithSolveFactoryBase<double> {
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

  /** @name Overridden from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
    ,LinearOpWithSolveBase<double>                             *Op
    ) const;

  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
    ,LinearOpWithSolveBase<double>                             *Op
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<double>                             *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >        *fwdOp
    ) const;

  //@}

  /** \name Overridden from PreconditionedLinearOpWithSolveBase */
  //@{

  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<double> >    &precOp
    ,const EPreconditionerInputType                             precOpType
    ,LinearOpWithSolveBase<double>                              *Op
    ) const;

  /** \brief . */
  void uninitializePreconditionedOp(
    LinearOpWithSolveBase<double>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *fwdOp
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *precOp
    ,EPreconditionerInputType                           *precOpType
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RefCountPtr<Teuchos::ParameterList>     fwdSolveParamlist_;
  bool                                             fwd_cerr_warning_if_unused_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>     adjSolveParamlist_;
  bool                                             adj_cerr_warning_if_unused_;

  // /////////////////////////
  // Private member functions

  void initializeOp_impl(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<double> >    &precOp
    ,const EPreconditionerInputType                             precOpType
    ,const bool                                                 reusePrec
    ,LinearOpWithSolveBase<double>                              *Op
    ) const;

  void uninitializeOp_impl(
    LinearOpWithSolveBase<double>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *fwdOp
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *precOp
    ,EPreconditionerInputType                           *precOpType
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
