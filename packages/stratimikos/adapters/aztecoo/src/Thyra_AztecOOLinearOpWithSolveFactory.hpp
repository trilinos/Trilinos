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
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
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

  /** \brief Construct uninitialized. */
   AztecOOLinearOpWithSolveFactory(
     Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList = Teuchos::null
     );
  
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS( EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{
  /** \brief Returns true . */
  bool acceptsPreconditionerFactory() const;
  /** \brief . */
  void setPreconditionerFactory(
    const Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >  &precFactory,
    const std::string  &precFactoryName
    );
  /** \brief . */
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> > getPreconditionerFactory() const;
  /** \brief . */
  void unsetPreconditionerFactory(
    Teuchos::RefCountPtr<PreconditionerFactoryBase<double> > *precFactory,
    std::string *precFactoryName
    );
  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > createOp() const;
  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op
    ) const;
  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<double> *Op,
    Teuchos::RefCountPtr<const LinearOpSourceBase<double> > *fwdOpSrc,
    Teuchos::RefCountPtr<const PreconditionerBase<double> > *prec,
    Teuchos::RefCountPtr<const LinearOpSourceBase<double> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;
  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;
  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RefCountPtr<const PreconditionerBase<double> > &prec,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void initializeApproxPreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &approxFwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> > precFactory_;
  std::string precFactoryName_;
  Teuchos::RefCountPtr<Teuchos::ParameterList> thisValidParamList_;
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList_;

  int defaultFwdMaxIterations_;
  double defaultFwdTolerance_;
  int defaultAdjMaxIterations_;
  double defaultAdjTolerance_;
  bool outputEveryRhs_;

  bool useAztecPrec_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RefCountPtr<const Teuchos::ParameterList> generateAndGetValidParameters();
  void updateThisValidParamList();

  void initializeOp_impl(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RefCountPtr<const PreconditionerBase<double> > &prec,
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &approxFwdOpSrc,
    const bool reusePrec,
    LinearOpWithSolveBase<double> *Op
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
