// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
 * \ingroup AztecOO_Thyra_adapters_grp
 */
class AztecOOLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct uninitialized. */
   AztecOOLinearOpWithSolveFactory(
     Teuchos::RCP<Teuchos::ParameterList> const& paramList = Teuchos::null
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
    const Teuchos::RCP<PreconditionerFactoryBase<double> >  &precFactory,
    const std::string  &precFactoryName
    );
  /** \brief . */
  Teuchos::RCP<PreconditionerFactoryBase<double> > getPreconditionerFactory() const;
  /** \brief . */
  void unsetPreconditionerFactory(
    Teuchos::RCP<PreconditionerFactoryBase<double> > *precFactory,
    std::string *precFactoryName
    );
  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;
  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<double> > createOp() const;
  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op
    ) const;
  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<double> *Op,
    Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
    Teuchos::RCP<const PreconditionerBase<double> > *prec,
    Teuchos::RCP<const LinearOpSourceBase<double> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;
  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;
  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RCP<const PreconditionerBase<double> > &prec,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void initializeApproxPreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RCP<const LinearOpSourceBase<double> > &approxFwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;
  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RCP<PreconditionerFactoryBase<double> > precFactory_;
  std::string precFactoryName_;
  Teuchos::RCP<Teuchos::ParameterList> thisValidParamList_;
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  int defaultFwdMaxIterations_;
  double defaultFwdTolerance_;
  int defaultAdjMaxIterations_;
  double defaultAdjTolerance_;
  bool outputEveryRhs_;

  bool useAztecPrec_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RCP<const Teuchos::ParameterList> generateAndGetValidParameters();
  void updateThisValidParamList();

  void initializeOp_impl(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const Teuchos::RCP<const PreconditionerBase<double> > &prec,
    const Teuchos::RCP<const LinearOpSourceBase<double> > &approxFwdOpSrc,
    const bool reusePrec,
    LinearOpWithSolveBase<double> *Op
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
