// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_StratimikosFactory_hpp__
#define __Teko_StratimikosFactory_hpp__

#include <vector>

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

#include "Teko_RequestHandler.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

#ifdef TEKO_HAVE_EPETRA
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Epetra_Operator.h"
#endif

namespace Teko {

/** \brief Concrete preconditioner factory subclass based on ML.
 *
 * ToDo: Finish documentation!
 */
class StratimikosFactory : public Thyra::PreconditionerFactoryBase<double> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  StratimikosFactory();

  StratimikosFactory(const Teuchos::RCP<Teko::RequestHandler> &rh);
  StratimikosFactory(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> &builder,
                     const Teuchos::RCP<Teko::RequestHandler> &rh);

#ifdef TEKO_HAVE_EPETRA
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS(Thyra::EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor);
#endif

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOp) const;
  /** \brief . */
  bool applySupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  Teuchos::RCP<Thyra::PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOp,
                      Thyra::PreconditionerBase<double> *prec,
                      const Thyra::ESupportSolveUse supportSolveUse) const;
  /** \brief . */
  void uninitializePrec(Thyra::PreconditionerBase<double> *prec,
                        Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > *fwdOp,
                        Thyra::ESupportSolveUse *supportSolveUse) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  // ToDo: Add an override of describe(...) to give more detail!

  //@}

  /** Setup an thyra preconditioner (most likely blocked)
   */
  void initializePrec_Thyra(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOp,
                            Thyra::PreconditionerBase<double> *prec,
                            const Thyra::ESupportSolveUse supportSolveUse) const;

#ifdef TEKO_HAVE_EPETRA
  /** Setup an epetra preconditioner.
   */
  void initializePrec_Epetra(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOp,
                             Thyra::PreconditionerBase<double> *prec,
                             const Thyra::ESupportSolveUse supportSolveUse) const;
#endif

  /** Access to the application communication request handling mechnism
   */
  void setRequestHandler(const Teuchos::RCP<Teko::RequestHandler> &rh) { reqHandler_ = rh; }

  /** Access to the application communication request handling mechnism
   */
  Teuchos::RCP<Teko::RequestHandler> getRequestHandler() const { return reqHandler_; }

  //! Get the decomposition vector in use by this factory
  const std::vector<int> &getDecomposition() const { return decomp_; }

 private:
#ifdef TEKO_HAVE_EPETRA
  /** Build the segragated jacobian operator according to
   * the input parameter list.
   *
   * \param[in] Jac Epetra_CrsMatrix (assumed) to be decomposed.
   * \param[in] wrapInput RCP to an Epetra operator that was either previously
   *                      wrapped, or it is <code>Teuchos::null</code>
   * \param[in] out Stream for use when reporting testing results
   *
   * \returns Blocked operator that was requested, if <code>wrapInput</code> is not
   *          null, then the returned pointer will match (i.e. wrapInput will be over-
   *          written).
   */
  Teuchos::RCP<Epetra_Operator> buildWrappedEpetraOperator(
      const Teuchos::RCP<const Epetra_Operator> &Jac,
      const Teuchos::RCP<Epetra_Operator> &wrapInput, std::ostream &out) const;

  /** Build strided vectors using the operator range map and
   * the decomposition vectors.
   *
   * \param[in] Jac Epetra_CrsMatrix (assumed) to be decomposed.
   * \param[in] decomp Decomposition vector.
   * \param[in,out] vars Vector of vectors of global ids specifying
   *                     how the operator is to be blocked.
   */
  void buildStridedVectors(const Epetra_Operator &Jac, const std::vector<int> &decomp,
                           std::vector<std::vector<int> > &vars) const;
#endif  // TEKO_HAVE_EPETRA

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  mutable Teuchos::RCP<Teko::InverseLibrary> invLib_;
  mutable Teuchos::RCP<Teko::InverseFactory> invFactory_;
  Teuchos::RCP<Teko::RequestHandler> reqHandler_;
  mutable std::vector<int> decomp_;
  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>
      builder_;  // builder to use for default solvers
};

/** Inject Teko into a solver builder with a specified name. Note that the builder
 * is used to define the default solvers in Teko. Therefore any dynamically specified
 * solver/preconditioner, will be inherited by Teko. This relationship does not persist,
 * so no solvers added to the builder after the call to addTekoToStratimikosBuilder will
 * be inherited. Finally, to avoid circular references Teko will not include itself in
 * that list of solvers/precondtioners.
 */
void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder &builder,
                                 const std::string &stratName = "Teko");

/** Inject Teko into a solver builder with a specified name. Note that the builder
 * is used to define the default solvers in Teko. Therefore any dynamically specified
 * solver/preconditioner, will be inherited by Teko. This relationship does not persist,
 * so no solvers added to the builder after the call to addTekoToStratimikosBuilder will
 * be inherited. Finally, to avoid circular references Teko will not include itself in
 * that list of solvers/precondtioners.
 */
void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder &builder,
                                 const Teuchos::RCP<Teko::RequestHandler> &rh,
                                 const std::string &stratName = "Teko");

}  // namespace Teko

#endif
