#ifndef __Teko_StratimikosFactory_hpp__
#define __Teko_StratimikosFactory_hpp__

#include <vector>

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

#include "Teko_RequestHandler.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

#include "Epetra_Operator.h"

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

  StratimikosFactory(const Teuchos::RCP<Teko::RequestHandler> & rh);
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS(
    Thyra::EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor );

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const Thyra::LinearOpSourceBase<double> &fwdOp ) const;
  /** \brief . */
  bool applySupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(Thyra::EConj conj) const;
  /** \brief . */
  Teuchos::RCP<Thyra::PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(
    const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOp,
    Thyra::PreconditionerBase<double> *prec,
    const Thyra::ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void uninitializePrec(
    Thyra::PreconditionerBase<double> *prec
    ,Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > *fwdOp
    ,Thyra::ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const& paramList);
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

  /** Access to the application communication request handling mechnism
    */
  void setRequestHandler(const Teuchos::RCP<Teko::RequestHandler> & rh)
  { reqHandler_ = rh; }

  /** Access to the application communication request handling mechnism
    */
  Teuchos::RCP<Teko::RequestHandler> getRequestHandler() const 
  { return reqHandler_; }

private:

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
                                                     const Teuchos::RCP<const Epetra_Operator> & Jac,
                                                     const Teuchos::RCP<Epetra_Operator> & wrapInput,
                                                     std::ostream & out) const;

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  mutable Teuchos::RCP<Teko::InverseLibrary> invLib_;
  mutable Teuchos::RCP<Teko::InverseFactory> invFactory_;
  Teuchos::RCP<Teko::RequestHandler> reqHandler_;
  mutable std::vector<int> decomp_;
};

void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder & builder,
                               const std::string & stratName="Teko");

void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder & builder,
                               const Teuchos::RCP<Teko::RequestHandler> & rh,
                               const std::string & stratName="Teko");

} // namespace Teko

#endif 
