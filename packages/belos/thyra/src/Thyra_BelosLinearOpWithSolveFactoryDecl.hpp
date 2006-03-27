
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief <tt>LinearOpWithSolveFactoryBase</tt> subclass implemented in terms
 * of <tt>Belos</tt>.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Belos_Thyra_adapters_grp
 */
template<class Scalar>
class BelosLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<Scalar> {
public:

  /** \name Parameter names for Paramter List */
  //@{

  /** \brief . */
  static const std::string SolverType_name;
  /** \brief . */
  static const std::string MaxIters_name;
  /** \brief . */
  static const std::string MaxRestarts_name;
  /** \brief . */
  static const std::string BlockSize_name;
  /** \brief . */
  static const std::string DefaultRelResNorm_name;
  /** \brief . */
  static const std::string GMRES_name;
  /** \brief . */
  static const std::string GMRES_Length_name;
  /** \brief . */
  static const std::string GMRES_Variant_name;
  /** \brief . */
  static const std::string Outputter_name;
  /** \brief . */
  static const std::string Outputter_OutputFrequency_name;
  /** \brief . */
  static const std::string Outputter_OutputMaxResOnly_name;
  /** \brief . */
  static const std::string Preconditioner_name;

  //@}

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct without preconditioner factory. */
  BelosLinearOpWithSolveFactory();

  /** \brief Calls <tt>this->setPreconditionerFactory(precFactory)</tt.  . */
  BelosLinearOpWithSolveFactory(
    const Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  &precFactory
    );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{
  /** \brief Returns true . */
  bool acceptsPreconditionerFactory() const;
  /** \brief . */
  void setPreconditionerFactory(
    const Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  &precFactory
    ,const std::string                                              &precFactoryName
    );
  /** \brief . */
  Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory() const;
  /** \brief . */
  void unsetPreconditionerFactory(
    Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  *precFactory
    ,std::string                                              *precFactoryName
    );
  /** \brief . */
  bool isCompatible( const LinearOpBase<Scalar> &fwdOp ) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > createOp() const;
  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
    ,LinearOpWithSolveBase<Scalar>                             *Op
    ,const ESupportSolveUse                                    supportSolveUse
    ) const;
  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
    ,LinearOpWithSolveBase<Scalar>                             *Op
    ) const;
  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;
  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >             &fwdOp
    ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >      &prec
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;
  /** \brief . */
  void initializeApproxPreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >             &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >            &approxFwdOp
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;
  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar>                               *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *fwdOp
    ,Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >    *prec
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *approxFwdOp
    ,ESupportSolveUse                                           *supportSolveUse
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

  Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  precFactory_;
  std::string                                               precFactoryName_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>              thisValidParamList_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>              paramList_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RefCountPtr<const Teuchos::ParameterList> generateAndGetValidParameters();
  void updateThisValidParamList();

};

//@}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
