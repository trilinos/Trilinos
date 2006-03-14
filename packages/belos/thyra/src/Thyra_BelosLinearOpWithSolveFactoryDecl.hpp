
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

  /** @name Constructors/initializers/accessors */
  //@{


  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{

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
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &precOp
    ,const EPreconditionerInputType                             precOpType
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *fwdOp
    ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *precOp
    ,EPreconditionerInputType                           *precOpType
    ,ESupportSolveUse                                   *supportSolveUse
    ) const;

  //@}


private:

  // /////////////////////////
  // Private data members

  // /////////////////////////
  // Private member functions

};

//@}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
