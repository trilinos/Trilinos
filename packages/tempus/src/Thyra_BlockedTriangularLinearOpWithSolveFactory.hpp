//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_BlockedTriangularLinearOpWithSolveFactory_hpp
#define Thyra_BlockedTriangularLinearOpWithSolveFactory_hpp

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Teuchos_Array.hpp"

namespace Thyra {

/** \brief Implicit subclass that takes a blocked triangular LOWB object and
 * turns it into a LOWSB object.
 *
 * This class takes any upper or lower triangular
 * <tt>PhysicallyBlockedLinearOpBase</tt> object and compatible
 * <tt>LinearOpWithSolveFactoryBase</tt> object(s) and creates a LOWSB version
 * by creating LOWSB objects along the diagonal.
 *
 *
 * For example, consider the lower block triangular linear operator:

 \verbatim

       [ M(0,0)                   ]
   M = [ M(1,0)   M(1,1)          ]
       [ M(2,0)   M(2,1)   M(2,2) ]

 \endverbatim

 * This class object will then create a new LOWSB object (of type
 * <tt>DefaultBlockedTriangularLinearOpWithSolve</tt>) that looks like:

 \verbatim

          [ invM(0,0)                       ]
   invM = [ M(1,0)     invM(1,1)            ]
          [ M(2,0)     M(2,1)     invM(2,2) ]

 \endverbatim

 * where <tt>invM(k,k)</tt> are LOWSB objects created from the LOB objects
 * <tt>M(k,k)</tt> given a LOWSFB object.
 *
 * This class is not very compliciated, see the function
 * <tt>initializeOp()</tt> see what this class actually does!
 *
 * Note, this is basically the same as
 * DefaultBlockedTriangularLinearOpWithSolveFactory except this version allows
 * you to set a different LOWSF for each diagonal block.
 */
template <class Scalar>
class BlockedTriangularLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar> {
 public:
  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  /** \brief Create given an array of non-const LOWSFB objects.
   *
   * \param lowsf [in,persisting] The LOWSFB objects that will be used to
   * create the LOWSB objects for the diagonal blocks.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  BlockedTriangularLinearOpWithSolveFactory(
      const Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > > &lowsf);

  /** \brief Create given an array of const LOWSFB objects.
   *
   * \param lowsf [in,persisting] The LOWSFB objects that will be used to
   * create the LOWSB objects for the diagonal blocks.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  BlockedTriangularLinearOpWithSolveFactory(
      const Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > > &lowsf);

  Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > > getUnderlyingLOWSF();

  Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > > getUnderlyingLOWSF()
      const;

  //@}

  /** \name Overridden from Teuchos::Describable. */
  //@{

  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions)
   */
  //@{

  void setParameterList(RCP<ParameterList> const &paramList);
  RCP<ParameterList> getNonconstParameterList();
  RCP<ParameterList> unsetParameterList();
  RCP<const ParameterList> getParameterList() const;
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from LinearOpWithSolveFactoyBase */
  //@{

  /** \brief returns false. */
  virtual bool acceptsPreconditionerFactory() const;

  /** \brief Throws exception. */
  virtual void setPreconditionerFactory(
      const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
      const std::string &precFactoryName);

  /** \brief Returns null . */
  virtual RCP<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory()
      const;

  /** \brief Throws exception. */
  virtual void unsetPreconditionerFactory(
      RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
      std::string *precFactoryName);

  virtual bool isCompatible(const LinearOpSourceBase<Scalar> &fwdOpSrc) const;

  virtual RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  virtual void initializeOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  virtual void initializeAndReuseOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op) const;

  virtual void uninitializeOp(
      LinearOpWithSolveBase<Scalar> *Op,
      RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
      RCP<const PreconditionerBase<Scalar> > *prec,
      RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
      ESupportSolveUse *supportSolveUse) const;

  virtual bool supportsPreconditionerInputType(
      const EPreconditionerInputType precOpType) const;

  virtual void initializePreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      const RCP<const PreconditionerBase<Scalar> > &prec,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  virtual void initializeApproxPreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  //@}

 protected:
  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  void informUpdatedVerbosityState() const;

  //@}

 private:
  typedef Teuchos::ConstNonconstObjectContainer<
      LinearOpWithSolveFactoryBase<Scalar> >
      LOWSF_t;

  Array<LOWSF_t> lowsf_;

  // Not defined and not to be called
  BlockedTriangularLinearOpWithSolveFactory();
};

/** \brief Nonmember constructor.
 *
 * \relates BlockedTriangularLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<BlockedTriangularLinearOpWithSolveFactory<Scalar> >
blockedTriangularLinearOpWithSolveFactory(
    const Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > > &lowsf)
{
  return Teuchos::rcp(
      new BlockedTriangularLinearOpWithSolveFactory<Scalar>(lowsf));
}

/** \brief Nonmember constructor.
 *
 * \relates BlockedTriangularLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<BlockedTriangularLinearOpWithSolveFactory<Scalar> >
blockedTriangularLinearOpWithSolveFactory(
    const Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > > &lowsf)
{
  return Teuchos::rcp(
      new BlockedTriangularLinearOpWithSolveFactory<Scalar>(lowsf));
}

// Overridden from Constructors/Initializers/Accessors

template <class Scalar>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    BlockedTriangularLinearOpWithSolveFactory(
        const Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > > &lowsf)
  : lowsf_(lowsf.size())
{
  for (Ordinal i = 0; i < lowsf.size(); ++i) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf[i]));
#endif
    lowsf_[i].initialize(lowsf[i]);
  }
}

template <class Scalar>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    BlockedTriangularLinearOpWithSolveFactory(
        const Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > > &lowsf)
  : lowsf_(lowsf.size())
{
  for (Ordinal i = 0; i < lowsf.size(); ++i) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf[i]));
#endif
    lowsf_[i].initialize(lowsf[i]);
  }
}

template <class Scalar>
Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > >
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  Array<RCP<LinearOpWithSolveFactoryBase<Scalar> > > lowsf(lowsf_.size());
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    lowsf[i] = lowsf_[i].getNonconstObj();
  }
  return lowsf;
}

template <class Scalar>
Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > >
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  Array<RCP<const LinearOpWithSolveFactoryBase<Scalar> > > lowsf(lowsf_.size());
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    lowsf[i] = lowsf_[i].getConstObj();
  }
  return lowsf;
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string BlockedTriangularLinearOpWithSolveFactory<Scalar>::description()
    const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description() << "{";
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    oss << "lowsf=";
    if (!is_null(lowsf_[i].getConstObj()))
      oss << lowsf_[i].getConstObj()->description();
    else
      oss << "NULL";
  }
  oss << "}";
  return oss.str();
}

// Overridden from ParameterListAcceptor

// Note, we should probably do something smarter with the parameter lists

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::setParameterList(
    RCP<ParameterList> const &paramList)
{
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    lowsf_[i].getNonconstObj()->setParameterList(paramList);
  }
}

template <class Scalar>
RCP<ParameterList>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_[0].getNonconstObj()->getNonconstParameterList();
}

template <class Scalar>
RCP<ParameterList>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  RCP<ParameterList> pl;
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    pl = lowsf_[i].getNonconstObj()->unsetParameterList();
  }
  return pl;
}

template <class Scalar>
RCP<const ParameterList>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_[0].getConstObj()->getParameterList();
}

template <class Scalar>
RCP<const ParameterList>
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_[0].getConstObj()->getValidParameters();
}

// Overridden from LinearOpWithSolveFactoyBase

template <class Scalar>
bool BlockedTriangularLinearOpWithSolveFactory<
    Scalar>::acceptsPreconditionerFactory() const
{
  return false;
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    setPreconditionerFactory(
        const RCP<PreconditionerFactoryBase<Scalar> > & /* precFactory */,
        const std::string & /* precFactoryName */
    )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't support a preconditioner factory!");
}

template <class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
BlockedTriangularLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory()
    const
{
  return Teuchos::null;
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<
    Scalar>::unsetPreconditionerFactory(RCP<PreconditionerFactoryBase<Scalar> >
                                            * /* precFactory */,
                                        std::string * /* precFactoryName */
)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't support a preconditioner factory!");
}

template <class Scalar>
bool BlockedTriangularLinearOpWithSolveFactory<Scalar>::isCompatible(
    const LinearOpSourceBase<Scalar> & /* fwdOpSrc */
) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  TEUCHOS_UNREACHABLE_RETURN(false);
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
BlockedTriangularLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return defaultBlockedTriangularLinearOpWithSolve<Scalar>();
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse /* supportSolveUse */
) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    lowsf_[i].getConstObj()->setOStream(this->getOStream());
    lowsf_[i].getConstObj()->setVerbLevel(this->getVerbLevel());
  }

  // Get the block interface to get at the blocks
  typedef PhysicallyBlockedLinearOpBase<Scalar> PBLOB;
  const RCP<const PBLOB> blo =
      rcp_dynamic_cast<const PBLOB>(fwdOpSrc->getOp().assert_not_null());

  // Dynamic cast to get the DefaultBlockedTriangularLinearOpWithSolveBase
  // interface that we will fill.

  typedef DefaultBlockedTriangularLinearOpWithSolve<Scalar> DBTLOWS;
  DBTLOWS &btlows = dyn_cast<DBTLOWS>(*Op);

  // Determine if this is the first time through or if we have already
  // initialized before.  This will be needed to allow efficient reuse of the
  // LOWSB objects for the diagonal blocks.
  const bool firstTime = is_null(btlows.range());

  // If this is the first time through, we need to fill and create the block
  // structure
  if (firstTime)
    btlows.beginBlockFill(blo->productRange(), blo->productDomain());

  const int N = blo->productRange()->numBlocks();
  for (int k = 0; k < N; ++k) {
    const RCP<const LinearOpBase<Scalar> > fwdOp_k =
        blo->getBlock(k, k).assert_not_null();
    if (firstTime) {
      // This is the first time through so reate and initialize a new LOWSB
      // object for each block
      btlows.setNonconstLOWSBlock(
          k, k, linearOpWithSolve<Scalar>(*lowsf_[k].getConstObj(), fwdOp_k));
    }
    else {
      // This is not the first time through so we need to just reinitiallize
      // the object that is already created.  This allows us to efficiently
      // reuse precreated structure and storage.
      RCP<LinearOpWithSolveBase<Scalar> > invOp_k =
          btlows.getNonconstLOWSBlock(k, k).assert_not_null();
      Thyra::initializeOp<Scalar>(*lowsf_[k].getConstObj(), fwdOp_k,
                                  invOp_k.ptr());
    }
  }

  // If this is the first time through, then we need to finalize the block
  // structure.
  if (firstTime) btlows.endBlockFill();

  // After the block structure has been setup, set the off-diagonal blocks.
  // Note that this also sets the diagonal blocks but these are ignored since
  // the LOWSB blocks created above override these.
  btlows.setBlocks(blo);

  // Set the verbosity settings
  btlows.setOStream(this->getOStream());
  btlows.setVerbLevel(this->getVerbLevel());
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
    LinearOpWithSolveBase<Scalar> * /* Op */
) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse * /* supportSolveUse */
) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  typedef DefaultBlockedTriangularLinearOpWithSolve<Scalar> DBTLOWS;
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
  DBTLOWS &btlowsOp = dyn_cast<DBTLOWS>(*Op);
  if (fwdOpSrc) {
    const RCP<const LinearOpBase<Scalar> > fwdOp = btlowsOp.getBlocks();
    if (!is_null(fwdOp))
      *fwdOpSrc = defaultLinearOpSource<Scalar>(fwdOp);
    else
      *fwdOpSrc = Teuchos::null;
  }
  if (prec) *prec = Teuchos::null;
  if (approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null;
}

template <class Scalar>
bool BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    supportsPreconditionerInputType(
        const EPreconditionerInputType /* precOpType */
    ) const
{
  // We don't support any external preconditioners!
  return false;
  // 20071006: rabartl: Note: We could support external preconditioners but it
  // will take some work.  We would have to extract out the individual
  // preconditioners from each block.  This would be pretty easy to do but I
  // am not going to do this until we have to.
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    initializePreconditionedOp(
        const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
        const RCP<const PreconditionerBase<Scalar> > & /* prec */,
        LinearOpWithSolveBase<Scalar> * /* Op */,
        const ESupportSolveUse /* supportSolveUse */
    ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't support an external preconditioner!");
}

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<Scalar>::
    initializeApproxPreconditionedOp(
        const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
        const RCP<const LinearOpSourceBase<Scalar> > & /* approxFwdOpSrc */,
        LinearOpWithSolveBase<Scalar> * /* Op */,
        const ESupportSolveUse /* supportSolveUse */
    ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error, we don't support an external preconditioner!");
}

// protected

template <class Scalar>
void BlockedTriangularLinearOpWithSolveFactory<
    Scalar>::informUpdatedVerbosityState() const
{
  for (Ordinal i = 0; i < lowsf_.size(); ++i) {
    lowsf_[i].getConstObj()->setVerbLevel(this->getVerbLevel());
    lowsf_[i].getConstObj()->setOStream(this->getOStream());
  }
}

}  // namespace Thyra

#endif  // Thyra_BlockedTriangularLinearOpWithSolveFactory_hpp
