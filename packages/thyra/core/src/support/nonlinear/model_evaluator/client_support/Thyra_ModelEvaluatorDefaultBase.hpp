// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MODEL_EVALUATOR_DEFAULT_BASE_HPP
#define THYRA_MODEL_EVALUATOR_DEFAULT_BASE_HPP

#include "Thyra_VectorBase.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"


#ifdef HAVE_THYRA_ME_POLYNOMIAL


// Define the polynomial traits class specializtaion
// Teuchos::PolynomialTraits<VectorBase > before there is any chance of an
// instantiation requiring the definition of this class.  The Intel C++ 9.1
// compiler is having trouble compiling Thyra_EpetraModelEvaluator.cpp because
// Thyra_ModelEvaluatorBase_decl.hpp is only including
// "Teuchos_Polynomial.hpp" which does not contain the needed specialization.
// By including it here, all client code is likely to compile just fine.  I am
// going through all of the in order to avoid having to put
// Thyra_PolynomialVectorTraits.hpp in the interfaces directory.  The problem
// with putting Thyra_PolynomialVectorTraits.hpp the interfaces directory is
// that it contains calls to code in the support directory and creates an
// unacceptable circular dependency that will cause problems once we move to
// subpackages in the CMake build system.
#include "Thyra_PolynomialVectorTraits.hpp"


#endif // HAVE_THYRA_ME_POLYNOMIAL


namespace Thyra {


//
// Undocumentated (from the user's perspective) types that are used to
// implement ModelEvaluatorDefaultBase.  These types are defined outside of
// the templated class since they do nt depend on the scalar type.
//


namespace ModelEvaluatorDefaultBaseTypes {


// Type used to determine if the ModelEvaluatorDefaultBase implementation will
// provide a LinearOpBase wrapper for derivative object given in
// MultiVectorBase form.
class DefaultDerivLinearOpSupport {
public:
  DefaultDerivLinearOpSupport()
    :provideDefaultLinearOp_(false),
     mvImplOrientation_(ModelEvaluatorBase::DERIV_MV_BY_COL)
    {}
  DefaultDerivLinearOpSupport(
    const ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvImplOrientation_in
    )
    :provideDefaultLinearOp_(true),
     mvImplOrientation_(mvImplOrientation_in)
    {}
  bool provideDefaultLinearOp() const
    { return provideDefaultLinearOp_; }
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvImplOrientation() const
    { return mvImplOrientation_; }
private:
  bool provideDefaultLinearOp_;
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvImplOrientation_;
};


// Type used to determine if the ModelEvaluatorDefaultBase implementation will
// provide an explicit transpose copy of a derivative object given in
// MultiVectorBase form.
class DefaultDerivMvAdjointSupport {
public:
  DefaultDerivMvAdjointSupport()
    :provideDefaultAdjoint_(false),
     mvAdjointCopyOrientation_(ModelEvaluatorBase::DERIV_MV_BY_COL)
    {}
  DefaultDerivMvAdjointSupport(
    const ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvAdjointCopyOrientation_in
    )
    :provideDefaultAdjoint_(true),
     mvAdjointCopyOrientation_(mvAdjointCopyOrientation_in)
    {}
  bool provideDefaultAdjoint() const
    { return provideDefaultAdjoint_; }
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvAdjointCopyOrientation() const
    { return mvAdjointCopyOrientation_; }
private:
  bool provideDefaultAdjoint_;
  ModelEvaluatorBase::EDerivativeMultiVectorOrientation mvAdjointCopyOrientation_;
};


// Type used to remember a pair of transposed multi-vectors to implement a
// adjoint copy.
template<class Scalar>
struct MultiVectorAdjointPair {
  MultiVectorAdjointPair()
    {}
  MultiVectorAdjointPair(
    const RCP<MultiVectorBase<Scalar> > &in_mvOuter,
    const RCP<const MultiVectorBase<Scalar> > &in_mvImplAdjoint
    )
    : mvOuter(in_mvOuter),
      mvImplAdjoint(in_mvImplAdjoint)
    {}
  RCP<MultiVectorBase<Scalar> > mvOuter;
  RCP<const MultiVectorBase<Scalar> > mvImplAdjoint;
private:
};




} // namespace ModelEvaluatorDefaultBaseTypes


/** \brief Default base class for concrete model evaluators.
 *
 * The primary purposes of this base class are to:
 *
 * <ul>
 *
 * <li> Provide default implementations for other forms of derivatives in
 * <tt>ModelEvaluatorBase::Derivative</tt>.  For example, if a multi-vector
 * form of a derivative (i.e. my column <tt>DhDy</tt> or by row
 * <tt>DhDy^T</tt>) is provided then the other form will be provided assuming
 * the range space has in-core vectors.  Also, if any multi-vector form of a
 * general derivative is provided, a <tt>LinearOpBase</tt> version is
 * automatically supported.
 *
 * <li> Provide a default implementation for computing the <tt>LOWSB</tt> from
 * <tt>W</tt> given the <tt>LOB</tt>-only form <tt>W_op</tt> given a
 * <tt>LOWSFB</tt> object <tt>W_factory</tt> supplied by the subclass.  If the
 * subclass wants to take this over, then it should override
 * <tt>create_W()</tt>.
 *
 * <li> Assert (in debug mode) that the underlying model has been set up
 * correctly.
 *
 * </ul>
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class ModelEvaluatorDefaultBase : virtual public ModelEvaluator<Scalar>
{
public:

  /** \name Overridden from ModelEvaluator */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DfDp_op(int l) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDx_dot_op(int j) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_DgDp_op(int j, int l) const;
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

protected:

  /** \name Setup functions called by subclasses */
  //@{

  /** \brief Function called by subclasses to fully initialize this object on
   * any important change.
   *
   * Note: This class will automatically call this function the first time to
   * set things up and does not need to be called by the client the first
   * time.  However, if the state of the object changes, then this function
   * should be called to reset the state of this object's implemention!
   */
  void initializeDefaultBase();

  //@}

private:

  /** \name Private functions with default implementaton to be overridden by subclasses */
  //@{

  /** \brief . */
  virtual RCP<LinearOpBase<Scalar> > create_DfDp_op_impl(int l) const;

  /** \brief . */
  virtual RCP<LinearOpBase<Scalar> > create_DgDx_dot_op_impl(int j) const;

  /** \brief . */
  virtual RCP<LinearOpBase<Scalar> > create_DgDx_op_impl(int j) const;

  /** \brief . */
  virtual RCP<LinearOpBase<Scalar> > create_DgDp_op_impl(int j, int l) const;
  
  //@}

  /** \name Private pure virtual functions that must be overridden by subclasses */
  //@{

  /** \brief . */
  virtual ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const = 0;

  /** \brief . */
  virtual void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const = 0;

  //@}

protected:

  /** \brief . */
  ModelEvaluatorDefaultBase();

private:
  
  // //////////////////////////////
  // Private tpyes

  typedef ModelEvaluatorDefaultBaseTypes::DefaultDerivLinearOpSupport
  DefaultDerivLinearOpSupport;

  typedef ModelEvaluatorDefaultBaseTypes::DefaultDerivMvAdjointSupport
  DefaultDerivMvAdjointSupport;

  typedef ModelEvaluatorDefaultBaseTypes::MultiVectorAdjointPair<Scalar>
  MultiVectorAdjointPair;
  
  // //////////////////////////////
  // Private data members

  bool isInitialized_;

  Array<DefaultDerivLinearOpSupport> DfDp_default_op_support_;

  Array<DefaultDerivLinearOpSupport> DgDx_dot_default_op_support_;

  Array<DefaultDerivLinearOpSupport> DgDx_default_op_support_;

  Array<Array<DefaultDerivLinearOpSupport> > DgDp_default_op_support_;
  Array<Array<DefaultDerivMvAdjointSupport> > DgDp_default_mv_support_;

  bool default_W_support_;

  ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  // //////////////////////////////
  // Private member functions

  void lazyInitializeDefaultBase() const;

  void assert_l(const int l) const;

  void assert_j(const int j) const;

  // //////////////////////////////
  // Private static functions

  static DefaultDerivLinearOpSupport
  determineDefaultDerivLinearOpSupport(
    const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl
    );

  static RCP<LinearOpBase<Scalar> >
  createDefaultLinearOp(
    const DefaultDerivLinearOpSupport &defaultLinearOpSupport,
    const RCP<const VectorSpaceBase<Scalar> > &fnc_space,
    const RCP<const VectorSpaceBase<Scalar> > &var_space
    );

  static ModelEvaluatorBase::DerivativeSupport
  updateDefaultLinearOpSupport(
    const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
    const DefaultDerivLinearOpSupport &defaultLinearOpSupport
    );

  static ModelEvaluatorBase::Derivative<Scalar>
  getOutArgImplForDefaultLinearOpSupport(
    const ModelEvaluatorBase::Derivative<Scalar> &deriv,
    const DefaultDerivLinearOpSupport &defaultLinearOpSupport
    );

  static DefaultDerivMvAdjointSupport
  determineDefaultDerivMvAdjointSupport(
    const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
    const VectorSpaceBase<Scalar> &fnc_space,
    const VectorSpaceBase<Scalar> &var_space
    );

  static ModelEvaluatorBase::DerivativeSupport
  updateDefaultDerivMvAdjointSupport(
    const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
    const DefaultDerivMvAdjointSupport &defaultMvAdjointSupport
    );
  
};


} // namespace Thyra


//
// Implementations
//


#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Overridden from ModelEvaluator


template<class Scalar>
int ModelEvaluatorDefaultBase<Scalar>::Np() const
{
  lazyInitializeDefaultBase();
  return prototypeOutArgs_.Np();
}


template<class Scalar>
int ModelEvaluatorDefaultBase<Scalar>::Ng() const
{
  lazyInitializeDefaultBase();
  return prototypeOutArgs_.Ng();
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_W() const
{
  lazyInitializeDefaultBase();
  if (default_W_support_)
    return this->get_W_factory()->createOp();
  return Teuchos::null;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DfDp_op(int l) const
{
  lazyInitializeDefaultBase();
#ifdef TEUCHOS_DEBUG
  assert_l(l);
#endif
  const DefaultDerivLinearOpSupport
    defaultLinearOpSupport = DfDp_default_op_support_[l];
  if (defaultLinearOpSupport.provideDefaultLinearOp()) {
    return createDefaultLinearOp(
      defaultLinearOpSupport,
      this->get_f_space(),
      this->get_p_space(l)
      );
  }
  return this->create_DfDp_op_impl(l);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDx_dot_op(int j) const
{
  lazyInitializeDefaultBase();
#ifdef TEUCHOS_DEBUG
  assert_j(j);
#endif
  const DefaultDerivLinearOpSupport
    defaultLinearOpSupport = DgDx_dot_default_op_support_[j];
  if (defaultLinearOpSupport.provideDefaultLinearOp()) {
    return createDefaultLinearOp(
      defaultLinearOpSupport,
      this->get_g_space(j),
      this->get_x_space()
      );
  }
  return this->create_DgDx_dot_op_impl(j);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDx_op(int j) const
{
  lazyInitializeDefaultBase();
#ifdef TEUCHOS_DEBUG
  assert_j(j);
#endif
  const DefaultDerivLinearOpSupport
    defaultLinearOpSupport = DgDx_default_op_support_[j];
  if (defaultLinearOpSupport.provideDefaultLinearOp()) {
    return createDefaultLinearOp(
      defaultLinearOpSupport,
      this->get_g_space(j),
      this->get_x_space()
      );
  }
  return this->create_DgDx_op_impl(j);
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDp_op(int j, int l) const
{
  lazyInitializeDefaultBase();
#ifdef TEUCHOS_DEBUG
  assert_j(j);
  assert_l(l);
#endif
  const DefaultDerivLinearOpSupport
    defaultLinearOpSupport = DgDp_default_op_support_[j][l];
  if (defaultLinearOpSupport.provideDefaultLinearOp()) {
    return createDefaultLinearOp(
      defaultLinearOpSupport,
      this->get_g_space(j),
      this->get_p_space(l)
      );
  }
  return this->create_DgDp_op_impl(j,l);
}


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluatorDefaultBase<Scalar>::createOutArgs() const
{
  lazyInitializeDefaultBase();
  return prototypeOutArgs_;
}


template<class Scalar>
void ModelEvaluatorDefaultBase<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::outArg;
  typedef ModelEvaluatorBase MEB;

  lazyInitializeDefaultBase();

  const int l_Np = outArgs.Np();
  const int l_Ng = outArgs.Ng();

  //
  // A) Assert that the inArgs and outArgs object match this class!
  //

#ifdef TEUCHOS_DEBUG
  assertInArgsEvalObjects(*this,inArgs);
  assertOutArgsEvalObjects(*this,outArgs,&inArgs);
#endif  
  
  //
  // B) Setup the OutArgs object for the underlying implementation's
  // evalModelImpl(...) function
  //

  MEB::OutArgs<Scalar> outArgsImpl = this->createOutArgsImpl();
  Array<MultiVectorAdjointPair> DgDp_temp_adjoint_copies;

  {

    outArgsImpl.setArgs(outArgs,true);

    // DfDp(l)
    if (outArgsImpl.supports(MEB::OUT_ARG_f)) {
      for ( int l = 0; l < l_Np; ++l ) {
        const DefaultDerivLinearOpSupport defaultLinearOpSupport =
          DfDp_default_op_support_[l];
        if (defaultLinearOpSupport.provideDefaultLinearOp()) {
          outArgsImpl.set_DfDp( l,
            getOutArgImplForDefaultLinearOpSupport(
              outArgs.get_DfDp(l), defaultLinearOpSupport
              )
            );
        }
        else {
          // DfDp(l) already set by outArgsImpl.setArgs(...)!
        }
      }
    }

    // DgDx_dot(j)
    for ( int j = 0; j < l_Ng; ++j ) {
      const DefaultDerivLinearOpSupport defaultLinearOpSupport =
        DgDx_dot_default_op_support_[j];
      if (defaultLinearOpSupport.provideDefaultLinearOp()) {
        outArgsImpl.set_DgDx_dot( j,
          getOutArgImplForDefaultLinearOpSupport(
            outArgs.get_DgDx_dot(j), defaultLinearOpSupport
            )
          );
      }
      else {
        // DgDx_dot(j) already set by outArgsImpl.setArgs(...)!
      }
    }

    // DgDx(j)
    for ( int j = 0; j < l_Ng; ++j ) {
      const DefaultDerivLinearOpSupport defaultLinearOpSupport =
        DgDx_default_op_support_[j];
      if (defaultLinearOpSupport.provideDefaultLinearOp()) {
        outArgsImpl.set_DgDx( j,
          getOutArgImplForDefaultLinearOpSupport(
            outArgs.get_DgDx(j), defaultLinearOpSupport
            )
          );
      }
      else {
        // DgDx(j) already set by outArgsImpl.setArgs(...)!
      }
    }

    // DgDp(j,l)
    for ( int j = 0; j < l_Ng; ++j ) {
      const Array<DefaultDerivLinearOpSupport> &DgDp_default_op_support_j =
        DgDp_default_op_support_[j];
      const Array<DefaultDerivMvAdjointSupport> &DgDp_default_mv_support_j =
        DgDp_default_mv_support_[j];
      for ( int l = 0; l < l_Np; ++l ) {
        const DefaultDerivLinearOpSupport defaultLinearOpSupport =
          DgDp_default_op_support_j[l];
        const DefaultDerivMvAdjointSupport defaultMvAdjointSupport =
          DgDp_default_mv_support_j[l];
        MEB::Derivative<Scalar> DgDp_j_l;
        if (!outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none())
          DgDp_j_l = outArgs.get_DgDp(j,l);
        if (
          defaultLinearOpSupport.provideDefaultLinearOp()
          && !is_null(DgDp_j_l.getLinearOp())
          )
        {
          outArgsImpl.set_DgDp( j, l,
            getOutArgImplForDefaultLinearOpSupport(
              DgDp_j_l, defaultLinearOpSupport
              )
            );
        }
        else if (
          defaultMvAdjointSupport.provideDefaultAdjoint()
          && !is_null(DgDp_j_l.getMultiVector())
          )
        {
          const RCP<MultiVectorBase<Scalar> > DgDp_j_l_mv = 
            DgDp_j_l.getMultiVector();
          if (
            defaultMvAdjointSupport.mvAdjointCopyOrientation()
            ==
            DgDp_j_l.getMultiVectorOrientation()
            )
          {
            // The orientation of the multi-vector is different so we need to
            // create a temporary copy to pass to evalModelImpl(...) and then
            // copy it back again!
            const RCP<MultiVectorBase<Scalar> > DgDp_j_l_mv_adj =
              createMembers(DgDp_j_l_mv->domain(), DgDp_j_l_mv->range()->dim());
            outArgsImpl.set_DgDp( j, l,
              MEB::Derivative<Scalar>(
                DgDp_j_l_mv_adj,
                getOtherDerivativeMultiVectorOrientation(
                  defaultMvAdjointSupport.mvAdjointCopyOrientation()
                  )
                )
              );
            // Remember these multi-vectors so that we can do the transpose copy
            // back after the evaluation!
            DgDp_temp_adjoint_copies.push_back(
              MultiVectorAdjointPair(DgDp_j_l_mv, DgDp_j_l_mv_adj)
              );
          }
          else {
            // The form of the multi-vector is supported by evalModelImpl(..)
            // and is already set on the outArgsImpl object.
          }
        }
        else {
          // DgDp(j,l) already set by outArgsImpl.setArgs(...)!
        }
      }
    }

    // W
    {
      RCP<LinearOpWithSolveBase<Scalar> > W;
      if ( default_W_support_ && !is_null(W=outArgs.get_W()) ) {
        const RCP<const LinearOpWithSolveFactoryBase<Scalar> >
          W_factory = this->get_W_factory();
        // Extract the underlying W_op object (if it exists)
        RCP<const LinearOpBase<Scalar> > W_op_const;
        uninitializeOp<Scalar>(*W_factory, W.ptr(), outArg(W_op_const));
        RCP<LinearOpBase<Scalar> > W_op;
        if (!is_null(W_op_const)) {
          // Here we remove the const.  This is perfectly safe since
          // w.r.t. this class, we put a non-const object in there and we can
          // expect to change that object after the fact.  That is our
          // prerogative.
          W_op = Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(W_op_const);
        }
        else {
          // The W_op object has not been initialized yet so create it.  The
          // next time through, we should not have to do this!
          W_op = this->create_W_op();
        }
        outArgsImpl.set_W_op(W_op);
      }
    }
    
  }

  //
  // C) Evaluate the underlying model implementation!
  //

  this->evalModelImpl( inArgs, outArgsImpl );

  //
  // D) Post-process the output arguments
  //

  // Do explicit transposes for DgDp(j,l) if needed
  const int numMvAdjointCopies = DgDp_temp_adjoint_copies.size();
  for ( int adj_copy_i = 0; adj_copy_i < numMvAdjointCopies; ++adj_copy_i ) {
    const MultiVectorAdjointPair adjPair =
      DgDp_temp_adjoint_copies[adj_copy_i];
    doExplicitMultiVectorAdjoint( *adjPair.mvImplAdjoint, &*adjPair.mvOuter );
  }
  
  // Update W given W_op and W_factory
  {
    RCP<LinearOpWithSolveBase<Scalar> > W;
    if ( default_W_support_ && !is_null(W=outArgs.get_W()) ) {
      const RCP<const LinearOpWithSolveFactoryBase<Scalar> >
        W_factory = this->get_W_factory();
      W_factory->setOStream(this->getOStream());
      W_factory->setVerbLevel(this->getVerbLevel());
      initializeOp<Scalar>(*W_factory, outArgsImpl.get_W_op().getConst(), W.ptr());
    }
  }
  
}


// protected


// Setup functions called by subclasses

template<class Scalar>
void ModelEvaluatorDefaultBase<Scalar>::initializeDefaultBase()
{

  typedef ModelEvaluatorBase MEB;

  // In case we throw half way thorugh, set to uninitialized
  isInitialized_ = false;
  default_W_support_ = false;

  //
  // A) Get the InArgs and OutArgs from the subclass
  //
  
  const MEB::InArgs<Scalar> inArgs = this->createInArgs();
  const MEB::OutArgs<Scalar> outArgsImpl = this->createOutArgsImpl();

  //
  // B) Validate the subclasses InArgs and OutArgs
  //

#ifdef TEUCHOS_DEBUG
  assertInArgsOutArgsSetup( this->description(), inArgs, outArgsImpl );
#endif // TEUCHOS_DEBUG

  //
  // C) Set up support for default derivative objects and prototype OutArgs
  //

  const int l_Ng = outArgsImpl.Ng();
  const int l_Np = outArgsImpl.Np();

  // Set support for all outputs supported in the underly implementation
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(l_Np,l_Ng);
  outArgs.setSupports(outArgsImpl);

  // DfDp
  DfDp_default_op_support_.clear();
  if (outArgs.supports(MEB::OUT_ARG_f)) {
    for ( int l = 0; l < l_Np; ++l ) {
      const MEB::DerivativeSupport DfDp_l_impl_support =
        outArgsImpl.supports(MEB::OUT_ARG_DfDp,l);
      const DefaultDerivLinearOpSupport DfDp_l_op_support =
        determineDefaultDerivLinearOpSupport(DfDp_l_impl_support);
      DfDp_default_op_support_.push_back(DfDp_l_op_support);
      outArgs.setSupports(
        MEB::OUT_ARG_DfDp, l,
        updateDefaultLinearOpSupport(
          DfDp_l_impl_support, DfDp_l_op_support
          )
        );
    }
  }

  // DgDx_dot
  DgDx_dot_default_op_support_.clear();
  for ( int j = 0; j < l_Ng; ++j ) {
    const MEB::DerivativeSupport DgDx_dot_j_impl_support =
      outArgsImpl.supports(MEB::OUT_ARG_DgDx_dot,j);
    const DefaultDerivLinearOpSupport DgDx_dot_j_op_support =
      determineDefaultDerivLinearOpSupport(DgDx_dot_j_impl_support);
    DgDx_dot_default_op_support_.push_back(DgDx_dot_j_op_support);
    outArgs.setSupports(
      MEB::OUT_ARG_DgDx_dot, j,
      updateDefaultLinearOpSupport(
        DgDx_dot_j_impl_support, DgDx_dot_j_op_support
        )
      );
  }

  // DgDx
  DgDx_default_op_support_.clear();
  for ( int j = 0; j < l_Ng; ++j ) {
    const MEB::DerivativeSupport DgDx_j_impl_support =
      outArgsImpl.supports(MEB::OUT_ARG_DgDx,j);
    const DefaultDerivLinearOpSupport DgDx_j_op_support =
      determineDefaultDerivLinearOpSupport(DgDx_j_impl_support);
    DgDx_default_op_support_.push_back(DgDx_j_op_support);
    outArgs.setSupports(
      MEB::OUT_ARG_DgDx, j,
      updateDefaultLinearOpSupport(
        DgDx_j_impl_support, DgDx_j_op_support
        )
      );
  }

  // DgDp
  DgDp_default_op_support_.clear();
  DgDp_default_mv_support_.clear();
  for ( int j = 0; j < l_Ng; ++j ) {
    DgDp_default_op_support_.push_back(Array<DefaultDerivLinearOpSupport>());
    DgDp_default_mv_support_.push_back(Array<DefaultDerivMvAdjointSupport>());
    for ( int l = 0; l < l_Np; ++l ) {
      const MEB::DerivativeSupport DgDp_j_l_impl_support =
        outArgsImpl.supports(MEB::OUT_ARG_DgDp,j,l);
      // LinearOpBase support
      const DefaultDerivLinearOpSupport DgDp_j_l_op_support =
        determineDefaultDerivLinearOpSupport(DgDp_j_l_impl_support);
      DgDp_default_op_support_[j].push_back(DgDp_j_l_op_support);
      outArgs.setSupports(
        MEB::OUT_ARG_DgDp, j, l,
        updateDefaultLinearOpSupport(
          DgDp_j_l_impl_support, DgDp_j_l_op_support
          )
        );
      // MultiVectorBase
      const DefaultDerivMvAdjointSupport DgDp_j_l_mv_support = 
        determineDefaultDerivMvAdjointSupport(
          DgDp_j_l_impl_support, *this->get_g_space(j), *this->get_p_space(l)
          );
      DgDp_default_mv_support_[j].push_back(DgDp_j_l_mv_support);
      outArgs.setSupports(
        MEB::OUT_ARG_DgDp, j, l,
        updateDefaultDerivMvAdjointSupport(
          outArgs.supports(MEB::OUT_ARG_DgDp, j, l),
          DgDp_j_l_mv_support
          )
        );
    }
  }
  // 2007/09/09: rabart: ToDo: Move the above code into a private helper
  // function!
  
  // W (given W_op and W_factory)
  default_W_support_ = false;
  if ( outArgsImpl.supports(MEB::OUT_ARG_W_op) && !is_null(this->get_W_factory())
    && !outArgsImpl.supports(MEB::OUT_ARG_W) )
  {
    default_W_support_ = true;
    outArgs.setSupports(MEB::OUT_ARG_W);
    outArgs.set_W_properties(outArgsImpl.get_W_properties());
  }
  
  //
  // D) All done!
  //

  prototypeOutArgs_ = outArgs;
  isInitialized_ = true;

}


// Private functions with default implementaton to be overridden by subclasses


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DfDp_op_impl(int l) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgsImpl();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DfDp,l).supports(MEB::DERIV_LINEAR_OP),
    std::logic_error,
    "Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DfDp("<<l<<") (as determined from its"
    " OutArgs object created by createOutArgsImpl())"
    " but this function create_DfDp_op_impl(...) has not been overriden"
    " to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDx_dot_op_impl(int j) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgsImpl();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DgDx_dot,j).supports(MEB::DERIV_LINEAR_OP),
    std::logic_error,
    "Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DgDx_dot("<<j<<") (as determined from"
    " its OutArgs object created by createOutArgsImpl())"
    " but this function create_DgDx_dot_op_impl(...) has not been overriden"
    " to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDx_op_impl(int j) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgsImpl();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DgDx,j).supports(MEB::DERIV_LINEAR_OP),
    std::logic_error,
    "Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DgDx("<<j<<") (as determined from"
    " its OutArgs object created by createOutArgsImpl())"
    " but this function create_DgDx_op_impl(...) has not been overriden"
    " to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgsImpl();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DgDp,j,l).supports(MEB::DERIV_LINEAR_OP),
    std::logic_error,
    "Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DgDp("<<j<<","<<l<<")"
    " (as determined from its OutArgs object created by createOutArgsImpl())"
    " but this function create_DgDp_op_impl(...) has not been overriden"
    " to create such an object!"
    );
  return Teuchos::null;
}


// protected


template<class Scalar>
ModelEvaluatorDefaultBase<Scalar>::ModelEvaluatorDefaultBase()
  :isInitialized_(false), default_W_support_(false)
{}


// private


template<class Scalar>
void ModelEvaluatorDefaultBase<Scalar>::lazyInitializeDefaultBase() const
{
  if (!isInitialized_)
    const_cast<ModelEvaluatorDefaultBase<Scalar>*>(this)->initializeDefaultBase();
}


template<class Scalar>
void ModelEvaluatorDefaultBase<Scalar>::assert_l(const int l) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l,0,this->Np());
}


template<class Scalar>
void ModelEvaluatorDefaultBase<Scalar>::assert_j(const int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j,0,this->Ng());
}


// Private static functions


template<class Scalar>
ModelEvaluatorDefaultBaseTypes::DefaultDerivLinearOpSupport
ModelEvaluatorDefaultBase<Scalar>::determineDefaultDerivLinearOpSupport(
  const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl
  )
{
  typedef ModelEvaluatorBase MEB;
  if (
    (
      derivSupportImpl.supports(MEB::DERIV_MV_BY_COL)
      ||
      derivSupportImpl.supports(MEB::DERIV_TRANS_MV_BY_ROW)
      )
    &&
    !derivSupportImpl.supports(MEB::DERIV_LINEAR_OP)
    )
  {
    return DefaultDerivLinearOpSupport(
      derivSupportImpl.supports(MEB::DERIV_MV_BY_COL)
      ? MEB::DERIV_MV_BY_COL
      : MEB::DERIV_TRANS_MV_BY_ROW
      );
  }
  return DefaultDerivLinearOpSupport();
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorDefaultBase<Scalar>::createDefaultLinearOp(
  const DefaultDerivLinearOpSupport &defaultLinearOpSupport,
  const RCP<const VectorSpaceBase<Scalar> > &fnc_space,
  const RCP<const VectorSpaceBase<Scalar> > &var_space
  )
{
  using Teuchos::rcp_implicit_cast;
  typedef LinearOpBase<Scalar> LOB;
  typedef ModelEvaluatorBase MEB;
  switch(defaultLinearOpSupport.mvImplOrientation()) {
    case MEB::DERIV_MV_BY_COL:
      // The MultiVector will do just fine as the LinearOpBase
      return createMembers(fnc_space, var_space->dim());
    case MEB::DERIV_TRANS_MV_BY_ROW:
      // We will have to implicitly transpose the underlying MultiVector
      return nonconstAdjoint<Scalar>(
        rcp_implicit_cast<LOB>(createMembers(var_space, fnc_space->dim()))
        );
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
  }
  return Teuchos::null; // Will never be called!
}


template<class Scalar>
ModelEvaluatorBase::DerivativeSupport
ModelEvaluatorDefaultBase<Scalar>::updateDefaultLinearOpSupport(
  const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
  const DefaultDerivLinearOpSupport &defaultLinearOpSupport
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::DerivativeSupport derivSupport = derivSupportImpl;
  if (defaultLinearOpSupport.provideDefaultLinearOp())
    derivSupport.plus(MEB::DERIV_LINEAR_OP);
  return derivSupport;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorDefaultBase<Scalar>::getOutArgImplForDefaultLinearOpSupport(
  const ModelEvaluatorBase::Derivative<Scalar> &deriv,
  const DefaultDerivLinearOpSupport &defaultLinearOpSupport
  )
{

  using Teuchos::rcp_dynamic_cast;
  typedef ModelEvaluatorBase MEB;
  typedef MultiVectorBase<Scalar> MVB;
  typedef ScaledAdjointLinearOpBase<Scalar> SALOB;

  // If the derivative is not a LinearOpBase then it it either empty or is a
  // MultiVectorBase object, and in either case we just return it since the
  // underlying evalModelImpl(...) functions should be able to handle it!
  if (is_null(deriv.getLinearOp()))
    return deriv;

  // The derivative is LinearOpBase so get out the underlying MultiVectorBase
  // object and return its derivative multi-vector form.
  switch(defaultLinearOpSupport.mvImplOrientation()) {
    case MEB::DERIV_MV_BY_COL: {
      return MEB::Derivative<Scalar>(
        rcp_dynamic_cast<MVB>(deriv.getLinearOp(),true),
        MEB::DERIV_MV_BY_COL
        );
    }
    case MEB::DERIV_TRANS_MV_BY_ROW: {
      return MEB::Derivative<Scalar>(
        rcp_dynamic_cast<MVB>(
          rcp_dynamic_cast<SALOB>(deriv.getLinearOp(),true)->getNonconstOrigOp()
          ),
        MEB::DERIV_TRANS_MV_BY_ROW
        );
    }
#ifdef TEUCHOS_DEBUG
    default:
      TEST_FOR_EXCEPT(true);
#endif
  }

  return ModelEvaluatorBase::Derivative<Scalar>(); // Should never get here!

}


template<class Scalar>
ModelEvaluatorDefaultBaseTypes::DefaultDerivMvAdjointSupport
ModelEvaluatorDefaultBase<Scalar>::determineDefaultDerivMvAdjointSupport(
  const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
  const VectorSpaceBase<Scalar> &fnc_space,
  const VectorSpaceBase<Scalar> &var_space
  )
{
  typedef ModelEvaluatorBase MEB;
  // Here we will support the adjoint copy for of a multi-vector if both
  // spaces give rise to in-core vectors.
  const bool implSupportsMv =
    ( derivSupportImpl.supports(MEB::DERIV_MV_BY_COL)
      || derivSupportImpl.supports(MEB::DERIV_TRANS_MV_BY_ROW) );
  const bool bothSpacesHaveInCoreViews =
    ( fnc_space.hasInCoreView() && var_space.hasInCoreView() );
  if ( implSupportsMv && bothSpacesHaveInCoreViews ) {
    return DefaultDerivMvAdjointSupport(
      derivSupportImpl.supports(MEB::DERIV_MV_BY_COL)
      ? MEB::DERIV_TRANS_MV_BY_ROW
      : MEB::DERIV_MV_BY_COL
      );
  }
  // We can't provide an adjoint copy or such a copy is not needed!
  return DefaultDerivMvAdjointSupport();
}


template<class Scalar>
ModelEvaluatorBase::DerivativeSupport
ModelEvaluatorDefaultBase<Scalar>::updateDefaultDerivMvAdjointSupport(
  const ModelEvaluatorBase::DerivativeSupport &derivSupportImpl,
  const DefaultDerivMvAdjointSupport &defaultMvAdjointSupport
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::DerivativeSupport derivSupport = derivSupportImpl;
  if (defaultMvAdjointSupport.provideDefaultAdjoint())
    derivSupport.plus(defaultMvAdjointSupport.mvAdjointCopyOrientation());
  return derivSupport;
}


} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_DEFAULT_BASE_HPP
