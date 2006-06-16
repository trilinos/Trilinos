// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef THYRA_DEFAULT_MULTI_PERIOD_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_MULTI_PERIOD_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluator.hpp"

namespace Thyra {

/** \brief Composite subclass that takes a list of individual
 * <tt>ModelEvaluator</tt> objects and represents them as a single aggregate
 * multi-period <tt>ModelEvalator</tt> object.
 *
 * This class represents fairly classic "composite" subclass that takes
 * <tt>N</tt> <tt>ModelEvaluator</tt> objects and makes them appear to be just
 * one.
 *
 * Mathematically, this subclass is used to form a multi-period problem of the form:
 
 \verbatim

   (x_bar,{p_l}\p_period) -> f_bar

   (x_bar,{p_l}\p_period) -> g_bar

 \endverbatim

 * where

 \verbatim

   x_bar = [ x[0], x[1], ..., x[N-1]]

   {p_l}\p_period = { p_0, p_1, ..., p_{period-1}, p_{period+1}, ..., p_Np

   f_bar(...) = [ f[0](x[0],{p_l}), f[1](x[1],{p_l}), ..., f[N-1](x[N-1],{p_l}) ]

   g_bar(...) = sum( w[i]*g(x[0],{p_l}\p_period,p_peroid[i]), i = 0,...,N-1 )

 \endverbatim

 * Above, the notation <tt>{p_l}\p_period</tt> is meant to represent the set
 * of all parameter subvectors in each of the constituent models minus the
 * parameter used to define period data.
 *
 * Also, above we use <tt>g(...)</tt> and <tt>g_bar(..)</tt> to represent the
 * zeroth axillary response function unless otherwise noted.
 *
 * This class could be made much more general but for now ...
 */
template<class Scalar>
class DefaultMultiPeriodModelEvaluator : virtual public ModelEvaluator<Scalar> {
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  DefaultMultiPeriodModelEvaluator();

  /** \brief . */
  DefaultMultiPeriodModelEvaluator(
    const int                                                     N
    ,const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >          models[]
    ,const Scalar                                                 weights[]
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        z[]
    ,const int                                                    z_index
    ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &x_bar_space = Teuchos::null
    ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &f_bar_space = Teuchos::null
    //,const AbstractFactory<BlockLinearOpWithSolveBase<Scalar> >   &blockLOWSFF = Teuchos::null
    );

  /** \brief Initialize.
   *
   * \param  N
   *           [in] The number of periods.
   * \param  models
   *           [in] Array (length <tt>N</tt>) of the constituent period models 
   * \param  weights
   *           [in] Array (length <tt>N</tt>) of the weights for the auxiliary
   *           response functions in the summation for <tt>g_bar(...)</tt>
   *           shown above.
   * \param  z
   *           [in] Array (length <tt>N</tt>) of the period defining auxiliary
   *           paramter subvectors.  If <tt>z==NULL</tt> then these are
   *           assumed to not be needed to define the period's model
   *           <tt>models[i]</tt>.  Likewise, if <tt>z!=NULL</tt> but
   *           <tt>z[i].get()==NULL</tt> then it is assumed that
   *           <tt>models[i]</tt> does not need a parameter subvector to
   *           define its state.
   * \param  z_index
   *           [in] Index of the model's parameter subvector that represents
   *           <tt>z[i]</tt>.  In other words <tt>p_{z_index} == z[i]</tt>.
   *           This assumes that each model <tt>models[i]</tt> will have the
   *           same parameter subvector structure.  <tt>z==NULL</tt> then this
   *           argument is ignored.
   * \param  x_bar_space
   *           [in] The product vector space that represents the space for
   *           <tt>x_bar</tt> as defined above.  If
   *           <tt>x_bar_space.get()==NULL</tt> then this will be created
   *           internally.
   * \param  f_bar_space
   *           [in] The product vector space that represents the space for
   *           <tt>f_bar</tt> as defined above.  If
   *           <tt>f_bar_space.get()==NULL</tt> then this will be created
   *           internally.
   * \param  blockLOWSFF
   *           [in] Factory object that is used to create the block LOWS
   *           object that will be used to represent the block diagonal
   *           object. <tt>W_bar</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>N > 0</tt>
   * <li> <tt>models!=NULL</tt> and <tt>models[i].get()!=NULL</tt>,
   *      for <tt>i=0...N-1</tt>.
   * <li> [<li>z_index < 0</tt>] <tt>z==NULL</tt>
   * <li> [<tt>z_index >= 0</tt>] <tt>z!=NULL</tt> and <tt>z[i].get()!=NULL</tt>,
   *      for <tt>i=0...N-1</tt>.
   * </ul>
   * 
   */
  void initialize(
    const int                                                     N
    ,const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >          models[]
    ,const Scalar                                                 weights[]
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        z[]
    ,const int                                                    z_index
    ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &x_bar_space = Teuchos::null
    ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &f_bar_space = Teuchos::null
    //,const AbstractFactory<BlockLinearOpWithSolveBase<Scalar> >   &blockLOWSFF = Teuchos::null
    );
  // ToDo: To be more general we would have to define a z_index for each
  // *models[i] object and we would need parameter index maps to map from the
  // parameter subvectors in *models[i] to the a "global" list of parameter
  // subvectors.  This could be added as a different initialization function
  // for the more general use case.
  
  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief. */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \breif . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DfDp_op(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>       &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>     &outArgs
    ) const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

private:

  // //////////////////////////
  // Private types

  typedef std::vector<Teuchos::RefCountPtr<ModelEvaluator<Scalar> > >    models_t;
  typedef std::vector<Scalar>                                            weights_t;
  typedef std::vector<Teuchos::RefCountPtr<const VectorBase<Scalar> > >  z_t;

  // /////////////////////////
  // Private data members

  models_t                                               models_;  // size == N
  weights_t                                              weights_; // size == N
  z_t                                                    z_;       // size == N
  int                                                    z_index_;
  Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  x_bar_space_;
  Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  f_bar_space_;
  
};

// /////////////////////////////////
// Implementations

// Constructors/Intializers/Accessors

template<class Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::DefaultMultiPeriodModelEvaluator()
  :z_index_(-1)
{}

template<class Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::DefaultMultiPeriodModelEvaluator(
  const int                                                     N
  ,const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >          models[]
  ,const Scalar                                                 weights[]
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        z[]
  ,const int                                                    z_index
  ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &x_bar_space
  ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &f_bar_space
  //,const AbstractFactory<BlockLinearOpWithSolveBase<Scalar> >   &blockLOWSFF
  )
{
  initialize(N,models,weights,z,z_index,x_bar_space,f_bar_space);
}

template<class Scalar>
void DefaultMultiPeriodModelEvaluator<Scalar>::initialize(
  const int                                                     N
  ,const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >          models[]
  ,const Scalar                                                 weights[]
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        z[]
  ,const int                                                    z_index
  ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &x_bar_space
  ,const Teuchos::RefCountPtr<ProductVectorSpaceBase<Scalar> >  &f_bar_space
  //,const AbstractFactory<BlockLinearOpWithSolveBase<Scalar> >   &blockLOWSFF
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(!(N>0));
  TEST_FOR_EXCEPT(!(models));
  TEST_FOR_EXCEPT(!((z_index>=0&&z!=NULL)||(z_index==0&&z==NULL)));
  models_.resize(N);
  weights_.resize(N);
  z_.resize(N);
  z_index_ = z_index;
  for( int i = 0; i < N; ++i ) {
    models_[i] = models[i].assert_not_null();
    weights_[i] = ( weights ? weights[i] : ST::one() );
    if( z && z[i].get() )
      z_[i] = z[i];
    else
      z_[i] = Teuchos::null;
  }
  if( x_bar_space.get() ) {
    TEST_FOR_EXCEPT(!(x_bar_space->numBlocks()==N));
    // ToDo: Check the constituent spaces more carefully against models[]->get_x_space().
    x_bar_space_ = x_bar_space;
  }
  else {
    TEST_FOR_EXCEPT(true); // ToDo: Create our own product space!
  }
  if( f_bar_space.get() ) {
    TEST_FOR_EXCEPT(!(f_bar_space->numBlocks()==N));
    // ToDo: Check the constituent spaces more carefully against models[]->get_f_space().
    f_bar_space_ = f_bar_space;
  }
  else {
    TEST_FOR_EXCEPT(true); // ToDo: Create our own product space!
  }
}

// Public functions overridden from ModelEvaulator

template<class Scalar>
int DefaultMultiPeriodModelEvaluator<Scalar>::Np() const
{
  TEST_FOR_EXCEPT(true);
  return 0;
}

template<class Scalar>
int DefaultMultiPeriodModelEvaluator<Scalar>::Ng() const
{
  TEST_FOR_EXCEPT(true);
  return 0;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::get_f_space() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::get_p_space(int l) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::get_g_space(int j) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::getNominalValues() const
{
  TEST_FOR_EXCEPT(true);
  return ModelEvaluatorBase::InArgs<Scalar>();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::getLowerBounds() const
{
  TEST_FOR_EXCEPT(true);
  return ModelEvaluatorBase::InArgs<Scalar>();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::getUpperBounds() const
{
  TEST_FOR_EXCEPT(true);
  return ModelEvaluatorBase::InArgs<Scalar>();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::create_W() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null; // Should never be called!
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null; // Should never be called!
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::create_DfDp_op(int l) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::create_DgDx_op(int j) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultMultiPeriodModelEvaluator<Scalar>::create_DgDp_op( int j, int l ) const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::createInArgs() const
{
  ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  TEST_FOR_EXCEPT(true);
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultMultiPeriodModelEvaluator<Scalar>::createOutArgs() const
{
  ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  TEST_FOR_EXCEPT(true);
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}

template<class Scalar>
void DefaultMultiPeriodModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>       &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>     &outArgs
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void DefaultMultiPeriodModelEvaluator<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  TEST_FOR_EXCEPT(true);
}

} // namespace Thyra

#endif // THYRA_DEFAULT_MULTI_PERIOD_MODEL_EVALUATOR_HPP
