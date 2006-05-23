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

#ifndef THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP
#define THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief This is a base class that delegetes almost all function to a
 * wrapped model evaluator object.
 *
 * This class makes it easy to write many different types of decorator
 * subclasses by only requiring those subclasses to override just the behavior
 * that they need to overide and leave the rest alone.  Note that whenever the
 * definition of <tt>ModelEvaluator</tt> changes, this class will have to be
 * updated.  However, the decorator subclasses that derive from this base
 * class might be able to remain unchanged.  Note that changes in only the
 * <tt>ModelEvaluatorBase::InArgs</tt> and
 * <tt>ModelEvaluatorBase::OutArgs</tt> classes should not require any changes
 * here.
 * 
 * The only functions that a client must override in order to create a
 * concrete subcalss is the <tt>evalModel()</tt> function.  All other
 * functions have implementations here that simply delegate to the model
 * evaluator object returned from <tt>getUnderlyingModel()</tt>.  However,
 * most decorator classes will need to override at least one other function.
 *
 * This class provides a default implemntation of the machinary to store and
 * access the wrapped model evaluator object.  A subclass can choose to ignore
 * this and override the functions <tt>isUnderlyingModelConst()<tt>,
 * <tt>getConstUnderlyingModel()</tt>, and <tt>getUnderlyingModel()</tt>.
 */
template<class Scalar>
class ModelEvaluatorDelegatorBase : virtual public ModelEvaluator<Scalar> {
public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Constructs to uninitialized. */
  ModelEvaluatorDelegatorBase();

  /** \brief Calls <tt>initialize()</tt>. */
  ModelEvaluatorDelegatorBase(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &model
    );

  /** \brief Calls <tt>initialize()</tt>. */
  ModelEvaluatorDelegatorBase(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >   &model
    );

  /** \brief Initialize given a non-const model evaluator. */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &model
    );

  /** \brief Initialize given a const model evaluator. */
  void initialize(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >   &model
    );

  /** \brief Uninitialize. */
  void uninitialize();

  //@}

  /** \name Virtual functions that can overriden */
  //@{

  /** \brief . */
  virtual bool isUnderlyingModelConst() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<ModelEvaluator<Scalar> > getNonconstUnderlyingModel();

  /** \brief . */
  virtual Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > getUnderlyingModel() const;

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
  /** \brief . */
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
  /** \brief . */
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
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

private:

  Teuchos::ConstNonconstObjectContainer<ModelEvaluator<Scalar> > model_;
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers

template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase()
{}

template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &model
  )
{
  this->initialize(model);
}

template<class Scalar>
ModelEvaluatorDelegatorBase<Scalar>::ModelEvaluatorDelegatorBase(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >   &model
  )
{
  this->initialize(model);
}

template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &model
  )
{
  model_.initialize(model);
}

template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >   &model
  )
{
  model_.initialize(model);
}

template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::uninitialize()
{
  model_.uninitialize();
}

// Virtual functions that can overriden

template<class Scalar>
bool ModelEvaluatorDelegatorBase<Scalar>::isUnderlyingModelConst() const
{
  return model_.isConst();
}

template<class Scalar>
Teuchos::RefCountPtr<ModelEvaluator<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::getNonconstUnderlyingModel()
{
  return model_.getNonconstObj();
}

template<class Scalar>
Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::getUnderlyingModel() const
{
  return model_.getConstObj();
}

// Overridden from ModelEvaulator.

template<class Scalar>
int ModelEvaluatorDelegatorBase<Scalar>::Np() const
{
  return getUnderlyingModel()->Np();
}

template<class Scalar>
int ModelEvaluatorDelegatorBase<Scalar>::Ng() const
{
  return getUnderlyingModel()->Ng();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_x_space() const
{
  return getUnderlyingModel()->get_x_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_f_space() const
{
  return getUnderlyingModel()->get_f_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_p_space(int l) const
{
  return getUnderlyingModel()->get_p_space(l);
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::get_g_space(int j) const
{
  return getUnderlyingModel()->get_g_space(j);
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getNominalValues() const
{
  return getUnderlyingModel()->getNominalValues();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getLowerBounds() const
{
  return getUnderlyingModel()->getLowerBounds();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::getUpperBounds() const
{
  return getUnderlyingModel()->getUpperBounds();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_W() const
{
  return getUnderlyingModel()->create_W();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_W_op() const
{
  return getUnderlyingModel()->create_W_op();
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DfDp_op(int l) const
{
  return getUnderlyingModel()->create_DfDp_op(l);
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DgDx_op(int j) const
{
  return getUnderlyingModel()->create_DgDx_op(j);
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
ModelEvaluatorDelegatorBase<Scalar>::create_DgDp_op( int j, int l ) const
{
  return getUnderlyingModel()->create_DgDp_op(j,l);
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::createInArgs() const
{
  ModelEvaluatorBase::InArgsSetup<Scalar> inArgs = getUnderlyingModel()->createInArgs();
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluatorDelegatorBase<Scalar>::createOutArgs() const
{
  ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs = getUnderlyingModel()->createOutArgs();
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}

template<class Scalar>
void ModelEvaluatorDelegatorBase<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  getNonconstUnderlyingModel()->reportFinalPoint(finalPoint,wasSolved);
}

} // namespace Thyra

#endif // THYRA_DEFAULT_MODEL_EVALUATOR_DELEGETOR_BASE_HPP
