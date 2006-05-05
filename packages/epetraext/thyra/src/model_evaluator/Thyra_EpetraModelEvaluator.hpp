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

#ifndef THYRA_EPETRA_MODEL_EVALUATOR_HPP
#define THYRA_EPETRA_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"

namespace Thyra {

/** \brief . */
class EpetraModelEvaluator : public ModelEvaluator<double> {
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  EpetraModelEvaluator();

  /** \brief . */
  EpetraModelEvaluator(
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>         &epetraModel
    ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >  &W_factory
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>         &epetraModel
    ,const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >  &W_factory
    );

  /** \brief . */
  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> getEpetraModel() const;

  /** \brief . */
  void setInitialGuess( const ModelEvaluatorBase::InArgs<double>& initialGuess );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>         *epetraModel = NULL
    ,Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >  *W_factory   = NULL
    );

  /** \brief . */
  const ModelEvaluatorBase::InArgs<double>& getFinalPoint() const;

  /** \brief . */
  bool finalPointWasSolved() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_f_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getUpperBounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_W_op() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DfDp_op(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DgDx_op(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<double> > create_DgDp_op( int j, int l ) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<double> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<double>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<double>  &outArgs
    ) const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<double>      &finalPoint
    ,const bool                                   wasSolved
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // ////////////////////
  // Private types

  typedef std::vector<Teuchos::RefCountPtr<const Epetra_Map> > p_map_t;
  typedef std::vector<Teuchos::RefCountPtr<const Epetra_Map> > g_map_t;

  typedef std::vector<Teuchos::RefCountPtr<const MPIVectorSpaceDefaultBase<double> > > p_space_t;
  typedef std::vector<Teuchos::RefCountPtr<const MPIVectorSpaceDefaultBase<double> > > g_space_t;

  // ////////////////////
  // Private data mebers

  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>              epetraModel_;
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >        W_factory_;
  Teuchos::RefCountPtr<const Epetra_Map>                             x_map_;
  p_map_t                                                            p_map_;
  g_map_t                                                            g_map_;
  Teuchos::RefCountPtr<const Epetra_Map>                             f_map_;
  Teuchos::RefCountPtr<const MPIVectorSpaceDefaultBase<double> >     x_space_;
  p_space_t                                                          p_space_;
  Teuchos::RefCountPtr<const MPIVectorSpaceDefaultBase<double> >     f_space_;
  g_space_t                                                          g_space_;
  ModelEvaluatorBase::InArgs<double>                                 initialGuess_;
  ModelEvaluatorBase::InArgs<double>                                 lowerBounds_;
  ModelEvaluatorBase::InArgs<double>                                 upperBounds_;
  ModelEvaluatorBase::InArgs<double>                                 finalPoint_;
  bool                                                               finalPointWasSolved_;
  
};

//
// Utility functions
//

/** \brief . */
ModelEvaluatorBase::EDerivativeMultiVectorOrientation
convert( const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation &mvOrientation );

/** \brief . */
EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation
convert( const ModelEvaluatorBase::EDerivativeMultiVectorOrientation &mvOrientation );

/** \brief . */
ModelEvaluatorBase::DerivativeProperties
convert( const EpetraExt::ModelEvaluator::DerivativeProperties &derivativeProperties );

/** \brief . */
ModelEvaluatorBase::DerivativeSupport
convert( const EpetraExt::ModelEvaluator::DerivativeSupport &derivativeSupport );

/** \brief . */
EpetraExt::ModelEvaluator::Derivative
convert(
  const ModelEvaluatorBase::Derivative<double>        &derivative
  ,const Teuchos::RefCountPtr<const Epetra_Map>       &fnc_map
  ,const Teuchos::RefCountPtr<const Epetra_Map>       &var_map
  );

} // namespace Thyra

#endif // THYRA_EPETRA_MODEL_EVALUATOR_HPP
