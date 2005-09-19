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
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>               &epetraModel
    ,const Teuchos::RefCountPtr<const LinearOpWithSolveFactoryBase<double> >  &W_factory
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>               &epetraModel
    ,const Teuchos::RefCountPtr<const LinearOpWithSolveFactoryBase<double> >  &W_factory
    );

  /** \brief . */
  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> getEpetraModel() const;

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>               *epetraModel = NULL
    ,Teuchos::RefCountPtr<const LinearOpWithSolveFactoryBase<double> >  *W_factory   = NULL
    );

  //@}

  /** \name Overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_x_space() const;

  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<double> > get_f_space() const;

  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<double> > get_x_init() const;

  /** \brief . */
  double get_t_init() const;

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > create_W() const;

  /** \brief . */
  ModelEvaluatorBase::InArgs<double> createInArgs() const;

  /** \brief . */
  ModelEvaluatorBase::OutArgs<double> createOutArgs() const;

  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<double>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<double>  &outArgs
    ) const;

  //@}

private:

  Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>              epetraModel_;
  Teuchos::RefCountPtr<const LinearOpWithSolveFactoryBase<double> >  W_factory_;
  Teuchos::RefCountPtr<const Epetra_Map>                    x_map_;
  Teuchos::RefCountPtr<const Epetra_Map>                    f_map_;
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >   x_space_;
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >   f_space_;
  
};

} // namespace Thyra

#endif // THYRA_EPETRA_MODEL_EVALUATOR_HPP
