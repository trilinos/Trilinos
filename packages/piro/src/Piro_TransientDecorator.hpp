// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TRANSIENTDECORATOR_H
#define PIRO_TRANSIENTDECORATOR_H

#include <iostream>

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace Piro {

template <typename Scalar>
class TransientDecorator
    : public Thyra::ModelEvaluatorDefaultBase<Scalar> {

  public:

  /** \name Constructors/initializers */
  //@{

  TransientDecorator(){}

  //@}

  Teuchos::RCP<Thyra::VectorBase<Scalar> > get_x_dotdot() const { return xDotDot; }

  void set_x_dotdot(const Teuchos::RCP<Thyra::VectorBase<Scalar> > &xdd) const {  xDotDot = xdd; }
  void set_x_dotdot_data(const Teuchos::RCP<Thyra::VectorBase<Scalar> > &xdd) const {  assign(xDotDot.ptr(), *xdd); }

  Scalar get_omega() const { return omega; };

  void set_omega(const Scalar oo) const {omega = oo; }

  void printVec(const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &vec) const {
       std::cout << "GAH Vector is: " << std::endl;
       for(int ii=0; ii < vec->space()->dim(); ii++)
          std::cout << get_ele(*vec, ii) << std::endl;
  }

  Thyra::ModelEvaluatorBase::Derivative<Scalar> get_DgDx_dotdot(int j) const { return DgDx_dotdot_[j]; }

  /** \brief Precondition: <tt>supports(OUT_ARG_DgDx_dot,j)==true</tt>.  */
  void set_DgDx_dotdot(int j, const Thyra::ModelEvaluatorBase::Derivative<Scalar> &DgDx_dotdot_j){
          DgDx_dotdot_[j] = DgDx_dotdot_j;
       }

  protected:

   mutable Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot;
   mutable Scalar omega;
   mutable Teuchos::Array<Thyra::ModelEvaluatorBase::Derivative<Scalar> > DgDx_dotdot_;

};


}

#endif
