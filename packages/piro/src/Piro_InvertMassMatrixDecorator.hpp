// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PIRO_INVERTMASSMATRIXDECORATOR_H
#define PIRO_INVERTMASSMATRIXDECORATOR_H

#include <iostream>

#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Piro_config.hpp"


namespace Piro {

template<typename Scalar>
class InvertMassMatrixDecorator
    : public Thyra::ModelEvaluatorDefaultBase<Scalar>
{

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  InvertMassMatrixDecorator( 
                Teuchos::RCP<Teuchos::ParameterList> stratParams,
                Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
                bool massMatrixIsConstant=true,
                bool lumpMassMatrix=false,
                bool massMatrixIsCoeffOfSecondDeriv=false
                );

  //@}

  ~InvertMassMatrixDecorator();


  /** \name Overridden from Thyra::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP< Thyra::LinearOpBase< Scalar > > create_W_op () const;
  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > create_W_prec() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgsImpl() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
 
  void reportFinalPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint, const bool wasSolved);  

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;

  Teuchos::RCP<const Teuchos::Array<std::string> >
    get_p_names(int l) const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;


  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot;

   Teuchos::RCP<Thyra::LinearOpBase<Scalar> > massMatrix;
   Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;

   bool massMatrixIsConstant; // User Setting
   bool lumpMassMatrix; // User Setting to rowSum Matrix
   bool massMatrixIsCoeffOfSecondDeriv; // Set to true for x_dotdot acceleration problems
   Teuchos::RCP<Thyra::VectorBase<Scalar> > invDiag;

   // The following get modified in evalModel and so are mutable
   mutable Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
   mutable bool calcMassMatrix; //Internal flag
   mutable Teuchos::RCP<const Thyra::LinearOpBase<double> > A;
};
}
/** \class Piro::RythmosSolver
 *  \ingroup Piro_Thyra_solver_grp
 * */

#include "Piro_InvertMassMatrixDecorator_Def.hpp"

#endif
