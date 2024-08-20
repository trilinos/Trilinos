// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_INVERSEMASSMATRIXDECORATOR_H
#define PIRO_EPETRA_INVERSEMASSMATRIXDECORATOR_H

#include <iostream>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"

#include "EpetraExt_ModelEvaluator.h"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"

namespace Piro {
namespace Epetra {

class InvertMassMatrixDecorator
    : public EpetraExt::ModelEvaluator
{

  typedef double Scalar;

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  InvertMassMatrixDecorator( 
                Teuchos::RCP<Teuchos::ParameterList> stratParams,
                Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
                bool massMatrixIsConstant=true,
                bool lumpMassMatrix=false,
                bool massMatrixIsCoeffOfSecondDeriv=false
                );

  //@}

  ~InvertMassMatrixDecorator();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<Epetra_Vector> x_dot;

   //Teuchos::RCP<Epetra_CrsMatrix> massMatrix;
   Teuchos::RCP<Epetra_Operator> massMatrix;
   Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;

   bool massMatrixIsConstant; // User Setting
   bool lumpMassMatrix; // User Setting to rowSum Matrix
   bool massMatrixIsCoeffOfSecondDeriv; // Set to true for x_dotdot acceleration problems
   Teuchos::RCP<Epetra_Vector> invDiag;

   // The following get modified in evalModel and so are mutable
   mutable Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
   mutable bool calcMassMatrix; //Internal flag
   mutable Teuchos::RCP<const Thyra::LinearOpBase<double> > A;
};
}
}
#endif
