#ifndef PIRO_EPETRA_INVERSEMASSMATRIXDECORATOR_H
#define PIRO_EPETRA_INVERSEMASSMATRIXDECORATOR_H

#include <iostream>

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"

#include "EpetraExt_ModelEvaluator.h"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

/** \brief  Decorator class that creates a W (Jacobian) operator
 *          using Matrix-Free directional derivatives. 
 *
 * This class takes a model evaluator that supports a residual
 * calculation, and adds support for the W matrix as an
 * epetra_operator.
 *
 * This class is the ModelEvaluator version of NOX_Epetra_MatrixFree.
 * One difference, is that it supports time-dependent problems
 * when x_dot != null, and using alpha and beta. 
 * 
 */

namespace Piro {
namespace Epetra {

class MatrixFreeDecorator
    : public EpetraExt::ModelEvaluator,
      public Epetra_Operator
{

  typedef double Scalar;

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MatrixFreeDecorator( 
                Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
                double lambda_ = 1.0e-6
                );

  //@}

  ~MatrixFreeDecorator();


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
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

  //@}
  /** \name Overridden from Epetra_Operator . */
  //@{
  
  virtual int SetUseTranspose(bool UseTranspose);
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  virtual double NormInf() const;
  virtual const char* Label () const;
  virtual bool UseTranspose() const;
  virtual bool HasNormInf() const;
  virtual const Epetra_Comm & Comm() const;
  virtual const Epetra_Map& OperatorDomainMap () const;
  virtual const Epetra_Map& OperatorRangeMap () const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   // Storage for Base resid vector for MatrixFree 
   Teuchos::RCP<Epetra_Vector> fBase;
   // Storage for perturbed resid vector for MatrixFree 
   Teuchos::RCP<Epetra_Vector> fPert;

   // Constant used in formulas to pick perturbation, typically 1.0e-6
   double lambda;

};
}
}
#endif
