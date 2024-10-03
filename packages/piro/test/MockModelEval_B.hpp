// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MOCKMODELEVAL_B_H
#define MOCKMODELEVAL_B_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
typedef int MPI_Comm;
#include "Epetra_SerialComm.h"
#endif 

#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_ModelEvaluator.h"


/** \brief Epetra-based Model Evaluator subclass for Charon!
 *
 * This class will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Aristos to
 * solve different types of problems with Charon.
 */

class MockModelEval_B
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_B(const MPI_Comm appComm);

  //@}

  ~MockModelEval_B();


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
  /** \brief Function that does regression testing. */

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
  Teuchos::RCP<const Epetra_Vector> get_x_dotdot_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Epetra_Map> x_map;
   Teuchos::RCP<Epetra_Vector> x_vec;
   Teuchos::RCP<Epetra_Vector> x_dot_vec;
   Teuchos::RCP<Epetra_Vector> x_dotdot_vec;
   Teuchos::RCP<Epetra_Map> p_map;
   Teuchos::RCP<Epetra_Vector> p_init;
   Teuchos::RCP<Epetra_Map> g_map;
   Teuchos::RCP<Epetra_Comm> Comm;
   Teuchos::RCP<Epetra_CrsMatrix> W;

};

#endif // SIMPLE_MODELEVAL_H
