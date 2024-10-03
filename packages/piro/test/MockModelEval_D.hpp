// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MOCKMODELEVAL_D_H
#define MOCKMODELEVAL_D_H

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
typedef int MPI_Comm;
#include "Epetra_SerialComm.h"
#endif 

#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"

#include "EpetraExt_ModelEvaluator.h"


/** \brief Simple model evaluator for testing SG response derivatives
 */
class MockModelEval_D : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_D(const Teuchos::RCP<const Epetra_Comm>& comm);

  //@}

  ~MockModelEval_D();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
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
  
  //@}

private:

  //These are set in the constructor and used in evalModel
  Teuchos::RCP<const Epetra_Comm> comm;
  Teuchos::RCP<Epetra_Map> x_map;
  Teuchos::RCP<Epetra_Vector> x_init;
  Teuchos::RCP<Epetra_Map> p1_map;
  Teuchos::RCP<Epetra_Map> p2_map;
  Teuchos::RCP<Epetra_Vector> p1_init;
  Teuchos::RCP<Epetra_Vector> p2_init;
  Teuchos::RCP<Epetra_Map> g_map;
  Teuchos::RCP<Epetra_CrsGraph> graph;
  
};

#endif // MOCKMODELEVAL_C_H
