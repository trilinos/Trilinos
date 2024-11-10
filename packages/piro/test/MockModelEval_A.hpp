// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MOCKMODELEVAL_A_H
#define MOCKMODELEVAL_A_H

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


/** \brief Epetra-based Model Evaluator
 *
 * Concrete model evaluator for the solution of the following PDE-Constrained problem:
 *
 * find (p_0,p_1) that minimizes
 * g = 0.5*(Sum(x)-Sum(p)-12)^2 + 0.5*(p0-1)^2
 * subject to:
 * f_0 = (x_0)^2 - p_0 = 0
 * f_i = x_i^2 - (i+p_1)^2 (for i != 0), for i = 1,2,3,4
 *
 * solution is p = (1,3).
 */

class MockModelEval_A
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  MockModelEval_A(const MPI_Comm appComm);

  //@}

  ~MockModelEval_A();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_lower_bounds(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_upper_bounds(int l) const;
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
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Epetra_Map> x_map;
   Teuchos::RCP<Epetra_Vector> x_vec;
   Teuchos::RCP<Epetra_Vector> x_dot_vec;
   Teuchos::RCP<Epetra_Map> p_map;
   Teuchos::RCP<Epetra_Vector> p_init;
   Teuchos::RCP<Epetra_Vector> p_lo;
   Teuchos::RCP<Epetra_Vector> p_up;
   Teuchos::RCP<Epetra_Map> g_map;
   Teuchos::RCP<Epetra_Comm> Comm;
   Teuchos::RCP<Epetra_CrsGraph> jacGraph;

};

#endif // SIMPLE_MODELEVAL_H
