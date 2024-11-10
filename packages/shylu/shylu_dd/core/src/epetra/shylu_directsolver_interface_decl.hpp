// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_directsolver_interface_decl.hpp

    \brief Eperta/Tpetra templated interface for call Amesos and Amesos2

    \author Joshua Dennis Booth
*/
#ifndef SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP
#define SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP

#include "ShyLU_DDCore_config.h"

//ShyLU
#include <shylu.h>
#include <shylu_util.h>

//Epetra
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>

//Isorropia
#include <Isorropia_config.h>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraProber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

//Tperta
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

//#include <Zoltan2_config.h>
#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE)
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#include <Teuchos_XMLParameterListHelpers.hpp>

//Amesos
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#ifdef HAVE_SHYLU_DDCORE_AMESOS2
#include <Amesos2.hpp>
#endif

namespace ShyLU{

  /** \brief DirectSolverInterface class templated on Epetra/Tpetra Matrix and Vector
   *
   * This class acts as an interface that will allow Shylu to call either Amesos/Amesos2 
   * without having to address if the matrix/submatrix is either of Epetra or Tpetra form.
   * Currently:  Only supporting limited solvers in both Amesos and Amesos2, will be updated
   */
template <class Matrix, class Vector>
class DirectSolverInterface
{
public:
  ~DirectSolverInterface();

  /** \brief Main constructor of class
   *
   * This constructor requires a Teuchos ParameterList that provides information on solver.
   * It assumes that if Tpetra matrix is given then Amesos2 ParameterList must be given.
   * Likewise, if Epetra matrix is given then Amesos ParameterList must be given.
   */

  DirectSolverInterface();
  DirectSolverInterface(Matrix *inA, Teuchos::ParameterList* pList);
  int init_matrix(Matrix *inA, Teuchos::ParameterList* pList);
  int factor();
  int solve(Vector* b, Vector* x);
private:
  //Private Functions
  int factorAmesos();
  int solveAmesos(Vector* b, Vector *x );
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
  int factorAmesos2();
  int solveAmesos2(Vector* b, Vector *x);
#endif

  //Private Variable
  Teuchos::ParameterList *pList;
  Matrix *A;
  int maxproc;
  
  //amesos
  Epetra_LinearProblem  problem_amesos;
  Amesos_BaseSolver*    solver_amesos;

  //amesos2
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
  Teuchos::RCP<Amesos2::Solver<Matrix,Vector> > solver_amesos2;
#endif


}; //end class

}// end namespace
#endif //end header if

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

