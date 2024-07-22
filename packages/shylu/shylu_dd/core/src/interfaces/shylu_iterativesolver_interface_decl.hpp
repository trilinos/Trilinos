// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLU_ITERATIVESOLVER_INTERFACE_DECL_HPP
#define SHYLU_ITERATIVESOLVER_INTERFACE_DECL_HPP

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
#include <BelosTpetraAdapter.hpp>
#endif

//#include <Zoltan2_config.h>
#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE)
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#include <Teuchos_XMLParameterListHelpers.hpp>

// Aztec
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

// Belos //May have to go back and add more
//#ifdef HAVE_SHYLU_DDCORE_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
//#endif

namespace ShyLU
{

  template <class Matrix, class Vector>
  class IterativeSolverInterface
  {
    //scalar type
    /*
#ifdef HAVE_SHYLU_DDCORE_TPETRA
    typedef typename Matrix::scalar_type           scalar_type;   
    typedef typename Matrix::local_ordinal_type    local_ordinal_type;
    typedef typename Matrix::global_ordinal_type   global_ordinal_type;
    typedef typename Matrix::node_type             node_type;
    typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>                                     OP;
#else
    //Nded to get the scaler type for the matrix
    typedef double                                 scalar_type; 
    typedef Teuchos::ScalarTraits<scalar_type>     scalar_type_traits;
    typedef scalar_type_traits::magnitudeType      magnitude_type;
    typedef Epetra_Operator                        OP;
#endif
    */


  public:
    ~IterativeSolverInterface();
    
    //doxygen
    IterativeSolverInterface();
    IterativeSolverInterface(Matrix *inA, Teuchos::ParameterList *pList);
    int init_matrix(Matrix *inA, Teuchos::ParameterList* pList);
    int solve(Vector *b, Vector *x);

    ///....

  private:
    //private functions ...
    int solveAztec(Vector *b, Vector *x);
    int solveBelos(Vector *b, Vector *x);


    //private variables ..
    Teuchos::ParameterList *pList;
    Matrix *A;
    int maxproc;


    Epetra_LinearProblem e_problem;
    AztecOO              *solver_aztec; 


    //add ifdefs around this
    //Belos::LinearProblem<scalar_type,Vector,OP> belos_problem;
    //Belos::SolverManager<scalar_type,Vector,OP> *solver_belos;
  
  };

}//end namespace ShyLU

#endif

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

