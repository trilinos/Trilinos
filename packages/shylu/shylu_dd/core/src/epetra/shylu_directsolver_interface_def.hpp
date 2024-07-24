// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_directsolver_interface_def.hpp

    \brief Eperta/Tpetra templated interface for call Amesos and Amesos2

    \author Joshua Dennis Booth
*/

#ifndef SHYLU_DIRECTSOLVER_INTERFACE_DEF_HPP
#define SHYLU_DIRECTSOLVER_INTERFACE_DEF_HPP

#include "ShyLU_DDCore_config.h"
#include "shylu_directsolver_interface_decl.hpp"

#include <Teuchos_XMLParameterListHelpers.hpp>


namespace ShyLU{

  template<class Matrix, class Vector>
  DirectSolverInterface<Matrix,Vector>::DirectSolverInterface()
  {}

  template<class Matrix, class Vector>
  DirectSolverInterface<Matrix, Vector>::DirectSolverInterface(Matrix *inA, Teuchos::ParameterList *inpList)
  {
    A=inA;
    pList = inpList;
  }

  template<class Matrix, class Vector>
  DirectSolverInterface<Matrix, Vector>::~DirectSolverInterface()
  { }

  template <class Matrix, class Vector>
  int
  DirectSolverInterface<Matrix, Vector>::init_matrix(Matrix *inA, Teuchos::ParameterList *inpList)
  {
    A=inA;
    pList = inpList;
    return 0;
  }




  template <class Matrix, class Vector>
  int
  DirectSolverInterface<Matrix,Vector>::factorAmesos()
  {
    std::cout << "**Error**: Amesos is only supported for Epetra Matrices \n";
    exit(1);
  }


  //Note: come back and update
  template <>
  int
  DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::factorAmesos()
  {

    Teuchos::ParameterList subList = pList->sublist("Amesos Input");
    std::string solvertype = Teuchos::getParameter<std::string>(subList, "Solver");
    Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");

    problem_amesos.SetOperator(A);
    //future reference
    //setLHS(x)
    //setRHS(B)

    Amesos Factory;
    solver_amesos  = Factory.Create(solvertype, problem_amesos);
    assert (solver_amesos != 0);
    solver_amesos->SetParameters(subsubList);
    //Add error checking
    solver_amesos->SymbolicFactorization();
    solver_amesos->NumericFactorization();

    return 0;
  }//end factorAmesos

#ifdef HAVE_SHYLU_DDCORE_AMESOS2
  template <class Matrix, class Vector>
  int DirectSolverInterface<Matrix,Vector>::factorAmesos2()
  {
    //#pragma message("solve amesos2 compiled")
    //std::cout << "odd call";
    Teuchos::ParameterList subList = pList->sublist("Amesos2 Input");
    std::string solvertype = Teuchos::getParameter<std::string>(subList, "Solver");
    Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");
    solver_amesos2 = Amesos2::create<Matrix, Vector>
      (solvertype, Teuchos::rcp(A,false) ); //right now only use default
    solver_amesos2->symbolicFactorization().numericFactorization();

    return 0;
  }//end factor_amesos2
#endif


   template <class Matrix, class Vector>
  int DirectSolverInterface<Matrix,Vector>::factor()
  {
    #ifdef HAVE_SHYLU_DDCORE_AMESOS2
    return factorAmesos2();
    #else
    return 1;
    #endif
  }//end factor()


  template < >
  int
  DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::factor()
  {
    std::string solverpackage = Teuchos::getParameter<std::string>(*pList,"Direct Solver Package");
    int returnvalue = 1;
    if(solverpackage.compare("Amesos")==0)
      {
        returnvalue = factorAmesos();
      }
    else if(solverpackage.compare("Amesos2")==0)
      {
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
        returnvalue = factorAmesos2();
#else
        std::cout << "Amesos2 is not installed \n";
        exit(1);
#endif
      }
    else
      {
        std::cout << "No Direct Solver Package Found";
      }
    return returnvalue;
  }//end factor<epetra>






template <class Matrix, class Vector>
int
DirectSolverInterface<Matrix, Vector>::solveAmesos(Vector* b, Vector *x)
{
  std::cout << "**Error**: Amesos is only supported for Epetra Matrices \n";
  exit(1);
}
template <>
int
DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::solveAmesos(Epetra_MultiVector *b, Epetra_MultiVector *x)
{

  Teuchos::ParameterList subList = pList->sublist("Amesos Input");
  std::string solvertype = Teuchos::getParameter<std::string>(subList, "Solver");
  Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");

  //Epetra_LinearProblem Problem(A, x, b);
  //Amesos Factory;
  // Amesos_BaseSolver* Solver = Factory.Create(solvertype, Problem);
  //assert (Solver != 0);
  problem_amesos.SetRHS(x);
  problem_amesos.SetLHS(b);

  //solver_amesos->SetParameters(subsubList);
  //Add error checking
  //Solver->SymbolicFactorization();
  //Solver->NumericFactorization();
  solver_amesos->Solve();

  return 0;
}
template <class Matrix, class Vector>
int
DirectSolverInterface<Matrix, Vector>::solve(Vector* b, Vector* x)
{
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
  return solveAmesos2(b, x);
#else
  return 1;
#endif
}

template < >
int
DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::solve(Epetra_MultiVector*b ,Epetra_MultiVector* x)
{
  std::string solverpackage = Teuchos::getParameter<std::string>(*pList,"Direct Solver Package");
  int returnvalue = 1;
  if(solverpackage.compare("Amesos")==0)
    {
      returnvalue = solveAmesos(b,x);
    }
  else if(solverpackage.compare("Amesos2")==0)
    {
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
        returnvalue = solveAmesos2(b,x);
#else
        std::cout << "Amesos2 is not installed \n";
        exit(1);
#endif
    }
    else
      {
        std::cout << "No Direct Solver Package Found";
      }
    return returnvalue;
}
#ifdef HAVE_SHYLU_DDCORE_AMESOS2
template <class Matrix, class Vector>
int DirectSolverInterface<Matrix,Vector>::solveAmesos2(Vector* b, Vector* x)
{
  //#pragma message("solve amesos2 compiled")
  //std::cout << "odd call";
  Teuchos::ParameterList subList = pList->sublist("Amesos2 Input");
  std::string solvertype = Teuchos::getParameter<std::string>(subList, "Solver");
  Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");
  solver_amesos2->solve(x, b);

  return 0;
}
#endif




}// end namespace

#endif // end header if

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

