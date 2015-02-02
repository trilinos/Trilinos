#ifndef SHYLU_DIRECTSOLVER_INTERFACE_DEF_HPP
#define SHYLU_DIRECTSOLVER_INTERFACE_DEF_HPP

#include "ShyLUCore_config.h"
#include "shylu_directsolver_interface_decl.hpp"

#include <Teuchos_XMLParameterListHelpers.hpp>


namespace ShyLU{

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
DirectSolverInterface<Matrix, Vector>::solveAmesos(Vector* b, Vector *x)
{
  cout << "**Error**: Amesos is only supported for Epetra Matrices \n";
  exit(1);
}
template <>
int
DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::solveAmesos(Epetra_MultiVector *b, Epetra_MultiVector *x)
{
  
  Teuchos::ParameterList subList = pList->sublist("Amesos Input");
  string solvertype = Teuchos::getParameter<string>(subList, "Solver");
  Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");
 
  Epetra_LinearProblem Problem(A, x, b);
  Amesos Factory;
  Amesos_BaseSolver* Solver = Factory.Create(solvertype, Problem);
  assert (Solver != 0);
  Solver->SetParameters(subsubList);
  //Add error checking
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();
  delete Solver;

  return 0;
}
template <class Matrix, class Vector>
int
DirectSolverInterface<Matrix, Vector>::solve(Vector* b, Vector* x)
{
#ifdef HAVE_SHYLU_AMESOS2 
  return solveAmesos2(b, x);
#else
  return 1;
#endif
}

template < >
int 
DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector>::solve(Epetra_MultiVector*b ,Epetra_MultiVector* x)
{
  string solverpackage = Teuchos::getParameter<string>(*pList,"Direct Solver Package");
  int returnvalue = 1;
  if(solverpackage.compare("Amesos")==0)
    {
      returnvalue = solveAmesos(b,x);
    }
  else if(solverpackage.compare("Amesos2")==0)
    {
#ifdef HAVE_SHYLU_AMESOS2
      	returnvalue = solveAmesos2(b,x);
#else
	cout << "Amesos2 is not installed \n";
	exit(1);
#endif
    }
    else
      {
	cout << "No Direct Solver Package Found";
      }
    return returnvalue;
}
#ifdef HAVE_SHYLU_AMESOS2
template <class Matrix, class Vector>
int DirectSolverInterface<Matrix,Vector>::solveAmesos2(Vector* b, Vector* x)
{

  //cout << "odd call";
  Teuchos::ParameterList subList = pList->sublist("Amesos2 Input");
  string solvertype = Teuchos::getParameter<string>(subList, "Solver");
  Teuchos::ParameterList subsubList = subList.sublist(solvertype + " Input");
  Teuchos::RCP<Amesos2::Solver<Matrix, Vector> > solver;
  solver = Amesos2::create<Matrix, Vector>
    (solvertype, Teuchos::rcp(A,false), Teuchos::rcp(x,false), Teuchos::rcp(b,false));

  solver->symbolicFactorization().numericFactorization().solve();

  return 0;
}
#endif

}// end namespace

#endif // end header if
