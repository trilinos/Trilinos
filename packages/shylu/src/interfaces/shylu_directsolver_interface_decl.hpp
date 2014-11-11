#ifndef SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP
#define SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP

#include "ShyLU_config.h"

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
#ifdef HAVE_SHYLU_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

//#include <Zoltan2_config.h>
#ifdef HAVE_SHYLU_ZOLTAN2
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#include <Teuchos_XMLParameterListHelpers.hpp>

//Amesos
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#ifdef HAVE_SHYLU_AMESOS2
#include <Amesos2.hpp>
#endif

namespace ShyLU{

template <class Matrix, class Vector>
class DirectSolverInterface
{
public:
  ~DirectSolverInterface();
  DirectSolverInterface(Matrix *inA, Teuchos::ParameterList* pList);
  int solve(Vector* b, Vector* x);
private:
  int solveAmesos(Vector* b, Vector *x );
  Teuchos::ParameterList *pList;
  Matrix *A;
  int maxproc;
#ifdef HAVE_SHYLU_AMESOS2
  int solveAmesos2(Vector* b, Vector *x);
#endif
}; //end class

}// end namespace

#endif //end header if
