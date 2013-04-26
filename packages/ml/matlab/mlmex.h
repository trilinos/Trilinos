/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#ifndef MLMEX_H
#define MLMEX_H

#include "ml_config.h"
#include "ml_common.h"

#ifdef HAVE_ML_MLAPI
/* MLAPI Headers */
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"
#endif

/* ML_Epetra Headers */
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

#include "ml_RefMaxwell.h"

#ifdef HAVE_ML_MATLAB
#include "mex.h"

/* The big classes of data */
class ml_data_pack;
class ml_data_pack{
public:
  ml_data_pack();
  virtual ~ml_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  virtual int setup(int N,int* rowind,int* colptr, double* vals)=0;

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  virtual int status()=0;

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  virtual int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x, int &iters)=0;  
public:
  int id;
  Teuchos::ParameterList *List;
  double operator_complexity;  
  ml_data_pack *next;
};

class mlapi_data_pack:public ml_data_pack{
public:
  mlapi_data_pack();
  ~mlapi_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  int status();


  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O] (NOT IMPLEMENTED)
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x, int &iters);  

private:
  MLAPI::Space * FineSpace;
  MLAPI::DistributedMatrix * A;
  MLAPI::MultiLevelAdaptiveSA *Prec;
};


class ml_epetra_data_pack:public ml_data_pack{
public:
  ml_epetra_data_pack();
  ~ml_epetra_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  int status();

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x,int &iters);  

  /* GetPreconditioner - returns a pointer to the preconditioner */     
  ML_Epetra::MultiLevelPreconditioner* GetPreconditioner(){return Prec;}
  
  
private:
  Epetra_CrsMatrix * A;
  ML_Epetra::MultiLevelPreconditioner *Prec;
};


class ml_maxwell_data_pack:public ml_data_pack{
public:
  ml_maxwell_data_pack();
  ~ml_maxwell_data_pack();
  
  /* setup - sets up an ml_data_pack object
     Parameters:
     name    - Name of matrix to create [I]
     mxa     - mxArray containing data
     rewrap_ints - if the ints don't line up correctly
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup_matrix(const char * name, const mxArray * mxa, bool rewrap_ints);
  /*
    Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup_preconditioner();

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  int status();

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x,int &iters);  

  /* GetPreconditioner - returns a pointer to the preconditioner */     
  ML_Epetra::RefMaxwellPreconditioner* GetPreconditioner(){return Prec;}
  
  
private:
  // Disabled
  int setup(int N,int* rowind,int* colptr, double* vals){return 0;}


  Epetra_CrsMatrix *EdgeMatrix, *GradMatrix, *NodeMatrix;
  //  ML_Epetra::MultiLevelPreconditioner *Prec;
  ML_Epetra::RefMaxwellPreconditioner *Prec;

};

#endif
#endif
