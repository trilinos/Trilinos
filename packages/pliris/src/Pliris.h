
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _PLIRIS_H_
#define _PLIRIS_H_


class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"

//! Pliris: An Obect-Oriented Interface to a Dense LU Solver.
/*! The Pliris class : Provides the functionality to interface to a dense LU
*/  


class Pliris {

 public:

  Pliris(Epetra_Vector * A, Epetra_MultiVector * X, Epetra_MultiVector * B);

  //! Pliris Default constructor.

  Pliris ();


  //! Pliris LHS Set
  /*! Associates an already defined Epetra_MultiVector (or Epetra_Vector) as the initial guess
      and location where the solution will be return.
   */
  int SetLHS(Epetra_MultiVector * X);


  //! Pliris RHS Set
  /*! Associates an already defined Epetra_MultiVector (or Epetra_Vector) as the right-hand-side of 
      the linear system.
   */

  int SetRHS(Epetra_MultiVector * B);

  //! Pliris Matrix Set
  /*! Associates an already defined  Epetra_Vector as the matrix (column ordered) of 
      the linear system.
   */
  int SetMatrix(Epetra_Vector * A);
 
  //! Pliris Matrix Set
  /*! Associates an already defined  Epetra_SerialDenseVector as the matrix (column ordered) of
      the linear system.
   */
  int SetMatrix(Epetra_SerialDenseVector * A); 

  //! Pliris GetDistirbution
  /*! Gives the distribution information that is required by the dense solver
   */

    /*! 
     \param nprocs_row(In)   - number of processors for a row
     \param number_of_unknowns(In)  - order of the dense matrix
     \param nrhs(In)  - number of right hand sides
     \param my_rows(Out)  - number of rows of the matrix on this processor
     \param my_cols (Out)  - number of columns of the matrix on this processor
     \param  my_first_row(Out)  - first (global) row number on this processor (array starts at index 1)
     \param  my_first_col (Out)  - first (global) column number on this processor (array starts at index 1)
     \param  my_rhs(Out)  - number of Right hand sides on this processor
     \param my_row(Out)  - row number in processor mesh, 0 to the  number of processors for a column -1
     \param my_col(Out)  - column  number in processor mesh, 0 to the  number of processors for a row -1
    */

  int GetDistribution( int * nprocs_row, 
                             int * number_of_unknowns,
                             int * nrhs,
                             int * my_rows,
                             int * my_cols,
                             int * my_first_row,
                             int * my_first_col,
                             int * my_rhs,
                             int * my_row,
                             int * my_col );

  //! Pliris FactorSolve
  /*! Factors and solves the dense matrix
   */

    /*! 
     \param A(In) --  Epetra Vector that has the matrix and rhs packed
     \param my_rows(In) --  number of rows of the matrix on this processor
     \param my_cols(In) --  number of columns of the matrix on this processor  
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param num_rhs(In) --  number of right hand sides 
     \param secs(Out) -- factor and solve time in seconds 
    */

   int FactorSolve( Epetra_Vector * A,
                    int my_rows,
                    int my_cols,
                    int* matrix_size,
                    int* num_procsr,
                    int* num_rhs,
                    double* secs);

  //! Pliris FactorSolve
  /*! Factors and solves the dense matrix
   */

    /*! 
     \param AA(In) --  Epetra Serial Dense Vector that has the matrix and rhs packed
     \param my_rows(In) --  number of rows of the matrix on this processor
     \param my_cols(In) --  number of columns of the matrix on this processor 
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param num_rhs(In) --  number of right hand sides 
     \param secs(Out) -- factor and solve time in seconds 
    */

   int FactorSolve( Epetra_SerialDenseVector * AA,
                    int my_rows,
                    int my_cols,
                    int* matrix_size,
                    int* num_procsr,
                    int* num_rhs,
                    double* secs);

  //! Pliris Factor
  /*! Factors the dense matrix
   */

    /*! 
     \param A(In) --  Epetra Vector that has the matrix packed
     \param matrix_size(In) -- order of the dense matrix
     \param num_procsr(In) --  number of processors for a row
     \param permute(In) -- permutation matrix
     \param secs(Out) -- factor and solve time in seconds 
    */


   int Factor( Epetra_Vector* A,
               int* matrix_size,
               int* num_procsr,
               int* permute,
	       double* secs);


  //! Pliris Default Destructor.
  /*! Completely deletes a Pliris object.  
  */
  virtual ~Pliris(void); 
  //@


  protected:

  double *x_;
  double *a_;
  double *b_;
  int x_LDA_;
  int b_LDA_;
  
    bool inConstructor_;
    Epetra_MultiVector * X_;
    Epetra_MultiVector * B_;
    Epetra_Vector * A_;
    Epetra_SerialDenseVector * AA_;

};



#endif /* _PLIRIS_H_ */
