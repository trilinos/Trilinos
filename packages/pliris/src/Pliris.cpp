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

#include "Pliris.h"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "distribute.h"
#include "xlu_solve.h"
#include "permute.h"
#include "x_factor.h"



//=============================================================================

Pliris::Pliris(Epetra_Vector * A, 
		Epetra_MultiVector * X,
	        Epetra_MultiVector * B)  {



  inConstructor_ = true;

  SetMatrix(A);

  SetLHS(X);

  SetRHS(B);

  inConstructor_ = false;
}

Pliris::Pliris( ) {

}


//=============================================================================
Pliris::~Pliris(void) {
 

}

//=============================================================================
 int Pliris::GetDistribution( int*  nprocs_row, 
                               int* number_of_unknowns,
			       int* nrhs,
                               int* my_rows,
                               int* my_cols,
                               int* my_first_row,
                               int* my_first_col,
                               int* my_rhs,
                               int* my_row,
                               int* my_col ) {

 

  // This function echoes the multiprocessor distribution of the matrix
 
  distmat_ ( nprocs_row,
             number_of_unknowns,
             nrhs,
             my_rows,
             my_cols,
             my_first_row,
             my_first_col,
             my_rhs,
             my_row,
             my_col);


  return(0);

} 
 
//=============================================================================
int Pliris::FactorSolve( Epetra_Vector* A,
                         int my_rows,
                         int my_cols,
                         int* matrix_size,
                         int* num_procsr,
                         int* num_rhs,
                         double* secs){

    

      SetMatrix(A);
    

      dlusolve_(a_, 
                matrix_size, 
                num_procsr, 
                (a_ + my_rows*my_cols),
                num_rhs,
                secs);

     return(0);

}

//=============================================================================
int Pliris::FactorSolve( Epetra_SerialDenseVector* A,
                         int my_rows,
                         int my_cols,
                         int* matrix_size,
                         int* num_procsr,
                         int* num_rhs,
                         double* secs){

     // Associate the matrix with a_

     SetMatrix(A);
    
     //  Factor and solve the matrix equation

     dlusolve_  (a_, 
                 matrix_size, 
                 num_procsr, 
                 (a_ + my_rows*my_cols),
                 num_rhs,
                 secs);

     return(0);

}

//  Factor the matrix for later solve

//=============================================================================

int Pliris::Factor( Epetra_Vector* A,
                     int* matrix_size,
                     int* num_procsr,
                     int* permute,
                     double* secs){


     SetMatrix(A);

     // Factor the Matrix

     dfactor_(a_,
              matrix_size,
              num_procsr,
              permute,
              secs);

      
     // Permute the lower triangular matrix

     dpermute_(a_, permute); 

     return(0);
} 




//=============================================================================
int Pliris::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int Pliris::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}
//=============================================================================
int Pliris::SetMatrix(Epetra_Vector * A) {

  if (A == 0 && inConstructor_ == true) return(0);
  if (A == 0) EPETRA_CHK_ERR(-1);
  A_ = A;
  A_->ExtractView(&a_);

 return(0);
}

//=============================================================================
int Pliris::SetMatrix(Epetra_SerialDenseVector * AA_) {

  a_ = AA_-> Values();

 return(0);

}  
