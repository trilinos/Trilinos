/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/



#include <iostream>
#include "Pliris.h"

#include "Epetra_MpiComm.h"
#include "mpi.h"

//  #include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_SerialDenseVector.h"


int main(int argc, char *argv[])
{


  int  my_rows;
  int  my_cols;
  int  my_first_row;
  int  my_first_col;
  int  my_rhs;
  int  my_row;
  int  my_col;
  int matrix_size;
  int nprocs_row;

  double mflops;

  MPI_Comm rowcomm;


  static int buf[4];
  int nrhs;

  int i, m;
  int k;
  int l;

  int mlen;   // Message length for input data

  int ierror;

  int num_global_length;

  unsigned int seed= 10;

  int num_my_length;

  double secs;

  double eps;

  double othird;

  double four_thirds = 4./3.;

  double tempc;

  double solinf;

  double ainf;

  double axbinf;

  double xh;

  double one =1.;


  // Enroll into MPI

  MPI_Init(&argc,&argv);

  Epetra_MpiComm comm(MPI_COMM_WORLD);


  // Pliris example using the Epetra_Vector interface


  // Initialize Input buffer

  for(i=0;i<4;i++) buf[i]=-1;

  std::cout << "proc " << comm.MyPID() << " is alive of   " << comm.NumProc()<< " Processors" << std::endl ;

  if( comm.MyPID() == 0 ) {


    // Check for commandline input

     if (argc > 1) {
       // argv[1] should be size of matrix
       buf[0] = atoi(argv[1]);
       if (argc > 2) {
          // argv[2] should be #procs per row
          buf[1] = atoi(argv[2]);
       }
       else
          // default is 1, but sqrt(p) would be better
          buf[1] = 1;
     }
     else {

    // Input Data about matrix and distribution

   	if (buf[0] < 0) {
	  std::cout << "Enter size of matrix " << std::endl;
	  std::cin >> buf[0];
	}
	if (buf[1] < 0) {
	  std::cout << "Enter number of processors to which each row is assigned "  << std::endl;
	  std::cin >> buf[1];
	}

     }

	std::cout << " Matrix Size "<< buf[0] <<"\n";

        std::cout << " Processors in a row  " <<  buf[1] << "\n";

  }

   mlen = 4* sizeof(int);

    /* Send the initilization data to each processor    */

    // Using Epetra Communicator

    comm.Broadcast(buf,mlen,0);


     // Set the values where needed

   matrix_size = buf[0];

   nprocs_row = buf[1];

   // Example for 1 RHS

   nrhs = 1;


   if( comm.MyPID() == 0) {
    std::cout << " ---- Building Pliris solver ----" << std::endl;
   }

   // Instantiate the Solver


   Pliris solver;

   // Get Info to build the matrix on a processor

   solver.GetDistribution( &nprocs_row,
                            &matrix_size,
			    &nrhs,
                            &my_rows,
                            &my_cols,
                            &my_first_row,
                            &my_first_col,
                            &my_rhs,
                            &my_row,
			    &my_col );


   //   Define a new communicator

   MPI_Comm_split(MPI_COMM_WORLD,my_row,my_col,&rowcomm);

   //if( comm.MyPID() == 0 ){
   std::cout << " ------ PARALLEL Distribution Info for : ---------" <<std::endl;

   std::cout << "   Processor  " << comm.MyPID() << std::endl
        << "    my rows  " << my_rows << std::endl
        << "    my cols  " << my_cols << std::endl
        << "    my rhs  " << my_rhs << std::endl
        << "    my first col  " << my_first_col  << std::endl
        << "    my first row  " << my_first_row << std::endl
        << "    my_row  " << my_row << std::endl
        << "    num procs row   " << nprocs_row << std::endl
        << "    my_col  " << my_col << std::endl;

   //}

   //  Local size -- my_rows  * (my_cols + my_rhs)


   // Some temp arrays

     double* temp = new double[my_rows];

     double* temp2 = new double[my_rows];

     double* rhs = new double[matrix_size];

     double* temp4 = new double[matrix_size];


     num_my_length = my_rows *(my_rhs + my_cols+6) ;


     num_global_length = -1 ;


     // Define Epetra_map


     Epetra_Map map(num_global_length,num_my_length,
  			 0, comm);


      Epetra_Vector A(map);


     // Set Random values

     if( comm.MyPID() == 0 )
              std::cout << " ****   Setting Random Matrix    ****" << std::endl;


     ierror = A.SetSeed(seed+comm.MyPID() );

     ierror = A.Random();


     // Now Create the RHS


      if( comm.MyPID() == 0 )
              std::cout << " ****   Creating RHS   ****" << std::endl;

     // Sum the portion of the row that I have


      for (k= 0; k < my_rows; k++) {
       temp[k] = 0.;
       for (m=0; m < my_cols; m++) {
        temp[k] = temp[k] + A[m*my_rows+k];
       }
     }

    // Sum to Processor 0


     MPI_Allreduce(temp,temp2,my_rows,MPI_DOUBLE,MPI_SUM,rowcomm);

      if( comm.MyPID() == 0 )
            std::cout << " ****   Packing RHS in Matrix   ****" << std::endl;

    // Now put the RHS in the Appropriate position

      for (k= 0; k < matrix_size; k++) {

        temp4[k]=0.;
        rhs[k]=0. ;
      }

    if( my_rhs > 0 ) {
     for (k= 0; k < my_rows; k++) {

      A[my_rows*my_cols + k] = temp2[k];

      rhs[k+ my_first_row - 1]= temp2[k];

     }
    }

     // Globally Sum the RHS needed for testing later


     MPI_Allreduce(rhs,temp4,matrix_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

     // Pack back into RHS

     for (k=0 ; k < matrix_size; k++) {

       rhs[k] = temp4[k];

     }

     delete [ ] temp4;


    // Now Solve the Problem

    if( comm.MyPID() == 0 )
             std::cout << " ****   Beginning Matrix Solve   ****" << std::endl;

     solver.FactorSolve(&A,
		        my_rows,
		        my_cols,
                        &matrix_size,
                        &nprocs_row,
                        &nrhs,
		        &secs);



      if( comm.MyPID() == 0)   {
         std::cout << " ----  Solution time  ----   "
	 << secs << "  in secs. " << std::endl;


	 mflops = 2./3.*pow(matrix_size,3.)/secs/1000000.;

	   std::cout << " *****   MFLOPS   *****  " << mflops << std::endl;
      }
     // Now Check the Solution

     // Delete temp variables they need to be resized

     delete [] temp;

     delete [] temp2;


     temp = new double[matrix_size];

     temp2 = new double[matrix_size];


     // Intialize the Matrices

      for (k= 0; k < matrix_size; k++) {

	temp[k]=0. ;
        temp2[k]=0. ;
      }


      // Pack the Answer into the apropriate position

     if ( my_rhs > 0) {

      for (k= 0; k < my_rows; k++) {

       temp[k+ my_first_row -1] = A[my_rows*my_cols + k];

      }
     }


     // All processors get the answer

      MPI_Allreduce(temp,temp2,matrix_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      solinf = fabs(temp2[0]);

      xh = fabs(temp2[0] -one);

      for (k= 0; k < matrix_size; k++) {


       if ( fabs (temp2[k]) > solinf ) solinf = fabs(temp2[k]);

       if ( fabs (temp2[k]- one) > xh ) xh = fabs(temp2[k]-one);

      }

     // Reset the matrix Random values


     ierror = A.SetSeed(seed + comm.MyPID() );


     ierror = A.Random();


     ierror = A.NormInf(&ainf);

    // perform the Matrix vector product


     for (k= 0; k < my_rows; k++) {
      temp[k+my_first_row-1] = 0.;
     for (l=0; l < my_cols; l++) {
       temp[k+my_first_row-1] = temp[k+my_first_row-1] + A[l*my_rows+k]*temp2[l+my_first_col-1];
     }
    }


     double * temp3 = new double[matrix_size];

      for (k= 0; k < matrix_size; k++) {

        temp3[k]=0. ;
      }


   //   Epetra_Map map(numGlobalEquations, numLocalEquations,
   //			update, 0, comm);

     MPI_Allreduce(temp,temp3,matrix_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


     if( comm.MyPID() == 0)   {

       std::cout <<  "======================================" << std::endl;
        std::cout << " ---- Error Calculation ----" << std::endl;

        axbinf = fabs(temp3[0]-rhs[0]);

        for (k= 0; k < matrix_size; k++) {
	
           if ( fabs (temp3[k]-rhs[k]) > axbinf ) axbinf = fabs(temp3[k]-rhs[k]);

         }

     }


    // Machine epsilon Calculation

      othird = four_thirds - 1.;

      tempc = othird + othird + othird;

      eps = fabs(tempc-1.0);

      if ( comm.MyPID() == 0 ) {

	std::cout << "   Machine eps  " << eps  << std::endl;

      }



     if ( comm.MyPID() == 0 ) {



       std::cout << "   ||Ax - b||_oo = " << axbinf << std::endl;

       std::cout << "   ||A||_oo ~  " << ainf*matrix_size << std::endl;

       std::cout << "   ||X||_oo = " << solinf << std::endl;

       std::cout << "   ||Ax - b||_oo / ||A||_oo/||X||_oo  = " << axbinf/(ainf*matrix_size)/solinf  << std::endl;

       std::cout << "  matrix_size*eps*1.  " << matrix_size*eps*1. << std::endl;

      if ( axbinf/(ainf*matrix_size)/solinf  > (matrix_size*eps*1.))

        std::cout << " ****    Solution Fails   ****" <<  std::endl;

      else

	std::cout << " ****   Solution Passes   ****" << std::endl;

      std::cout <<  "======================================" << std::endl;


     }


     //  Delete all the vectors used


     delete [] temp2;

     delete [] temp;

     delete [] rhs;

     delete [] temp3;



  MPI_Finalize() ;

return (0);
}
