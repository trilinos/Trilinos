
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

// Trilinos Tutorial
// -----------------
// Using IFPACK factorizations as Aztec's preconditioners
//
// NOTE: this example implemenets minor modifications to one of the
// examples included in the AztecOO package. Please give a look
// to file ${TRILINOS_HOME}/packages/aztecoo/examples/IfpackAztecOO/cxx_main.cpp
// for more details.
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 19-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

// function for fancy output

string toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

string toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// main driver

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  // B E G I N   O F   M A T R I X   C O N S T R U C T I O N
  
  /* this example creates a tridiagonal matrix of type
   *
   *     |  2   0  -1              |
   *     |  0   2  0   -1          | 
   * A = |  -1  0  2   0   -1  ... |
   *     |             ... ... ... |
   *     |        -1    0   2   0  |        
   *     |             -1   0   2  |
   */
  
  // set global dimension to 5, could be any number > 5
  int NumGlobalElements = 100;
  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create a Epetra_Matrix
  // create a CSR matrix

  Epetra_CrsMatrix A(Copy,Map,3);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
  int NumEntries;

  for( int i=0 ; i<NumMyElements ; ++i ) {
    // off-diagonal terms
    if( MyGlobalElements[i]==0 ) {
      Indices[0] = 2;
      NumEntries = 1;
    } else if(  MyGlobalElements[i]==1 ) {
      Indices[0] = 3;
      NumEntries = 1;
    } else if( MyGlobalElements[i]==NumGlobalElements-2 ) {
      Indices[0] = NumGlobalElements-4;
      NumEntries = 1;
    } else if(  MyGlobalElements[i]==NumGlobalElements-1 ) {
      Indices[0] = NumGlobalElements-3;
      NumEntries = 1;
    } else {
      Indices[0] = MyGlobalElements[i]-2;
      Indices[1] = MyGlobalElements[i]+2;
      NumEntries = 2;
    }
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i)==0);
  }

  // Finish up
  assert(A.TransformToLocal()==0);

  // E N D   O F   M A T R I X   C O N S T R U C T I O N  

  // ============================ //
  // Construct ILU preconditioner //
  // ---------------------------- //

  //  modify those parameters 
  int    LevelFill = 0;
  int    Overlap = 0;
  double Athresh = 0.0;
  double Rthresh = 1.0;

  Ifpack_IlukGraph * IlukGraph = 0;
  Ifpack_CrsRiluk * ILUK = 0;

  IlukGraph = new Ifpack_IlukGraph(A.Graph(), LevelFill, Overlap);
  assert(IlukGraph->ConstructFilledGraph()==0);

  Epetra_Flops fact_counter;

  Epetra_Time timer(Comm);
  
  double elapsed_time = timer.ElapsedTime();
  ILUK = new Ifpack_CrsRiluk(*IlukGraph);
  ILUK->SetFlopCounter(fact_counter);
  ILUK->SetAbsoluteThreshold(Athresh);
  ILUK->SetRelativeThreshold(Rthresh);

  int initerr = ILUK->InitValues(A);
  if (initerr!=0) cout << Comm << "InitValues error = " << initerr;

  assert(ILUK->Factor()==0);
  elapsed_time = timer.ElapsedTime() - elapsed_time;
  double total_flops = ILUK->Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "Time to compute preconditioner values = " 
		    << elapsed_time << endl
		    << "MFLOPS for Factorization = " << MFLOPs << endl;
  // try to uncomment the following for some output
  // cout << *ILUK << endl;
  
  bool transA = true;
  double Condest;
  ILUK->Condest(transA, Condest);

  if (verbose) cout << "Condition number estimate for this preconditioner = "
		    << Condest << endl;

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRiluk Preconditioner: LevelFill = " + toString(LevelFill) + 
                                                 " Overlap = " + toString(Overlap) + 
                                                 " Athresh = " + toString(Athresh) + 
                                                 " Rthresh = " + toString(Rthresh); 
  ILUK->SetLabel(label.c_str());

  //int Maxiter = 500;
  //double Tolerance = 1.0E-14;


  Epetra_Vector x(Map);
  x.PutScalar(0.0);
  Epetra_Vector b(Map);
  b.Random();
  
  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  x.SetFlopCounter(A);
  b.SetFlopCounter(A);
  ILUK->SetFlopCounter(A);

  // Here we create an AztecOO object
  AztecOO solver;
  solver.SetUserMatrix(&A);
  solver.SetLHS(&x);
  solver.SetRHS(&b);

  // Here we set the IFPACK preconditioner and specify few parameters
  
  solver.SetPrecOperator(ILUK);

  int Niters = 1200;
  solver.SetAztecOption(AZ_kspace, Niters); 
  solver.Iterate(Niters, 5.0e-10);

  // number of flops
  elapsed_time = timer.ElapsedTime() - elapsed_time;
  total_flops = counter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "Time to compute solution = " 
		    << elapsed_time << endl
		    << "Number of operations in solve = " << total_flops << endl
		    << "MFLOPS for Solve = " << MFLOPs<< endl << endl;

  if (ILUK!=0) delete ILUK;
  if (IlukGraph!=0) delete IlukGraph;
				       
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:ifpack]> mpirun -np 2 ./ex1.exe
Time to compute preconditioner values = 0.00021
MFLOPS for Factorization = 1.84762
Condition number estimate for this preconditioner = 84.5

                *******************************************************
                ***** Preconditioned GMRES solution
                ***** Ifpack_CrsRiluk Preconditioner: LevelFill = 0 Overlap = 0 Athresh = 0 Rthresh = 1
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 5.731834e-01
                iter:    2           residual = 5.401640e-01
                iter:    3           residual = 1.289039e-14


                Solution time: 0.003409 (sec.)
                total iterations: 3
Time to compute solution = 0.00533505
Number of operations in solve = 5192
MFLOPS for Solve = 0.973186

*/
