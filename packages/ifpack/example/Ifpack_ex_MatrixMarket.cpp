// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AZTECOO)  && defined(HAVE_IFPACK_TEUCHOS) && defined(HAVE_IFPACK_TRIUTILS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AztecOO.h"
#include "Ifpack.h"

#include "Ifpack_Utils.h"
// This example reads a matrix in MatrixMarket format on process 0,
// then it redistributed this matrix, using all the available
// processes, build an IFPACK preconditioner, and uses AztecOO to
// solve a linear problem with this matrix, and solution corresponding
// to a random RHS.
//
// Parameters for IFPACK preconditioners can be tuned at the 
// command line. Type 
// 
// $ ./Ifpack_ex_MatrixMarket.exe --help
// 
// for a list of available options. For example, the user amy type:
//
// $ ./Ifpack_ex_MatrixMarket --matrix=e05r0500.mtx --prec=ILUT --fill=2
//
// to solve using ILUT preconditioner, and a level-of-fill of 2.

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // process the command line. 

  Teuchos::CommandLineProcessor CLP;
  // amount of overlap among processes
  int Overlap = 0;
  CLP.setOption("overlap", &Overlap, "Overlap among processes");
  // amount of overlap among local blocks for Ifpack_BlockJacobi
  int BlockOverlap = 0;
  CLP.setOption("block-overlap", &BlockOverlap, "Overlap among blocks");
  // chosen preconditioner
  string PrecType = "ILUT";
  CLP.setOption("prec", &PrecType, "IFPACK preconditioner");
  // level-of-fill for IC and ILU
  int LevelOfFill = 0;
  CLP.setOption("fill", &LevelOfFill, "level-of-fill for IC and ILU");
  // relaxation value for ICT and ILUT
  double Relax = 0.0;
  CLP.setOption("relax", &Relax, "relaxation value for ICT and ILUT");
  double AddToDiag = 1e-5;
  CLP.setOption("add-diag", &AddToDiag, "value to be added to diagonals");
  // dropping value 
  double DropTol = 0.0;
  CLP.setOption("drop", &DropTol, "drop all elements below this value");
  // number of local blocks for Ifpack_Jacobi, GaussSeidel and 
  // symmetric GaussSeidel
  int LocalParts = 4;
  CLP.setOption("local-parts", &LocalParts, "number of local blocks");
  // matrix name
  string FileName = "not-set";
  CLP.setOption("matrix", &FileName, "Name of file containing the MTX matrix");
  // is true, the matrix contains only half of the matrix
  int SymFormat = false;
  CLP.setOption("sym-matrix", &SymFormat, "Set to non-zero value if the matrix is stored in symmetric format (that is, only the nonzeros of the lower or upper part are stored)");

  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  if (FileName == "not-set") {
    cerr << "You must at least specify the name of the file" << endl;
    cerr << "containing the matrix, using --matrix=my_file.mtx" << endl;
    cerr << "Run this example with option --help for more details" << endl;
    // return success not to break tests
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }
      
  int NumRows;
  int NumCols;
  int Offset = 1;
  int NumElements;
  std::ifstream data_file;

  if (Comm.MyPID() == 0) {

    string Title;

    data_file.open(FileName.c_str());

    if (!data_file.good()) {
      cerr << "Error opening file `" << FileName << "'" << endl;
      exit(EXIT_FAILURE);
    }

    char title[100];
    data_file.getline(title,100);
    
    data_file >> NumRows;
    data_file >> NumCols;
    data_file >> NumElements;

    assert (NumCols == NumRows);

    cout << "Label = `" << title << "'" << endl;
    cout << "Number of rows             = " << NumRows << endl;
    cout << "Number of nonzero elements = " << NumElements << endl;
    cout << "Offset                     = " << Offset << endl;
  }
  else
    NumRows = 0;
  
  // creates a map with all elements on proc 0
  Epetra_Map* SerialMap = new Epetra_Map(-1,NumRows,0,Comm);
  Epetra_CrsMatrix* SerialMatrix;
  SerialMatrix = new Epetra_CrsMatrix(Copy,*SerialMap,0);

  if (Comm.MyPID() == 0) {

    // now proc 0 read the actual matrix, element by element
    for (int i = 0 ; i < NumElements ; ++i) {
      int row;
      int col;
      double val;
      data_file >> row;
      data_file >> col;
      data_file >> val;
      row -= Offset;
      col -= Offset;
      IFPACK_CHK_ERR(SerialMatrix->InsertGlobalValues(row,1,&val,&col));
      if (col != row && SymFormat)
        IFPACK_CHK_ERR(SerialMatrix->InsertGlobalValues(col,1,&val,&row));
    }

    // want to be sure that the diagonal element is set,
    // AztecOO may suffer in parallel.
    for (int i = 0 ; i < NumRows ; ++i) {
      double val = 0.0;
      SerialMatrix->InsertGlobalValues(i,1,&val,&i);
    }

  }
  SerialMatrix->FillComplete();

  // need to create the distributed map, this
  // is for simplicity linear
  Comm.Broadcast(&NumRows,1,0);
  Epetra_Map DistributedMap(NumRows, 0, Comm);

  Epetra_CrsMatrix A(Copy, DistributedMap,0);

  // creates the import 
  Epetra_Import Importer(DistributedMap,*SerialMap);

  IFPACK_CHK_ERR(A.Import(*SerialMatrix,
			  Importer, Insert));
  
  IFPACK_CHK_ERR(A.FillComplete());
  
  // can delete serial objects, no longer needed
  delete SerialMap;
  delete SerialMatrix;

  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // allocates an IFPACK factory. No data is associated 
  // to this object.
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &A, Overlap);
  assert(Prec != 0);

  // specify parameters
  List.set("fact: level-of-fill", LevelOfFill);
  List.set("fact: relax value", Relax);

  // other parameters
  List.set("partitioner: overlap", BlockOverlap);
#ifdef HAVE_IFPACK_METIS
  List.set("partitioner: type", "METIS");
#else
  List.set("partitioner: type", "greedy");
#endif
  List.set("partitioner: local parts", LocalParts);

  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(List));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Builds the preconditioners, by looking for the values of 
  // the matrix.
  IFPACK_CHK_ERR(Prec->Compute());

  // One might want to check the estimated condition number.
  // This is a very crude estimate, but if this value is very
  // large, probably one should re-build the preconditioner.

  double Estimated = Prec->Condest();
  if (Comm.MyPID() == 0)
    cout << "Estimated condition number = " << Estimated << endl;

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //

  // At this point, we need some additional objects
  // to define and solve the linear system.

  // defines LHS and RHS
  Epetra_Vector LHS(A.OperatorDomainMap());
  Epetra_Vector RHS(A.OperatorDomainMap());

  // solution is constant
  LHS.PutScalar(1.0);
  // now build corresponding RHS
  A.Apply(LHS,RHS);

  // now randomize the solution
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(&A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify solver and output.
  Solver.SetAztecOption(AZ_solver,AZ_gmres);
  Solver.SetAztecOption(AZ_output,32);

  // HERE WE SET THE IFPACK PRECONDITIONER
  Solver.SetPrecOperator(Prec);

  // .. and here we solve
  Solver.Iterate(1550,1e-5);

  // delete the preconditioner
  delete Prec;
  
#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --enable-aztecoo");
  puts("--enable-triutils --enable-teuchos");
  puts("to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
