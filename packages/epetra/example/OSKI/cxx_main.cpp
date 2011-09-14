/*
This example compares OSKI matrix operations to native Epetra operations.  To invoke this example, use

  cxx_main.exe filename

where filename is a Matrix Market formatted data file.  The first two lines *must* be of the form

%%MatrixMarket matrix coordinate real general
XX YY ZZ

where the triplet XX YY ZZ contains the number of global rows, columns, and nonzeros, respectively.

To compile this example, use a makefile similar to that below:

### start of makefile ####

include /home/ikarlin/Trilinos/build_mpi/packages/epetraext/Makefile.export.epetraext

ROOT=cxx_main

FILES=/home/ikarlin/OSKI/install-debug/lib/oski/liboski.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboskilt.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_CSR_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_CSC_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_BCSR_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_MBCSR_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_GCSR_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_CB_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_VBR_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_mat_DENSE_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_heur_regprof_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_heur_symmrb_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski_heur_mregblock_Tid.a \
/home/ikarlin/OSKI/install-debug/lib/oski/liboski.a


all: ${ROOT}.exe

${ROOT}.exe: ${ROOT}.o
        mpicxx -g -O0 -o ${ROOT}.exe ${ROOT}.o ${EPETRAEXT_INCLUDES} ${EPETRAEXT_LIBS} ${FILES} -ldl -lltdl
#       mpicxx -o ${ROOT}.exe ${ROOT}.o ${EPETRAEXT_INCLUDES} ${EPETRAEXT_LIBS} -Wl,--whole-archive `/bin/cat /lib/oski/site-modules-static.txt` -Wl,--no-whole-archive -ldl -lltdl

${ROOT}.o: ${ROOT}.cpp
        mpicxx -g -O0 -c -I/include -DHAVE_CONFIG_H ${EPETRAEXT_INCLUDES} ${ROOT}.cpp

### end of makefile

*/

//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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


#include "Epetra_Time.h"

#ifdef HAVE_OSKI
#ifdef HAVE_EPETRA_TEUCHOS
#include "Epetra_OskiMatrix.h"
#include "Epetra_OskiVector.h"
#include "Epetra_OskiUtils.h"
#include "Epetra_OskiPermutation.h"
#include <unistd.h>

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"

int nonzero;
using namespace Teuchos;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int mypid = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
  int mypid = 0;
#endif

  ParameterList List;
  ParameterList List2;

  const char *datafile;
  if (argc > 1) datafile = argv[1];
  else          datafile = "A.dat";

  // ===================================================== //
  // READ IN MATRICES FROM FILE                            //
  // ===================================================== //

  Epetra_CrsMatrix *Amat=NULL;
  int errCode=0;

  const int lineLength = 1025;
  char line[lineLength];
  int Nrows,Ncols,NZ;
  FILE *handle = fopen(datafile,"r");
  if (handle == 0) {
    if (mypid==0) {
      printf("Cannot open file \"%s\" for reading.\n",datafile);
      printf("usage: cxx_main.exe <filename>, where filename defaults to \"A.dat\".\n");
    }
#   ifdef HAVE_MPI
    MPI_Finalize();
#   endif
    exit(EXIT_FAILURE);
  }
    
  // Strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) {if (handle!=0) fclose(handle);}
  } while (line[0] == '%');
  // Get problem dimensions: #global rows, #global cols, #global nonzeros
  if(sscanf(line,"%d %d %d", &Nrows, &Ncols, &NZ)==0) {if (handle!=0) fclose(handle);}
  fclose(handle);
  
  Epetra_Map* rowmap;
  Epetra_Map* colmap;

  rowmap = new Epetra_Map (Nrows, 0, Comm);
  colmap = new Epetra_Map (Ncols, 0, Comm);

  if (mypid==0) printf("Reading matrix with %d rows, %d columns, %d nonzeros from file \"%s\".\n",Nrows,Ncols,NZ,datafile);
  errCode=EpetraExt::MatrixMarketFileToCrsMatrix(datafile, *rowmap, *colmap, Amat);
  Amat->OptimizeStorage();

  //begin epetra code
  Epetra_Vector x(Amat->DomainMap()); x.Random();
  Epetra_Vector w(Amat->DomainMap()); w.Random();
  Epetra_Vector z(Amat->RowMap()); z.Random();
  Epetra_Vector y(Amat->RowMap()); y.Random();
  Epetra_MultiVector X(Amat->DomainMap(), 5); X.Random();
  Epetra_MultiVector Y(Amat->RowMap(), 5); Y.Random();

  //begin example oski code
  
   
  Epetra_OskiUtils object; //Create an Oski utility object
  object.Init();  //Initialize Oski now we can make Oski matrices and vectors

  //Parameters to create a matrix based on.
  List.set("zerobased", true); //We are zero based in Epetra
  List.set("unique", true); //All entries are unique because optimize storage was called

  Epetra_OskiMatrix* OskiAmat = new Epetra_OskiMatrix(*Amat, List); //Create an OskiMatrix from an Epetra Matrix
  OskiAmat->Multiply(false, x, y);  //Perform y = A*x using OSKI multiply with Epetra Vectors
  OskiAmat->Multiply(true, y, x, 2, 1);  //Perform x = 2*A^T*y + x using OSKI multiply with Epetra Vectors
  
  //Create OskiVectors
  Epetra_OskiVector Oskix(x);
  Epetra_OskiVector Oskiy(y);
  Epetra_OskiVector Oskiw(w);
  Epetra_OskiVector Oskiz(z);

  OskiAmat->Multiply(false, Oskix, Oskiy);  //Perform y = A*x using OSKI multiply with Oski Vectors
  OskiAmat->Multiply(true, Oskiy, Oskix, 2, 1);  //Perform x = 2*A^T*y + x using OSKI multiply with Oski Vectors
  
  //Create OskiMultiVectors
  Epetra_OskiMultiVector OskiX(X);
  Epetra_OskiMultiVector OskiY(Y);

  OskiAmat->Multiply(false, OskiX, OskiY);  //Perform Y = A*X using OSKI multiply with Oski Vectors
  OskiAmat->Multiply(true, OskiY, OskiX, 2.0, 1.0);  //Perform X = 2*A^T*Y + X using OSKI multiply with Oski Vectors

  //Tune Multiply aggressively
  List2.set("singleblocksize", true); //Set machine specific hints here for blocks
  List2.set("row", 3); 
  List2.set("col", 3); 
  List2.set("alignedblocks", true); 
  OskiAmat->SetHintMultiply(false, 1.0, Oskix, 0.0, Oskiy, ALWAYS_TUNE_AGGRESSIVELY, List); //Pass routine specific hints
  OskiAmat->SetHint(List2); //Pass matrix specific hints
  OskiAmat->TuneMatrix();  //Tune matrix
  char* trans;  
  trans = OskiAmat->GetMatrixTransforms();  //Get and print out transforms performed
  std::cout << "Aggressive transforms performed are: " << trans << "\n";
  OskiAmat->Multiply(false, Oskix, Oskiy); //Perform the tuned multiply

  //Done for demonstration purposes
  delete OskiAmat;
  OskiAmat = new Epetra_OskiMatrix(*Amat, List); //Create an OskiMatrix from an Epetra Matrix

  //Tune MultiVec moderately
  //need tuning list params set here
  OskiAmat->SetHintMultiply(true, 2.0, OskiX, 1.0, OskiY, ALWAYS_TUNE, List); //Pass routine specific hints
  OskiAmat->SetHint(List2); //Pass matrix specific hints
  OskiAmat->TuneMatrix();  //Tune matrix
  trans = OskiAmat->GetMatrixTransforms();  //Get and print out transforms performed
  std::cout << "Moderate transforms performed are: " << trans << "\n";
  OskiAmat->Multiply(true, OskiX, OskiY, 2, 1); //Perform the tuned multiply

  //Done for demonstration purposes
  delete OskiAmat;
  OskiAmat = new Epetra_OskiMatrix(*Amat, List); //Create an OskiMatrix from an Epetra Matrix

  //Tune MultiVec based on calls
  OskiAmat->SetHintMultiply(true, 2.0, OskiX, 1.0, OskiY, 10, List); //Pass routine specific hints for 10 calls
  OskiAmat->SetHint(List2); //Pass matrix specific hints
  OskiAmat->TuneMatrix();  //Tune matrix
  trans = OskiAmat->GetMatrixTransforms();  //Get and print out transforms performed
  std::cout << "Moderate transforms performed are: " << trans << "\n";
  OskiAmat->Multiply(true, OskiX, OskiY, 2, 1); //Perform the tuned multiply

  //Composed multiplies
  //ATA
  OskiAmat->MatTransMatMultiply(true, Oskix, Oskiw, NULL);  //Perform the multi not saving the intermediate result.  This will not be tuned or composed by OSKI.
  List.set("invalidinter", true);  //Tell the tune function we're not storing the intermediate.
  OskiAmat->SetHintMatTransMatMultiply(true, 1.0, Oskix, 0.0, Oskiw, Oskiy, ALWAYS_TUNE, List); //Set our tuning hints.
  OskiAmat->TuneMatrix();  //Tune the matrix.
  OskiAmat->MatTransMatMultiply(true, Oskix, Oskiw, NULL); //Call the tuned matrix-vector multiply which will be composed.

  //AAT  In parallel AAT will never be composed and will instead always be two OSKI matvecs since AAT cannot be performed as an atomic operation
  //     without communication in parallel.
  OskiAmat->MatTransMatMultiply(false, Oskiy, Oskiz, NULL); //Same as above
  OskiAmat->SetHintMatTransMatMultiply(false, 1.0, Oskiy, 0.0, Oskiz, Oskiy, ALWAYS_TUNE, List); //This time lets store the intermediate value.
  OskiAmat->TuneMatrix();  //Call the tune function
  OskiAmat->MatTransMatMultiply(false, Oskiy, Oskiz, &Oskix);  //Store the intermediate.  

  //2Mult
  OskiAmat->MultiplyAndMatTransMultiply(true, Oskix, Oskiy, Oskiz, Oskiw);  //Perform the two multiply routine.
  OskiAmat->SetHintMatTransMatMultiply(true, 1.0, Oskix, 0.0, Oskiy, Oskiz, ALWAYS_TUNE, List); //Tune the routine so we use the composed routine.
  OskiAmat->TuneMatrix();  //Tune the Matrix
  OskiAmat->MultiplyAndMatTransMultiply(true, Oskix, Oskiy, Oskiz, Oskiw); //Call the composed routine.

  OskiAmat->MultiplyAndMatTransMultiply(false, Oskix, Oskiy, Oskiw, Oskiz); //Don't transpose the second calculation.

  delete OskiAmat; 
  object.Close(); //close the OSKI object allows OSKI to do any garbage collection or freeing it needs.
  delete rowmap;
  free(trans);
  delete colmap;
  delete Amat;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} //main
#endif
#endif
