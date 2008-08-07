/*
This comment includes my makefile I used to compile this example.  It is ugly but this is what I had to do to get things to work.  Good luck if you try.  IK 08-06-2008

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
*/

//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_Time.h"

#ifdef HAVE_OSKI
#ifdef HAVE_EPETRA_TEUCHOS
#ifdef HAVE_EPETRAEXT
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

  // Read XML input deck
  ParameterList masterList;
  if (argc > 1) {
    if (strncmp("-h",argv[1],2) == 0) {
      cout << "help" << endl;
    }
    else {
      int i=0,j;
      FILE* fid = fopen(argv[1],"r");
      if (fid) {
        i++;
        fclose(fid);
      }
      Comm.SumAll(&i, &j, 1);
      if (j!=Comm.NumProc()) {
        cout << "Could not open input file." << endl;
      }
      FileInputSource fileSrc(argv[1]);
      XMLObject fileXML = fileSrc.getObject();
      XMLParameterListReader ListReader;
      masterList = ListReader.toParameterList(fileXML);
    }
  } else {
    cout << "No input file specified." << endl;
}

  char hostname[100];
  char buf[100];
  char go;  


  ParameterList *fileList, *tests, *tunings;
  ParameterList List;
  ParameterList List2;
  fileList = &(masterList.sublist("data files",true));
  tests = &(masterList.sublist("tests"));
  tunings = &(masterList.sublist("tunings"));

  
  string matrixfile = fileList->get("matrix input file","A.dat");
  const char *datafile = matrixfile.c_str();
  int MultiVecs =  tests->get("multivecs", 5);
  int Power = tests->get("power", 2);
  bool Multi = tests->get("multi", false);
  bool MatVec = tests->get("matvec", false);
  bool Solve = tests->get("solve", false);
  bool ATA = tests->get("ata", false);
  bool AAT = tests->get("aat", false);
  bool APowX = tests->get("apowx", false);
  bool TwoMult = tests->get("twomult", false);
  bool TwoTrans = tests->get("twotrans", false);
  bool Trans = tests->get("trans", false);
  bool singleblock = tunings->get("singleblock", false);
  if(singleblock)
    List2.set("singleblocksize", true); 
  int rows = tunings->get("row", 0);
  if(rows)
    List2.set("row", 3); 
  int cols = tunings->get("col", 0);
  if(cols)
    List2.set("col", 3); 
  bool alligned = tunings->get("alligned", false);
  if(alligned)
    List2.set("alignedblocks", true); 
  bool correlated = tunings->get("correlated", false);
  if(correlated)
    List2.set("correlatedpattern", true); 
  int numGlobalRows;

  // ===================================================== //
  // READ IN MATRICES FROM FILE                            //
  // ===================================================== //

  if (!mypid) printf("reading %s\n",datafile); fflush(stdout);
  Epetra_CrsMatrix *Amat=NULL;
  //Epetra_Map *RowMap=NULL;
  int errCode=0;
  
  int CRAP[1];
  int N[1];
  int NZ[1];
  int NZ2[1];

  N[0] = 21879;
  NZ[0] = 21879;
  NZ2[0] = 1571926;

  MPI_Barrier(MPI_COMM_WORLD);

  std::ifstream inFile("A.dat", std::ios::in | std::ios::binary);

  if (inFile == NULL)
    EPETRA_CHK_ERR(-1); // file not found

  inFile.read((char*)N, sizeof(N));
  inFile.read((char*)NZ, sizeof(NZ));
  inFile.read((char*)NZ2, sizeof(NZ2));
  
  std::cerr << N[0] << " " << NZ[0] << " " << NZ2[0] << "\n";
  inFile.close(); 
  Epetra_Map* rowmap;
  Epetra_Map* colmap;

  if(Comm.NumProc() - 1 != mypid) {
    rowmap = new Epetra_Map (N[0], 0, Comm);
    colmap = new Epetra_Map (NZ[0], 0, Comm);
  }
  else {
    rowmap = new Epetra_Map (N[0], 0, Comm);
    colmap = new Epetra_Map (NZ[0], 0, Comm);
  }

  errCode=EpetraExt::BinaryFileToCrsMatrix(datafile, *rowmap, *colmap, Amat);
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
  //need tuning list params set here
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
  delete [] trans;
  delete colmap;
  delete Amat;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} //main
#endif
#endif
