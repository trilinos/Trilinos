/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_Version.h"
#ifdef EPETRA_MPI
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
#include "Epetra_Export.h"

#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_VectorIn.h"
// include the header files again to check if guards are in place
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_VectorIn.h"
#include <string>
#include <vector>
#include "Poisson2dOperator.h"
// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
    return(1); 
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}
int runTests(Epetra_Map & map, Epetra_CrsMatrix & A, Epetra_Vector & x, Epetra_Vector & b, Epetra_Vector & xexact, bool verbose);
int runOperatorTests(Epetra_Operator & A, bool verbose);
int generateHyprePrintOut(const char* filename, const Epetra_Comm &comm);
int runHypreTest(Epetra_CrsMatrix &A);

int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int MyPID = comm.MyPID();

  bool verbose = false;
  bool verbose1 = false; 
  // Check if we should print results to standard out
  if (argc > 1) {
    if ((argv[1][0] == '-') && (argv[1][1] == 'v')) {
      verbose1 = true;
      if (MyPID==0) verbose = true;
    }
  }
  if (verbose)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  if (verbose1) cout << comm << endl;


  // Uncomment the next three lines to debug in mpi mode
  //int tmp;
  //if (MyPID==0) cin >> tmp;
  //comm.Barrier();

  Epetra_Map * map;
  Epetra_CrsMatrix * A; 
  Epetra_Vector * x; 
  Epetra_Vector * b;
  Epetra_Vector * xexact;

  int nx = 20*comm.NumProc();
  int ny = 30;
  int npoints = 7;
  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

   
  int ierr = 0;
  // Call routine to read in HB problem 0-base
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  ierr += runTests(*map, *A, *x, *b, *xexact, verbose);

  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

  // Call routine to read in HB problem 1-base
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact, 1);

  ierr += runTests(*map, *A, *x, *b, *xexact, verbose);

  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

  // Call routine to read in HB problem -1-base
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact, -1);

  ierr += runTests(*map, *A, *x, *b, *xexact, verbose);

  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

  int nx1 = 5;
  int ny1 = 4;
  Poisson2dOperator Op(nx1, ny1, comm);
  ierr += runOperatorTests(Op, verbose);

  generateHyprePrintOut("MyMatrixFile", comm);

  EPETRA_CHK_ERR(EpetraExt::HypreFileToCrsMatrix("MyMatrixFile", comm, A));
  
  runHypreTest(*A);
  delete A;

  #ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(ierr);
}

int runHypreTest(Epetra_CrsMatrix &A){
  
  Epetra_Vector X(A.RowMap());
  EPETRA_CHK_ERR(X.Random());
  Epetra_Vector Y(A.RowMap());
  EPETRA_CHK_ERR(A.Multiply(false, X, Y));
  
  return 0;
}

int runTests(Epetra_Map & map, Epetra_CrsMatrix & A, Epetra_Vector & x, Epetra_Vector & b, Epetra_Vector & xexact, bool verbose) {

  int ierr = 0;

  // Create MultiVectors and put x, b, xexact in both columns of X, B, and Xexact, respectively.
  Epetra_MultiVector X( map, 2, false );
  Epetra_MultiVector B( map, 2, false );
  Epetra_MultiVector Xexact( map, 2, false );

  for (int i=0; i<X.NumVectors(); ++i) {
    *X(i) = x;
    *B(i) = b;
    *Xexact(i) = xexact;
  }
  double residual;
  std::vector<double> residualmv(2);
  residual = A.NormInf(); double rAInf = residual;
  if (verbose) cout << "Inf Norm of A                                                     = " << residual << endl;
  residual = A.NormOne(); double rAOne = residual;
  if (verbose) cout << "One Norm of A                                                     = " << residual << endl;
  xexact.Norm2(&residual); double rxx = residual;	
  Xexact.Norm2(&residualmv[0]); std::vector<double> rXX( residualmv );	
  if (verbose) cout << "Norm of xexact                                                    = " << residual << endl;
  if (verbose) cout << "Norm of Xexact                                                    = (" << residualmv[0] << ", " <<residualmv[1] <<")"<< endl;
  Epetra_Vector tmp1(map);
  Epetra_MultiVector tmp1mv(map,2,false);
  A.Multiply(false, xexact, tmp1);
  A.Multiply(false, Xexact, tmp1mv);
  tmp1.Norm2(&residual); double rAx = residual;
  tmp1mv.Norm2(&residualmv[0]); std::vector<double> rAX( residualmv );
  if (verbose) cout << "Norm of Ax                                                        = " << residual << endl;
  if (verbose) cout << "Norm of AX                                                        = (" << residualmv[0] << ", " << residualmv[1] <<")"<< endl;
  b.Norm2(&residual); double rb = residual;
  B.Norm2(&residualmv[0]); std::vector<double> rB( residualmv );
  if (verbose) cout << "Norm of b (should equal norm of Ax)                               = " << residual << endl;
  if (verbose) cout << "Norm of B (should equal norm of AX)                               = (" << residualmv[0] << ", " << residualmv[1] <<")"<< endl;
  tmp1.Update(1.0, b, -1.0);
  tmp1mv.Update(1.0, B, -1.0);
  tmp1.Norm2(&residual);
  tmp1mv.Norm2(&residualmv[0]);
  if (verbose) cout << "Norm of difference between compute Ax and Ax from file            = " << residual << endl;
  if (verbose) cout << "Norm of difference between compute AX and AX from file            = (" << residualmv[0] << ", " << residualmv[1] <<")"<< endl;
  map.Comm().Barrier();

  EPETRA_CHK_ERR(EpetraExt::BlockMapToMatrixMarketFile("Test_map.mm", map, "Official EpetraExt test map", 
						       "This is the official EpetraExt test map generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::RowMatrixToMatrixMarketFile("Test_A.mm", A, "Official EpetraExt test matrix", 
							"This is the official EpetraExt test matrix generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_x.mm", x, "Official EpetraExt test initial guess", 
						     "This is the official EpetraExt test initial guess generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::MultiVectorToMatrixMarketFile("Test_mvX.mm", X, "Official EpetraExt test initial guess", 
					      	          "This is the official EpetraExt test initial guess generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_xexact.mm", xexact, "Official EpetraExt test exact solution", 
						     "This is the official EpetraExt test exact solution generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::MultiVectorToMatrixMarketFile("Test_mvXexact.mm", Xexact, "Official EpetraExt test exact solution", 
						          "This is the official EpetraExt test exact solution generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_b.mm", b, "Official EpetraExt test right hand side", 
						     "This is the official EpetraExt test right hand side generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::MultiVectorToMatrixMarketFile("Test_mvB.mm", B, "Official EpetraExt test right hand side", 
						          "This is the official EpetraExt test right hand side generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::MultiVectorToMatlabFile("Test_mvB.mat", B));
 
  EPETRA_CHK_ERR(EpetraExt::RowMatrixToMatlabFile("Test_A.dat", A));

  Epetra_Map * map1;
  Epetra_CrsMatrix * A1; 
  Epetra_CrsMatrix * A2; 
  Epetra_CrsMatrix * A3; 
  Epetra_Vector * x1; 
  Epetra_Vector * b1;
  Epetra_Vector * xexact1;
  Epetra_MultiVector * X1; 
  Epetra_MultiVector * B1;
  Epetra_MultiVector * Xexact1;

  EpetraExt::MatrixMarketFileToMap("Test_map.mm", map.Comm(), map1);

  if (map.SameAs(*map1)) {
    if (verbose) cout << "Maps are equal.  In/Out works." << endl;
  }
  else {
    if (verbose) cout << "Maps are not equal.  In/Out fails." << endl;
    ierr += 1;
  }
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToCrsMatrix("Test_A.mm", *map1, A1));
  // If map is zero-based, then we can compare to the convenient reading versions
  if (map1->IndexBase()==0) EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToCrsMatrix("Test_A.mm", map1->Comm(), A2));
  if (map1->IndexBase()==0) EPETRA_CHK_ERR(EpetraExt::MatlabFileToCrsMatrix("Test_A.dat", map1->Comm(), A3));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToVector("Test_x.mm", *map1, x1));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToVector("Test_xexact.mm", *map1, xexact1));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToVector("Test_b.mm", *map1, b1));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector("Test_mvX.mm", *map1, X1));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector("Test_mvXexact.mm", *map1, Xexact1));
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector("Test_mvB.mm", *map1, B1));

  residual = A1->NormInf(); double rA1Inf = residual;
  if (verbose) cout << "Inf Norm of A1                                                    = " << residual << endl;
  ierr += checkValues(rA1Inf,rAInf,"Inf Norm of A", verbose);

  residual = A1->NormOne(); double rA1One = residual;
  if (verbose) cout << "One Norm of A1                                                    = " << residual << endl;
  ierr += checkValues(rA1One,rAOne,"One Norm of A", verbose);

  xexact1->Norm2(&residual); double rxx1 = residual;
  if (verbose) cout << "Norm of xexact1                                                   = " << residual << endl;
  ierr += checkValues(rxx1,rxx,"Norm of xexact", verbose);

  Xexact1->Norm2(&residualmv[0]); std::vector<double> rXX1(residualmv);
  if (verbose) cout << "Norm of Xexact1                                                   = (" << residualmv[0] <<", " <<residualmv[1]<<")"<< endl;
  ierr += checkValues(rXX1[0],rXX[0],"Norm of Xexact", verbose);
  ierr += checkValues(rXX1[1],rXX[1],"Norm of Xexact", verbose);

  Epetra_Vector tmp11(*map1);
  A1->Multiply(false, *xexact1, tmp11);

  Epetra_MultiVector tmp11mv(*map1,2,false);
  A1->Multiply(false, *Xexact1, tmp11mv);

  tmp11.Norm2(&residual); double rAx1 = residual;
  if (verbose) cout << "Norm of A1*x1                                                     = " << residual << endl;
  ierr += checkValues(rAx1,rAx,"Norm of A1*x", verbose);

  tmp11mv.Norm2(&residualmv[0]); std::vector<double> rAX1(residualmv);
  if (verbose) cout << "Norm of A1*X1                                                     = (" << residualmv[0] <<", "<<residualmv[1]<<")"<< endl;
  ierr += checkValues(rAX1[0],rAX[0],"Norm of A1*X", verbose);
  ierr += checkValues(rAX1[1],rAX[1],"Norm of A1*X", verbose);

  if (map1->IndexBase()==0) {
    Epetra_Vector tmp12(*map1);
    A2->Multiply(false, *xexact1, tmp12);
    
    tmp12.Norm2(&residual); double rAx2 = residual;
    if (verbose) cout << "Norm of A2*x1                                                     = " << residual << endl;
    ierr += checkValues(rAx2,rAx,"Norm of A2*x", verbose);

    Epetra_Vector tmp13(*map1);
    A3->Multiply(false, *xexact1, tmp13);
    
    tmp13.Norm2(&residual); double rAx3 = residual;
    if (verbose) cout << "Norm of A3*x1                                                     = " << residual << endl;
    ierr += checkValues(rAx3,rAx,"Norm of A3*x", verbose);
  }
  b1->Norm2(&residual); double rb1 = residual;
  if (verbose) cout << "Norm of b1 (should equal norm of Ax)                              = " << residual << endl;
  ierr += checkValues(rb1,rb,"Norm of b", verbose);

  B1->Norm2(&residualmv[0]); std::vector<double> rB1(residualmv);
  if (verbose) cout << "Norm of B1 (should equal norm of AX)                              = (" << residualmv[0] <<", "<<residualmv[1]<<")"<< endl;
  ierr += checkValues(rB1[0],rB[0],"Norm of B", verbose);
  ierr += checkValues(rB1[1],rB[1],"Norm of B", verbose);

  tmp11.Update(1.0, *b1, -1.0);
  tmp11.Norm2(&residual);
  if (verbose) cout << "Norm of difference between computed A1x1 and A1x1 from file        = " << residual << endl;
  ierr += checkValues(residual,0.0,"Norm of difference between computed A1x1 and A1x1 from file", verbose);

  tmp11mv.Update(1.0, *B1, -1.0);
  tmp11mv.Norm2(&residualmv[0]);
  if (verbose) cout << "Norm of difference between computed A1X1 and A1X1 from file        = (" << residualmv[0] << ", "<<residualmv[1]<<")"<< endl;
  ierr += checkValues(residualmv[0],0.0,"Norm of difference between computed A1X1 and A1X1 from file", verbose);
  ierr += checkValues(residualmv[1],0.0,"Norm of difference between computed A1X1 and A1X1 from file", verbose);

  if (map1->IndexBase()==0) {delete A2; delete A3;}
  delete A1;
  delete x1;
  delete b1;
  delete xexact1;
  delete X1;
  delete B1;
  delete Xexact1;
  delete map1;

  return(ierr);
}

int runOperatorTests(Epetra_Operator & A, bool verbose) {

  int ierr = 0;


  double residual;
  EPETRA_CHK_ERR(EpetraExt::OperatorToMatrixMarketFile("Test_A1.mm", A, "Official EpetraExt test operator", 
							"This is the official EpetraExt test operator generated by the EpetraExt regression tests"));
  EPETRA_CHK_ERR(EpetraExt::OperatorToMatlabFile("Test_A1.dat", A));

  A.OperatorRangeMap().Comm().Barrier();
  A.OperatorRangeMap().Comm().Barrier();
  Epetra_CrsMatrix * A1; 
  Epetra_CrsMatrix * A2; 
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToCrsMatrix("Test_A1.mm", A.OperatorRangeMap(), A1));
  EPETRA_CHK_ERR(EpetraExt::MatlabFileToCrsMatrix("Test_A1.dat", A.OperatorRangeMap().Comm(), A2));


  residual = A.NormInf(); double rAInf = residual;
  if (verbose) cout << "Inf Norm of Operator A                                            = " << residual << endl;
  residual = A1->NormInf(); double rA1Inf = residual;
  if (verbose) cout << "Inf Norm of Matrix A1                                             = " << residual << endl;
  ierr += checkValues(rA1Inf,rAInf,"Inf Norm of A", verbose);


  Epetra_Vector x(A.OperatorDomainMap()); x.Random();
  Epetra_Vector y1(A.OperatorRangeMap());
  Epetra_Vector y2(A.OperatorRangeMap());
  Epetra_Vector y3(A.OperatorRangeMap());
  A.Apply(x,y1);
  A1->Multiply(false, x, y2);
  A2->Multiply(false, x, y3);

  y1.Norm2(&residual); double rAx1 = residual;
  if (verbose) cout << "Norm of A*x                                                       = " << residual << endl;

  y2.Norm2(&residual); double rAx2 = residual;
  if (verbose) cout << "Norm of A1*x                                                      = " << residual << endl;
  ierr += checkValues(rAx1,rAx2,"Norm of A1*x", verbose);

  y3.Norm2(&residual); double rAx3 = residual;
  if (verbose) cout << "Norm of A2*x                                                      = " << residual << endl;
  ierr += checkValues(rAx1,rAx3,"Norm of A2*x", verbose);

  delete A1;
  delete A2;

  return(ierr);
}

int generateHyprePrintOut(const char *filename, const Epetra_Comm &comm){
  int MyPID = comm.MyPID();
  int NumProc = comm.NumProc();

  int N = 100;
  int ilower = MyPID * N;
  int iupper = (MyPID+1)*N-1;

  double filePID = (double)MyPID/(double)100000;
  std::ostringstream stream;
  // Using setprecision() puts it in the string
  stream << std::setiosflags(std::ios::fixed) << std::setprecision(5) << filePID;
  // Then just ignore the first character
  std::string fileName(filename);
  fileName += stream.str().substr(1,7);

  std::ofstream myfile(fileName.c_str());

  if(myfile.is_open()){
    myfile << ilower << " " << iupper << " " << ilower << " " << iupper << endl;
    for(int i = ilower; i <= iupper; i++){
      for(int j=i-5; j <= i+5; j++){
        if(j >= 0 && j < N*NumProc)
          myfile << i << " " << j << " " << (double)rand()/(double)RAND_MAX << endl;
      }
    }
    myfile.close();
    return 0;
  } else {
    cout << "\nERROR:\nCouldn't open file.\n";
    return -1;
  }
}
