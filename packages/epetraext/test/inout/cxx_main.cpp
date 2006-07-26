/*@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
// include the header files again to check if guards are in place
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include <string>
#include "Poisson2dOperator.h"
// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    return(1); 
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}
int runTests(Epetra_Map & map, Epetra_CrsMatrix & A, Epetra_Vector & x, Epetra_Vector & b, Epetra_Vector & xexact, bool verbose);
int runOperatorTests(Epetra_Operator & A, bool verbose);

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

  #ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(ierr);
}

int runTests(Epetra_Map & map, Epetra_CrsMatrix & A, Epetra_Vector & x, Epetra_Vector & b, Epetra_Vector & xexact, bool verbose) {

  int ierr = 0;

  double residual;
  residual = A.NormInf(); double rAInf = residual;
  if (verbose) cout << "Inf Norm of A                                                     = " << residual << endl;
  residual = A.NormOne(); double rAOne = residual;
  if (verbose) cout << "One Norm of A                                                     = " << residual << endl;
  xexact.Norm2(&residual); double rxx = residual;
  if (verbose) cout << "Norm of xexact                                                    = " << residual << endl;
  Epetra_Vector tmp1(map);
  A.Multiply(false, xexact, tmp1);
  tmp1.Norm2(&residual); double rAx = residual;
  if (verbose) cout << "Norm of Ax                                                        = " << residual << endl;
  b.Norm2(&residual); double rb = residual;
  if (verbose) cout << "Norm of b (should equal norm of Ax)                               = " << residual << endl;
  tmp1.Update(1.0, b, -1.0);
  tmp1.Norm2(&residual);
  if (verbose) cout << "Norm of difference between compute Ax and Ax from file            = " << residual << endl;
  map.Comm().Barrier();

  EPETRA_CHK_ERR(EpetraExt::BlockMapToMatrixMarketFile("Test_map.mm", map, "Official EpetraExt test map", 
						       "This is the official EpetraExt test map generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::RowMatrixToMatrixMarketFile("Test_A.mm", A, "Official EpetraExt test matrix", 
							"This is the official EpetraExt test matrix generated by the EpetraExt regression tests"));

  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_x.mm", x, "Official EpetraExt test initial guess", 
						     "This is the official EpetraExt test initial guess generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_xexact.mm", xexact, "Official EpetraExt test exact solution", 
						     "This is the official EpetraExt test exact solution generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::VectorToMatrixMarketFile("Test_b.mm", b, "Official EpetraExt test right hand side", 
						     "This is the official EpetraExt test right hand side generated by the EpetraExt regression tests"));
				       
  EPETRA_CHK_ERR(EpetraExt::RowMatrixToMatlabFile("Test_A.dat", A));

  Epetra_Map * map1;
  Epetra_CrsMatrix * A1; 
  Epetra_CrsMatrix * A2; 
  Epetra_CrsMatrix * A3; 
  Epetra_Vector * x1; 
  Epetra_Vector * b1;
  Epetra_Vector * xexact1;

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


  residual = A1->NormInf(); double rA1Inf = residual;
  if (verbose) cout << "Inf Norm of A1                                                    = " << residual << endl;
  ierr += checkValues(rA1Inf,rAInf,"Inf Norm of A", verbose);

  residual = A1->NormOne(); double rA1One = residual;
  if (verbose) cout << "One Norm of A1                                                    = " << residual << endl;
  ierr += checkValues(rA1One,rAOne,"One Norm of A", verbose);

  xexact1->Norm2(&residual); double rxx1 = residual;
  if (verbose) cout << "Norm of xexact1                                                   = " << residual << endl;
  ierr += checkValues(rxx1,rxx,"Norm of xexact", verbose);

  Epetra_Vector tmp11(*map1);
  A1->Multiply(false, *xexact1, tmp11);

  tmp11.Norm2(&residual); double rAx1 = residual;
  if (verbose) cout << "Norm of A1*x1                                                     = " << residual << endl;
  ierr += checkValues(rAx1,rAx,"Norm of A1*x", verbose);

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

  tmp11.Update(1.0, *b1, -1.0);
  tmp11.Norm2(&residual);
  if (verbose) cout << "Norm of difference between computed A1x1 and A1x1 from file        = " << residual << endl;
  ierr += checkValues(residual,0.0,"Norm of difference between computed A1x1 and A1x1 from file", verbose);

  if (map1->IndexBase()==0) {delete A2; delete A3;}
  delete A1;
  delete x1;
  delete b1;
  delete xexact1;
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
