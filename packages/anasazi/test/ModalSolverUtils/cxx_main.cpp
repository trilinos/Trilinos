// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the internal utilities that are used by Anasazi solvers.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziSolverUtils.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziBasicSort.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MySDMHelpers.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

typedef Teuchos::ScalarTraits<double> SCT;

int main(int argc, char *argv[]) 
{
  
#ifdef EPETRA_MPI
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  Epetra_SerialComm Comm;

#endif
  
  int MyPID = Comm.MyPID();
  
  bool verbose = false;
  if (argc>1) if (strncmp("-v",argv[1],2) == 0) verbose = true;
  
  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }
  
  int numberFailedTests = 0;

  //  Create SolverUtils object
  typedef Anasazi::SolverUtils<double, Epetra_MultiVector, Epetra_Operator> Utils;
  
  //  Dimension of the multivector
  int NumGlobalElements = 99;
  int NumColumns = 7;
  
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);
  
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double,MV> MVT;

  // declare an orthomanager for use below
  Anasazi::BasicOrthoManager<double,MV,OP> orthman;

  //--------------------------------------------------------------------------
  //  test Householder code
  //--------------------------------------------------------------------------
  {
    Teuchos::LAPACK<int,double> lapack;
    double err;

    if (verbose && MyPID == 0) {
      cout << endl << "************* Householder Apply Test *************" << endl << endl;
    }


    // generate random multivector V and orthonormalize it
    Epetra_MultiVector V(Map,NumColumns), VQ(Map,NumColumns);
    MVT::MvRandom(V);
    orthman.normalize(V,Teuchos::null);

    // generate random orthogonal matrix
    int info;
    std::vector<double> qrwork(NumColumns), tau(NumColumns);
    Teuchos::SerialDenseMatrix<int,double> Q(NumColumns,NumColumns), H(NumColumns,NumColumns);
    Anasazi::randomSDM(H);
    // QR of random matrix H: this puts Householder reflectors in H,tau
    lapack.GEQRF(H.numRows(),H.numCols(),H.values(),H.stride(),&tau[0],&qrwork[0],(int)qrwork.size(),&info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,"error in LAPACK::GEQRF()");
    // generate Q of the QR=H
    Q.assign(H);
    lapack.ORGQR(Q.numRows(),Q.numCols(),Q.numCols(),Q.values(),Q.stride(),&tau[0],&qrwork[0],(int)qrwork.size(),&info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,"error in LAPACK::ORGQR()");

    // test orthonormality of V
    err = orthman.orthonormError(V);
    if (verbose && MyPID == 0) {
      cout << "             orthonorm error of V: " << err << endl;
    }

    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        cout<< "ERROR:  orthonormalization failed." << endl;
      }
    }

    // generate explicit V*Q
    MVT::MvTimesMatAddMv(1.0,V,Q,0.0,VQ);

    // test orthonormality of V*Q
    err = orthman.orthonormError(VQ);
    if (verbose && MyPID == 0) {
      cout << "            orthonorm error of VQ: " << err << endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        cout<< "ERROR:  V*Q failed." << endl;
      }
    }

    // apply house(V,H,tau)
    Utils::applyHouse(H.numCols(),V,H,tau);
    err = orthman.orthonormError(V);
    if (verbose && MyPID == 0) {
      cout << "    orthonorm error of applyHouse: " << err << endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        cout<< "ERROR:  applyHouse failed." << endl;
      }
    }

    // test house(V,H,tau) == V*Q
    err = Utils::errorEquality(V,VQ);
    if (verbose && MyPID == 0) {
      cout << "        error(VQ - house(V,H,tau): " << err << endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        cout<< "ERROR:  applyHouse failed." << endl;
      }
    }

  }

  //--------------------------------------------------------------------------
  //  test directSolver, permuteVectors
  //--------------------------------------------------------------------------
  {
    if (verbose && MyPID == 0) {
      cout << endl << "************* DirectSolver Test *************" << endl << endl;
    }

    int size = 11;
    int nev  = 7;
    Anasazi::BasicSort<double> sorter("SR");

    // form random eigenvalues
    std::vector<double> lambda1(nev);
    SerialDenseMatrix<int,double> Lambda(nev,1)
    Anasazi::randomSDM(Lambda);
    for (int i=0; i<nev; ++i) {
      lambda1[i] = Lambda(i,0);
    }
    // this will order the eigenvalues and give us a random permutation 
    // to use below
    std::vector<int> rperm(nev);
    sorter.sort(lambda1,Teuchos::rcp(&rperm,false),nev);

    // step one: eigenvalues of diag(k) are k
    {
      Teuchos::SerialDenseMatrix<int,double> K(size,size), Q(nev,nev);
      std::vector<double> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::null,Q,lambda2,rank,10);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(10) returned error " << info << endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(10) didn't return all eigenpairs" << endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev; i++) {
          if (SCT::magnitude(lambda2[i] - lambda1[i]) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              cout << "ERROR: directSolve(diag(lambda)) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << lambda1[i] << " vs. " << lambda2[i] << endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            cout << "pass: directSolve(diag(lambda)) correct." << endl;
          }
        }
      }
    }

    // step two: eigenvalues of diag(k),diag(m) are k./m
    {
      Teuchos::SerialDenseMatrix<int,double> K(size,size), M(size,size), Q(nev,nev);
      std::vector<double> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
        M(i,i) = 2.0;
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::rcp(&M,false),Q,lambda2,rank,1);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(1) returned error " << info << endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(10) didn't return all eigenpairs" << endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev; i++) {
          if (SCT::magnitude(lambda2[i] - K(i,i)/M(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            cout << "pass: directSolve(diag(lambda),2I) correct." << endl;
          }
        }
      }
    }


    // step three: directsolve of diag(k),diag([m 0]) fails appropriately
    {
      Teuchos::SerialDenseMatrix<int,double> K(size,size), M(size,size), Q(nev,nev);
      std::vector<double> lambda2(nev);
      // KK,MM have only rank nev-2
      for (int i=0; i<nev-2; i++) {
        K(i,i) = lambda1[i];
        M(i,i) = 2.0;
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::rcp(&M,false),Q,lambda2,rank,0);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(0) returned error " << info << endl;
        }
      }
      else if (rank != nev-2) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(10) didn't return all eigenpairs" << endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev-2; i++) {
          if (SCT::magnitude(lambda2[i] - K(i,i)/M(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            cout << "pass: directSolve(diag(lambda),2I) correct." << endl;
          }
        }
      }
    }

    // step four: 1) solve K = Q*L*Q'
    //            2) permute columns of Q
    //            3) shows that Q'*K*Q gives permuted L
    // this tests the eigenvectors and the permutation routine
    {
      Teuchos::SerialDenseMatrix<int,double> K(size,size), M(size,size), Q(nev,nev), 
                                             T1(nev,nev), TK(nev,nev), TM(nev,nev);
      std::vector<double> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
        M(i,i) = 2.0;
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::rcp(&M,false),Q,lambda2,rank,0);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(0) returned error " << info << endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          cout << "ERROR: directSolve(10) didn't return all eigenpairs" << endl;
        }
      }
      else {
        Teuchos::SerialDenseMatrix<int,double> KK(Teuchos::View,K,nev,nev), 
                                               MM(Teuchos::View,M,nev,nev);
        // permute Q
        Utils::permuteVectors(rperm,Q);
        // compute Q'*K*Q
        info = T1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,KK,Q,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Erroneous call to Teuchos::SDM::multiply.");
        info = TK.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,1.0,Q,T1,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Erroneous call to Teuchos::SDM::multiply.");
        // compute Q'*M*Q
        info = T1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,MM,Q,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Erroneous call to Teuchos::SDM::multiply.");
        info = TM.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,1.0,Q,T1,0.0);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Erroneous call to Teuchos::SDM::multiply.");
        bool testfailed = false;
        for (int i=0; i<nev; i++) {
          // check lambda2
          if (SCT::magnitude(lambda2[i] - K(i,i)/M(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " << endl
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << endl;
            }
            break;
          }
          // check permuted Q'*K*Q
          if (SCT::magnitude(lambda2[rperm[i]] - TK(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              cout << "ERROR: Q'*K*Q diagonals don't match lambdas" << endl
                   << "i: " << i << "   " << TK(i,i) << " vs. " << lambda2[rperm[i]] << endl;
            }
            break;
          }
        }
        // check Q'*K*Q == L
        for (int i=0; i<nev; i++) {
          TK(i,i) = 0.0;
        }
        if (TK.normFrobenius() > 1e-14) {
          testfailed = true;
          numberFailedTests++;
          if (verbose && MyPID==0) {
            cout << "ERROR: permuted directSolve(diag(lambda),2I),  produced non-K-orthogonal Ritz vectors: " << endl
              << "| Q'*K*Q - L |: " <<  TK.normFrobenius() << endl;
          }
        }
        // check Q'*M*Q == I
        for (int i=0; i<nev; i++) {
          TM(i,i) -= 1.0;
        }
        if (TM.normFrobenius() > 1e-14) {
          numberFailedTests++;
          if (verbose && MyPID==0) {
            cout << "ERROR: permuted directSolve(diag(lambda),2I) produced non-M-orthonormal Ritz vectors: " << endl
              << "| Q'*M*Q - I |: " <<  TM.normFrobenius() << endl;
          }
        }

        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            cout << "pass: directSolve(diag(lambda),2I) with permute correct." << endl;
          }
        }
      }
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

 if (numberFailedTests) {
    if (verbose && MyPID==0) {
      cout << endl << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << endl << "End Result: TEST PASSED" << endl;
  }
  return 0;

}
