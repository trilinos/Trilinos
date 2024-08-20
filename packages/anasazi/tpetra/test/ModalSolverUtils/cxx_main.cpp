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

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_Vector_fwd.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

// Anasazi
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicOrthoManager.hpp"

#include "../../test/MVOPTester/MySDMHelpers.hpp"

template<typename ScalarType>
int run(int argc, char *argv[]) {
  using ST  = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO  = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO  = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT  = typename Tpetra::MultiVector<>::node_type;
  using SCT = typename Teuchos::ScalarTraits<ScalarType>;

  using OP = Tpetra::Operator<ST,LO,GO,NT>;
  using MV = Tpetra::MultiVector<ST,LO,GO,NT>;

  using tmap_t = Tpetra::Map<LO,GO,NT>;

  using MVT   = Anasazi::MultiVecTraits<ST,MV>;
  using Utils = Anasazi::SolverUtils<ST,MV,OP>;

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  const auto comm = Tpetra::getDefaultComm();
  const int MyPID = comm->getRank();

  bool verbose = false;
  if (argc>1) if (strncmp("-v",argv[1],2) == 0) verbose = true;

  if (verbose && MyPID == 0) {
    std::cout << Anasazi::Anasazi_Version() << std::endl << std::endl;
  }

  int numberFailedTests = 0;

  //  Dimension of the multivector
  int NumGlobalElements = 99;
  const int NumColumns = 7;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  RCP<tmap_t> Map = rcp(new tmap_t(NumGlobalElements, 0, comm));

  //--------------------------------------------------------------------------
  //  test Householder code
  //--------------------------------------------------------------------------
  {
    Teuchos::LAPACK<int,ST> lapack;
    ST err;

    if (verbose && MyPID == 0) {
      std::cout << std::endl << "************* Householder Apply Test *************" << std::endl << std::endl;
    }

    // initialize an OrthoManager
    Anasazi::BasicOrthoManager<ST,MV,OP> orthman;

    // generate random multivector V and orthonormalize it
    MV V(Map,NumColumns);
    MV VQ(Map,NumColumns);
    MVT::MvRandom(V);
    orthman.normalize(V,Teuchos::null);

    // generate random orthogonal matrix
    int info;
    std::vector<ST> qrwork(NumColumns), tau(NumColumns);
    Teuchos::SerialDenseMatrix<int,ST> Q(NumColumns,NumColumns), H(NumColumns,NumColumns);
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
      std::cout << "             orthonorm error of V: " << err << std::endl;
    }

    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        std::cout<< "ERROR:  orthonormalization failed." << std::endl;
      }
    }

    // generate explicit V*Q
    MVT::MvTimesMatAddMv(1.0,V,Q,0.0,VQ);

    // test orthonormality of V*Q
    err = orthman.orthonormError(VQ);
    if (verbose && MyPID == 0) {
      std::cout << "            orthonorm error of VQ: " << err << std::endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        std::cout<< "ERROR:  V*Q failed." << std::endl;
      }
    }

    // apply house(V,H,tau)
    Utils::applyHouse(H.numCols(),V,H,tau);
    err = orthman.orthonormError(V);
    if (verbose && MyPID == 0) {
      std::cout << "    orthonorm error of applyHouse: " << err << std::endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        std::cout<< "ERROR:  applyHouse failed." << std::endl;
      }
    }

    // test house(V,H,tau) == V*Q
    err = Utils::errorEquality(V,VQ);
    if (verbose && MyPID == 0) {
      std::cout << "        error(VQ - house(V,H,tau): " << err << std::endl;
    }
    if (err > 1e-14) {
      numberFailedTests++;
      if (verbose && MyPID == 0) {
        std::cout<< "ERROR:  applyHouse failed." << std::endl;
      }
    }

  }

  //--------------------------------------------------------------------------
  //  test directSolver, permuteVectors
  //--------------------------------------------------------------------------
  {
    if (verbose && MyPID == 0) {
      std::cout << std::endl << "************* DirectSolver Test *************" << std::endl << std::endl;
    }

    int size = 11;
    int nev  = 7;
    Anasazi::BasicSort<ST> sorter("SR");

    // form random eigenvalues
    std::vector<ST> lambda1(nev);
    Teuchos::SerialDenseMatrix<int,ST> Lambda(nev,1);
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
      Teuchos::SerialDenseMatrix<int,ST> K(size,size), Q(nev,nev);
      std::vector<ST> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::null,Q,lambda2,rank,10);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(10) returned error " << info << std::endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(10) didn't return all eigenpairs" << std::endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev; i++) {
          if (SCT::magnitude(lambda2[i] - lambda1[i]) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              std::cout << "ERROR: directSolve(diag(lambda)) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << lambda1[i] << " vs. " << lambda2[i] << std::endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            std::cout << "pass: directSolve(diag(lambda)) correct." << std::endl;
          }
        }
      }
    }

    // step two: eigenvalues of diag(k),diag(m) are k./m
    {
      Teuchos::SerialDenseMatrix<int,ST> K(size,size), M(size,size), Q(nev,nev);
      std::vector<ST> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
        M(i,i) = 2.0;
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::rcp(&M,false),Q,lambda2,rank,1);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(1) returned error " << info << std::endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(10) didn't return all eigenpairs" << std::endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev; i++) {
          if (SCT::magnitude(lambda2[i] - K(i,i)/M(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              std::cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << std::endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            std::cout << "pass: directSolve(diag(lambda),2I) correct." << std::endl;
          }
        }
      }
    }


    // step three: directsolve of diag(k),diag([m 0]) fails appropriately
    {
      Teuchos::SerialDenseMatrix<int,ST> K(size,size), M(size,size), Q(nev,nev);
      std::vector<ST> lambda2(nev);
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
          std::cout << "ERROR: directSolve(0) returned error " << info << std::endl;
        }
      }
      else if (rank != nev-2) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(10) didn't return all eigenpairs" << std::endl;
        }
      }
      else {
        bool testfailed = false;
        for (int i=0; i<nev-2; i++) {
          if (SCT::magnitude(lambda2[i] - K(i,i)/M(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              std::cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " 
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << std::endl;
            }
            break;
          }
        }
        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            std::cout << "pass: directSolve(diag(lambda),2I) correct." << std::endl;
          }
        }
      }
    }

    // step four: 1) solve K = Q*L*Q'
    //            2) permute columns of Q
    //            3) shows that Q'*K*Q gives permuted L
    // this tests the eigenvectors and the permutation routine
    {
      Teuchos::SerialDenseMatrix<int,ST> K(size,size), M(size,size), Q(nev,nev),
                                             T1(nev,nev), TK(nev,nev), TM(nev,nev);
      std::vector<ST> lambda2(nev);
      for (int i=0; i<nev; i++) {
        K(i,i) = lambda1[i];
        M(i,i) = 2.0;
      }
      int rank = nev;
      int info = Utils::directSolver(nev,K,Teuchos::rcp(&M,false),Q,lambda2,rank,0);
      if (info != 0) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(0) returned error " << info << std::endl;
        }
      }
      else if (rank != nev) {
        numberFailedTests++;
        if (verbose && MyPID == 0) {
          std::cout << "ERROR: directSolve(10) didn't return all eigenpairs" << std::endl;
        }
      }
      else {
        Teuchos::SerialDenseMatrix<int,ST> KK(Teuchos::View,K,nev,nev),
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
              std::cout << "ERROR: directSolve(diag(lambda),2I) produced wrong eigenvalues: " << std::endl
                   << "i: " << i << "   " << K(i,i)/M(i,i) << " vs. " << lambda2[i] << std::endl;
            }
            break;
          }
          // check permuted Q'*K*Q
          if (SCT::magnitude(lambda2[rperm[i]] - TK(i,i)) > 1e-14) {
            testfailed = true;
            numberFailedTests++;
            if (verbose && MyPID==0) {
              std::cout << "ERROR: Q'*K*Q diagonals don't match lambdas" << std::endl
                   << "i: " << i << "   " << TK(i,i) << " vs. " << lambda2[rperm[i]] << std::endl;
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
            std::cout << "ERROR: permuted directSolve(diag(lambda),2I),  produced non-K-orthogonal Ritz vectors: " << std::endl
              << "| Q'*K*Q - L |: " <<  TK.normFrobenius() << std::endl;
          }
        }
        // check Q'*M*Q == I
        for (int i=0; i<nev; i++) {
          TM(i,i) -= 1.0;
        }
        if (TM.normFrobenius() > 1e-14) {
          numberFailedTests++;
          if (verbose && MyPID==0) {
            std::cout << "ERROR: permuted directSolve(diag(lambda),2I) produced non-M-orthonormal Ritz vectors: " << std::endl
              << "| Q'*M*Q - I |: " <<  TM.normFrobenius() << std::endl;
          }
        }

        if (testfailed == false) {
          if (verbose && MyPID == 0) {
            std::cout << "pass: directSolve(diag(lambda),2I) with permute correct." << std::endl;
          }
        }
      }
    }
  }

  if (numberFailedTests) {
    if (verbose && MyPID==0) {
      std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}

int main(int argc, char *argv[]) {
  run<double>(argc,argv);
  // run<float>(argc,argv);
  return 0;
}
