/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#include <iostream>

#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "KokkosBlas.hpp"
#include "KokkosSparse.hpp"

#include <vector>

class Test {

public:

  Test(const char* filename, int nIter) {
    testTpetra(filename, nIter);
    testEpetra(filename, nIter);
    Teuchos::TimeMonitor::summarize();
  }


private:
  
  ////////////////////////////////////////////////////////////////////////
  void testEpetra(const char* filename, int nIter)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    using crsmatrix_t = Epetra_CrsMatrix;
    using multivector_t = Epetra_MultiVector;

    std::cout << "testEpetra:  begin file " << filename 
              << "; " << nIter << " iterations" << std::endl;

    // Read in the matrix
    crsmatrix_t *Amat;
    EpetraExt::MatrixMarketFileToCrsMatrix(filename, comm, Amat);

    std::cout << comm.MyPID() << " Epetra Matrix:  "
              << Amat->NumMyRows() << " rows; " 
              << Amat->NumMyCols() << " cols; " 
              << Amat->NumMyNonzeros() << " nnz; " 
              << Amat->NumGlobalRows() << " global rows; "
              << Amat->NumGlobalNonzeros() << " global nnz"
              << std::endl;

    crsmatrix_t Bmat(*Amat);
    crsmatrix_t Cmat(*Amat);

    // Create some new matrix values to be used in tests
    int maxValuesPerRow = Amat->MaxNumEntries();
    double *newMatValues = new double[maxValuesPerRow];
    for (int i = 0; i < maxValuesPerRow; i++) newMatValues[i] = i;

    // Create multivectors for SpMV and LinearCombination
    std::vector<multivector_t> Y;
    std::vector<multivector_t> X;
    std::vector<multivector_t> W;
    std::vector<multivector_t> Z;
    for (int it = 0; it < nIter; it++) {
      multivector_t y(Amat->RangeMap(), 2);
      multivector_t x(Amat->DomainMap(), 2);
      x.Random();

      multivector_t w(x);
      multivector_t z(x);

      Y.push_back(y);
      X.push_back(x);
      W.push_back(w);
      Z.push_back(z);
    }

    // For several iterations, measure time for replacing values
    using ttm = Teuchos::TimeMonitor;

    //////////////////////////////////////////////////////////////
    {
      // Grab and hold pointers into matrix values (ick!)
      // Replace matrix values 
      int nrows;
      double **matValues;
      int *nEntries;
      {
        nrows = Amat->RowMap().NumMyElements();
        matValues = new double*[nrows];
        nEntries = new int[nrows];

        auto tm = ttm::getNewTimer("CRS Setup Epetra Xyce Held pointers");
        tm->start();

        for (int it = 0; it < nIter; it++) {
          for (int i = 0; i < nrows; i++) {
            int nEnt;
            double *mVals;
            int *indices;
            Amat->ExtractMyRowView(i, nEnt, mVals, indices);
            matValues[i] = mVals;
            nEntries[i] = nEnt;
          }
        }

        tm->stop();
      }

      {
        auto tm = ttm::getNewTimer("CRS Replace Epetra Xyce Held pointers");
        tm->start();

        for (int it = 0; it < nIter; it++) {
          for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < nEntries[i]; j++)
              matValues[i][j] = newMatValues[j];
          }
        }

        tm->stop();
      }

      delete [] nEntries;
      delete [] matValues;
    }

    // Replace matrix values using view of values
    {
      auto tm = ttm::getNewTimer("CRS Replace Epetra using views");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        int nEntries;
        double *matValues;
        int *indices;
        int nrows = Amat->RowMap().NumMyElements();
        for (int i = 0; i < nrows; i++) {
          Amat->ExtractMyRowView(i, nEntries, matValues, indices);
          for (int j = 0; j < nEntries; j++)
            matValues[j] = newMatValues[j];
        }
      }

      tm->stop();
    }

    // Replace matrix values using ReplaceMyValues (as in Xyce)
    {
      auto tm = ttm::getNewTimer("CRS Replace Epetra ReplaceLocalValues");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        int nEntries;
        double *matValues;
        int *indices;
        int nrows = Amat->RowMap().NumMyElements();
        for (int i = 0; i < nrows; i++) {
          Amat->ExtractMyRowView(i, nEntries, matValues, indices);
          Amat->ReplaceMyValues(i, nEntries, newMatValues, indices);
        }
      }

      tm->stop();
    }

    // SpMV
    // Warm-up -- don't time it
    Amat->Multiply(false, X[nIter-1], Y[nIter-1]);
    {
      auto tm = ttm::getNewTimer("CRS Apply Epetra");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        Amat->Multiply(false, X[it], Y[it]);
      }

      tm->stop();
    }
  
    // LinearCombination
    {
      double a = 10.;
      double b = -10;

      auto tm = ttm::getNewTimer("CRS LinearCombo Epetra "
                                 "using bracket operator");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        double *avalues, *bvalues;
        int *aind, *bind;
        int anEntries, bnEntries;

        int nrows = Amat->RowMap().NumMyElements();
        for (int i = 0; i < nrows; i++) {
          Amat->ExtractMyRowView(i, anEntries, avalues, aind);
          Bmat.ExtractMyRowView(i, bnEntries, bvalues, bind);
          for (int j = 0; j < anEntries; j++) 
            Cmat[i][j] = a * avalues[j] + b * bvalues[j];
        }
      }

      tm->stop();
    }

    {
      double a = 10.;
      double b = -10;

      auto tm = ttm::getNewTimer("CRS LinearCombo Epetra using views");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        double *avalues, *bvalues, *cvalues;
        int *aind, *bind, *cind;
        int anEntries, bnEntries, cnEntries;

        int nrows = Amat->RowMap().NumMyElements();
        for (int i = 0; i < nrows; i++) {
          Amat->ExtractMyRowView(i, anEntries, avalues, aind);
          Bmat.ExtractMyRowView(i, bnEntries, bvalues, bind);
          Cmat.ExtractMyRowView(i, cnEntries, cvalues, cind);
          for (int j = 0; j < anEntries; j++) 
            cvalues[j] = a * avalues[j] + b * bvalues[j];
        }
      }

      tm->stop();
    }

    // Set all matrix values to a single value
    {
      auto tm = ttm::getNewTimer("CRS Init Epetra set all to scalar ");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        Amat->PutScalar(1.0);
      }

      tm->stop();
    }

    //////////////////////////////////////////////////////////////
    // Replace vector values using view (not sure how Xyce does it)
    {
      auto tm = ttm::getNewTimer("MV Replace Epetra using views");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        double **xvalues;
        int len = X[it].MyLength();
        X[it].ExtractView(&xvalues);
        for (int k = 0; k < X[it].NumVectors(); k++) {
          for (int i = 0; i < len; i++) {
            xvalues[k][i] *= 3.;
          }
        }
      }

      tm->stop();
    }

    // Replace vector values using bracket operator 
    {
      auto tm = ttm::getNewTimer("MV Replace Epetra BracketOperator");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        int len = X[it].MyLength();
        for (int k = 0; k < X[it].NumVectors(); k++) {
          for (int i = 0; i < len; i++) {
            X[it][k][i] *= 3.;
          }
        }
      }

      tm->stop();
    }

    // Vector linear combination
    {
      double a = 1.2;
      double b = 2.1;

      auto tm = ttm::getNewTimer("MV Update Epetra a*A+b*this");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        W[it].Update(a, X[it], b);
      }

      tm->stop();
    }

    {
      double a = 1.2;
      double b = 2.1;
      double c = 3.2;

      auto tm = ttm::getNewTimer("MV Update Epetra a*A+b*B+c*this");
      tm->start();

      for (int it = 0; it < nIter; it++) {
        W[it].Update(a, X[it], b, Z[it], c);
      }

      tm->stop();
    }

    delete [] newMatValues;
    std::cout << "testEpetra:  done " << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  void testTpetra(const char* filename, int nIter)
  {
    std::cout << "testTpetra:  begin file " << filename
              << "; " << nIter << " iterations" << std::endl;

    auto comm = Tpetra::getDefaultComm();

    // Read in the matrix
    using crsmatrix_t = Tpetra::CrsMatrix<>;
    using multivector_t = Tpetra::MultiVector<>;

    using lno_t = typename multivector_t::local_ordinal_type;
    using scalar_t = typename multivector_t::scalar_type;

    using kind_t = typename crsmatrix_t::local_inds_host_view_type;
    using kval_t = typename crsmatrix_t::values_host_view_type;
    using nc_kval_t = typename crsmatrix_t::nonconst_values_host_view_type;

    using reader_t = Tpetra::MatrixMarket::Reader<crsmatrix_t>;

    Teuchos::ParameterList params;
    Teuchos::RCP<crsmatrix_t> Amat = 
      reader_t::readSparseFile(filename, comm, params);

    std::cout << comm->getRank() << " Tpetra Matrix:  "
              << Amat->getNodeNumRows() << " rows; " 
              << Amat->getNodeNumCols() << " cols; " 
              << Amat->getNodeNumEntries() << " nnz; " 
              << Amat->getGlobalNumRows() << " global rows; "
              << Amat->getGlobalNumEntries() << " global nnz"
              << std::endl;

    crsmatrix_t Bmat(*Amat);
    crsmatrix_t Cmat(*Amat);

    // Create some new matrix values to be used in tests
    lno_t maxValuesPerRow = Amat->getNodeMaxNumRowEntries();
    nc_kval_t newMatValues("newMatValues", maxValuesPerRow);
    for (lno_t i = 0; i < maxValuesPerRow; i++) newMatValues[i] = scalar_t(i);

    // Create multivectors for SpMV
    std::vector<multivector_t> Y;
    std::vector<multivector_t> X;
    std::vector<multivector_t> W;
    std::vector<multivector_t> Z;
    for (int it = 0; it < nIter; it++) {
      multivector_t y(Amat->getRangeMap(), 2);
      multivector_t x(Amat->getDomainMap(), 2);
      x.randomize();

      multivector_t w(x);
      multivector_t z(x);

      Y.push_back(y);
      X.push_back(x);
      W.push_back(w);
      Z.push_back(z);
    }

    // For several iterations, measure time for replacing values
    using ttm = Teuchos::TimeMonitor;

    //////////////////////////////////////////////////////////////
    // Grab and hold pointers into matrix values (ick!)
    // Replace matrix values 
    {
      lno_t nrows;
      scalar_t **matValues;
      lno_t *nEntries;
      {
        nrows = Amat->getRowMap()->getNodeNumElements();
        matValues = new scalar_t*[nrows];
        nEntries = new lno_t[nrows];

        auto tm = ttm::getNewTimer("CRS Setup Tpetra Xyce Held pointers");
        tm->start();

        for (int it = 0; it < nIter; it++) {
          auto lclMat = Amat->getLocalMatrixHost();
          for (lno_t i = 0; i < nrows; i++) {
            nc_kval_t mVals;
            kind_t indices;

            // Bad! Bad! BAD!  Don't store pointers to data.
            // This code mimics what Xyce does for its devices.
            // But this code is completely unsafe and incorrect when Kokkos 
            // uses separate host and device spaces.
            // Since we don't check the results of the matrix operations
            // in this test, we get away with it, but a normal application
            // would see inconsistent data between host and device.
            // Please, never ever store pointers this way in a real application!
            matValues[i] = lclMat.values.data() + lclMat.graph.row_map[i]; 
            // Bad! Bad! BAD!  Don't store pointers to data

            nEntries[i] = lclMat.graph.row_map[i+1] - lclMat.graph.row_map[i];
          }
        }
        tm->stop();
      }

      {
        auto tm = ttm::getNewTimer("CRS Replace Tpetra Xyce Held pointers");
        tm->start();
        for (int it = 0; it < nIter; it++) {
          for (lno_t i = 0; i < nrows; i++) {
            for (lno_t j = 0; j < nEntries[i]; j++)
              matValues[i][j] = newMatValues[j];
          }
        }
        tm->stop();
      }

      delete [] nEntries;
      delete [] matValues;
    }

    // Replace matrix values using local matrix
    {
      auto tm = ttm::getNewTimer("CRS Replace Tpetra using views");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        lno_t nrows = Amat->getRowMap()->getNodeNumElements();
        auto lclMat = Amat->getLocalMatrixHost();
        auto offsets = lclMat.graph.row_map;
        auto matValues = lclMat.values;
        for (lno_t i = 0; i < nrows; i++) {
          lno_t nEntries = offsets[i+1] - offsets[i];
          for (lno_t j = offsets[i], k = 0; k < nEntries; j++, k++) {
            matValues[j] = newMatValues[k];
          }
        }
      }
      tm->stop();
    }

    // Replace matrix values using ReplaceMyValues (as in Xyce)
    {
      auto tm = ttm::getNewTimer("CRS Replace Tpetra ReplaceLocalValues");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        kval_t matValues;
        kind_t indices;
        lno_t nrows = Amat->getRowMap()->getNodeNumElements();
        for (lno_t i = 0; i < nrows; i++) {
          Amat->getLocalRowView(i, indices, matValues);
          Amat->replaceLocalValues(i, indices, newMatValues);
        }
      }
      tm->stop();
    }

    // SpMV
    // Warm-up -- don't time it
    Amat->apply(X[nIter-1], Y[nIter-1]);
    {
      auto tm = ttm::getNewTimer("CRS Apply Tpetra");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        Amat->apply(X[it], Y[it]);
      }
      tm->stop();
    }
  
    if (comm->getSize() == 1)
    {
      auto lclmat = Amat->getLocalMatrixHost();
      using mvdata_t = typename multivector_t::dual_view_type::t_host;
      std::vector<typename mvdata_t::const_type> xvec;
      std::vector<mvdata_t> yvec;
      for (int it = 0; it < nIter; it++) {
        xvec.push_back(X[it].getLocalViewHost(Tpetra::Access::ReadOnly));
        yvec.push_back(Y[it].getLocalViewHost(Tpetra::Access::ReadWrite));
      }

      auto tm = ttm::getNewTimer("CRS Apply Kokkos "
                                 "-- local only, KokkosSparse::spmv");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        KokkosSparse::spmv(KokkosSparse::NoTranspose, 1., lclmat, xvec[it],
                                                      0., yvec[it]);
      }
      tm->stop();
    }
  
    {
      scalar_t a(10.);
      scalar_t b(-10);

      auto tm = ttm::getNewTimer("CRS LinearCombo Tpetra using views");
      tm->start();
      for (int it = 0; it < nIter; it++) {

        auto avalues = Amat->getLocalMatrixHost().values;
        auto bvalues = Bmat.getLocalMatrixHost().values;
        auto cvalues = Cmat.getLocalMatrixHost().values;

        int nEntries = avalues.extent(0);

        for (int i = 0; i < nEntries; i++) {
          cvalues[i] = a * avalues[i] + b * bvalues[i];
        }
      }
      tm->stop();
    }

    // Set all matrix values to a single value
    {
      auto tm = ttm::getNewTimer("CRS Init Tpetra set all to scalar ");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        Amat->setAllToScalar(1.0);
      }
      tm->stop();
    }

    {
      auto vals = Amat->getLocalMatrixHost().values;

      auto tm = ttm::getNewTimer("CRS Init Kokkos set all to scalar "
                                 " -- Kokkos::deep_copy");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        Kokkos::deep_copy(vals, 1.0);
      }
      tm->stop();
    }

    //////////////////////////////////////////////////////////////
    // Replace vector values using view (not sure how Xyce does it)
    {
      auto tm = ttm::getNewTimer("MV Replace Tpetra using views");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        auto xvalues = X[it].getLocalViewHost(Tpetra::Access::ReadWrite);
        lno_t len = X[it].getLocalLength();
        int nVec = X[it].getNumVectors();
        for (int k = 0; k < nVec; k++) {
          for (int i = 0; i < len; i++) {
            xvalues(i,k) *= 3.;
          }
        }
      }
      tm->stop();
    }

    // Vector linear combination
    {
      double a = 1.2;
      double b = 2.1;

      auto tm = ttm::getNewTimer("MV Update Tpetra a*A+b*this");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        W[it].update(a, X[it], b);
      }
      tm->stop();
    }

    {
      double a = 1.2;
      double b = 2.1;

      using mvdata_t = typename multivector_t::dual_view_type::t_host;
      std::vector<typename mvdata_t::const_type> xvec;
      std::vector<mvdata_t> wvec;
      for (int it = 0; it < nIter; it++) {
        xvec.push_back(X[it].getLocalViewHost(Tpetra::Access::ReadOnly));
        wvec.push_back(W[it].getLocalViewHost(Tpetra::Access::ReadWrite));
      }

      auto tm = ttm::getNewTimer("MV Update Kokkos a*A+b*this"
                                 " -- local only, KokkosBlas:axpby");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        KokkosBlas::axpby (a, xvec[it], b, wvec[it]);
      }
      tm->stop();
    }

    {
      double a = 1.2;
      double b = 2.1;
      double c = 3.2;

      auto tm = ttm::getNewTimer("MV Update Tpetra a*A+b*B+c*this");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        W[it].update(a, X[it], b, Z[it], c);
      }
      tm->stop();
    }

    {
      double a = 1.2;
      double b = 2.1;
      double c = 3.2;

      using mvdata_t = typename multivector_t::dual_view_type::t_host;
      std::vector<typename mvdata_t::const_type> xvec;
      std::vector<typename mvdata_t::const_type> zvec;
      std::vector<mvdata_t> wvec;

      for (int it = 0; it < nIter; it++) {
        xvec.push_back(X[it].getLocalViewHost(Tpetra::Access::ReadOnly));
        zvec.push_back(Z[it].getLocalViewHost(Tpetra::Access::ReadOnly));
        wvec.push_back(W[it].getLocalViewHost(Tpetra::Access::ReadWrite));
      }

      auto tm = ttm::getNewTimer("MV Update Kokkos a*A+b*B+c*this"
                                 " -- local only, KokkosBlas:update");
      tm->start();
      for (int it = 0; it < nIter; it++) {
        KokkosBlas::update(a, xvec[it], b, zvec[it], c, wvec[it]);
      }
      tm->stop();
    }
    std::cout << "testTpetra:  done " << std::endl;
  }

private:
};

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard scope(&narg, &arg);

  if (narg < 2) {
    std::cout << "Usage:  a.out filename.mtx niter" << std::endl;
    return -1;
  }

  auto comm = Tpetra::getDefaultComm();

  int niter = 10;
  if (narg == 3) niter = std::atoi(arg[2]);

  Test test(arg[1], niter);

  if (comm->getRank() == 0) std::cout << "TEST PASSED" << std::endl;
  return 0;
}

