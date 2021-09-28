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

class Test {

public:

  Test(const char* filename, int nIter) {
    testEpetra(filename, nIter);
    testTpetra(filename, nIter);
    Teuchos::TimeMonitor::summarize();
  }


private:
  
  ////////////////////////////////////////////////////////////////////////
  void testEpetra(const char* filename, int nIter)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    using crsmatrix_t = Epetra_CrsMatrix;
    using multivector_t = Epetra_MultiVector;

    std::cout << "testEpetra:  begin file " << filename << std::endl;

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
    multivector_t y(Amat->RangeMap(), 2);
    multivector_t x(Amat->DomainMap(), 2);
    x.Random();

    multivector_t w(x);
    multivector_t z(x);

    // For several iterations, measure time for replacing values
    using ttm = Teuchos::TimeMonitor;
    for (int it = 0; it < nIter; it++) {

      //////////////////////////////////////////////////////////////
      // Grab and hold pointers into matrix values (ick!)
      // Replace matrix values 
      {
        int nrows;
        double **matValues;
        int *nEntries;
        {
          ttm t(*ttm::getNewTimer("CRS Setup Epetra Xyce Held pointers"));
          nrows = Amat->RowMap().NumMyElements();
          matValues = new double*[nrows];
          nEntries = new int[nrows];

          for (int i = 0; i < nrows; i++) {
            int nEnt;
            double *mVals;
            int *indices;
            Amat->ExtractMyRowView(i, nEnt, mVals, indices);
            matValues[i] = mVals;
            nEntries[i] = nEnt;
          }
        }

        {
          ttm t(*ttm::getNewTimer("CRS Replace Epetra Xyce Held pointers"));
          for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < nEntries[i]; j++)
              matValues[i][j] = newMatValues[j];
          }
        }

        delete [] nEntries;
        delete [] matValues;
      }

      // Replace matrix values using view of values
      {
        ttm t(*ttm::getNewTimer("CRS Replace Epetra using views"));
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

      // Replace matrix values using ReplaceMyValues (as in Xyce)
      {
        ttm t(*ttm::getNewTimer("CRS Replace Epetra ReplaceLocalValues"));
        int nEntries;
        double *matValues;
        int *indices;
        int nrows = Amat->RowMap().NumMyElements();
        for (int i = 0; i < nrows; i++) {
          Amat->ExtractMyRowView(i, nEntries, matValues, indices);
          Amat->ReplaceMyValues(i, nEntries, newMatValues, indices);
        }
      }

      // SpMV
      {
        ttm t(*ttm::getNewTimer("CRS Apply Epetra"));
        Amat->Multiply(false, x, y);
      }
  
      // LinearCombination
      {
        ttm t(*ttm::getNewTimer("CRS LinearCombo Epetra "
                              "using bracket operator as in Xyce"));

        double a = 10.;
        double b = -10;
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

      {
        ttm t(*ttm::getNewTimer("CRS LinearCombo Epetra using views"));

        double a = 10.;
        double b = -10;
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

      // Set all matrix values to a single value
      {
        ttm t(*ttm::getNewTimer("CRS Init Epetra set all to scalar "));
        Amat->PutScalar(1.0);
      }

      //////////////////////////////////////////////////////////////
      // Replace vector values using view (not sure how Xyce does it)
      {
        ttm t(*ttm::getNewTimer("MV Replace Epetra using views"));
        double **xvalues;
        int len = x.MyLength();
        x.ExtractView(&xvalues);
        for (int k = 0; k < x.NumVectors(); k++) {
          for (int i = 0; i < len; i++) {
            xvalues[k][i] *= 3.;
          }
        }
      }

      // Replace vector values using bracket operator 
      {
        ttm t(*ttm::getNewTimer("MV Replace Epetra BracketOperator"));
        int len = x.MyLength();
        for (int k = 0; k < x.NumVectors(); k++) {
          for (int i = 0; i < len; i++) {
            x[k][i] *= 3.;
          }
        }
      }

      // Vector linear combination
      {
        ttm t(*ttm::getNewTimer("MV Update Epetra a*A+b*this"));
        double a = 1.2;
        double b = 2.1;
        w.Update(a, x, b);
      }

      {
        ttm t(*ttm::getNewTimer("MV Update Epetra a*A+b*B+c*this"));
        double a = 1.2;
        double b = 2.1;
        double c = 3.2;
        w.Update(a, x, b, z, c);
      }
    }
    delete [] newMatValues;
    std::cout << "testEpetra:  done " << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////
  void testTpetra(const char* filename, int nIter)
  {
    std::cout << "testTpetra:  begin file " << filename << std::endl;

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
    multivector_t y(Amat->getRangeMap(), 2);
    multivector_t x(Amat->getDomainMap(), 2);
    x.randomize();
    multivector_t w(x);
    multivector_t z(x);

    // For several iterations, measure time for replacing values
    using ttm = Teuchos::TimeMonitor;
    for (int it = 0; it < nIter; it++) {

      //////////////////////////////////////////////////////////////
      // Grab and hold pointers into matrix values (ick!)
      // Replace matrix values 
      {
        lno_t nrows;
        scalar_t **matValues;
        lno_t *nEntries;
        {
          ttm t(*ttm::getNewTimer("CRS Setup Tpetra Xyce Held pointers"));
          nrows = Amat->getRowMap()->getNodeNumElements();
          matValues = new scalar_t*[nrows];
          nEntries = new lno_t[nrows];

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

        {
          ttm t(*ttm::getNewTimer("CRS Replace Tpetra Xyce Held pointers"));
          for (lno_t i = 0; i < nrows; i++) {
            for (lno_t j = 0; j < nEntries[i]; j++)
              matValues[i][j] = newMatValues[j];
          }
        }

        delete [] nEntries;
        delete [] matValues;
      }


      // Replace matrix values using local matrix
      {
        ttm t(*ttm::getNewTimer("CRS Replace Tpetra using views"));
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

      // Replace matrix values using ReplaceMyValues (as in Xyce)
      {
        ttm t(*ttm::getNewTimer("CRS Replace Tpetra ReplaceLocalValues"));
        kval_t matValues;
        kind_t indices;
        lno_t nrows = Amat->getRowMap()->getNodeNumElements();
        for (lno_t i = 0; i < nrows; i++) {
          Amat->getLocalRowView(i, indices, matValues);
          Amat->replaceLocalValues(i, indices, newMatValues);
        }
      }

      // SpMV
      {
        ttm t(*ttm::getNewTimer("CRS Apply Tpetra"));
        Amat->apply(x, y);
      }
  
      {
        ttm t(*ttm::getNewTimer("CRS LinearCombo Tpetra using views"));

        scalar_t a(10.);
        scalar_t b(-10);

        auto avalues = Amat->getLocalMatrixHost().values;
        auto bvalues = Bmat.getLocalMatrixHost().values;
        auto cvalues = Cmat.getLocalMatrixHost().values;

        int nEntries = avalues.extent(0);

        for (int i = 0; i < nEntries; i++) {
          cvalues[i] = a * avalues[i] + b * bvalues[i];
        }
      }

      // Set all matrix values to a single value
      if (!(it % 10) && (comm->getRank() == 0)) 
        std::cout << it << " of " << nIter 
                  << ": CRS Init Tpetra set all to scalar" << std::endl;

      {
        ttm t(*ttm::getNewTimer("CRS Init Tpetra set all to scalar "));
        Amat->setAllToScalar(1.0);
      }

      //////////////////////////////////////////////////////////////
      // Replace vector values using view (not sure how Xyce does it)
      {
        ttm t(*ttm::getNewTimer("MV Replace Tpetra using views"));
        auto xvalues = x.getLocalViewHost(Tpetra::Access::ReadWrite);
        lno_t len = x.getLocalLength();
        int nVec = x.getNumVectors();
        for (int k = 0; k < nVec; k++) {
          for (int i = 0; i < len; i++) {
            xvalues(i,k) *= 3.;
          }
        }
      }

      // Vector linear combination
      {
        ttm t(*ttm::getNewTimer("MV Update Tpetra a*A+b*this"));
        double a = 1.2;
        double b = 2.1;
        w.update(a, x, b);
      }

      {
        ttm t(*ttm::getNewTimer("MV Update Tpetra a*A+b*B+c*this"));
        double a = 1.2;
        double b = 2.1;
        double c = 3.2;
        w.update(a, x, b, z, c);
      }
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

