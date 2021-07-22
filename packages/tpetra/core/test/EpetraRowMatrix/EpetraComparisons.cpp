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
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

class Test {

public:

  Test(const char* filename, int nIter) {
    testEpetra(filename, nIter);
    testTpetra(filename, nIter);
    Teuchos::TimeMonitor::summarize();
  }


private:
  
  void testEpetra(const char* filename, int nIter)
  {
    Epetra_SerialComm Comm;

    // Read in the matrix
    Epetra_CrsMatrix *Amat;
    EpetraExt::MatrixMarketFileToCrsMatrix(filename, Comm, Amat);

    Epetra_CrsMatrix Bmat(*Amat);
    Epetra_CrsMatrix Cmat(*Amat);

    // Create some new matrix values to be used in tests
    int maxValuesPerRow = Amat->MaxNumEntries();
    double *newMatValues = new double[maxValuesPerRow];
    for (int i = 0; i < maxValuesPerRow; i++) newMatValues[i] = i;

    // Create multivectors for SpMV
    Epetra_MultiVector y(Amat->RangeMap(), 2);
    Epetra_MultiVector x(Amat->DomainMap(), 2);
    x.Random();
    Epetra_MultiVector w(x);
    Epetra_MultiVector z(x);

    // For several iterations, measure time for replacing values
    using tm = Teuchos::TimeMonitor;
    for (int it = 0; it < nIter; it++) {

      //////////////////////////////////////////////////////////////
      // Grab and hold pointers into matrix values (ick!)
      // Replace matrix values 
      {
        int nrows;
        double **matValues;
        int *nEntries;
        {
          tm t(*tm::getNewTimer("Epetra_CrsMatrix Xyce GrabPointers Setup"));
          nrows = Amat->RowMap().NumMyElements();
          matValues = new double*[nrows];
          nEntries = new int[nrows];

          for (int i = 0; i < nrows; i++) {
            int nEnt;
            double *mVals;
            int *indices;
            Amat->ExtractMyRowView(i, nEnt, mVals, indices);
            matValues[i] = new double[nEnt];
            nEntries[i] = nEnt;
          }
        }

        {
          tm t(*tm::getNewTimer("Epetra_CrsMatrix Xyce GrabPointers"));
          for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < nEntries[i]; j++)
              matValues[i][j] = newMatValues[j];
          }
        }

        delete [] nEntries;
        for (int i = 0; i < nrows; i++) delete [] matValues[i];
        delete [] matValues;
      }

      // Replace matrix values using view of values
      {
        tm t(*tm::getNewTimer("Epetra_CrsMatrix ExtractRowView"));
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
        tm t(*tm::getNewTimer("Epetra_CrsMatrix ReplaceMyValues"));
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
        tm t(*tm::getNewTimer("Epetra_Apply"));
        Amat->Multiply(false, x, y);
      }
  
      // LinearCombination
      {
        tm t(*tm::getNewTimer("Epetra_CrsMatrix LinearCombo "
                              "(with bracket operator as in Xyce)"));

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
        tm t(*tm::getNewTimer("Epetra_CrsMatrix LinearCombo (better way)"));

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
        tm t(*tm::getNewTimer("Epetra_CrsMatrix PutScalar "));
        Amat->PutScalar(1.0);
      }

      //////////////////////////////////////////////////////////////
      // Replace vector values using view (not sure how Xyce does it)
      {
        tm t(*tm::getNewTimer("Epetra_MultiVector  ExtractView"));
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
        tm t(*tm::getNewTimer("Epetra_MultiVector  BracketOperator"));
        int len = x.MyLength();
        for (int k = 0; k < x.NumVectors(); k++) {
          for (int i = 0; i < len; i++) {
            x[k][i] *= 3.;
          }
        }
      }

      // Vector linear combination
      {
        tm t(*tm::getNewTimer("Epetra_MultiVector  Update a*A+b*this"));
        double a = 1.2;
        double b = 2.1;
        w.Update(a, x, b);
      }
      {
        tm t(*tm::getNewTimer("Epetra_MultiVector  Update a*A+b*B+c*this"));
        double a = 1.2;
        double b = 2.1;
        double c = 3.2;
        w.Update(a, x, b, z, c);
      }
    }
    delete [] newMatValues;
  }

  void testTpetra(const char* filename, int nIter)
  {
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
  if (comm->getSize() > 1) {
    std::cout << "Usage:  run on one processor only" << std::endl;
    return 0;
  }

  int niter = 10;
  if (narg == 3) niter = std::atoi(arg[2]);

  Test test(arg[1], niter);
  return 0;
}

