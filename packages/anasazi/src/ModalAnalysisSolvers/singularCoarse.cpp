// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER
#include "singularCoarse.h"

namespace singularCoarse {

  void setNullSpace(double *V, int row, int col, double *VtV, const Epetra_Comm *_Comm) {

    Qcoarse = V;
    rowQcoarse = row;
    colQcoarse = col;
    QcoarseTQcoarse = VtV; 
    commCoarse = _Comm;

  }

  void projection(double *z, int *options, int *proc_config, double *params,
                  AZ_MATRIX_STRUCT *Amat, AZ_PREC_STRUCT *prec) {

    Epetra_BLAS callBLAS;
    Epetra_LAPACK callLAPACK;

//    // Do one Jacobi step
//    if (Amat->data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
//      if (commCoarse->MyPID() == 0)
//        cout << " Do ONE JACOBI STEP " << endl;
//      for (int i = 0; i < rowQcoarse; ++i)
//        z[i] /= Amat->val[i]; 
//    }

    double *tmp = new double[2*colQcoarse];
    memset(tmp, 0, 2*colQcoarse*sizeof(double));

    int info = 0;

    for (int i = 0; i < 2; ++i) {
      if (rowQcoarse > 0) {
        callBLAS.GEMV('T',rowQcoarse,colQcoarse,1.0,Qcoarse,rowQcoarse,z,0.0,tmp+colQcoarse);
      }
      commCoarse->SumAll(tmp + colQcoarse, tmp, colQcoarse);
      if (rowQcoarse > 0) {
        callLAPACK.POTRS('U', colQcoarse, 1, QcoarseTQcoarse, colQcoarse, tmp, colQcoarse, &info); 
        callBLAS.GEMV('N', rowQcoarse, colQcoarse, -1.0, Qcoarse, rowQcoarse, tmp, 1.0, z);
      }
    }

    delete[] tmp;

  }

}


