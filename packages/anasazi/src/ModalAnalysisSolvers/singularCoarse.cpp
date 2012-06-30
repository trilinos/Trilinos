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


