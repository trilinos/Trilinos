#ifndef MLAPI_MATLABSTREAM_H
#define MLAPI_MATLABSTREAM_H

#include "MLAPI_Operator.h"

namespace MLAPI {

class MATLABStream
{
public:

  MATLABStream(const string& FileName, bool UseSparse = true) 
  {
    MyPID_ = GetMLComm()->ML_mypid;
    NumProc_ = GetMLComm()->ML_nprocs;
    FileName_ = FileName;
    if (MyPID_ == 0) {
      fp_ = fopen(FileName_.c_str(),"w");
      fclose(fp_);
    }
  }

  ~MATLABStream()
  { }

  MATLABStream& operator << (const int i)
  {
    if (MyPID_ == 0)
      fprintf(fp_,"%d",i);
    return(*this);
  }

  MATLABStream& operator << (const string s)
  {
    if (MyPID_ == 0)
      fprintf(fp_,"%s",s.c_str());
    return(*this);
  }

  MATLABStream& operator << (const Operator& Op)
  {
    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = Op.GetOperator();

    if (matrix->getrow == NULL) 
      throw("getrow not set");

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    int NumGlobalRows = Op.DomainSpace().NumGlobalElements();
    int NumGlobalCols = Op.RangeSpace().NumGlobalElements();

    if (MyPID_ == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d,%d);\n",
              Op.Name().c_str(), NumGlobalRows, NumGlobalCols);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < NumProc_ ; ++iproc) {

      if (MyPID_ == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = Op.DomainSpace()(i) + 1;
            int GlobalCol = Op.ColumnSpace()(bindx[j]) + 1;
            fprintf(fp_,"%s(%d,%d) = %e;\n",
                    Op.Name().c_str(), GlobalRow, GlobalCol, val[j]);
          }
        }
        fclose(fp_);
      }
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    ML_free(val);
    ML_free(bindx);
    return (*this);
  }

  MATLABStream& operator << (const DoubleVector& v)
  {
    int NumMyRows = v.VectorSpace().NumMyElements();
    int NumGlobalRows = v.VectorSpace().NumGlobalElements();

    if (MyPID_ == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d);\n",
              v.Name().c_str(), NumGlobalRows);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < NumProc_ ; ++iproc) {

      if (MyPID_ == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < NumMyRows ; ++i)
          fprintf(fp_,"%s(%d) = %e;\n",
                  v.Name().c_str(), v.VectorSpace()(i) + 1, v(i));
        fclose(fp_);
      }
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    return (*this);
  }

private:

  int MyPID_;
  int NumProc_;
  string FileName_;
  FILE* fp_;

};

}

#endif
