#ifndef MLAPI_MATLABSTREAM_H
#define MLAPI_MATLABSTREAM_H

#include "MLAPI_Operator.h"

namespace MLAPI {

/*!
\class MATLABStream

\brief Basic stream to save in a MATLAB-compatible file MLAPI objects.

Class MATLABStream allows the user to save in a file MLAPI objects, 
in a format that MATLAB can understand. Both serial and distributed objects 
can be save in just one file. MATLABStream behaves like a "normal" stream. 
Objects can be saved in file by using the operator << as in the following example:
\code
MLAPI::Space MySpace(5);
MLAPI::DoubleVector V(MySpace);
V = 2;
MLAPI::Operator A(MySpace,MySpace,EpetraMatrix);

MLAPI::MATLABStream Stream;

Stream << "% MySpace" << endl;
Stream << MySpace;
Stream << "% vector" << endl;
Stream << V;
Stream << "% operator" << endl;
Stream << A;
\endcode

It is important to give a name to each object, in order to make the MATLAB file more readable. This can be done as follows:
\code
MySpace = "space";
V = "vector";
A = "A";
\endcode
\e before writing a given object to MATLABStream. Names must be valid MATLAB namnes (for instance, without spaces);

\author Marzio Sala, SNL 9214

\date Last updated on 07-Jan-05
*/

class MATLABStream
{
public:

  //! Opens the specified file for writing.
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

  //! Destructor.
  ~MATLABStream()
  { }

  //! Writes on file input integer on process 0 only.
  MATLABStream& operator << (const int i)
  {
    if (MyPID_ == 0)
      fprintf(fp_,"%d",i);
    return(*this);
  }

  //! Writes on file input string on process 0 only.
  MATLABStream& operator << (const string s)
  {
    if (MyPID_ == 0)
      fprintf(fp_,"%s",s.c_str());
    return(*this);
  }

  //! Writes on file input Operator, one process at-a-time, using global ordering.
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

  //! Writes on file input DoubleVector, one process at-a-time.
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

  //! Writes on file input Space, one process at-a-time.
  MATLABStream& operator << (const Space& S)
  {
    int NumMyRows = S.NumMyElements();
    int NumGlobalRows = S.NumGlobalElements();

    if (MyPID_ == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d);\n",
              S.Name().c_str(), NumGlobalRows);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < NumProc_ ; ++iproc) {

      if (MyPID_ == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < NumMyRows ; ++i)
          fprintf(fp_,"%s(%d) = %d;\n",
                  S.Name().c_str(), i, S(i) + 1);
        fclose(fp_);
      }
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    return (*this);
  }

private:

  //! ID of calling process.
  int MyPID_;
  //! Total number of available processes.
  int NumProc_;
  //! Name of output file.
  string FileName_;
  //! FILE pointer.
  FILE* fp_;

};

} // namespace MLAPI

#endif
