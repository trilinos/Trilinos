#ifndef MLAPI_MATLABSTREAM_H
#define MLAPI_MATLABSTREAM_H

#include "MLAPI_Operator.h"

namespace MLAPI {

/*!
\class MATLABStream

\brief Basic stream to save in a MATLAB-compatible file MLAPI objects.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.
*/

class MATLABStream
{
public:

  // @{ Constructos and destructors
  
  //! Opens the specified file for writing.
  MATLABStream(const string& FileName, bool UseSparse = true) 
  {
    FileName_ = FileName;
    if (GetMyPID() == 0) {
      fp_ = fopen(FileName_.c_str(),"w");
      fclose(fp_);
    }
  }

  //! Destructor.
  ~MATLABStream()
  { }

  // @}
  // @{ Overloaded operators

  //! Writes on file input integer on process 0 only.
  MATLABStream& operator << (const int i)
  {
    if (GetMyPID() == 0)
      fprintf(fp_,"%d",i);
    return(*this);
  }

  //! Writes on file input string on process 0 only.
  MATLABStream& operator << (const string s)
  {
    if (GetMyPID() == 0)
      fprintf(fp_,"%s",s.c_str());
    return(*this);
  }

  //! Writes on file input Operator, one process at-a-time, using global ordering.
  MATLABStream& operator << (const Operator& Op)
  {
    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = Op.GetML_Operator();

    if (matrix->getrow == NULL) {
      cerr << "ERROR: In MATLABStream::operator<<" << endl;
      cerr << "ERROR: (file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      cerr << "ERROR: getrow() not set!" << endl;
      throw(-1);
    }

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    int NumGlobalRows = Op.GetDomainSpace().GetNumGlobalElements();
    int NumGlobalCols = Op.GetRangeSpace().GetNumGlobalElements();

    if (GetMyPID() == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d,%d);\n",
              Op.GetLabel().c_str(), NumGlobalRows, NumGlobalCols);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = Op.GetDomainSpace()(i) + 1;
            int GlobalCol = Op.GetColumnSpace()(bindx[j]) + 1;
            fprintf(fp_,"%s(%d,%d) = %e;\n",
                    Op.GetLabel().c_str(), GlobalRow, GlobalCol, val[j]);
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

  //! Writes on file input MultiVector, one process at-a-time.
  MATLABStream& operator << (const MultiVector& v)
  {
    int NumMyRows = v.GetVectorSpace().GetNumMyElements();
    int NumGlobalRows = v.GetVectorSpace().GetNumGlobalElements();

    if (GetMyPID() == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d);\n",
              v.GetLabel().c_str(), NumGlobalRows);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < NumMyRows ; ++i)
          fprintf(fp_,"%s(%d) = %e;\n",
                  v.GetLabel().c_str(), v.GetVectorSpace()(i) + 1, v(i));
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
    int NumMyRows = S.GetNumMyElements();
    int NumGlobalRows = S.GetNumGlobalElements();

    if (GetMyPID() == 0) {
      fp_ = fopen(FileName_.c_str(),"a");
      fprintf(fp_,"%s = zeros(%d);\n",
              S.GetLabel().c_str(), NumGlobalRows);
      fclose(fp_);
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        fp_ = fopen(FileName_.c_str(),"a");

        for (int i = 0 ; i < NumMyRows ; ++i)
          fprintf(fp_,"%s(%d) = %d;\n",
                  S.GetLabel().c_str(), i, S(i) + 1);
        fclose(fp_);
      }
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    return (*this);
  }

private:

  // @}
  // @{ Internal data
  
  //! Name of output file.
  string FileName_;
  //! FILE pointer.
  FILE* fp_;

  // @}
  
}; // class MATLABStream

} // namespace MLAPI

#endif
