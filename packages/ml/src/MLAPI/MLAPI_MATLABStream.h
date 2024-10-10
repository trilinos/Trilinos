#ifndef MLAPI_MATLABSTREAM_H
#define MLAPI_MATLABSTREAM_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_MATLABStream.h

\brief Basic stream to save in a MATLAB-compatible file MLAPI objects.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"

#include "MLAPI_Error.h"
#include "MLAPI_Operator.h"

namespace MLAPI {

/*!
\class MATLABStream

\brief Basic stream to save in a MATLAB-compatible file MLAPI objects.

For an example of usage, see \ref ml_blackboard_cpp

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.
*/

class MATLABStream
{
public:

  // @{ \name Constructors and destructors.

  //! Opens the specified file for writing.
  MATLABStream(const std::string& FileName, bool UseSparse = true)
  {
    FileName_ = FileName;
    SetUseSparse(UseSparse);

    // Creates a new stream, deletes old stream with the same name.
    if (GetMyPID() == 0) {
      Open(true);
      fprintf(fp_, "%% beginning of MLAPI::MATLABStream\n");
      Close();
    }
  }

  //! Finally closes the output file.
  ~MATLABStream()
  {
    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_, "%% end of MLAPI::MATLABStream\n");
      Close();
    }
  }

  // @}
  // @{ \name Overloaded operators

  //! Writes on file the specified integer (on process 0 only).
  MATLABStream& operator << (const int obj)
  {
    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_,"%d",obj);
      Close();
    }
    return(*this);
  }

  //! Writes on file the specified double (on process 0 only).
  MATLABStream& operator << (const double obj)
  {
    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_,"%f",obj);
      Close();
    }
    return(*this);
  }

  //! Writes on file the specified std::string on process 0 only.
  MATLABStream& operator << (const std::string obj)
  {
    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_,"%s",obj.c_str());
      Close();
    }
    return(*this);
  }

  //! Writes on file input Operator, one process at-a-time, using global ordering.
  MATLABStream& operator << (const Operator& obj)
  {
    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = obj.GetML_Operator();

    if (matrix->getrow == NULL)
      ML_THROW("getrow() not set", -1);

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    int NumGlobalRows = obj.GetDomainSpace().GetNumGlobalElements();
    int NumGlobalCols = obj.GetRangeSpace().GetNumGlobalElements();

    if (GetMyPID() == 0) {
      Open();
      if (GetUseSparse())
        fprintf(fp_,"%s = sparse(%d,%d);\n",
                obj.GetLabel().c_str(), NumGlobalRows, NumGlobalCols);
      else
        fprintf(fp_,"%s = zeros(%d,%d);\n",
                obj.GetLabel().c_str(), NumGlobalRows, NumGlobalCols);

      Close();
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        Open();

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = obj.GetGRID(i) + 1;
            int GlobalCol = obj.GetGCID(bindx[j]) + 1;
            fprintf(fp_,"%s(%d,%d) = %e;\n",
                    obj.GetLabel().c_str(), GlobalRow, GlobalCol, val[j]);
          }
        }

        Close();
      }
      Barrier();
    }

    ML_free(val);
    ML_free(bindx);
    return (*this);
  }

  //! Writes on file the input MultiVector, one process at-a-time.
  MATLABStream& operator << (const MultiVector& obj)
  {
    int NumMyRows = obj.GetVectorSpace().GetNumMyElements();
    int NumGlobalRows = obj.GetVectorSpace().GetNumGlobalElements();
    int NumVectors = obj.GetNumVectors();

    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_,"%s = zeros(%d, %d);\n",
              obj.GetLabel().c_str(), NumGlobalRows, NumVectors);
      Close();
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        Open();

        for (int i = 0 ; i < NumMyRows ; ++i) {
          for (int j = 0 ; j < NumVectors ; ++j) {
            fprintf(fp_,"%s(%d, %d) = %e;\n",
                    obj.GetLabel().c_str(), obj.GetVectorSpace()(i) + 1, j + 1, obj(i,j));
          }
        }
        Close();
      }
      Barrier();
    }

    return (*this);
  }

  //! Writes on file input Space, one process at-a-time.
  MATLABStream& operator << (const Space& obj)
  {
    int NumMyRows = obj.GetNumMyElements();
    int NumGlobalRows = obj.GetNumGlobalElements();

    if (GetMyPID() == 0) {
      Open();
      fprintf(fp_,"%s = zeros(%d);\n",
              obj.GetLabel().c_str(), NumGlobalRows);
      Close();
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        Open();

        for (int i = 0 ; i < NumMyRows ; ++i)
          fprintf(fp_,"%s(%d) = %d;\n",
                  obj.GetLabel().c_str(), i, obj(i) + 1);
        Close();
      }
      Barrier();
    }

    return (*this);
  }

  // @}
  // @{ \name Sets and gets methods

  //! Returns \c true if the stream uses sparse MATLAB format.
  bool GetUseSparse() const
  {
    return(UseSparse_);
  }

  //! Toggles the use of sparse MATLAB formats.
  void SetUseSparse(const bool UseSparse)
  {
    UseSparse_ = UseSparse;
  }

  //! Returns the name of the output file.
  inline std::string GetFileName() const
  {
    return(FileName_);
  }

  //@}

private:

  //! Opens the file stream in append mode, or in write more if \c FirstTime == \c true.
  void Open(const bool FirstTime = false)
  {
    if (FirstTime)
      fp_ = fopen(FileName_.c_str(),"w");
    else
      fp_ = fopen(FileName_.c_str(),"a");
  }

  //! Closes the file stream.
  void Close()
  {
    fclose(fp_);
  }

  //! Name of output file.
  std::string FileName_;
  //! If \c true, prints out using sparse MATLAB commands.
  bool UseSparse_;
  //! FILE pointer.
  FILE* fp_;

}; // class MATLABStream

} // namespace MLAPI

#endif
