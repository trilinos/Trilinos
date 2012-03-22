
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
//@HEADER

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Util.h>

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#endif

namespace epetra_test {

bool global_check_for_flag_on_proc_0(const char* flag,
                                     int numargs,
                                     char** strargs,
                                     const Epetra_Comm& comm)
{
  int mypid = comm.MyPID();
  int numprocs = comm.NumProc();

  int flag_found = 0;
  if (mypid==0) {
    for(int i=0; i<numargs; ++i) {
      if (strargs[i]==0) continue;

      if (strcmp(flag, strargs[i]) == 0) {
        flag_found = 1; break;
      }
    }
  }

  if (numprocs > 1) {
    comm.Broadcast(&flag_found, 1, 0);
  }

  bool return_value = flag_found==1 ? true : false;

  return( return_value );
}

Epetra_Comm* create_comm(int argc, char** argv)
{
#ifdef EPETRA_MPI
  MPI_Init(&argc, &argv);

  return( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  return( new Epetra_SerialComm );
#endif
}

bool compare_matrices(const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B)
{
  const Epetra_Map& Amap = A.RowMap();
  const Epetra_Map& Bmap = B.RowMap();

  if (!Amap.PointSameAs(Bmap)) {
    return(false);
  }

  int numRows = Amap.NumMyElements();

  if(Amap.GlobalIndicesInt())
  {
    int* rows = Amap.MyGlobalElements();
 
    Epetra_Util util;
 
    for(int i=0; i<numRows; ++i) {
      int row = rows[i];
      int rowLen = A.NumGlobalEntries(row);
      if (rowLen != B.NumGlobalEntries(row)) {
        return(false);
      }
 
      int* indices = new int[rowLen*2];
      int* Bindices = indices+rowLen;
 
      double* values = new double[rowLen*2];
      double* Bvalues = values+rowLen;
 
      A.ExtractGlobalRowCopy(row, rowLen, rowLen, values, indices);
      B.ExtractGlobalRowCopy(row, rowLen, rowLen, Bvalues, Bindices);
 
      util.Sort(true, rowLen, indices, 1, &values, 0, 0, 0, 0);
      util.Sort(true, rowLen, Bindices, 1, &Bvalues, 0, 0, 0, 0);
 
      bool same = true;
      for(int j=0; j<rowLen; ++j) {
        if (indices[j] != Bindices[j]) {
          same = false; break;
        }
        if (values[j] != Bvalues[j]) {
          same = false; break;
        }
      }
 
      delete [] indices;
      delete [] values;
 
      if (!same) {
        return(false);
      }
    }
  }
  else if(Amap.GlobalIndicesLongLong()) {
    long long* rows = Amap.MyGlobalElements_LL();
 
    Epetra_Util util;
 
    for(int i=0; i<numRows; ++i) {
      long long row = rows[i];
      int rowLen = A.NumGlobalEntries(row);
      if (rowLen != B.NumGlobalEntries(row)) {
        return(false);
      }
 
      long long* indices = new long long[rowLen*2];
      long long* Bindices = indices+rowLen;
 
      double* values = new double[rowLen*2];
      double* Bvalues = values+rowLen;
 
      A.ExtractGlobalRowCopy(row, rowLen, rowLen, values, indices);
      B.ExtractGlobalRowCopy(row, rowLen, rowLen, Bvalues, Bindices);
 
      util.Sort(true, rowLen, indices, 1, &values, 0, 0, 0, 0);
      util.Sort(true, rowLen, Bindices, 1, &Bvalues, 0, 0, 0, 0);
 
      bool same = true;
      for(int j=0; j<rowLen; ++j) {
        if (indices[j] != Bindices[j]) {
          same = false; break;
        }
        if (values[j] != Bvalues[j]) {
          same = false; break;
        }
      }
 
      delete [] indices;
      delete [] values;
 
      if (!same) {
        return(false);
      }
    }
  }
  else {
    return(false);
  }
  

  return(true);
}

bool compare_matrices_LL(const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B)
{
  const Epetra_Map& Amap = A.RowMap();
  const Epetra_Map& Bmap = B.RowMap();

  if (!Amap.PointSameAs(Bmap)) {
    return(false);
  }

  int numRows = Amap.NumMyElements();
  return(true);
}

}//namespace epetra_test

