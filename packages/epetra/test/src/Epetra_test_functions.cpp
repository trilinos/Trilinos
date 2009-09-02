
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

    util.Sort(true, rowLen, indices, 1, &values, 0, 0);
    util.Sort(true, rowLen, Bindices, 1, &Bvalues, 0, 0);

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

  return(true);
}

}//namespace epetra_test

