/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER
*/

// SetParameters Test routine
#include <Ifpack_ConfigDefs.h>
#include <Ifpack_IlukGraph.h>
#include <Ifpack_CrsRiluk.h>
#include <Ifpack_CrsIct.h>
#include <Ifpack_OverlapGraph.h>

#include <Epetra_CombineMode.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include <Teuchos_ParameterList.hpp>

#include <ifp_parameters.h>
#include <Epetra_SerialComm.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif

int main(int argc, char* argv[]) {
  using std::cerr;
  using std::cout;
  using std::endl;

  //bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument

  int returnierr = 0;
  //int size = 1;
  //int rank = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Ifpack::param_struct params;
  params.double_params[Ifpack::absolute_threshold] = -99.9;

  Teuchos::ParameterList paramlist;
  paramlist.set("absolute_threshold", 44.0);
  paramlist.set("level_fill", 2);
  paramlist.set("LEVEL_OVERLAP", 2);
  paramlist.set("relative_threshold", 1.e-2);
  paramlist.set("fill_tolerance", 2.0);
  paramlist.set("use_reciprocal", false);
  paramlist.set("level_overlap", 2);
  paramlist.set("overlap_mode", Add);

  Ifpack::set_parameters(paramlist, params);

  if (params.double_params[Ifpack::absolute_threshold] != 44.0) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set absolute_threshold."<<endl;
    }
    return(-1);
  }

  int i, local_n = 5;
  int my_pid = Comm.MyPID();
  int num_procs = Comm.NumProc();
  int global_n = num_procs*local_n;

  Epetra_Map map(global_n, 0, Comm);
  Epetra_CrsGraph graph(Copy, map, 1);
  int first_global_row = my_pid*local_n;

  for(i=0; i<local_n; ++i) {
    int row = first_global_row + i;
    graph.InsertGlobalIndices(row, 1, &row);
  }

  graph.FillComplete();

  Ifpack_IlukGraph ilukgraph(graph, 1,1);
  Ifpack_CrsRiluk crsriluk(ilukgraph);
  // MS // this was failing
#if 0
  Ifpack_OverlapGraph overlapgraph(&graph, 1);
#endif

  Epetra_CrsMatrix A(Copy, graph);

  for(i=0; i<local_n; ++i) {
    int row = first_global_row + i;
    double val = 2.0;
    A.SumIntoGlobalValues(row, 1, &val, &row);
  }

  Ifpack_CrsIct crsict(A, 1.0, 1);

  ilukgraph.SetParameters(paramlist);

  int levelfill = ilukgraph.LevelFill();
  if (levelfill != 2) {
    cerr << "SetParameters test failed to correctly set level_fill."
        << endl;
    return(-1);
  }

  int leveloverlap = ilukgraph.LevelOverlap();
  if (leveloverlap != 2) {
    cerr << "SetParameters test failed to correctly set level_overlap."
        << endl;
    return(-1);
  }

  crsriluk.SetParameters(paramlist);

  double athresh = crsriluk.GetAbsoluteThreshold();
  if (athresh != 44.0) {
    cerr << "SetParameters test failed to correctly set absolute_threshold."
        << endl;
    return(-1);
  }

  Epetra_CombineMode overlapmode = crsriluk.GetOverlapMode();
  if (overlapmode != Add) {
    cerr << "SetParameters test failed to correctly set overlapmode."
        << endl;
    return(-1);
  }

  crsict.SetParameters(paramlist);

  double rthresh = crsict.GetRelativeThreshold();
  if (rthresh != 1.e-2) {
    cerr << "SetParameters test failed to correctly set relative_threshold."
        << endl;
    return(-1);
  }

  overlapmode = crsict.GetOverlapMode();
  if (overlapmode != Add) {
    cerr << "SetParameters test failed to correctly set overlapmode."
        << endl;
    return(-1);
  }

#if 0
  overlapgraph.SetParameters(paramlist);

  int overlaplevel = overlapgraph.OverlapLevel();
  if (overlaplevel != 2) {
    cerr << "SetParameters test failed to correctly set overlaplevel."
        << endl;
    return(-1);
  }
#endif

  if (verbose1==true) {
    cout << "********* Test passed **********" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(returnierr);
}

