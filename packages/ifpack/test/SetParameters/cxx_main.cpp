//@HEADER
// ************************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

// SetParameters Test routine
#include <Ifpack_ConfigDefs.h>
#include <Ifpack_IlukGraph.h>
#include <Ifpack_CrsRiluk.h>
#include <Ifpack_CrsIct.h>
#include <Ifpack_OverlapGraph.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#ifdef HAVE_IFPACK_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#endif

#include <ifp_parameters.h>
#include <Epetra_SerialComm.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif

int main(int argc, char* argv[]) {
  bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;

  int err;
  int returnierr = 0;
  int size = 1;
  int rank = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (verbose1) {
  }

#ifdef HAVE_IFPACK_TEUCHOS
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

  Ifpack::set_parameters(paramlist, params);

  if (params.double_params[Ifpack::absolute_threshold] != 44.0) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set absolute_threshold."<<endl;
    }
    return(-1);
  }
#endif

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
  Ifpack_CrsRiluk crsiluk(ilukgraph);
  Ifpack_OverlapGraph overlapgraph(&graph, 1);

  Epetra_CrsMatrix A(Copy, graph);

  for(i=0; i<local_n; ++i) {
    int row = first_global_row + i;
    double val = 2.0;
    A.SumIntoGlobalValues(row, 1, &val, &row);
  }

  Ifpack_CrsIct crsict(A, 1.0, 1);

#ifdef HAVE_IFPACK_TEUCHOS
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
#endif

#ifdef HAVE_IFPACK_TEUCHOS
  crsiluk.SetParameters(paramlist);

  double athresh = crsiluk.GetAbsoluteThreshold();
  if (athresh != 44.0) {
    cerr << "SetParameters test failed to correctly set absolute_threshold."
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

