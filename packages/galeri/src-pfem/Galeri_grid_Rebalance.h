/*
// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
*/

#ifndef GALERI_GRID_LOADBALANCE_H
#define GALERI_GRID_LOADBALANCE_H

#include <set>
#include <map>

#include "parmetis.h"

#include "Epetra_MpiComm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"

namespace Galeri {
namespace grid {

class Rebalance
{
public:
  static Epetra_Map*
  balanceMap(const Epetra_CrsGraph& graph)
  {
    assert (graph.Filled());

    assert (graph.RowMap().LinearMap());

    assert (graph.NumMyRows() == graph.NumMyBlockRows());

    const Epetra_MpiComm& comm = dynamic_cast<const Epetra_MpiComm&>(graph.Comm());

    const Epetra_BlockMap& blockRowMap = graph.RowMap();
    Epetra_Map rowMap(-1, blockRowMap.NumMyElements(), blockRowMap.MyGlobalElements(), 0, comm);

    //build datastructures for parmetis
    int numIndices;
    vector<int> xadj(graph.NumMyBlockRows() * 10 + 1);
    vector<int> adjacy(graph.NumMyNonzeros() * 10);

    int numMyBlockRows = graph.NumMyBlockRows();

    vector<int> vtxdst(comm.NumProc()+1 + 10);  
    vtxdst[0] = 0;
    int myMaxGid;
    comm.ScanSum(&numMyBlockRows, &myMaxGid, 1);
    comm.GatherAll(&myMaxGid, &(vtxdst[1]), 1);

    Epetra_Map linearMap(-1, numMyBlockRows, 0, comm);
    Epetra_IntVector linearRowVector(linearMap);
    for (int i = 0; i < numMyBlockRows; ++i)
      linearRowVector[i] = vtxdst[comm.MyPID()] + i;

    Epetra_IntVector linearColVector(graph.ColMap());
    if (comm.NumProc() > 1)
      linearColVector.Import(linearRowVector, *graph.Importer(), Insert);
    else
      linearColVector = linearRowVector;

    xadj[0] = 0;
    int length = graph.MaxNumIndices();
    vector<int> indices(length * 10);

    for (int i = 0; i < numMyBlockRows; ++i) 
    {
      int ierr;
      ierr = graph.ExtractMyRowCopy(i, length, numIndices, &indices[0]);

      TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::logic_error,
                         "graph.ExtractMyRowCopy() returned a negative value, "
                         << ierr);
      
      int count = 0;
      for (int j = 0; j < numIndices; ++j)
      {
        if (indices[j] != i)
          adjacy[xadj[i] + (count++)] = linearColVector[indices[j]];
      }

      xadj[i+1] = xadj[i] + count;
    }

    int zero = 0;
    int ncon = 1;
    int nparts = comm.NumProc();
    int edgecut;
    int options[3] = {1, 0, 7};
    vector<int> part(graph.NumMyBlockRows());
    vector<float> tpwgts(nparts);
    for (int i = 0; i < nparts; ++i)
      tpwgts[i] = 1.0 / nparts;
    float ubvec = 1.05;

    MPI_Comm MpiComm = comm.Comm();

    ParMETIS_V3_PartKway(&vtxdst[0], &xadj[0], &adjacy[0], NULL, NULL,
                         &zero, &zero, &ncon, &nparts, &tpwgts[0], &ubvec,
                         options, &edgecut, &part[0], &MpiComm);

    Epetra_Map procMap(-1, 1, 0, comm);
    Epetra_FECrsGraph procGraph(Copy, procMap, 10); // FIXME

    for (int i = 0; i < numMyBlockRows; ++i)
    {
      int owner = part[i];
      int GID = rowMap.GID(i);
      procGraph.InsertGlobalIndices(1, &owner, 1, &GID);
    }

    procGraph.GlobalAssemble(rowMap, procMap, false);

    int numMyElements;
    int* myGlobalElements;

    procGraph.ExtractGlobalRowView(procMap.GID(0), numMyElements, myGlobalElements);

    Epetra_Map* balancedMap = new Epetra_Map(-1, numMyElements, myGlobalElements, 0, comm);

    return(balancedMap);
  }

  static Epetra_FECrsGraph*
  balanceGraph(const Epetra_CrsGraph& graph)
  {
    Epetra_Map* newMap = balanceMap(graph);
    Epetra_Import importer(*newMap, graph.Map());
    Epetra_FECrsGraph* newGraph = new Epetra_FECrsGraph(Copy, *newMap, 0);

    newGraph->Import(graph, importer, Insert);

    newGraph->FillComplete();

    return(newGraph);
  }

private:

}; // class Rebalance

} // namespace grid
} // namespace Galeri

#endif
