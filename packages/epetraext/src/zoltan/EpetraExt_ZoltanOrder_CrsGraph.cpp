//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#ifdef ZOLTAN_ORDER

#include <EpetraExt_ZoltanOrder_CrsGraph.h>

#include <EpetraExt_Transpose_CrsGraph.h>

#include <EpetraExt_ZoltanQuery.h>
#include <Zoltan_LoadBalance.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

#include <Epetra_MpiComm.h>

#include <vector>
#include <set>

using std::vector;
using std::set;

namespace EpetraExt {

EPETRAEXT_DEPRECATED
ZoltanOrder_CrsGraph::
~ZoltanOrder_CrsGraph()
{
  if( newObj_ ) delete newObj_;
  if( NewRowMap_ ) delete NewRowMap_;
}

EPETRAEXT_DEPRECATED
ZoltanOrder_CrsGraph::NewTypeRef
ZoltanOrder_CrsGraph::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  //Setup Load Balance Object
  float version;
  char * dummy = 0;
  Zoltan::LoadBalance LB( 0, &dummy, &version );
  err = LB.Create( dynamic_cast<const Epetra_MpiComm&>(orig.Comm()).Comm() );
  LB.Set_Param( "ORDER_METHOD", "METIS" );
  LB.Set_Param( "ORDER_TYPE", "LOCAL" );

  //Setup Query Object
  CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & TransGraph = transposeTransform( orig );
  ZoltanQuery Query( orig, &TransGraph, true );
  if( err == ZOLTAN_OK ) err = LB.Set_QueryObject( &Query );

  if( err != ZOLTAN_OK )
  { cout << "Setup of Zoltan Load Balancing Objects FAILED!\n"; exit(0); }

  //Generate Reorder
  int num_elements = orig.RowMap().NumMyElements();
  int num_gid_entries, num_lid_entries;
  vector<ZOLTAN_ID_TYPE> global_ids(num_elements,0);
  vector<ZOLTAN_ID_TYPE> local_ids(num_elements,0);
  vector<int> rank(num_elements);
  vector<int> iperm(num_elements);

  orig.Comm().Barrier();
  err = LB.Order( &num_gid_entries,
                   &num_lid_entries,
                   num_elements,
                   &global_ids[0],
                   &local_ids[0],
                   &rank[0],
                   &iperm[0] );
  orig.Comm().Barrier();

#ifdef EDT_ZOLTANORDER_DEBUG
cout << "------------------------------\n";
cout << "#GIDs: " << num_gid_entries << endl;
cout << "#LIDs: " << num_lid_entries << endl;
cout << "GIDs and LIDs\n";
for( int i = 0; i < num_elements; ++i ) cout << global_ids[i] << " " << local_ids[i] << endl;
cout << "Rank and Perm\n";
for( int i = 0; i < num_elements; ++i ) cout << rank[i] << " " << iperm[i] << endl;
cout << "------------------------------\n";
#endif

  //Generate New Row Map
  vector<int> gids(num_elements);
  for( int i = 0; i < num_elements; ++i )
    gids[ (num_elements-1) - rank[i] ] = global_ids[i];
  NewRowMap_ = new Epetra_Map( orig.RowMap().NumGlobalElements(),
                                           num_elements,
                                           &gids[0],
                                           orig.RowMap().IndexBase(),
                                           orig.RowMap().Comm() );

  //Create Importer
  Epetra_Import Importer( *NewRowMap_, orig.RowMap() );

  //Create New Graph
  Epetra_CrsGraph * NewGraph( new Epetra_CrsGraph( Copy, *NewRowMap_, 0 ) );
  NewGraph->Import( orig, Importer, Insert );
  NewGraph->FillComplete( true );
  
  newObj_ = NewGraph;

  return NewGraph;
}

} // namespace EpetraExt

#endif //ZOLTAN_ORDER
