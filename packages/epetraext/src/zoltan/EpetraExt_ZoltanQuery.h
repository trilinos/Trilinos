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

#ifndef EPETRAEXT_ZOLTANQUERY_H
#define EPETRAEXT_ZOLTANQUERY_H

#include "EpetraExt_ConfigDefs.h"

#include <Zoltan_QueryObject.h>

#include <vector>

class Epetra_CrsGraph;

namespace EpetraExt {

///
/** Query helper object to be used form Zoltan partitioning/ordering.
 * This object allows Zoltan to query an Epetra_CrsGraph object for it's
 * partitioning/ordering algorithms
 */
class EPETRAEXT_DEPRECATED ZoltanQuery : public Zoltan::QueryObject
{

  const Epetra_CrsGraph & graph_;
  const Epetra_CrsGraph * tgraph_;

  std::vector< std::vector<int> > LBProc_;
  std::vector< std::vector<int> > LBProc_Trans_;

  const bool localEdgesOnly_;

 public:

  ///
  /** Constructor
   */
  ZoltanQuery( const Epetra_CrsGraph & graph,
               const Epetra_CrsGraph * tgraph = 0,
               bool localEdgesOnly = false );

  //General Functions
  int Number_Objects  ( void * data,
                        int * ierr );

  void Object_List    ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_ids,
                        ZOLTAN_ID_PTR local_ids,
                        int weight_dim,
                        float * object_weights,
                        int * ierr );

  //Graph Based Functions
  int Number_Edges    ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        int * ierr );
  
  void Edge_List      ( void * data,
                        int num_gid_entries,
                        int num_lid_entries,
                        ZOLTAN_ID_PTR global_id,
                        ZOLTAN_ID_PTR local_id,
                        ZOLTAN_ID_PTR neighbor_global_ids,
                        int * neighbor_procs,
                        int weight_dim,
                        float * edge_weights,
                        int * ierr );
  
  //HyperGraph Functions
  int Number_HG_Edges ( void * data,
    		        int * ierr );

  int Number_HG_Pins ( void * data,
		        int * ierr );

  int HG_Edge_List   ( void * data,
                       int num_gid_entries,
                       int ewgt_dim,
                       int nedge,
                       int maxsize,
                       int * edge_sizes,
                       ZOLTAN_ID_PTR edge_verts,
                       int * edge_procs,
                       float * edge_weights );

};

} //namespace EpetraExt

#endif //EPETRAEXT_ZOLTANQUERY_H

