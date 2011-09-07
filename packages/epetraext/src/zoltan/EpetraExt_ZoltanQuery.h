//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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

