/*
  This file will contain all the application-supplied query functions
  that zoltan needs. These functions actually make calls to methods on
  ports. These pointers to ports are used here but initialized elsewhere in
  ::ZoltanSpace::PartitionerFactory_JR::init() 

  This file also has the following methods

  a) ::ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesEntityInfoPort() 
  b) ::ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesGeomInfoPort()
  c) ::ZoltanSpace::PartitionerFactory_JR::SetZoltanQueryEdgeInfoPort()
  d) ::ZoltanSpace::PartitionerFactory_JR::SetZoltanQueryTreeInfoPort()
  which must be called by the different load-balancers depending upon
  which load-balancer needs what.

  e) ::ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesDataMigrationPort();
  which is called by the partitionerfactory itself to set the query functions
  for Zoltan that do  the pre-, in- and post-migration processing, as well
  as packing and unpacking of elements.
*/

#include "include/zoltan.h"
#include "PartitionerFactory.h"

namespace ZoltanSpace
{
  namespace Globals
  {
    ::LoadPartitionerSpace::EntityInfo *pEntityInfo = 0;
    ::LoadPartitionerSpace::GeomInfo *pGeomInfo = 0;
    ::LoadPartitionerSpace::EdgeInfo *pEdgeInfo = 0;
    ::LoadPartitionerSpace::TreeInfo *pTreeInfo = 0;
    ::LoadPartitionerSpace::DataMigrationPort *pDataMigPort = 0 ;
  } ;
} ;
    
int ZC_JR_Num_Obj_Fn(void *data, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pEntityInfo->GetNumberOfEntities() ; 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return (res) ;
}

void ZC_JR_Obj_List_Fn(void *data, int num_gid_entries, int num_lid_entries,
		      ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim,
		      float *wgts, int *ierr)
{
  bool wgt_reqd = (wgt_dim == 1) ? true : false ;
  int err = ::ZoltanSpace::Globals::pEntityInfo->EntitiesOnThisProc(
 					  num_gid_entries, num_lid_entries, 
					  global_ids,local_ids, wgt_reqd, wgts) ;
  *ierr = (err < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;
}

int ZC_JR_First_Obj_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR first_gid,
		       ZOLTAN_ID_PTR first_lid, int wgt_dim, float *wgt, int *ierr)
{
  bool wgt_reqd = (wgt_dim == 1) ? true : false ;
  int err = ::ZoltanSpace::Globals::pEntityInfo->GetFirstEntity(
                                          num_gid, num_lid, first_gid, first_lid, 
					  wgt_reqd, wgt) ;
  *ierr = (err < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return( (*ierr == 0) ? 1 : 0 ) ;
}

int ZC_JR_Next_Obj_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR prev_gid,
		      ZOLTAN_ID_PTR prev_lid, ZOLTAN_ID_PTR next_gid,
		      ZOLTAN_ID_PTR next_lid, int wgt_dim, float *wgt, int *ierr)
{
  bool wgt_reqd = (wgt_dim == 1) ? true : false ;
  int err = ::ZoltanSpace::Globals::pEntityInfo->GetNextEntity(
                                          num_gid, num_lid, prev_gid, prev_lid, 
					  next_gid, next_lid, wgt_reqd, wgt) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return( (*ierr == 0) ? 1 : 0 ) ;
}

int ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesEntityInfoPort(
	  			            struct Zoltan_Struct *zz)
{
  if ( pEntityInfo == 0 ) return (-1) ;
  ::ZoltanSpace::Globals::pEntityInfo = pEntityInfo ;
  Zoltan_Set_Num_Obj_Fn(zz, &(ZC_JR_Num_Obj_Fn), NULL) ;

  if ( pIterate->value == false )
    Zoltan_Set_Obj_List_Fn(zz, &(ZC_JR_Obj_List_Fn), NULL) ;
  else
  {
    Zoltan_Set_First_Obj_Fn(zz, &(ZC_JR_First_Obj_Fn), NULL) ;
    Zoltan_Set_Next_Obj_Fn(zz, &(ZC_JR_Next_Obj_Fn), NULL) ;
  }
  return(0) ;
}
//-------------------------------------------------------------------

int ZC_JR_Num_Dimension_Fn(void *data, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pGeomInfo->GetNumberOfDims() ; 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return(res);
}

void ZC_JR_Coords_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
		    ZOLTAN_ID_PTR lid, double *vec, int *ierr)
{
  int err = ::ZoltanSpace::Globals::pGeomInfo->GetGeomDescription(
                                          num_gid, num_lid, gid, lid, vec);
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
}

int ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesGeomInfoPort(
					    struct Zoltan_Struct *zz)
{
  if ( pGeomInfo  == 0 ) return (-1) ;

  ::ZoltanSpace::Globals::pGeomInfo = pGeomInfo ;
  Zoltan_Set_Num_Geom_Fn(zz, &(ZC_JR_Num_Dimension_Fn), NULL) ;
  Zoltan_Set_Geom_Fn(zz, &(ZC_JR_Coords_Fn), NULL) ;

  return(0) ;
}
//---------------------------------------------------------------------------

int ZC_JR_Num_Edges_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
		       ZOLTAN_ID_PTR lid, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pEdgeInfo->GetNumberOfEdges( num_gid, 
							     num_lid, gid, lid);
  *ierr = ( res < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return(res) ;
}

void ZC_JR_Edge_List_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
			ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR nbor_gid, int *nbor_procs,
			int wgt_dim, float *wgt, int *ierr)
{
  bool wgt_reqd = (wgt_dim == 1) ? true : false ;
  int err = ::ZoltanSpace::Globals::pEdgeInfo->EdgesListForThisEntity(
                                          num_gid, num_lid, gid, lid, nbor_gid, 
					  nbor_procs, wgt_reqd, wgt) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
}

int ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesEdgeInfoPort(
                                            struct Zoltan_Struct *zz)
{
  if ( pEdgeInfo == 0) return (-1) ;

  ::ZoltanSpace::Globals::pEdgeInfo = pEdgeInfo ;

  Zoltan_Set_Num_Edges_Fn(zz, &(ZC_JR_Num_Edges_Fn), NULL);
  Zoltan_Set_Edge_List_Fn(zz, &(ZC_JR_Edge_List_Fn), NULL);

  return(0);
}
//-----------------------------------------------------------------------------
  
int ZC_JR_Num_Coarse_Elem_Fn(void *data, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pTreeInfo->GetNumberOfCoarseElements() ; 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return(res) ;
}

void ZC_JR_Coarse_Obj_List_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
			      ZOLTAN_ID_PTR lid, int *on_proc, int *num_vert, 
			      ZOLTAN_ID_TYPE *vertices, int *in_order, 
			      ZOLTAN_ID_TYPE *in_vertex, 
			      ZOLTAN_ID_TYPE *out_vertex, 
			      int *ierr)
{
  int num_obj = ::ZoltanSpace::Globals::pTreeInfo->GetNumberOfCoarseElements();
  bool *assigned = new bool [num_obj] ; 
  bool order ;
  int err = ::ZoltanSpace::Globals::pTreeInfo->CoarseMeshEntitiesList(
                                          num_gid, num_lid, gid, lid, assigned, 
					  num_vert, vertices, &order, in_vertex,
					  out_vertex) ;
  *ierr = (err < 0) ?  ZOLTAN_FATAL : ZOLTAN_OK ;
  for(int i = 0; i < num_obj; i++) 
    on_proc[i] = (assigned[i] == true) ? 1 : 0;
  *in_order = (order == true) ? 1 : 0 ;
  delete [] assigned ;
}

int ZC_JR_First_Coarse_Obj_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
			      ZOLTAN_ID_PTR lid, int *assigned, int *num_vertices,
			      ZOLTAN_ID_TYPE *vertices, int *in_order, 
			      ZOLTAN_ID_TYPE *in_vertex, ZOLTAN_ID_TYPE *out_vertex, 
			      int *ierr)
{
  bool assign, order ;
  int err = ::ZoltanSpace::Globals::pTreeInfo->GetFirstCoarseMeshEntity(
					  num_gid, num_lid, gid, lid, &assign, 
					  num_vertices, vertices, &order, in_vertex,
					  out_vertex) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  *assigned = (assign == true) ? 1 : 0 ;
  *in_order = (order == true) ? 1 : 0 ;

  return( (*ierr == 0) ? 1 : 0 ) ;
}

int ZC_JR_Next_Coarse_Obj_Fn(void *data, int num_gid, int num_lid, 
			     ZOLTAN_ID_PTR prev_gid, ZOLTAN_ID_PTR prev_lid, 
			     ZOLTAN_ID_PTR next_gid, ZOLTAN_ID_PTR next_lid,
			     int *assigned, int *num_vertices, 
			     ZOLTAN_ID_TYPE *vertex_list,
			     ZOLTAN_ID_TYPE *in_vertex, 
			     ZOLTAN_ID_TYPE *out_vertex, int *ierr)
{
  bool on_proc ;
  int err = ::ZoltanSpace::Globals::pTreeInfo->GetNextCoarseMeshEntity(
                                          num_gid, num_lid, prev_gid, prev_lid,
					  next_gid, next_lid, &on_proc, num_vertices,
					  vertex_list, in_vertex, out_vertex) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  *assigned = (on_proc == true) ?  1 : 0 ;
  return( (*ierr == 0) ? 1 : 0) ;
}

int ZC_JR_Num_Child_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
		       ZOLTAN_ID_PTR lid, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pTreeInfo->GetNumberOfChildren(
					    num_gid, num_lid, gid, lid) ;
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return( res ) ;
}

void ZC_JR_Child_List_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
			 ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR child_gid, 
			 ZOLTAN_ID_PTR child_lid, int *assigned, int *num_vertices, 
			 ZOLTAN_ID_TYPE *vertices, ZOLTAN_REF_TYPE *ref, 
			 ZOLTAN_ID_TYPE *in_vertex, ZOLTAN_ID_TYPE *out_vertex, 
			 int *ierr)
{
  int refinement, *vertex_list ;
  bool on_proc ; 
  int err = ::ZoltanSpace::Globals::pTreeInfo->GetChildrenListOfEntity(
                                          num_gid, num_lid, gid, lid, child_gid, 
					  child_lid, &on_proc, num_vertices,
					  &vertex_list, &refinement, in_vertex, 
					  out_vertex) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  *assigned = (on_proc == true ) ? 1 : 0;
  switch (refinement)
  {
    case (1) :
      *ref = ZOLTAN_TRI_BISECT ;   break ;
    case (2) :
      *ref = ZOLTAN_QUAD_QUAD ;    break ;
    case (3) :
      *ref = ZOLTAN_OTHER_REF ;    break ;
    case (4) :
      *ref = ZOLTAN_IN_ORDER ;     break ;
    default :
      *ref = ZOLTAN_IN_ORDER ;     *ierr = -1 ;
  }

  for(int i = 0; i < *num_vertices; i++) vertices[i] = vertex_list[i] ;
  ::ZoltanSpace::Globals::pTreeInfo->cleanup() ;
}

void ZC_JR_Child_Weight_Fn(void *data, int num_gid, int num_lid, ZOLTAN_ID_PTR gid,
			   ZOLTAN_ID_PTR lid, int wgt_dim, float *wgt, int *ierr)
{
  bool wgt_reqd  = (wgt_dim == 0) ? false : true ;
  int err = ::ZoltanSpace::Globals::pTreeInfo->GetChildrenWeight(num_gid,
					  num_lid, gid, lid, wgt_reqd, wgt) ;
  *ierr = (err < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
}

int ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesTreeInfoPort(
                                            struct Zoltan_Struct *zz)
{

  if ( pTreeInfo == 0 ) return (-1) ;

  ::ZoltanSpace::Globals::pTreeInfo = pTreeInfo ;

  Zoltan_Set_Num_Coarse_Obj_Fn( zz, &(ZC_JR_Num_Coarse_Elem_Fn), NULL) ;

  if ( pIterate->value == false )
    Zoltan_Set_Coarse_Obj_List_Fn( zz, &(ZC_JR_Coarse_Obj_List_Fn), NULL) ;
  else
  {
    Zoltan_Set_First_Coarse_Obj_Fn( zz, &(ZC_JR_First_Coarse_Obj_Fn), NULL) ;
    Zoltan_Set_Next_Coarse_Obj_Fn( zz, &(ZC_JR_Next_Coarse_Obj_Fn), NULL) ;
  }

  if ( pIterate->value == false )
    Zoltan_Set_Num_Child_Fn( zz, &(ZC_JR_Num_Child_Fn), NULL) ;
  else
  {
    Zoltan_Set_Child_List_Fn( zz, &(ZC_JR_Child_List_Fn), NULL) ;
    Zoltan_Set_Child_Weight_Fn( zz, &(ZC_JR_Child_Weight_Fn), NULL) ;
  }

  return(0) ;
}
//-------------------------------------------------------------------------

// methods that do packing and unpacking; things in DataMigrationPort

int ZC_JR_Obj_Size_Fn(void *data, int num_gid_entries, int num_lid_entries,
		      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pDataMigPort->DataBufferForThisEntity(
			  num_gid_entries, num_lid_entries, global_id, local_id) ;
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK ;
  return(res) ;
}

void ZC_JR_Pack_Obj_Fn(void *data, int num_gid_entries, int num_lid_entries,
		       ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
		       int dest_proc, int size, char *buf, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pDataMigPort->PackThisEntity(num_gid_entries,
			  num_lid_entries, global_id, local_id, dest_proc, size, buf) ;
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  
  return  ;
}

void ZC_JR_Unpack_Obj_Fn(void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, 
			 int size, char *buf, int *ierr)
{
  int res = ::ZoltanSpace::Globals::pDataMigPort->UnpackThisEntity(num_gid_entries,
			  global_id, size, buf) ;
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  
  return  ;
}

void ZC_JR_Pre_Migrate_Fn(void *data, int num_gid_entries, int num_lid_entries,
			  int num_import, ZOLTAN_ID_PTR import_glob_id, 
			  ZOLTAN_ID_PTR import_local_id, int *import_procs, 
			  int num_export, ZOLTAN_ID_PTR export_glob_id, 
			  ZOLTAN_ID_PTR export_local_id, int *export_procs,
			  int *ierr)
{
  ::ZoltanSpace::EntityListImpl outgoing(num_export), incoming(num_import) ;
  
  incoming.create_gid_list(num_import*num_gid_entries, import_glob_id) ;
  incoming.create_lid_list(num_import*num_lid_entries, import_local_id) ;
  incoming.create_proc_list(num_import, import_procs) ;

  outgoing.create_gid_list(num_export*num_gid_entries, export_glob_id) ;
  outgoing.create_lid_list(num_export*num_lid_entries, export_local_id) ;
  outgoing.create_proc_list(num_export, export_procs) ;
  
  incoming.addRef() ; outgoing.addRef() ;
  int res = ::ZoltanSpace::Globals::pDataMigPort->PrePackingProcessing( 
			   dynamic_cast< ::LoadPartitionerSpace::EntityList *> (&outgoing),
			   dynamic_cast< ::LoadPartitionerSpace::EntityList *>  (&incoming) );
  *ierr = (res < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;

  return ;
}

void ZC_JR_Mid_Migrate_Fn(void *data, int num_gid_entries, int num_lid_entries,
			  int num_import, ZOLTAN_ID_PTR import_glob_id, 
			  ZOLTAN_ID_PTR import_local_id, int *import_procs, 
			  int num_export, ZOLTAN_ID_PTR export_glob_id, 
			  ZOLTAN_ID_PTR export_local_id, int *export_procs,
			  int *ierr)
{
  ::ZoltanSpace::EntityListImpl outgoing(num_export), incoming(num_import) ;
  
  incoming.create_gid_list(num_import*num_gid_entries, import_glob_id) ;
  incoming.create_lid_list(num_import*num_lid_entries, import_local_id) ;
  incoming.create_proc_list(num_import, import_procs) ;

  outgoing.create_gid_list(num_export*num_gid_entries, export_glob_id) ;
  outgoing.create_lid_list(num_export*num_lid_entries, export_local_id) ;
  outgoing.create_proc_list(num_export, export_procs) ;

  incoming.addRef() ; outgoing.addRef() ;
  int res = ::ZoltanSpace::Globals::pDataMigPort->PostPackingProcessing( 
			dynamic_cast< ::LoadPartitionerSpace::EntityList *> (&outgoing),
		        dynamic_cast< ::LoadPartitionerSpace::EntityList *> (&incoming) ) ;
  *ierr = (res < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;

  return ;
}

void ZC_JR_Post_Migrate_Fn(void *data, int num_gid_entries, int num_lid_entries,
			   int num_import, ZOLTAN_ID_PTR import_glob_id, 
			   ZOLTAN_ID_PTR import_local_id, int *import_procs, 
			   int num_export, ZOLTAN_ID_PTR export_glob_id, 
			   ZOLTAN_ID_PTR export_local_id, int *export_procs,
			  int *ierr)
{
  ::ZoltanSpace::EntityListImpl outgoing(num_export), incoming(num_import) ;
  
  incoming.create_gid_list(num_import*num_gid_entries, import_glob_id) ;
  incoming.create_lid_list(num_import*num_lid_entries, import_local_id) ;
  incoming.create_proc_list(num_import, import_procs) ;

  outgoing.create_gid_list(num_export*num_gid_entries, export_glob_id) ;
  outgoing.create_lid_list(num_export*num_lid_entries, export_local_id) ;
  outgoing.create_proc_list(num_export, export_procs) ;

  incoming.addRef() ; outgoing.addRef() ;
  int res = ::ZoltanSpace::Globals::pDataMigPort->PostUnpackingProcessing( 
		      dynamic_cast< ::LoadPartitionerSpace::EntityList *> (&outgoing),	
		      dynamic_cast< ::LoadPartitionerSpace::EntityList *> (&incoming) ) ;

  *ierr = (res < 0 ) ? ZOLTAN_FATAL : ZOLTAN_OK ;

  return ;
}

int ZoltanSpace::PartitionerFactory_JR::SetZoltanQueriesDataMigrationPort(
                                           struct Zoltan_Struct *zz)
{
  if ( pDataMig == 0 ) return (-1) ;

  ::ZoltanSpace::Globals::pDataMigPort = pDataMig ;
  
  Zoltan_Set_Obj_Size_Fn(zz, &(ZC_JR_Obj_Size_Fn), NULL) ;
  Zoltan_Set_Pack_Obj_Fn(zz, &(ZC_JR_Pack_Obj_Fn), NULL) ;
  Zoltan_Set_Unpack_Obj_Fn(zz, &(ZC_JR_Unpack_Obj_Fn), NULL) ;
  Zoltan_Set_Pre_Migrate_Fn(zz, &(ZC_JR_Pre_Migrate_Fn), NULL) ;
  Zoltan_Set_Mid_Migrate_Fn(zz, &(ZC_JR_Mid_Migrate_Fn), NULL) ;
  Zoltan_Set_Post_Migrate_Fn(zz, &(ZC_JR_Post_Migrate_Fn), NULL) ;

  return(0) ;
}

// --------------------------------------------------------------------------------
