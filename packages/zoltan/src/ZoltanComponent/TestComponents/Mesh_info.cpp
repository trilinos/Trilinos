/*
  methods that give info about the mesh itself
*/

#include "Mesh.h"

extern "C"
{
#include "include/zoltan.h"

  int get_num_elements(void *data, int *ierr) ;
  
  int get_first_element(void *data, int num_gid_entries, int num_lid_entries,
			ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
			int wdim, float *wgt, int *ierr) ;
  
  int get_next_element(void *data, int num_gid_entries, int num_lid_entries,
		       ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
		       ZOLTAN_ID_PTR next_global_id, ZOLTAN_ID_PTR next_local_id, 
		       int wdim, float *next_wgt, int *ierr) ;

  int get_num_geom(void *data, int *ierr) ;

  void get_geom(void *data, int num_gid_entries, int num_lid_entries,
		ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
		double *coor, int *ierr) ;

  int get_num_edges(void *data, int num_gid_entries, int num_lid_entries,
		    ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) ;

  void get_edge_list (void *data, int num_gid_entries, int num_lid_entries, 
		      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
		      ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
		      int get_ewgts, float *nbor_ewgts, int *ierr) ;
}

int ZoltanTestSpace::Mesh::SetTypeOfEntity(char *name) 
{
  if (is_init == false) { init() ; is_init = true ; }
  /* There will only be one type of entity in this mesh */
  return(0) ;
}

int  ZoltanTestSpace::Mesh::GetNumberOfEntities()
{
  if (is_init == false) { init() ; is_init = true ; }
  int ierr = 0 ;
  return( get_num_elements( &mesh, &ierr) ) ;
}


int ZoltanTestSpace::Mesh::EntitiesOnThisProc(int num_gid_entries, 
					      int num_lid_entries,
					      unsigned int *gids, 
					      unsigned int *lids,
					      bool wgt_reqd, 
					      float *weights) 
{
  if (is_init == false) { init() ; is_init = true ; }
  return(-1) ; // not implemented; iterator access only
}

int ZoltanTestSpace::Mesh::GetFirstEntity(int num_gid_entries, 
					  int num_lid_entries,
					  unsigned int *first_gid, 
					  unsigned int *first_lid, 
					  bool wgt_reqd, 
					  float *first_weight) 
{
  if (is_init == false) { init() ; is_init = true ; }
  int wdim = ( wgt_reqd == true ? 1 : 0 ) ;
  int ierr = 0 ;
  get_first_element(&mesh, num_gid_entries, num_lid_entries, first_gid,
		    first_lid, wdim, first_weight, &ierr) ;
  if (ierr == ZOLTAN_OK)   return(0) ;   else     return(-1) ;

}

int ZoltanTestSpace::Mesh::GetNextEntity(int num_gid_entries, 
					 int num_lid_entries,
					 unsigned int *prev_gid, 
					 unsigned int *prev_lid,
					 unsigned int *next_gid, 
					 unsigned int *next_lid,
					 bool wgt_reqd, 
					 float *next_weight)
{
  if (is_init == false) { init() ; is_init = true ; }
  int wdim = (wgt_reqd == true) ? 1 : 0 ;
  int ierr = 0 ;

  get_next_element(&mesh, num_gid_entries, num_lid_entries, prev_gid, prev_lid,
		   next_gid, next_lid, wdim, next_weight, &ierr);

  if (ierr == ZOLTAN_OK)   return(0) ;   else     return(-1) ; 
}

int ZoltanTestSpace::Mesh::SetTypeOfGeom(char *name)
{
  if (is_init == false) { init() ; is_init = true ; }
  // Only one type of geom
  return(0) ;
}

int ZoltanTestSpace::Mesh::GetNumberOfDims()
{
  if (is_init == false) { init() ; is_init = true ; }
  int ierr = 0 ;
  return( get_num_geom(&mesh, &ierr) ) ;
}

int ZoltanTestSpace::Mesh::GetGeomDescription(int num_gid_entries, 
					      int num_lid_entries,
					      unsigned int *gids, 
					      unsigned int *lids,
					      double *geom_vec)
{
  if (is_init == false) { init() ; is_init = true ; }
  int ierr = 0 ;
  get_geom(&mesh, num_gid_entries, num_lid_entries, gids, lids, geom_vec,
	   &ierr) ;
  if ( ierr == ZOLTAN_OK) return (0) ;  else return(-1) ;
}

int ZoltanTestSpace::Mesh::SetTypeOfEdge(char *name)
{
  if (is_init == false) { init() ; is_init = true ; }
  return(0) ; // not implemented, just supplies the default.
}

int ZoltanTestSpace::Mesh::GetNumberOfEdges(int num_gid_entries, 
					    int num_lid_entries,
					    unsigned int *gids, 
					    unsigned int *lids)
{
  if (is_init == false) { init() ; is_init = true ; }
  int ierr = 0 ;
  return( get_num_edges(&mesh, num_gid_entries, num_lid_entries, gids, lids, &ierr) ) ;
}

int ZoltanTestSpace::Mesh::EdgesListForThisEntity(int num_gid_entries, 
						  int num_lid_entries,
						  unsigned int *gid, 
						  unsigned int *lid,
						  unsigned int *neighbor_gids, 
						  int *neighbor_procs, 
						  bool wgt_reqd,
						  float *edge_weights)
{
  if (is_init == false) { init() ; is_init = true ; }
  int ierr = 0 ;
  int wdim = (wgt_reqd == true) ? 1 : 0 ;

  get_edge_list( &mesh, num_gid_entries, num_lid_entries, gid, lid, 
		 neighbor_gids, neighbor_procs, wdim, edge_weights, &ierr) ;

  if (ierr == ZOLTAN_OK)   return(0) ;   else     return(-1) ;
}


