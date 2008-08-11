/*
  Methods that help with migrating entities

*/

#include "Mesh.h"

extern "C"
{
#include "include/zoltan.h"

  /*
    CAUTION !!

    These functions are defined in ZOLTAN_HOME/driver/ and can and do
    change their function signatures. Check that they haven't changed.
    They are implemented in ZOLTAN_HOME/driver/dr_migrate.c

    migrate_pre_process and migrate_post_process now take in 
    import_to_part[num_import] and export_to_part[num_export]. As of
    09/03/03, each processor can have > 1 partitions, unlike the old
    days when a processor could have just 1 partition. Now, an imported
    or exported element not only needs a processor to be specified, but
    also a partition number. This partition number is assigned by 
    Zoltan_LB_Partition().
    
    Jaideep Ray
  */
  void migrate_pre_process(void *data, int num_gid_entries, int num_lid_entries, 
			   int num_import, ZOLTAN_ID_PTR import_global_ids,
			   ZOLTAN_ID_PTR import_local_ids, int *import_procs,
			   int *import_to_part,
			   int num_export, ZOLTAN_ID_PTR export_global_ids,
			   ZOLTAN_ID_PTR export_local_ids, int *export_procs,
			   int *export_to_part,
			   int *ierr)  ;


  void migrate_post_process(void *data, int num_gid_entries, int num_lid_entries,
			    int num_import, ZOLTAN_ID_PTR import_global_ids,
			    ZOLTAN_ID_PTR import_local_ids, int *import_procs,
			    int *import_to_part,
			    int num_export, ZOLTAN_ID_PTR export_global_ids,
			    ZOLTAN_ID_PTR export_local_ids, int *export_procs,
			    int *export_to_part,
			    int *ierr) ;

  int migrate_elem_size(void *data, int num_gid_entries, int num_lid_entries,
			ZOLTAN_ID_PTR elem_gid, ZOLTAN_ID_PTR elem_lid, int *ierr) ;

  void migrate_pack_elem(void *data, int num_gid_entries, int num_lid_entries,
			 ZOLTAN_ID_PTR elem_gid, ZOLTAN_ID_PTR elem_lid,
			 int mig_proc, int elem_data_size, char *buf, int *ierr) ;

  void migrate_unpack_elem(void *data, int num_gid_entries, ZOLTAN_ID_PTR elem_gid, 
			   int elem_data_size, char *buf, int *ierr) ;
} 

int ZoltanTestSpace::Mesh::PrePackingProcessing( 
			 ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
			 ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) 
{
  if (is_init == false) { init() ; is_init = true ; }

  int num_gid_entries, num_lid_entries ;
  if (  MyEntitiesGoingOut->GetListLength() > 0 ) // Not an empty list
  {
    num_gid_entries = MyEntitiesGoingOut->GetGIDListLength() / 
                      MyEntitiesGoingOut->GetListLength() ; 
    num_lid_entries = MyEntitiesGoingOut->GetLIDListLength() / 
                      MyEntitiesGoingOut->GetListLength() ; 
  }
  else
  {
    // just put something here. it's a null list anyway
    num_gid_entries = num_lid_entries = 1;
  }
  int num_import = OthersEntitiesComingIn->GetListLength() ;
  unsigned int *import_global_ids = OthersEntitiesComingIn->GetAllGIDs() ;
  unsigned int *import_local_ids =  OthersEntitiesComingIn->GetAllLIDs() ;
  int *import_procs = OthersEntitiesComingIn->GetResidentProcsList() ;

  int num_export = MyEntitiesGoingOut->GetListLength() ;
  unsigned int *export_global_ids = MyEntitiesGoingOut->GetAllGIDs() ;
  unsigned int *export_local_ids  = MyEntitiesGoingOut->GetAllLIDs() ;
  int *export_procs = MyEntitiesGoingOut->GetResidentProcsList() ;

  int ierr = 0 ;
  /*
    JR CAUTION :
    I am going to assume that each processor has only 1 partition and its
    number is the same as the processor number. THIS MIGHT NOT AGREE WITH
    THE NUMBER RETURNED BY ZOLTAN_LB_PARTITION. However, ZoltanComponent will
    not call Zoltan_LB_Partition, just Zoltan_LB_Balance
  */
  migrate_pre_process( static_cast< void * > (&mesh), num_gid_entries, num_lid_entries,
		       num_import, import_global_ids, import_local_ids, 
		       import_procs, 
		       import_procs,
		       num_export, export_global_ids, 
		       export_local_ids, export_procs, 
		       export_procs,
		       &ierr) ;

  MyEntitiesGoingOut ->deleteRef() ; OthersEntitiesComingIn->deleteRef() ;

  if (ierr == ZOLTAN_OK) return (0) ; else return( ierr ) ;
}

int ZoltanTestSpace::Mesh::PostPackingProcessing( 
			 ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
			 ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) 
{
  if (is_init == false) { init() ; is_init = true ; }

  MyEntitiesGoingOut ->deleteRef() ; OthersEntitiesComingIn->deleteRef() ;
  // Nothing to do here 
  return(0) ;
}

int ZoltanTestSpace::Mesh::PostUnpackingProcessing( 
			 ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
			 ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) 
{
  if (is_init == false) { init() ; is_init = true ; }

  int num_gid_entries, num_lid_entries ;
  if ( MyEntitiesGoingOut->GetListLength() > 0 ) // Not an empty list
  {
    num_gid_entries = MyEntitiesGoingOut->GetGIDListLength() /
                      MyEntitiesGoingOut->GetListLength() ;
    num_lid_entries = MyEntitiesGoingOut->GetLIDListLength() / 
                      MyEntitiesGoingOut->GetListLength() ;
  }
  else
  {
    // just put something here. it's a null list anyway
    num_gid_entries = num_lid_entries = 1 ;
  }

  int num_import = OthersEntitiesComingIn->GetListLength() ;
  unsigned int *import_global_ids = OthersEntitiesComingIn->GetAllGIDs() ;
  unsigned int *import_local_ids =  OthersEntitiesComingIn->GetAllLIDs() ;
  int *import_procs = OthersEntitiesComingIn->GetResidentProcsList() ;

  int num_export = MyEntitiesGoingOut->GetListLength() ;
  unsigned int *export_global_ids = MyEntitiesGoingOut->GetAllGIDs() ;
  unsigned int *export_local_ids  = MyEntitiesGoingOut->GetAllLIDs() ;
  int *export_procs = MyEntitiesGoingOut->GetResidentProcsList() ;

  int ierr = 0 ;

  /*
    JR CAUTION :
    I am going to assume that each processor has only 1 partition and its
    number is the same as the processor number. THIS MIGHT NOT AGREE WITH
    THE NUMBER RETURNED BY ZOLTAN_LB_PARTITION. However, ZoltanComponent will
    not call Zoltan_LB_Partition, just Zoltan_LB_Balance
  */
  migrate_post_process( static_cast< void * > (&mesh), num_gid_entries, num_lid_entries,
			num_import, import_global_ids, import_local_ids, 
			import_procs, 
			import_procs,
			num_export, export_global_ids, export_local_ids, 
			export_procs, 
			export_procs,
			&ierr) ;

  MyEntitiesGoingOut ->deleteRef() ; OthersEntitiesComingIn->deleteRef() ;

  if (ierr == ZOLTAN_OK) return (0) ; else return( ierr ) ;
}

int ZoltanTestSpace::Mesh::DataBufferForThisEntity(int num_gid_entries, 
						   int num_lid_entries,
						   unsigned int *gid,
						   unsigned int *lid)
{
  if (is_init == false) { init() ; is_init = true ; }

  int ierr = 0 ;
  int size = migrate_elem_size( (void *) (&mesh), num_gid_entries, 
				num_lid_entries, gid, lid, &ierr) ;

  if (ierr != ZOLTAN_OK) 
    return ( -1 ) ; 
  else 
    return( static_cast<unsigned int>(size) ) ;
}

int ZoltanTestSpace::Mesh::PackThisEntity(int num_gid_entries, int num_lid_entries,
			  unsigned int *gid, unsigned int *lid, 
			  int dest_proc, int buffer_size, char *buffer)
{
  if (is_init == false) { init() ; is_init = true ; }

  int ierr = 0 ;
  migrate_pack_elem( (void *) (&mesh), num_gid_entries, num_lid_entries, gid, lid,
		     dest_proc, buffer_size, buffer, &ierr) ;

  if ( ierr == ZOLTAN_OK) return (0) ; else return( ierr ) ;
}

int ZoltanTestSpace::Mesh::UnpackThisEntity(int num_gid_entries, 
                             unsigned int *gid, int buffer_size, char *buffer)
{
  if (is_init == false) { init() ; is_init = true ; }

  int ierr = 0 ;
  migrate_unpack_elem( (void *) (&mesh) , num_gid_entries, gid, buffer_size,
		       buffer, &ierr) ;
  if ( ierr == ZOLTAN_OK ) return (0) ; else return ( ierr ) ;
}





