#include "PartitionerFactory.h"
#include "BaseLB.h"

int ZoltanSpace::PartitionerFactory_JR::ComputeDestinations(
               ::LoadPartitionerSpace::LoadBalancer *A, 
	       ::LoadPartitionerSpace::EntityList *OffProcEntitiesComingToMe,
               ::LoadPartitionerSpace::EntityList **MyEntitiesMigratingElsewhere)
{
  // Get the relevant info out of the incoming list of entities
  int num_known = OffProcEntitiesComingToMe->GetListLength() ;
  unsigned int *known_global_id = OffProcEntitiesComingToMe->GetAllGIDs() ;
  unsigned int *known_local_id = OffProcEntitiesComingToMe->GetAllLIDs() ;
  int *known_procs = OffProcEntitiesComingToMe->GetResidentProcsList() ;
  
  // Get ready for the outgoing list
  int num_export ;
  unsigned int *export_glob_id, *export_loc_id ;
  int *export_procs ;

  Zoltan_Struct *zz = ( dynamic_cast< ::ZoltanSpace::BaseLB * > (A) )->get_zz() ;

  Zoltan_Compute_Destinations( zz, num_known, known_global_id, known_local_id, 
			      known_procs, &num_export, &export_glob_id, &export_loc_id,
			      &export_procs) ;

  if (temp_outgoing != 0)
  {
    unsigned int *l1 = temp_outgoing->GetAllGIDs() ;
    unsigned int *l2 = temp_outgoing->GetAllLIDs() ;
    int *l3 = temp_outgoing->GetResidentProcsList() ;
    
    Zoltan_LB_Free_Data(NULL, NULL, NULL, &l1, &l2, &l3) ;

    delete (temp_outgoing) ; temp_outgoing = 0 ;
  }

  int num_gid_entries = OffProcEntitiesComingToMe->GetGIDListLength() /
                        OffProcEntitiesComingToMe->GetListLength() ;
  int num_lid_entries = OffProcEntitiesComingToMe->GetLIDListLength() /
                        OffProcEntitiesComingToMe->GetListLength() ;

  temp_outgoing = new EntityListImpl( num_export ) ;
  temp_outgoing->create_gid_list(num_gid_entries*num_export, export_glob_id) ;
  temp_outgoing->create_lid_list(num_lid_entries*num_export, export_loc_id) ;
  temp_outgoing->create_proc_list(num_export, export_procs) ;

  *MyEntitiesMigratingElsewhere = dynamic_cast< ::LoadPartitionerSpace::EntityList * >
                                               ( temp_outgoing ) ;

  OffProcEntitiesComingToMe->deleteRef() ; (*MyEntitiesMigratingElsewhere)->addRef() ;
  export_glob_id = 0 ; export_loc_id = 0 ; export_procs = 0 ;

  return(0) ;
}

int ZoltanSpace::PartitionerFactory_JR::ComputeSources(
			 ::LoadPartitionerSpace::LoadBalancer *A, 
			 ::LoadPartitionerSpace::EntityList *MyEntitiesMigratingElsewhere,
                         ::LoadPartitionerSpace::EntityList **OffProcEntitiesComingToMe)
{
  // Get the relevant info out of the incoming list of entities
  int num_known = MyEntitiesMigratingElsewhere->GetListLength() ;
  unsigned int *known_global_id = MyEntitiesMigratingElsewhere->GetAllGIDs() ;
  unsigned int *known_local_id = MyEntitiesMigratingElsewhere->GetAllLIDs() ;
  int *known_procs = MyEntitiesMigratingElsewhere->GetResidentProcsList() ;
  
  // Get ready for the outgoing list
  int num_import ;
  unsigned int *import_glob_id, *import_loc_id ;
  int *import_procs ;

  Zoltan_Struct *zz = ( dynamic_cast< ::ZoltanSpace::BaseLB * > (A) )->get_zz() ;

  Zoltan_Compute_Destinations( zz, num_known, known_global_id, known_local_id, 
			      known_procs, &num_import, &import_glob_id, &import_loc_id,
			      &import_procs) ;

  if (temp_incoming != 0)
  {
    unsigned int *l1 = temp_incoming->GetAllGIDs() ;
    unsigned int *l2 = temp_incoming->GetAllLIDs() ;
    int *l3 = temp_incoming->GetResidentProcsList() ;
    
    Zoltan_LB_Free_Data(NULL, NULL, NULL, &l1, &l2, &l3) ;

    delete (temp_incoming) ; temp_incoming = 0 ;
  }

  int num_gid_entries = MyEntitiesMigratingElsewhere->GetGIDListLength() /
                        MyEntitiesMigratingElsewhere->GetListLength() ;
  int num_lid_entries = MyEntitiesMigratingElsewhere->GetLIDListLength() /
                        MyEntitiesMigratingElsewhere->GetListLength() ;

  temp_incoming = new EntityListImpl( num_import ) ;
  temp_incoming->create_gid_list(num_gid_entries*num_import, import_glob_id) ;
  temp_incoming->create_lid_list(num_lid_entries*num_import, import_loc_id) ;
  temp_incoming->create_proc_list(num_import, import_procs) ;

  import_glob_id = 0 ; import_loc_id = 0 ; import_procs = 0 ;
  MyEntitiesMigratingElsewhere->deleteRef() ; (*OffProcEntitiesComingToMe)->addRef() ;

  return(0) ;
}

int ZoltanSpace::PartitionerFactory_JR::MigrateEntities(
	             ::LoadPartitionerSpace::LoadBalancer *A, 
                     ::LoadPartitionerSpace::EntityList *MyEntitiesMigratingElsewhere, 
                     ::LoadPartitionerSpace::EntityList *OffProcEntitiesComingToMe)
{
  ::ZoltanSpace::BaseLB *jj = dynamic_cast< ZoltanSpace::BaseLB * > (A) ;
  Zoltan_Struct *zz = jj->get_zz() ;
  int res = SetZoltanQueriesDataMigrationPort( jj->get_zz() ) ;

  if ( res != 0 ) 
  {
    cerr << " ZoltanSpace::PartitionerFactory_JR::MigrateEntities() "
	 << " could not set DataMigrationPort - doing no migration"
	 << endl ;
    return ( res ) ;
  }
  
  int num_import = OffProcEntitiesComingToMe->GetListLength() ;
  unsigned int *import_glob_ids =  OffProcEntitiesComingToMe->GetAllGIDs() ;
  unsigned int *import_local_ids = OffProcEntitiesComingToMe->GetAllLIDs() ;
  int *import_procs = OffProcEntitiesComingToMe->GetResidentProcsList() ;
  
  int num_export = MyEntitiesMigratingElsewhere->GetListLength() ;
  unsigned int *export_glob_ids = MyEntitiesMigratingElsewhere->GetAllGIDs() ;
  unsigned int *export_local_ids = MyEntitiesMigratingElsewhere->GetAllLIDs() ;
  int *export_procs = MyEntitiesMigratingElsewhere->GetResidentProcsList() ;

  res = Zoltan_Help_Migrate(zz, num_import, import_glob_ids, import_local_ids,
			    import_procs, num_export, export_glob_ids, export_local_ids,
			    export_procs) ;

  res = (res == ZOLTAN_OK) ? 0 : res ;

  MyEntitiesMigratingElsewhere->deleteRef() ;
  OffProcEntitiesComingToMe->deleteRef() ;

  return(res) ;
}

void ZoltanSpace::PartitionerFactory_JR::Cleanup( ::LoadPartitionerSpace::LoadBalancer *A )
{
  if (temp_incoming != 0)
  {
    unsigned int *l1 = temp_incoming->GetAllGIDs() ;
    unsigned int *l2 = temp_incoming->GetAllLIDs() ;
    int *l3 = temp_incoming->GetResidentProcsList() ;
    
    Zoltan_LB_Free_Data(NULL, NULL, NULL, &l1, &l2, &l3) ;

    delete (temp_incoming) ; temp_incoming = 0 ;
  }

  if (temp_outgoing != 0)
  {
    unsigned int *l1 = temp_outgoing->GetAllGIDs() ;
    unsigned int *l2 = temp_outgoing->GetAllLIDs() ;
    int *l3 = temp_outgoing->GetResidentProcsList() ;
    
    Zoltan_LB_Free_Data(NULL, NULL, NULL, &l1, &l2, &l3) ;

    delete (temp_outgoing) ; temp_outgoing = 0 ;
  }

  return ;
}
