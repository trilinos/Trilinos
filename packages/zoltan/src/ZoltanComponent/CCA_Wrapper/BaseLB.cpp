/*
  Implementation of the functionality common to all load balancers
  Jaideep Ray, SNL, 08/26/02
*/

#include "BaseLB.h"

int ZoltanSpace::BaseLB::CreateBalancedPartition(bool *did_it_change, 
			   ::LoadPartitionerSpace::EntityList **ElementsToImport,
			   ::LoadPartitionerSpace::EntityList **ElementsToExport)
{

  if (is_init == false ) init() ;

  int changes, num_gid_entries, num_lid_entries ;
  int num_import, *import_from_procs, num_export, *export_to_procs ;
  ZOLTAN_ID_PTR import_gid, import_lid, export_gid, export_lid ;
  
  int res = -1;
  res = Zoltan_LB_Balance(my_zz, &changes, &num_gid_entries, &num_lid_entries,
  			  &num_import, &import_gid, &import_lid, &import_from_procs,
  			  &num_export, &export_gid, &export_lid, &export_to_procs) ;

  cout << "ZoltanSpace::BaseLB::CreateBalancedPartition() done with LB_Balance res = "
       << res << endl ;
  if ( res == ZOLTAN_FATAL || res == ZOLTAN_MEMERR ) return (-1) ;

  if (incoming == 0) 
  {
    incoming = new EntityListImpl(num_import); assert(incoming != 0) ;
    if ( incoming == 0 ) return (-1) ; else incoming->addRef() ;
  }
  incoming->create_gid_list(num_import*num_gid_entries, import_gid) ;
  incoming->create_lid_list(num_import*num_lid_entries, import_lid);
  incoming->create_proc_list(num_import, import_from_procs);
  incoming->addRef() ;

  cout << " BaseLB::CreateBalancedPartition() finished incoming " << endl ;
  if (outgoing == 0) 
  {
    outgoing = new EntityListImpl(num_export); assert(incoming != 0) ;
    if ( outgoing == 0 ) return(-1) ; else outgoing->addRef() ;
  }
  outgoing->create_gid_list(num_export*num_gid_entries, export_gid) ;
  outgoing->create_lid_list(num_export*num_lid_entries, export_lid);
  outgoing->create_proc_list(num_export, export_to_procs);
  outgoing->addRef() ;
  cout << " BaseLB::CreateBalancedPartition() finished outgoing " << endl ;
  
  *did_it_change = (changes == 1) ? true : false ;
  *ElementsToImport = dynamic_cast< ::LoadPartitionerSpace::EntityList *> (incoming) ;
  *ElementsToExport = dynamic_cast< ::LoadPartitionerSpace::EntityList *> (outgoing) ;

  return(res) ;
}

int ZoltanSpace::BaseLB::EvaluateDecomposition(bool print_stats, int *nobjs_on_proc,
					       float *obj_wgt_on_proc, 
					       int *ncuts_on_proc, 
					       float *cut_wgts_on_proc, int *nbndry,
					       int *n_my_adjacent_procs)
{
  if (is_init == false ) init() ;

  int print = (print_stats == true) ? 1 : 0 ;
  int res = Zoltan_LB_Eval(my_zz, print, nobjs_on_proc, obj_wgt_on_proc, ncuts_on_proc,
			   cut_wgts_on_proc, nbndry, n_my_adjacent_procs) ;

  if ( res == ZOLTAN_FATAL || res == ZOLTAN_MEMERR ) return (-1) ;

  return(res) ;
}

int ZoltanSpace::BaseLB::cleanup()
{

  if ( incoming == 0 && outgoing == 0 ) return(0) ;

  ZOLTAN_ID_PTR igid = incoming->GetAllGIDs() ; 
  ZOLTAN_ID_PTR ilid = incoming->GetAllLIDs() ; 
  int *iproc = incoming->GetResidentProcsList() ;
  ZOLTAN_ID_PTR ogid = outgoing->GetAllGIDs() ; 
  ZOLTAN_ID_PTR olid = outgoing->GetAllLIDs() ; 
  int *oproc = outgoing->GetResidentProcsList(); 

  int res = ZOLTAN_OK;
  res = Zoltan_LB_Free_Data(&igid, &ilid, &iproc, &ogid, &olid, &oproc) ;
  if ( res == ZOLTAN_FATAL || res == ZOLTAN_MEMERR ) return (-1) ;

  incoming->zero_out_lists() ; incoming->deleteRef() ;  delete(incoming) ; incoming = 0 ;
  outgoing->zero_out_lists() ; outgoing->deleteRef() ;  delete(outgoing) ; outgoing = 0 ;

  return(0) ;
}


