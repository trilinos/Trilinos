/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <ctype.h>
#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 * This function gets a list of objects one way or the other,
 * i.e., by calling either Get_Obj_List or Get_First_Obj+Get_Next_Obj.
 */

void Zoltan_Get_Obj_List(
  ZZ *zz, 
  ZOLTAN_ID_PTR global_ids, 
  ZOLTAN_ID_PTR local_ids, 
  int wdim, 
  float *objwgts, 
  int *ierr
)
{
  int i, n;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_off, lid_off;
  ZOLTAN_ID_PTR lid, next_lid;  /* Temporary pointers to local IDs; used to pass 
                               NULL to query functions when 
                               NUM_LID_ENTRIES == 0. */

  *ierr = ZOLTAN_OK;
  if (zz->Get_Obj_List != NULL){
    /* Get object list directly */
    zz->Get_Obj_List(zz->Get_Obj_List_Data, 
                     num_gid_entries, num_lid_entries,
                     global_ids, local_ids, 
                     wdim, objwgts, ierr);
  }
  else if ((zz->Get_First_Obj != NULL) && (zz->Get_Next_Obj != NULL)){
    /* Use iterator functions to loop through object list */
    if (zz->Get_First_Obj(zz->Get_First_Obj_Data, 
                          num_gid_entries, num_lid_entries, 
                          global_ids, local_ids, 
                          wdim, objwgts, ierr)){
      /* Determine the number of objects since we don't trust the user
         to write the Get_Next_Obj query function in a safe way! */
      n = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, ierr);
      i = 0;
      while (!(*ierr) && (i<n-1)){ 
        gid_off = i * num_gid_entries;
        lid_off = i * num_lid_entries;
        lid = (num_lid_entries ? &(local_ids[lid_off]) : NULL);
        next_lid = (num_lid_entries ? &(local_ids[lid_off+num_lid_entries]) 
                                    : NULL);
        zz->Get_Next_Obj(zz->Get_Next_Obj_Data, 
                         num_gid_entries, num_lid_entries, 
                         &(global_ids[gid_off]), lid, 
                         &(global_ids[gid_off+num_gid_entries]),
                         next_lid,
                         wdim, &(objwgts[(i+1)*wdim]), ierr);
        i++;
      }
    }
  }
  else { /* No way to get objects */
    *ierr = ZOLTAN_FATAL;
  }
}
