/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
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
#include "lb_const.h"
#include "lb_util_const.h"

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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Zoltan_Hash is a hash function for Zoltan ids (local or global). 
 *
 * Input:
 *   key: a key to hash of type ZOLTAN_ID_PTR
 *   num_id_entries: the number of (ZOLTAN_ID_TYPE-sized) entries of the key to use
 *   n: the range of the hash function is 0..n-1
 *
 * Return value:
 *   the hash value, an unsigned integer between 0 and n-1
 *
 * Algorithm: 
 *   This hash function is based on Don Knuth's golden ratio
 *   multiplicative method. Bitwise xor is used for keys
 *   longer than an int. The method works well for keys
 *   of size one or two ints, which is typically the case.
 *
 *   This hash function should be replaced with a stronger method
 *   if good hashing of a large number of keys is important.
 *
 * Author: 
 *   Erik Boman, eboman@cs.sandia.gov (SNL 9226)
 */


unsigned int Zoltan_Hash(ZOLTAN_ID_PTR key, int num_id_entries, unsigned int n)
{
  unsigned int h, rest, *p;
  char *byteptr;
  int bytes;
  int num_bytes = num_id_entries * sizeof(ZOLTAN_ID_TYPE);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= sizeof(int); 
       bytes-=sizeof(int), p++){
    h = (h*2654435761U) ^ (*p);
  }

  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Merge the two parts */
  if (rest)
    h = (h*2654435761U) ^ rest;

  /* Return h mod n */
  return (h%n);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Remove leading & trailing white space and convert to upper case. */

int Zoltan_Clean_String(
char *string1,			/* original string */
char **pstring2) 		/* cleaned string to return */
{

    char     *string2;		/* cleaned up string */
    int       start, end;	/* indices bounding true string */
    int       length1;		/* length of string 1 */
    int       i;		/* loop counter */

    length1 = strlen(string1);
    start = 0;
    end = length1;
    while (start < length1 && isspace((int)(string1[start])))
	start++;
    while (end > start && isspace((int)(string1[end])))
	end--;

    string2 = (char *) ZOLTAN_MALLOC((end - start + 1) * sizeof(char));
    *pstring2 = string2;

    if (string2 == NULL)
	return (ZOLTAN_MEMERR);

    for (i = start; i < end; i++) {
	*string2++ = toupper(string1[i]);
    }
    *string2 = '\0';

    return (ZOLTAN_OK);
}

