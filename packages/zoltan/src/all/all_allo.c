// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include "zoltan_util.h"
#include "all_allo_const.h"
#include "params_const.h"

/* Fortran memory allocation callback functions */

static ZOLTAN_FORT_MALLOC_INT_FN *Zoltan_Fort_Malloc_int;
static ZOLTAN_FORT_FREE_INT_FN *Zoltan_Fort_Free_int;
static ZOLTAN_FORT_MALLOC_SET_STRUCT_FN *Zoltan_Fort_Malloc_Set_Struct;

int Zoltan_Set_Malloc_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    PARAM_VARS malloc_params[] = {
	{ "DEBUG_MEMORY", NULL, "INT", 0 },
	{ NULL, NULL, NULL, 0 }
    };

    status = Zoltan_Check_Param(name, val, malloc_params, &result, &index);
    if (status == 0 && index == 0) {
	Zoltan_Memory_Debug(result.ival);
	status = 3;
    }

    return(status);
}


/******************************************************************************
 * Special allocation for routines that allocate an array and return pointer.
 *
 * Zoltan_Special_Malloc allows the allocation to be done from either C or Fortran.
 *
 * Zoltan_Special_Free frees memory allocated by Zoltan_Special_Malloc
 *
 * Zoltan_Register_Fort_Malloc is called by the wrappers for the Fortran
 * interface to provide pointers to the Fortran allocation/free routines.
 *
 * int Zoltan_Special_Malloc(ZZ *zz, void **array, int size,
 *                       ZOLTAN_SPECIAL_MALLOC_TYPE type)
 *
 *   zz    -- the Zoltan structure in use
 *   array -- int**; returned as a
 *            pointer to the allocated space
 *   size  -- number of elements to be allocated in the array
 *   type  -- the type of array; ZOLTAN_SPECIAL_MALLOC_INT, ZOLTAN_SPECIAL_MALLOC_GID,
 *            or ZOLTAN_SPECIAL_MALLOC_LID
 *
 * The return value is 1 if the allocation succeeded, 0 if it failed.
 *
 * int Zoltan_Special_Free(ZZ *zz, void **array,
                       ZOLTAN_SPECIAL_MALLOC_TYPE type)
 *
 *****************************************************************************/

void Zoltan_Register_Fort_Malloc(ZOLTAN_FORT_MALLOC_INT_FN *fort_malloc_int,
                                 ZOLTAN_FORT_FREE_INT_FN *fort_free_int,
				 ZOLTAN_FORT_MALLOC_SET_STRUCT_FN *fort_malloc_set_struct)
{
   Zoltan_Fort_Malloc_int = fort_malloc_int;
   Zoltan_Fort_Free_int = fort_free_int;
   Zoltan_Fort_Malloc_Set_Struct = fort_malloc_set_struct;
}

int Zoltan_Special_Malloc(ZZ *zz, void **array, int size,
                      ZOLTAN_SPECIAL_MALLOC_TYPE type)
{
   int *ret_addr, success;
   char *yo = "Zoltan_Special_Malloc";

   success = 1;
   if (zz->Fortran) {

/* allocation from Fortran */

      switch(type) {
      case ZOLTAN_SPECIAL_MALLOC_INT:
#ifdef PGI /* special case for PGI Fortran compiler */
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU /* special case for Fujitsu and Lahey Fortran compilers */
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         if (ret_addr==0) success=0;
         break;
      case ZOLTAN_SPECIAL_MALLOC_GID:
         size *= zz->Num_GID;
#ifdef PGI
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         if (ret_addr==0) success=0;
         break;
      case ZOLTAN_SPECIAL_MALLOC_LID:
         size *= zz->Num_LID;
#ifdef PGI
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
         Zoltan_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         if (ret_addr==0) success=0;
         break;
      default:
	 ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Illegal value passed for type");
         success = 0;
      }
      if (success) {
         array[0] = ret_addr;
      }else{
         array[0] = NULL;
      }

   }else{

/* allocation from C */

      switch(type) {
      case ZOLTAN_SPECIAL_MALLOC_INT:
         *array = (int *) ZOLTAN_MALLOC(size*sizeof(int));
         if (*array==NULL) success=0;
         break;
      case ZOLTAN_SPECIAL_MALLOC_GID:
         *array = ZOLTAN_MALLOC_GID_ARRAY(zz, size);
         if (*array==NULL) success=0;
         break;
      case ZOLTAN_SPECIAL_MALLOC_LID:
         *array = ZOLTAN_MALLOC_LID_ARRAY(zz, size);
         if (zz->Num_LID > 0 && *array==NULL) success = 0;
         break;
      default:
	 ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Illegal value passed for type");
         *array = NULL;
         success = 0;
      }
   }
   return success;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Special_Free(ZZ *zz, void **array,
                    ZOLTAN_SPECIAL_MALLOC_TYPE type)
{
   int success;
   char *yo = "Zoltan_Special_Free";

   success = 1;
   if (zz->Fortran) {

/* deallocation from Fortran */

      switch(type) {
      case ZOLTAN_SPECIAL_MALLOC_INT:
#ifdef PGI /* special case for PGI Fortran compiler */
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU /* special case for Fujitsu and Lahey Fortran compilers */
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
         Zoltan_Fort_Free_int((int *)(array[1]));
#endif
#endif
         break;
      case ZOLTAN_SPECIAL_MALLOC_GID:
#ifdef PGI
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
         Zoltan_Fort_Free_int((int *)(array[1]));
#endif
#endif
         break;
      case ZOLTAN_SPECIAL_MALLOC_LID:
#ifdef PGI
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU
         Zoltan_Fort_Free_int((int *)(array[1]),array[2]);
#else
         Zoltan_Fort_Free_int((int *)(array[1]));
#endif
#endif
         break;
      default:
	 ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Illegal value passed for type");
         success = 0;
      }

   }else{

/* deallocation from C */

      ZOLTAN_FREE(array);
   }
   return success;
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Special_Fort_Malloc_Set_Struct(int *zz_addr_bytes, int **fort_zz) 
{
  Zoltan_Fort_Malloc_Set_Struct(zz_addr_bytes, fort_zz);
  return 1;
}

/*****************************************************************************/
/*                      END of all_allo.c                                     */
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
