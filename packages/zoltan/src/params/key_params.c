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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include "zz_util_const.h"
#include "key_params.h"
#include "params_const.h"
#include "timer_const.h"
#include "zz_rand.h"
#include "third_library_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static PARAM_VARS Key_params[] = {
  { "IMBALANCE_TOL", NULL, "FLOAT", 1 },
  { "AUTO_MIGRATE", NULL, "INT", 0 },
  { "OBJ_WEIGHT_DIM", NULL, "INT", 0 },
  { "EDGE_WEIGHT_DIM", NULL, "INT", 0 },
  { "DEBUG_LEVEL", NULL, "INT", 0 },
  { "DEBUG_PROCESSOR", NULL, "INT", 0 },
  { "DETERMINISTIC", NULL, "INT", 0 },
  { "TIMER", NULL, "STRING", 0 },
  { "NUM_GID_ENTRIES", NULL, "INT", 0 },
  { "NUM_LID_ENTRIES", NULL, "INT", 0 },
  { "RETURN_LISTS", NULL, "STRING", 0 },
  { "LB_METHOD", NULL, "STRING", 0 },
  { "TFLOPS_SPECIAL", NULL, "INT", 0 },
  { "COMM_WEIGHT_DIM", NULL, "INT", 0 }, /* For backward compatibility only. */
                                         /* Prefer use of EDGE_WEIGHT_DIM.   */
  { "NUM_GLOBAL_PARTS", NULL, "INT", 0 },
  { "NUM_GLOBAL_PARTITIONS", NULL, "INT", 0 }, /* Deprecated */
  { "NUM_LOCAL_PARTS", NULL, "INT", 0 },
  { "NUM_LOCAL_PARTITIONS", NULL, "INT", 0 },  /* Deprecated */
  { "MIGRATE_ONLY_PROC_CHANGES", NULL, "INT", 0 },
  { "REMAP", NULL, "INT", 0 },
  { "SEED", NULL, "INT", 0 },
  { "LB_APPROACH", NULL, "STRING", 0 },
  { NULL, NULL, NULL, 0 } };
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* 
 * Handle parameter changes for variables stored in Zoltan structure.
 */

int Zoltan_Set_Key_Param(
ZZ *zz,                         /* Zoltan structure */
const char *name,		/* name of variable */
const char *val,		/* value of variable */
int  idx 			/* index of vector param, -1 if scalar */
)
{
    char *yo = "Zoltan_Set_Key_Param";
    char msg[256];
    int status;			/* return code */
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    int tmp;
    int export, import;

    status = Zoltan_Check_Param(name, val, Key_params, &result, &index);

    if (status == 0) {

      switch (index) {

      case 0:  		/* Imbalance_Tol */
        if (result.def) 
            result.fval = ZOLTAN_LB_IMBALANCE_TOL_DEF;
	if (result.fval < 1.0) {
	    sprintf(msg, "Invalid Imbalance_Tol value (%g) "
		"being set to %g.", result.fval, ZOLTAN_LB_IMBALANCE_TOL_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.fval = ZOLTAN_LB_IMBALANCE_TOL_DEF;
	}
        if (idx > zz->Obj_Weight_Dim){
          sprintf(msg, "Imbalance_Tol index %d > Obj_Weight_Dim = %d\n",
            idx, zz->Obj_Weight_Dim);
          ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
        }
        else if (idx < -1){
          sprintf(msg, "Invalid Imbalance_Tol index %d\n", idx);
          ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
        }
        else if (idx == -1){
          /* Set all entries to the same value. */
          for (idx=0; idx<zz->LB.Imb_Tol_Len; idx++)
	    zz->LB.Imbalance_Tol[idx] = result.fval;
        }
        else
	  zz->LB.Imbalance_Tol[idx] = result.fval;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 1:		/* Help_Migrate */
        if (result.def)
            result.ival = ZOLTAN_AUTO_MIGRATE_DEF;
	zz->Migrate.Auto_Migrate = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 2:		/* Object weight dim.  */
        if (result.def)
            result.ival = ZOLTAN_OBJ_WEIGHT_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Obj_Weight_Dim value (%d) "
		"being set to %d.", result.ival, ZOLTAN_OBJ_WEIGHT_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_OBJ_WEIGHT_DEF;
	}
	zz->Obj_Weight_Dim = result.ival;
        if (zz->Obj_Weight_Dim > zz->LB.Imb_Tol_Len){
          /* Resize and reallocate Imb_Tol. */
          zz->LB.Imb_Tol_Len += 10;
          zz->LB.Imbalance_Tol = (float *) ZOLTAN_REALLOC(zz->LB.Imbalance_Tol, zz->LB.Imb_Tol_Len * sizeof(float));
        }
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 3: 		/* Edge weight dim.  */
      case 13:
        if (result.def)
            result.ival = ZOLTAN_EDGE_WEIGHT_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Edge_Weight_Dim value (%d) "
		"being set to %d.", result.ival, ZOLTAN_EDGE_WEIGHT_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_EDGE_WEIGHT_DEF;
	}
	zz->Edge_Weight_Dim = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 4: 		/* Debug level  */
        if (result.def)
            result.ival = ZOLTAN_DEBUG_LEVEL_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Debug_Level value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DEBUG_LEVEL_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_DEBUG_LEVEL_DEF;
	}
	zz->Debug_Level = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 5: 		/* Debug processor  */
        if (result.def)
            result.ival = ZOLTAN_DEBUG_PROC_DEF;
	if (result.ival < 0 || result.ival > zz->Num_Proc) {
	    sprintf(msg, "Invalid Debug_Processor value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DEBUG_PROC_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_DEBUG_PROC_DEF;
	}
	zz->Debug_Proc = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;
       
      case 6: 		/* Deterministic flag */
        if (result.def)
            result.ival = ZOLTAN_DETERMINISTIC_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Deterministic value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DETERMINISTIC_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_DETERMINISTIC_DEF;
	}
	zz->Deterministic = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 7: 		/* Timer */
	status = Zoltan_Set_Timer_Param(name, val, &tmp);
        zz->Timer = tmp;
        Zoltan_Timer_ChangeFlag(zz->ZTime, zz->Timer);

	if (status==0) status = 3;	/* Don't add to Params field of ZZ */
        break;

      case 8:           /* Num_GID_Entries */
        if (result.def)
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        if (result.ival < 1) {
	    sprintf(msg, "Invalid Num_GID_Entries value (%d); "
		"being set to %d.", result.ival, ZOLTAN_NUM_ID_ENTRIES_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        }
        zz->Num_GID = result.ival;
        status = 3;
        break;

      case 9:           /* Num_LID_Entries */
        if (result.def)
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        if (result.ival < 0) {
	    sprintf(msg, "Invalid Num_LID_Entries value (%d); "
		"being set to %d.", result.ival, ZOLTAN_NUM_ID_ENTRIES_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        }
        zz->Num_LID = result.ival;
        status = 3;
        break;

      case 10:          /* LB.Return_Lists */
        export = (strstr(result.sval, "EXPORT") != NULL);
        import = (strstr(result.sval, "IMPORT") != NULL);
        if ((export && import) || (strcmp(result.sval, "ALL") == 0)) {
          tmp = ZOLTAN_LB_ALL_LISTS;  /* export AND import lists */
          status = 3;
        }
        else if (import){
          tmp = ZOLTAN_LB_IMPORT_LISTS;  /* import lists */
          status = 3;
        }
        else if (export){
          tmp = ZOLTAN_LB_EXPORT_LISTS;  /* export lists */
          status = 3;
        }
        else if (strstr(result.sval, "PART")!=NULL) {
          /* list of every object's part assignment */
          tmp = ZOLTAN_LB_COMPLETE_EXPORT_LISTS; 
          status = 3;
        }
        else if (strcmp(result.sval, "NONE")==0) {
          tmp = ZOLTAN_LB_NO_LISTS;      /* no lists */
          status = 3;
        }
        else if (strcmp(result.sval, "CANDIDATE_LISTS")==0) {
          tmp = ZOLTAN_LB_CANDIDATE_LISTS;   /* candidates needed in matching */
          status = 3;
        }
        else{
          tmp = ZOLTAN_LB_RETURN_LISTS_DEF;
          sprintf(msg, "Unknown return_lists option %s.", result.sval);
          ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
          status = 2; /* Illegal parameter */
        }
	zz->LB.Return_Lists = tmp;
        break;

      case 11:          /* LB_Method */
        status = Zoltan_LB_Set_LB_Method(zz,result.sval);
        if (status == ZOLTAN_OK)
          status = 3;
        break;

      case 12: 		/* Tflops Special flag */
        if (result.def)
            result.ival = ZOLTAN_TFLOPS_SPECIAL_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Tflops Special value (%d) "
		"being set to %d.", result.ival, ZOLTAN_TFLOPS_SPECIAL_DEF);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	    result.ival = ZOLTAN_TFLOPS_SPECIAL_DEF;
	}
	zz->Tflops_Special = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 14:          /* Num_Global_Parts */
      case 15:
        if (result.def)
            result.ival = zz->Num_Proc;
        if (result.ival < 1) {
	    sprintf(msg, "Invalid Num_Global_Parts value (%d); "
		"being set to %d.", result.ival,zz->Num_Proc);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
            result.ival = zz->Num_Proc;
        }
        zz->LB.Num_Global_Parts_Param = result.ival;
        status = 3;
        break;

      case 16:          /* Num_Local_Parts */
      case 17:
        if (result.def)
            result.ival = -1;
        if (result.ival < -1) {
	    sprintf(msg, "Invalid Num_Local_Parts value (%d); "
		"being set to %d.", result.ival,-1);
            ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
            result.ival = -1;
        }
        zz->LB.Num_Local_Parts_Param = result.ival;
        status = 3;
        break;

      case 18:		/* Migrate_Only_Proc_Changes */
        if (result.def)
            result.ival = ZOLTAN_MIGRATE_ONLY_PROC_CHANGES_DEF;
	zz->Migrate.Only_Proc_Changes = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 19:		/* LB.Remap */
        if (result.def)
            result.ival = 0;
	zz->LB.Remap_Flag = result.ival;
	status = 3;		/* Don't add to Params field of ZZ */
        break;

      case 20:          /* Seed */
        if (result.def)
            result.ival = Zoltan_Seed();
        zz->Seed = result.ival;
        Zoltan_Srand(result.ival, NULL);
        status = 3;
        break;

      case 21:          /* LB_APPROACH */
        if (result.def)
          strcpy(result.sval, ZOLTAN_LB_APPROACH_DEF);
        strcpy(zz->LB.Approach, result.sval);
        status = 3;
        break;

      }  /* end switch (index) */
    }

    return(status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Print key parameters.
 */
void Zoltan_Print_Key_Params(ZZ const *zz)
{
  int i;
  for (i=0; i<(zz->Obj_Weight_Dim?zz->Obj_Weight_Dim:1); i++)
    printf("ZOLTAN Parameter %s[%1d] = %f\n", Key_params[0].name, 
         i, zz->LB.Imbalance_Tol[i]);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[1].name, 
         (zz->Migrate.Auto_Migrate ? "TRUE" : "FALSE"));
  printf("ZOLTAN Parameter %s = %d\n", Key_params[18].name, 
         zz->Migrate.Only_Proc_Changes);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[2].name, 
         zz->Obj_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[3].name, 
         zz->Edge_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[4].name, 
         zz->Debug_Level);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[5].name, 
         zz->Debug_Proc);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[6].name, 
         (zz->Deterministic ? "TRUE" : "FALSE"));
  printf("ZOLTAN Parameter %s = %d ", Key_params[7].name, zz->Timer);
  if (zz->Timer == ZOLTAN_TIME_WALL)
     printf("(wall)");
  else if (zz->Timer == ZOLTAN_TIME_CPU)
     printf("(cpu)");
  printf("\n");
  printf("ZOLTAN Parameter %s = %d\n", Key_params[8].name, 
         zz->Num_GID);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[9].name, 
         zz->Num_LID);
  printf("ZOLTAN Parameter %s = ", Key_params[10].name);
  switch (zz->LB.Return_Lists) {
  case ZOLTAN_LB_ALL_LISTS:
    printf("IMPORT AND EXPORT\n");
    break;
  case ZOLTAN_LB_IMPORT_LISTS:
    printf("IMPORT\n");
    break;
  case ZOLTAN_LB_EXPORT_LISTS:
    printf("EXPORT\n");
    break;
  case ZOLTAN_LB_COMPLETE_EXPORT_LISTS:
    printf("PARTITION ASSIGNMENTS\n");
    break;
  case ZOLTAN_LB_CANDIDATE_LISTS:
    printf("CANDIDATE LISTS\n");
    break;
  case ZOLTAN_LB_NO_LISTS:
    printf("NONE\n");
    break;
  }
  if (zz->Tflops_Special)   /* print only if set */
     printf("ZOLTAN Parameter %s = %s\n", Key_params[12].name, "TRUE");
  printf("ZOLTAN Parameter %s = %d\n", Key_params[14].name, 
         zz->LB.Num_Global_Parts_Param);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[16].name, 
         zz->LB.Num_Local_Parts_Param);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[19].name, 
         zz->LB.Remap_Flag);
  printf("ZOLTAN Parameter %s = %d (%u)\n", Key_params[20].name, 
         Zoltan_Seed(), Zoltan_Seed());
  printf("ZOLTAN Parameter %s = %s\n", Key_params[21].name, 
         zz->LB.Approach);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Print some compile time configuration information.
 */
void Zoltan_Print_Configuration(char *indent)
{
  if (indent == NULL){
    indent = "";
  }
  printf("\n%sZOLTAN_ID_TYPE: %s (%ld bytes)\n",
    indent, zoltan_id_datatype_name, sizeof(ZOLTAN_ID_TYPE));
  printf("%sZOLTAN_GNO_TYPE: %s, (%ld bytes)\n",
    indent, zoltan_gno_datatype_name, sizeof(ZOLTAN_GNO_TYPE));
  printf("%sMPI_Datatype for ZOLTAN_ID_TYPE: %s\n", indent,
    zoltan_mpi_id_datatype_name);
  printf("%sMPI_Datatype for ZOLTAN_GNO_TYPE: %s\n", indent, 
    Zoltan_mpi_gno_name());

  /* Metis and ParMetis have different version numbers.  Some
   * older versions do not define version numbers.
   */
#ifdef ZOLTAN_PARMETIS
  printf("%sThird party library: ParMetis ", indent);
  printf("version %d.%d", PARMETIS_MAJOR_VERSION, PARMETIS_MINOR_VERSION);
  #ifdef PARMETIS_SUBMINOR_VERSION
    printf(".%d", PARMETIS_SUBMINOR_VERSION);
  #endif
  printf("\n");
#endif

#ifdef ZOLTAN_METIS
  #ifdef METISTITLE
    printf("%sThird party library: %s\n",indent, METISTITLE);
  #else
    printf("%sThird party library: METIS\n",indent);
  #endif
#endif

  /* Scotch and PTScotch have the same version number.  Version
   * numbers are not defined in older versions.
   */

#ifdef ZOLTAN_PTSCOTCH
  printf("%sThird party library: PTScotch ", indent);
  #ifdef SCOTCH_VERSION
    printf("version %d.%d.%d\n", 
           SCOTCH_VERSION, SCOTCH_RELEASE, SCOTCH_PATCHLEVEL);
  #endif
#endif

#ifdef ZOLTAN_SCOTCH
  printf("%sThird party library: Scotch ", indent);
  #ifdef SCOTCH_VERSION
    printf("version %d.%d.%d\n", 
           SCOTCH_VERSION, SCOTCH_RELEASE, SCOTCH_PATCHLEVEL);
  #endif
#endif
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
