/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#include <stdio.h>
#include "lb_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static PARAM_VARS Key_params[] = {
	{ "IMBALANCE_TOL", NULL, "DOUBLE" },
	{ "AUTO_MIGRATE", NULL, "INT" },
	{ "OBJ_WEIGHT_DIM", NULL, "INT" },
	{ "COMM_WEIGHT_DIM", NULL, "INT" },
	{ "DEBUG_LEVEL", NULL, "INT" },
	{ "DETERMINISTIC", NULL, "INT" },
	{ NULL, NULL, NULL } };
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* 
 * Handle parameter changes for variables stored in LB structure.
 */

int LB_Set_Key_Param(
LB *lb,                         /* load balance structure */
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;			/* return code */
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */

    status = LB_Check_Param(name, val, Key_params, &result, &index);

    if (status == 0) {

      switch (index) {

      case 0:  		/* Imbalance_Tol */
	if (result.dval < 1.0) {
	    fprintf(stderr, "WARNING: Invalid Imbalance_Tol value (%g) "
		"being set to 1.0\n", result.dval);
	    result.dval = 1.0;
	}
	lb->Imbalance_Tol = result.dval;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 1:		/* Help_Migrate */
	lb->Migrate.Auto_Migrate = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 2:		/* Object weight dim.  */
	if (result.ival < 0) {
	    fprintf(stderr, "WARNING: Invalid Obj_Weight_Dim value (%d) "
		"being set to 0\n", result.ival);
	    result.ival = 0;
	}
	lb->Obj_Weight_Dim = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 3:		/* Communication weight dim.  */
	if (result.ival < 0) {
	    fprintf(stderr, "WARNING: Invalid Comm_Weight_Dim value (%d) "
		"being set to 0\n", result.ival);
	    result.ival = 0;
	}
	lb->Comm_Weight_Dim = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 4: 		/* Debug level  */
	if (result.ival < 0) {
	    fprintf(stderr, "WARNING: Invalid Debug_Level value (%d) "
		"being set to 0\n", result.ival);
	    result.ival = 0;
	}
	lb->Debug_Level = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;
       
      case  5: 		/* Deterministic flag */
	if (result.ival < 0) {
	    fprintf(stderr, "WARNING: Invalid Deterministic value (%d) "
		"being set to TRUE\n", result.ival);
	    result.ival = TRUE;
	}
	lb->Deterministic = result.ival;
	status = 3;		/* Don't add to Params field of LB */
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
void LB_Print_Key_Params(LB *lb)
{
  printf("ZOLTAN Parameter %s = %f\n", Key_params[0].name, 
         lb->Imbalance_Tol);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[1].name, 
         (lb->Migrate.Auto_Migrate ? "TRUE" : "FALSE"));
  printf("ZOLTAN Parameter %s = %d\n", Key_params[2].name, 
         lb->Obj_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[3].name, 
         lb->Comm_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[4].name, 
         lb->Debug_Level);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[5].name, 
         (lb->Deterministic ? "TRUE" : "FALSE"));
}
