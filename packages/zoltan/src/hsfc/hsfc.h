/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/



#ifndef ZOLTAN_HSFC_H
#define ZOLTAN_HSFC_H

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <values.h>
#include <limits.h>
#include <string.h>

#include "zz_const.h"
#include "timer_const.h"
#include "lbi_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "zoltan_util.h"  
#include "hsfc_hilbert_const.h"
#include "hsfc_const.h"


static const int ZOLTAN_HSFC_OK                =   0  ;
static const int ZOLTAN_HSFC_ILLEGAL_DIMENSION = -100 ;
static const int ZOLTAN_HSFC_POINT_NOT_FOUND   = -101 ;
static const int ZOLTAN_HSFC_FATAL_ERROR       = -102 ;

static const double  ZOLTAN_HSFC_EPSILON = 1.0e-7L ;       /* Andy's value */

static PARAM_VARS HSFC_params[] =
   {
   {"KEEP_CUTS", NULL, "INT"},   {NULL, NULL, NULL}
   } ;

static const int TWO = 2 ;
static const int N   = 8 ; /* number of "bins per processor", small positive integer */

#define ZOLTAN_HSFC_ERROR(error,str) {err = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto free ;}


typedef struct Partition
   {
   double r ;            /* rightmost boundary of partition part */
   double l ;            /* leftmost boundary of partition part  */
   int index ;           /* number of "owner" of data falling in partition */
   } Partition ;



typedef struct HSFC_Data
   {
   Partition *final_partition ;
   double     bbox_hi[3] ;       /* smallest bounding box, high point */
   double     bbox_lo[3] ;       /* smallest bounding box, low point */
   double     bbox_extent[3] ;   /* length of each side of bounding box */
   int        ndimension ;       /* number of dimensions in problem (2 or 3) */
   void     (*fhsfc)(double*, unsigned*, unsigned*) ;   /* space filling curve function */
   } HSFC_Data ;                 /* data preserved for point & box drop later */



typedef struct Dots
   {
   double fsfc ;         /* computed normalized SFC coordinate     */
   double x[3] ;         /* dots coordinates in problem domain     */
   float weight ;        /* scalar computed from weight vector     */
   int   proc ;          /* processor owning dots                  */
   } Dots ;              /* represents objects being load-balanced */


/* local function prototypes */
static int compare (const void *key, const void *arg) ;
static double convert_key (unsigned int *key) ;

#endif
