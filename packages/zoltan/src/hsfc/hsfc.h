/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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
#include <float.h>

#include "zz_const.h"
#include "timer_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "zoltan_util.h"  
#include "hsfc_hilbert_const.h"
#include "hsfc_const.h"



static const double  HSFC_EPSILON     = 1.0e-7L ;    /* Andy's value */
static const double  REFINEMENT_LIMIT = 10.0L * DBL_EPSILON ;   /* bin can't be divided */
static const double  DEFAULT_WEIGHT   = 1.0 ;        /* when dots have no weight */
static const int N = 8 ; /* number of "bins per processor", small positive integer */
static const int MAX_LOOPS = 16 ; /* refinement now hitting limit of double precision */

/* The following are used in determining the max/min hsfc coordinate in a bin, if bin empty */
static const double  DEFAULT_BIN_MAX = -2.0 ;
static const double  DEFAULT_BIN_MIN =  2.0 ;


#define ZOLTAN_HSFC_ERROR(error,str) {err = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto free ;}


typedef struct Partition
   {
   double r ;            /* rightmost boundary of partition interval */
   double l ;            /* leftmost boundary of partition interval  */
   int index ;           /* number of "owner" of data in interval */
   } Partition ;         /* interval is half open, [l,r) */



typedef struct HSFC_Data
   {
   Partition *final_partition ;
   double     bbox_hi[3] ;        /* smallest bounding box, high point */
   double     bbox_lo[3] ;        /* smallest bounding box, low point */
   double     bbox_extent[3] ;    /* length of each side of bounding box */
   int        nloops ;            /* number of loops for load balancing */
   int        ndimension ;        /* number of dimensions in problem (2 or 3) */
   double    (*fhsfc)(double*) ;  /* space filling curve function */
   } HSFC_Data ;                  /* data preserved for point & box drop later */



typedef struct Dots
   {
   double fsfc ;         /* computed normalized SFC coordinate     */
   double x[3] ;         /* dots coordinates in problem domain     */
   float weight ;        /* scalar computed from weight vector     */
   int   part ;          /* partition owning dot                   */
   } Dots ;              /* represents objects being load-balanced */


extern int  Zoltan_HSFC_compare (const void *key, const void *arg) ;

#endif
