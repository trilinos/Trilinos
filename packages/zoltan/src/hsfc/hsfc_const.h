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

#ifndef ZOLTAN_HSFC_CONST_H
#define ZOLTAN_HSFC_CONST_H

/* function prototypes */
int  Zoltan_HSFC_Point_Drop (ZZ *zz, double *x) ;
void Zoltan_HSFC_Box_Drop(ZZ *zz, int *array, double *lo, double *hi, int *n) ;
void Zoltan_HSFC_Free_Structure (ZZ *zz) ;
int  Zoltan_HSFC_Set_Param (char *name, char *val) ;

#endif
