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

#ifndef __HSFC_HILBERT_CONST_H
#define __HSFC_HILBERT_CONST_H

/* Bits per unsigned word */

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  \
  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

extern void Zoltan_HSFC_fhsfc2d(double coord[], int nkey, unsigned int key[]);
extern void Zoltan_HSFC_fhsfc3d(double coord[], int nkey, unsigned int key[]);

/* The following functions were added by rheaphy for use in hsfc */

extern double Zoltan_HSFC_IHilbert1d(double *coord); /* allows 1D problems */
extern double Zoltan_HSFC_IHilbert2d(double *coord);
extern double Zoltan_HSFC_IHilbert3d(double *coord);

#endif
