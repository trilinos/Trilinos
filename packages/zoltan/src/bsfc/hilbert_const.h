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

#ifndef __HILBERT_CONST_H
#define __HILBERT_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Bits per unsigned word */

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  \
  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

extern void Zoltan_BSFC_fhsfc2d(double coord[], unsigned *nkey, unsigned key[]);
extern void Zoltan_BSFC_fhsfc3d(double coord[], unsigned *nkey, unsigned key[]);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif
