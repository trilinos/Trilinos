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
#ifndef __HG_UTIL_H
#define __HG_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define EPS                     1e-6
    
#define MIN(A,B)                (((A) < (B)) ? (A) : (B))
#define MAX(A,B)                (((A) > (B)) ? (A) : (B))

extern unsigned int Zoltan_HG_Rand (void);
extern void         Zoltan_HG_Srand (unsigned int);
extern void         Zoltan_HG_Rand_Perm_Int (int*, int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_UTIL_H */
