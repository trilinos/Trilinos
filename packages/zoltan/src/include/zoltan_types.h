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

#ifndef __ZOLTAN_TYPES_H
#define __ZOLTAN_TYPES_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*
 *  Data type ZOLTAN_ID for identifiers used in Zoltan.
 */
/*****************************************************************************/

typedef unsigned int        ZOLTAN_ID_TYPE;
typedef ZOLTAN_ID_TYPE     *ZOLTAN_ID_PTR;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED

/*****************************************************************************/
/*
 * Error codes for Zoltan library
 *   ZOLTAN_OK     - no errors
 *   ZOLTAN_WARN   - some warning occurred in Zoltan library;
 *                   application should be able to continue running
 *   ZOLTAN_FATAL  - a fatal error occurred
 *   ZOLTAN_MEMERR - memory allocation failed; with this error, it could be
 *                   possible to try a different, more memory-friendly,
 *                   algorithm.
 */
/*****************************************************************************/
#define ZOLTAN_OK     0
#define ZOLTAN_WARN   1
#define ZOLTAN_FATAL  -1
#define ZOLTAN_MEMERR -2

/*****************************************************************************/
/* Hypergraph query function types 
 */
/*****************************************************************************/
#define ZOLTAN_COMPRESSED_EDGE   1
#define ZOLTAN_COMPRESSED_VERTEX 2

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__ZOLTAN_TYPES_H */
