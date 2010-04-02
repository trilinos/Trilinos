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
/* to get PRIdMAX, etc */
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Set at runtime in Zoltan_set_mpi_integer_data_types */
extern MPI_Datatype _mpi_int_32_type;
extern MPI_Datatype _mpi_int_max_type;
extern MPI_Datatype _mpi_uint_32_type;
extern MPI_Datatype _mpi_uint_max_type;

typedef int32_t  Z_INT;
typedef uint32_t Z_UINT;

#ifdef USE_32_BIT_ADDRESS_SPACE
typedef int32_t Z_INT_L;
typedef uint32_t Z_UINT_L;
#else
typedef intmax_t  Z_INT_L;
typedef uintmax_t Z_UINT_L;
#endif

#define Z_MPI_INT       _mpi_int_32_type
#define Z_MPI_UNSIGNED  _mpi_uint_32_type
#define Z_MPI_LONG          _mpi_int_max_type
#define Z_MPI_UNSIGNED_LONG _mpi_uint_max_type

/* string to specify these types in scanf, printf, etc as a decimal
    printf("%" Z_INT_L_SPECIFIER "\n",mygid);
 */

#ifdef USE_32_BIT_ADDRESS_SPACE
#define Z_INT_L_SPECIFIER    PRId32
#define Z_UINT_L_SPECIFIER   PRId32
#else
#define Z_INT_L_SPECIFIER    PRIdMAX 
#define Z_UINT_L_SPECIFIER   PRIdMAX
#endif

#define Z_INT_SPECIFIER      PRId32
#define Z_UINT_SPECIFIER     PRId32

/* A constant of type Z_INT_L can be created with the macro INTMAX_C(val)
 *            of type Z_UINT_L                             UINTMAX_C(val)
 *            of type Z_INT                                  INT32_C(val)
 *            of type Z_UINT                                UINT32_C(val)
 */

/*****************************************************************************/
/*
 *  Data type ZOLTAN_ID for identifiers used in Zoltan.
 */
/*****************************************************************************/

typedef Z_INT_L            ZOLTAN_ID_TYPE;
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
