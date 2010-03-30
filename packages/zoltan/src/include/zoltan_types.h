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
#include <inttypes.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

typedef int16_t Z_INT16;
typedef int32_t Z_INT32;
typedef uint16_t Z_UINT16;
typedef uint32_t Z_UINT32;

/* Data type macro definitions for 64 bit platform *************************/
#if __WORDSIZE==64

#ifdef USE_32_BIT_ADDRESS_SPACE       /* compile time option for faster code */
typedef int32_t Z_INT64;
typedef uint32_t Z_UINT64;
#else
typedef int64_t Z_INT64;
typedef uint64_t Z_UINT64;
#endif

/* Data type macro definitions for 32 bit platform *************************/
#else

typedef int32_t Z_INT64;
typedef uint32_t Z_UINT64;
 
#endif

/* We figure out which MPI datatypes to use at runtime in Zoltan_set_mpi_integer_data_types */

extern MPI_Datatype mpi_two_byte_int_type;
extern MPI_Datatype mpi_four_byte_int_type;
extern MPI_Datatype mpi_eight_byte_int_type;
extern MPI_Datatype mpi_two_byte_unsigned_int_type;
extern MPI_Datatype mpi_four_byte_unsigned_int_type;
extern MPI_Datatype mpi_eight_byte_unsigned_int_type;

#define Z_MPI_INT16  mpi_two_byte_int_type
#define Z_MPI_INT32  mpi_four_byte_int_type
#define Z_MPI_INT64  mpi_eight_byte_int_type
#define Z_MPI_UINT16  mpi_two_byte_unsigned_int_type
#define Z_MPI_UINT32  mpi_four_byte_unsigned_int_type
#define Z_MPI_UINT64  mpi_eight_byte_unsigned_int_type

/*****************************************************************************/
/*
 *  Data type ZOLTAN_ID for identifiers used in Zoltan.
 */
/*****************************************************************************/

typedef Z_INT64              ZOLTAN_ID_TYPE;
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
