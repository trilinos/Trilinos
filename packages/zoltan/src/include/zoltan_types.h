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
#include <unistd.h>
#include <limits.h>

/* int64_t is needed by 64-bit PT-Scotch header file */

#ifdef _WIN32
typedef __int64 int64_t;
#else
#include <stdint.h>
#endif

#ifndef ssize_t
typedef long ssize_t;
#endif

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "Zoltan_config.h"

/* The default ZOLTAN_ID_TYPE is "unsigned int" but this can be over-ridden on the compile command line.  
 *
 * The type of a Zoltan object global ID is ZOLTAN_ID_TYPE.  A pointer to it is ZOLTAN_ID_PTR.
 *
 * It's decimal type specifier: printf(ZOLTAN_ID_SPEC "\n", global_id);
 *
 * A constant of the same type:   ZOLTAN_ID_TYPE global_id = ZOLTAN_ID_CONSTANT(0);  (Do we need this?)
 *
 * The MPI_Datatype for a ZOLTAN_ID_TYPE is a ZOLTAN_ID_MPI_TYPE.
 *
 * We assume the local number of objects fits in a 32 bit integer, but the global number may require
 * the maximum integer width available on the machine.
 *
 * ZOLTAN_GNO_TYPE  is the global number/count type.
 *
 * The underlying type is: ssize_t (signed size_t).  This will be 32 or 64
 *   bits depending on whether the machine is a 32 or 64 bit machine.  (We use ssize_t
 *   instead of intmax_t because intmax_t may still be
 *   64 bits on a 32 bit machine because the compiler constructs a 64 bit int.)
 *
 * It's decimal type specifier is "z":    printf("%zd\n",globalNum);
 *
 * The MPI_Datatype for ssize_t is returned by Zoltan_mpi_gno_type().
 *
 * We don't assume a pointer is the same size as any size of int.  If we want to store
 * a pointer in an int we use types intptr_t or uintptr_t.
 */

#undef HAVE_LONG_LONG_INT

#ifdef ULLONG_MAX
#define HAVE_LONG_LONG_INT
#endif

#undef ZOLTAN_ID_MPI_TYPE
#undef ZOLTAN_ID_SPEC
#undef ZOLTAN_ID_CONSTANT

/* 
 * Autoconf build: --with-id-type={uint, ulong, ullong}
 *
 * CMake build:    -D Zoltan_ENABLE_UINT_IDS:Bool=ON
 *                 -D Zoltan_ENABLE_ULONG_IDS:Bool=ON
 *                 -D Zoltan_ENABLE_ULLONG_IDS:Bool=ON
 *
 */

#ifndef UNSIGNED_INT_GLOBAL_IDS
#ifndef UNSIGNED_LONG_GLOBAL_IDS
#ifndef UNSIGNED_LONG_LONG_GLOBAL_IDS

#define UNSIGNED_INT_GLOBAL_IDS

#endif
#endif
#endif

#ifdef UNSIGNED_LONG_LONG_GLOBAL_IDS

#ifndef HAVE_LONG_LONG_INT

#undef UNSIGNED_LONG_LONG_GLOBAL_IDS
#define UNSIGNED_LONG_GLOBAL_IDS
#warning Global ID type "unsigned long long int" was requested at configure time.
#warning The compiler does not fully support this data type.
#warning Instead ZOLTAN_ID_TYPE will be "unsigned long".

#endif

#endif

#ifdef UNSIGNED_LONG_GLOBAL_IDS

typedef unsigned long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED_LONG
#define zoltan_mpi_id_datatype_name "MPI_UNSIGNED_LONG"
#define zoltan_id_datatype_name "unsigned long"
#define ZOLTAN_ID_SPEC  "%lu"
#define ZOLTAN_ID_CONSTANT(z)  z ## L
#define ZOLTAN_ID_INVALID  ULONG_MAX

#endif

#ifdef UNSIGNED_LONG_LONG_GLOBAL_IDS

typedef unsigned long long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_LONG_LONG_INT
#define zoltan_mpi_id_datatype_name "MPI_LONG_LONG_INT"
#define zoltan_id_datatype_name "unsigned long long"
#define ZOLTAN_ID_SPEC  "%Lu"
#define ZOLTAN_ID_CONSTANT(z)  z ## LL
#define ZOLTAN_ID_INVALID  ULLONG_MAX

#endif

#ifdef UNSIGNED_INT_GLOBAL_IDS

typedef unsigned int ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED
#define zoltan_mpi_id_datatype_name "MPI_UNSIGNED"
#define zoltan_id_datatype_name "unsigned int"
#define ZOLTAN_ID_SPEC  "%u"
#define ZOLTAN_ID_CONSTANT(z)  z
#define ZOLTAN_ID_INVALID  UINT_MAX

#endif

typedef ZOLTAN_ID_TYPE     *ZOLTAN_ID_PTR;

/* 
 * The MPI_Datatype for ZOLTAN_GNO_TYPE is returned by Zoltan_mpi_gno_type().
 */

#define ZOLTAN_GNO_TYPE      ssize_t

/*
* 12/21/10: Trilinos compiles with warnings about c99 features (%zd).  
*           So we use %ld to for a ZOLTAN_GNO_TYPE.
*/

#define ZOLTAN_GNO_SPEC   "%ld"

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
