/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _PE_COMMON_H
#define _PE_COMMON_H

/*****************************************************************************/
/* This file contains defines that are common to both nem_spread, and        */
/* nem_join. Most of them were taken from rf_salsa.h which is unique         */
/* to nem_spread.                                                            */
/*****************************************************************************/

/*
 * Default value of the chunk size (in bytes) for use in I/O and message
 * passing
 */

#ifndef MAX_CHUNK_SIZE
#define MAX_CHUNK_SIZE 1073741824
/* Small message size for debugging purposes */
/* #define MAX_CHUNK_SIZE 16384 */
#endif

#define PEX_MAX(x, y) (((x) > (y)) ? (x) : (y)) /* max function */
#define PEX_MIN(x, y) (((x) < (y)) ? (x) : (y)) /* min function */

#endif /* _PE_COMMON_H */
