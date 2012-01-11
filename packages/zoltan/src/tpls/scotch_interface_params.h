/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2008 Sandia National Laboratories.                          *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __SCOTCH_INTERFACE_PARAMS_H
#define __SCOTCH_INTERFACE_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

  /**********  parameters structure for Scotch methods **********/
static PARAM_VARS Scotch_params[] = {
  { "SCOTCH_METHOD", NULL, "STRING", 0 },
  { "SCOTCH_TYPE", NULL, "STRING", 0 },
  { "SCOTCH_STRAT", NULL, "STRING", 0 },
  { "SCOTCH_STRAT_FILE", NULL, "STRING", 0 },
  { NULL, NULL, NULL, 0 } };

#ifdef __cplusplus
}
#endif

#endif
