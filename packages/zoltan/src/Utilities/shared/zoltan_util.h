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


#ifndef __ZOLTAN_UTIL_H
#define __ZOLTAN_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************   TRILINOS BUILD ENVIRONMENT  *******************/
/* This block should be executed only for an Autotools build. */
#ifndef TRILINOS_NO_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package 
 * and need to
 * be undef'd here to avoid warnings when this file is included from another 
 * package.
 * KL 11/25/02
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

/* This file passes values from configure to the source code. */
#include "Zoltan_config.h"

#ifdef HAVE_PARMETIS
#define ZOLTAN_PARMETIS
#endif

#ifdef HAVE_SCOTCH
#define ZOLTAN_SCOTCH
#endif

#ifdef HAVE_PATOH
#define ZOLTAN_PATOH
#endif

#ifdef HAVE_DRUM
#define ZOLTAN_DRUM
#endif

#ifdef HAVE_PARKWAY
#define ZOLTAN_PARKWAY
#endif

#ifdef HAVE_ZOLTAN_OCT
#define ZOLTAN_OCT
#endif

#endif /* TRILINOS_NO_CONFIG_H */

#define ZOLTAN_HIER
/*****************************************************************************/
/* 
 *  Macros and definitions that are common to all Zoltan modules and 
 *  utilities.
 */
/*****************************************************************************/

/*****************************************************************************/
/*
 *  Macros for consistently printing error and warning messages.
 */
/*****************************************************************************/

#define ZOLTAN_PRINT_ERROR(proc,yo,str) \
  fprintf(stderr, "[%d] Zoltan ERROR in %s (line %d of %s):  %s\n", \
          proc, yo, __LINE__, __FILE__, str);

#define ZOLTAN_PRINT_WARN(proc,yo,str) \
  fprintf(stderr, "[%d] Zoltan WARNING in %s (line %d of %s):  %s\n", \
          proc, yo, __LINE__, __FILE__, str);

#define ZOLTAN_TRACE(proc,where,yo,str) \
  printf("ZOLTAN (Processor %d) %s %s  %s\n", (proc), (where), (yo), \
         ((str) != NULL ? (char *)(str) : " "));

#define ZOLTAN_TRACE_IN(proc,yo,str) \
  ZOLTAN_TRACE((proc),"Entering",(yo),(str));

#define ZOLTAN_TRACE_OUT(proc,yo,str) \
  ZOLTAN_TRACE((proc),"Exiting",(yo),(str));

#define ZOLTAN_PRINT_INFO(proc,yo,str) \
  printf("ZOLTAN (Processor %d) %s: %s\n", (proc), (yo), \
         ((str) != NULL ? (char *)(str) : " "));


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__ZOLTAN_UTIL_H */
