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

#ifndef __ZOLTAN_UTIL_H
#define __ZOLTAN_UTIL_H

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

#define ZOLTAN_TRACE_ENTER(proc,yo,str) \
  printf("ZOLTAN (Processor %d) Entering %s  %s\n", (proc), (yo), \
         ((str) ? (str) : " "));

#define ZOLTAN_TRACE_EXIT(proc,yo,str) \
  printf("ZOLTAN (Processor %d) Leaving %s  %s\n", (proc), (yo), \
         ((str) ? (str) : " "));

#define ZOLTAN_PRINT_INFO(proc,yo,str) \
  printf("ZOLTAN (Processor %d) %s: %s\n", (proc), (yo), \
         ((str) ? (str) : " "));


#endif /* !__ZOLTAN_UTIL_H */
