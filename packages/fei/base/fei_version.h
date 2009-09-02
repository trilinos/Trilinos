/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_version_h_
#define _fei_version_h_

/*
  This 'fei_version_number' provides the version number of the
  fei implementation code. This number should be updated whenever
  an fei update is released.
*/
static const char fei_version_number[16] = {"2.24.01"};

/* IMPORTANT: Keep the version-number portion of the following macros
   synchronized with the above version number. These macros, which redefine
   prominent FEI symbols, guard against header-mismatch errors that arise
   if an application uses one version of FEI headers with another version
   of FEI libraries.
*/
#define fei_VERSION fei_2_24_01

#define FEI_MAJOR_VERSION 2
#define FEI_MINOR_VERSION 24
#define FEI_PATCH_VERSION 01

#define FEI_Implementation FEI_Implementation_2_24_01
#define FEI_create FEI_create_2_24_01

#endif // _fei_version_h_

