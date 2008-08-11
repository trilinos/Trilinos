#ifndef _fei_LinSysCore_struct_h_
#define _fei_LinSysCore_struct_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
   Define a "Linear System Core" struct. This is the object that
   handles all solver-library-specific functionality like sumIntoMatrix,
   launchSolver, etc., etc. The pointer 'lsc_' needs to hold an instance
   of an object which implements the C++ interface defined in
   LinearSystemCore.h. An implementation-specific
   function will be required to create one of these.

   e.g., Aztec_LinSysCore_create(LinSysCore** lsc, MPI_Comm comm);

   This function would be found in ../support-Trilinos/cfei_aztec.h, in the case
   of an Aztec implementation.
------------------------------------------------------------------------------*/

struct LinSysCore_struct {
   void* lsc_;
};
typedef struct LinSysCore_struct LinSysCore;


/* The following macro is an artifact of a previous bad design where the
   above LinSysCore declaration was located in a bad place, and could be
   multiply defined. The following macro prevents the multiple-definition
   in headers (such as cfei_prometheus.h in the Prometheus library) that
   haven't been revised.
*/
#ifndef CFEI_LinSysCore_DEFINED
#define CFEI_LinSysCore_DEFINED
#endif

#endif

