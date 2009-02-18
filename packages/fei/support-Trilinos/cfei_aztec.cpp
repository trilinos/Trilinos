/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <fei_SharedPtr.hpp>
#include <fei_Data.hpp>
#include <fei_LinearSystemCore.hpp>

#include <fei_mpi.h>

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif
#ifndef PLL
#define PLL PLL
#endif

#endif

#include <az_aztec.h>

#include "cfei_aztec.h"

#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_Vector.hpp>
#include <fei_AztecDMSR_Matrix.hpp>
#include <fei_Aztec_BlockMap.hpp>
#include <fei_AztecDVBR_Matrix.hpp>

#include <fei_Aztec_LinSysCore.hpp>

/*============================================================================*/
/* Create function for a Aztec_LinSysCore object.
*/
extern "C" int Aztec_LinSysCore_create(LinSysCore** lsc, MPI_Comm comm)
{
   fei::SharedPtr<LinearSystemCore>* linSys =
     new fei::SharedPtr<LinearSystemCore>(new fei_trilinos::Aztec_LinSysCore(comm));

   if (linSys->get() == NULL) return(1);

   *lsc = new LinSysCore;

   (*lsc)->lsc_ = (void*)linSys;

   return(0);
}

/*============================================================================*/
/* Destroy function, to de-allocate a Aztec_LinSysCore object.
*/
extern "C" int Aztec_LinSysCore_destroy(LinSysCore** lsc)
{
   if (*lsc == NULL) return(1);

   fei::SharedPtr<LinearSystemCore>* linSys =
     (fei::SharedPtr<LinearSystemCore>*)((*lsc)->lsc_);

   if (linSys->get() == NULL) return(1);

   delete linSys;

   delete *lsc;
   *lsc = NULL;

   return(0);
}

#endif
//HAVE_FEI_AZTECOO
