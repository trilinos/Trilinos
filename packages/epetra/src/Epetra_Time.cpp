
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Time.h"


//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Comm& Comm) 
  : Comm_(&Comm)
{
  StartTime_ = WallTime();
}
//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Time& Time) 
  : Comm_(Time.Comm_)
{
}
//=============================================================================
Epetra_Time::~Epetra_Time(void)  
{
}
//=============================================================================
double Epetra_Time::WallTime(void) const
{
#ifdef EPETRA_MPI

  return(MPI_Wtime());

#elif ICL

   clock_t start;
   //double duration;

   start = clock();
  return (double)( start ) / CLOCKS_PER_SEC;

#else

#ifndef MINGW
   struct timeval tp;
   static long start=0, startu;
   if (!start)
   {
      gettimeofday(&tp, NULL);
      start = tp.tv_sec;
      startu = tp.tv_usec;
      return(0.0);
   }
   gettimeofday(&tp, NULL);
   return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
#else
   return( (double) clock() / CLOCKS_PER_SEC );
#endif

#endif

}
//=============================================================================
void Epetra_Time::ResetStartTime(void)
{
  StartTime_ = WallTime();
  return;
}
//=============================================================================
double Epetra_Time::ElapsedTime(void) const
{
  return(WallTime()-StartTime_);
}
