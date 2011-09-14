
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Time.h"

#ifdef EPETRA_MPI
#include <time.h>
#endif
#ifdef EPETRA_HAVE_OMP
#include <omp.h>
#endif

//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Comm& Comm) 
  : StartTime_(0.0),
    Comm_(&Comm)
{
  StartTime_ = WallTime();
}
//=============================================================================
Epetra_Time::Epetra_Time(const Epetra_Time& Time) 
  : StartTime_(Time.StartTime_),
    Comm_(Time.Comm_)
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

	int mpiInitialized;
	MPI_Initialized(&mpiInitialized);

	if( mpiInitialized ) {

		return(MPI_Wtime());

	}
	else {

		clock_t start;

		start = clock();
		return( (double)( start ) / CLOCKS_PER_SEC );

	}

#else

#ifdef EPETRA_HAVE_OMP
       return(omp_get_wtime());
#else
#if ICL || defined(_WIN32)

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
#endif // MINGW

#endif // ICL || WIN32

#endif // EPETRA_HAVE_OMP
#endif // EPETRA_MPI

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
