/*
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
*/

#ifndef EPETRA_TIME_H
#define EPETRA_TIME_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



//! Epetra_Time:  The Epetra Timing Class.
/*! The Epetra_Time class is a wrapper that encapsulates the general
  information needed getting timing information.  Currently it return
  the elapsed time for each calling processor..
  A Epetra_Comm object is required for building all Epetra_Time objects.

  Epetra_Time support both serial execution and (via MPI) parallel
  distributed memory execution.  It is meant to insulate the user from
  the specifics of timing across a variety of platforms.
*/

#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#else
#if ICL || defined(_WIN32)
#include <time.h>
#else
#include <sys/time.h>
#ifndef MINGW
#include <sys/resource.h>
#endif
#endif
#endif

class EPETRA_LIB_DLL_EXPORT Epetra_Time: public Epetra_Object {

  public:
  //! Epetra_Time Constructor.
  /*! Creates a Epetra_Time instance. This instance can be queried for
      elapsed time on the calling processor.  StartTime is also set
      for use with the ElapsedTime function.
  */
  Epetra_Time(const Epetra_Comm & Comm);

  //! Epetra_Time Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Time instance.
  */
  Epetra_Time(const Epetra_Time& Time);

  //! Epetra_Time wall-clock time function.
  /*! Returns the wall-clock time in seconds.  A code section can be
      timed by putting it between two calls to WallTime and taking the
      difference of the times.
  */
  double WallTime(void) const;

  //! Epetra_Time function to reset the start time for a timer object.
  /*! Resets the start time for the timer object to the current time
      A code section can be
      timed by putting it between a call to ResetStartTime and ElapsedTime.
  */
  void ResetStartTime(void);

  //! Epetra_Time elapsed time function.
  /*! Returns the elapsed time in seconds since the timer object was
      constructed, or since the ResetStartTime function was called.
      A code section can be
      timed by putting it between the Epetra_Time constructor and a call to
      ElapsedTime, or between a call to ResetStartTime and ElapsedTime.
  */
  double ElapsedTime(void) const;

  //! Epetra_Time Destructor.
  /*! Completely deletes a Epetra_Time object.
  */
  virtual ~Epetra_Time(void);

  Epetra_Time& operator=(const Epetra_Time& src)
    {
      StartTime_ = src.StartTime_;
      Comm_ = src.Comm_;
      return( *this );
    }

 private:

  double StartTime_;
  const Epetra_Comm * Comm_;

};

#endif /* EPETRA_TIME_H */
