
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

#ifndef _EPETRA_TIME_H_
#define _EPETRA_TIME_H_

//! Epetra_Time:  The Epetra Timing Class.
/*! The Epetra_Time class is a wrapper that encapsulates the general
  information needed getting timing information.  Currently it return
  the elapsed time for each calling processor..
  A Epetra_Comm object is required for building all Epetra_Time objects.
  
  Epetra_Time support both serial execution and (via MPI) parallel 
  distributed memory execution.  It is meant to insulate the user from
  the specifics of timing across a variety of platforms.
*/

#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#else
#include <sys/time.h>
#ifndef MINGW
#include <sys/resource.h>
#else
#include <winsock.h>
#endif
#endif

class Epetra_Time: public Epetra_Object {
    
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

 private:

  double StartTime_;
  const Epetra_Comm * Comm_;
  
};

#endif /* _EPETRA_TIME_H_ */
