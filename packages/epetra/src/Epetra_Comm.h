
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

#ifndef _EPETRA_COMM_H_
#define _EPETRA_COMM_H_

#include "Epetra_Object.h"
class Epetra_Distributor;
class Epetra_Directory;
class Epetra_BlockMap;

//! Epetra_Comm:  The Epetra Communication Abstract Base Class.
/*! The Epetra_Comm class is an interface that encapsulates the general
  information and services needed for other Epetra classes to run on a parallel computer.
  An Epetra_Comm object is required for building all Epetra Map objects, which
  in turn are required for all other Epetra classes.
  
  Epetra_Comm has default implementations, via Epetra_SerialComm and 
  Epetra_MpiComm, for both serial execution and MPI
  distributed memory execution.  It is meant to insulate the user from
  the specifics of communication that are not required for normal
  manipulation of linear algebra objects.  Most Epetra_Comm interfaces are similar to MPI
  interfaces, except that the type of data is not required as an argument since C++ can bind
  to the appropriate interface based on argument typing.

  Any implementation of the Epetra_Comm interface is also responsible for generating an
  Epetra_Distributor and Epetra_Directory object.
*/

class Epetra_Comm {
    
  public:
  //@{ \name Destructor
  //! Epetra_Comm Destructor.
  virtual ~Epetra_Comm(){};
  //@}

  //@{ \name Barrier Methods
  //! Epetra_Comm Barrier function.
  /*! Each processor must wait at the point the barrier is called until all processors have arrived.
  */
  virtual void Barrier() const = 0;
  //@}

  //@{ \name Broadcast Methods
  //! Epetra_Comm Broadcast function.
  /*!Take list of input values from the root processor and sends to all other processors.
    \param MyVals InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  virtual int Broadcast(double * MyVals, int Count, int Root) const = 0;

  //! Epetra_Comm Broadcast function.
  /*!Take list of input values from the root processor and sends to all other processors.
    \param MyVals InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  virtual int Broadcast(int * MyVals, int Count, int Root) const = 0;
  //@}

  //@{ \name Gather Methods
  //! Epetra_Comm All Gather function.
  /*! Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must be of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  virtual int GatherAll(double * MyVals, double * AllVals, int Count) const = 0;

  //! Epetra_Comm All Gather function.
  /*!Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must be of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  virtual int GatherAll(int * MyVals, int * AllVals, int Count) const = 0;
  //@}

  //@{ \name Sum Methods
  //! Epetra_Comm Global Sum function.
  /*!Take list of input values from all processors in the communicator, computes the sum and returns the
    sum to all processors.
    \param PartialSums In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all processors.
    \param GlobalSums Out
           On exit, contains the list of values summed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */

  virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const = 0;

  //! Epetra_Comm Global Sum function.
  /*! Take list of input values from all processors in the communicator, computes the sum and returns the
    sum to all processors.
    \param PartialSums In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all processors.
    \param GlobalSums Out
           On exit, contains the list of values summed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const = 0;
  //@}

  //@{ \name Max/Min Methods
  //! Epetra_Comm Global Max function.
  /*! Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
	   using these Partial Maxs, the max across all processors will be computed.
    \param GlobalMaxs Out
           On exit, contains the list of maxs computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const = 0;


  //! Epetra_Comm Global Max function.
  /*! Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
	   using these Partial Maxs, the max across all processors will be computed.
    \param GlobalMaxs Out
           On exit, contains the list of maxs computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const = 0;

  //! Epetra_Comm Global Min function.
  /*! Take list of input values from all processors in the communicator, computes the min and returns the
    min to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
	   using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const = 0;


  //! Epetra_Comm Global Min function.
  /*! Take list of input values from all processors in the communicator, computes the min and returns the
    min to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
	   using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const = 0;
  //@}

  //@{ \name Parallel Prefix Methods
  //! Epetra_Comm Scan Sum function.
  /*! Take list of input values from all processors in the communicator, computes the scan sum and returns it 
    to all processors such that processor i contains the sum of values from processor 0 up to and including
    processor i.
    \param MyValss In
           On entry, contains the list of values to be summed across all processors.
    \param ScanSums Out
           On exit, contains the list of values summed across processors 0 through i.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int ScanSum(double * MyVals, double * ScanSums, int Count) const = 0;


  //! Epetra_Comm Scan Sum function.
  /*! Take list of input values from all processors in the communicator, computes the scan sum and returns it 
    to all processors such that processor i contains the sum of values from processor 0 up to and including
    processor i.
    \param MyVals In
           On entry, contains the list of values to be summed across all processors.
    \param ScanSums Out
           On exit, contains the list of values summed across processors 0 through i.
    \param Count In
           On entry, contains the length of the list of values.
  */
  virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const = 0;
  //@}

  //@{ \name Attribute Accessor Methods
    
  //! Return my process ID. 
  /*! In MPI mode returns the rank of the calling process.  In serial mode
    returns 0.
  */
  virtual int MyPID() const = 0;
  
  //! Returns total number of processes. 
  /*! In MPI mode returns the size of the MPI communicator.  In serial mode
    returns 1.
  */
  virtual int NumProc() const = 0;
  //@}

  //@{ \name Gather/Scatter and Directory Constructors
  //! Create a distributor object.
  virtual Epetra_Distributor * CreateDistributor() const = 0;
  //! Create a directory object for the given Epetra_BlockMap.
  virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const = 0;
  //@}

  //@{ \name I/O methods
  //! Print object to an output stream
  virtual void PrintInfo(ostream & os) const = 0;
  //@}
};
#endif /* _EPETRA_COMM_H_ */
