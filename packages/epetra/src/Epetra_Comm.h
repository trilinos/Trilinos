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

#ifndef EPETRA_COMM_H
#define EPETRA_COMM_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"

class Epetra_Distributor;
class Epetra_Directory;
class Epetra_BlockMap;

//! Epetra_Comm:  The Epetra Communication Abstract Base Class.
/*! The Epetra_Comm class is an interface that encapsulates the general
    information and services needed for other Epetra classes to run on a
    parallel computer. An Epetra_Comm object is required for building all
    Epetra Map objects, which in turn are required for all other Epetra
    classes.

    Epetra_Comm has default implementations, via Epetra_SerialComm and
    Epetra_MpiComm, for both serial execution and MPI distributed memory
    execution.  It is meant to insulate the user from the specifics of
    communication that are not required for normal manipulation of linear
    algebra objects.  Most Epetra_Comm interfaces are similar to MPI
    interfaces, except that the type of data is not required as an argument
    since C++ can bind to the appropriate interface based on argument typing.

    Any implementation of the Epetra_Comm interface is also responsible for
    generating an Epetra_Distributor and Epetra_Directory object.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_Comm {

  public:
    //! @name Constructor / Destructor
  //@{
	//! Epetra_Comm clone constructor.
	/*! The clone function will return a new heap-allocated Comm instance.
            It is the responsibility of the caller to ensure that this new instance
            is properly destroyed.
        */
	virtual Epetra_Comm * Clone() const = 0;
  //! Epetra_Comm Destructor.
  virtual ~Epetra_Comm() {};
  //@}

  //! @name Barrier Methods
  //@{
  //! Epetra_Comm Barrier function.
  /*! Each processor must wait at the point the barrier is called until all processors have arrived.
  */
  virtual void Barrier() const = 0;
  //@}

  //! @name Broadcast Methods
  //@{
  //! Epetra_Comm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
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
  /*! Take list of input values from the root processor and sends to all other processors.
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

  //! Epetra_Comm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param MyVals InOut
           On entry, the root processor contains the list of values.  On exit,
					 all processors will have the same list of values.  Note that values must be
					 allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  virtual int Broadcast(long * MyVals, int Count, int Root) const = 0;

  //! Epetra_Comm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param MyVals InOut
           On entry, the root processor contains the list of values.  On exit,
					 all processors will have the same list of values.  Note that values must be
					 allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  virtual int Broadcast(long long * MyVals, int Count, int Root) const = 0;

  //! Epetra_Comm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param MyVals InOut
           On entry, the root processor contains the list of values.  On exit,
					 all processors will have the same list of values.  Note that values must be
					 allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  virtual int Broadcast(char * MyVals, int Count, int Root) const = 0;

  //@}

  //! @name Gather Methods
  //@{
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
  /*! Take list of input values from all processors in the communicator and creates an ordered contiguous list of
      those values on each processor.
    \param MyVals In
           On entry, contains the list of values to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must be of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  virtual int GatherAll(int * MyVals, int * AllVals, int Count) const = 0;

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

  virtual int GatherAll(long * MyVals, long * AllVals, int Count) const = 0;

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

  virtual int GatherAll(long long * MyVals, long long * AllVals, int Count) const = 0;
  //@}

  //! @name Sum Methods
  //@{
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
  virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const = 0;

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
  virtual int SumAll(long long * PartialSums, long long * GlobalSums, int Count) const = 0;
  //@}

  //! @name Max/Min Methods
  //@{
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
  virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const = 0;

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
  virtual int MaxAll(long long * PartialMaxs, long long * GlobalMaxs, int Count) const = 0;

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
  virtual int MinAll(long * PartialMins, long * GlobalMins, int Count) const = 0;

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
  virtual int MinAll(long long * PartialMins, long long * GlobalMins, int Count) const = 0;
  //@}

  //! @name Parallel Prefix Methods
  //@{
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
  virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const = 0;

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
  virtual int ScanSum(long long * MyVals, long long * ScanSums, int Count) const = 0;
  //@}

  //! @name Attribute Accessor Methods
  //@{

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

  //! @name Gather/Scatter and Directory Constructors
  //@{
  //! Create a distributor object.
  virtual Epetra_Distributor * CreateDistributor() const = 0;
  //! Create a directory object for the given Epetra_BlockMap.
// CreateDirectory is defined in Winbase.h as a macro!
#ifdef CreateDirectory
#undef CreateDirectory
#endif
  virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const = 0;
  //@}

  //! @name I/O methods
  //@{
  //! Print object to an output stream
  virtual void PrintInfo(std::ostream & os) const = 0;
  //@}
};
#endif /* EPETRA_COMM_H */
