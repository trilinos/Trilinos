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

#ifndef EPETRA_MPICOMM_H
#define EPETRA_MPICOMM_H
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiDistributor.h"
class Epetra_Distributor;
#include "Epetra_BasicDirectory.h"
class Epetra_Directory;
class Epetra_BlockMap;
#include <mpi.h>
#include "Epetra_MpiCommData.h"
#ifdef EPETRA_HAVE_OMP
#include <omp.h>
#endif

//! Epetra_MpiComm:  The Epetra MPI Communication Class.
/*! The Epetra_MpiComm class is an implementation of Epetra_Comm that encapsulates the general
  information and services needed for other Epetra classes to run on a parallel computer using MPI.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_MpiComm: public Epetra_Object, public virtual Epetra_Comm {
    
  public:

    //! @name Constructor/Destructor Methods
  //@{ 
  //! Epetra_MpiComm MPI Constructor.
  /*! Creates a Epetra_MpiComm instance for use with MPI.  If no specialized
    MPI communicator is needed, this constuctor can be called with the
    argument MPI_COMM_WORLD.  
  */
  Epetra_MpiComm(MPI_Comm comm);


  //! Epetra_MpiComm Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_MpiComm instance.
  */
  Epetra_MpiComm(const Epetra_MpiComm & Comm);

  //! Clone method.
  Epetra_Comm * Clone() const
    {
      return(new Epetra_MpiComm(*this));
    };

  //! Epetra_MpiComm Destructor.
  /*! Completely deletes a Epetra_MpiComm object.  
    \warning Note:  All objects that depend
    on a Epetra_MpiComm instance should be destroyed prior to calling this
    function.
  */
  virtual ~Epetra_MpiComm();
  //@}

  //! @name Barrier Methods
  //@{ 
  //! Epetra_MpiComm Barrier function.
  /*!Causes each processor in the communicator to wait until all processors
    have arrived.
  */
  void Barrier() const;
  //@}

  //! @name Broadcast Methods
  //@{ 
  //! Epetra_MpiComm Broadcast function.
  /*!Takes list of input values from the root processor and sends to all other processors.
    \param Values InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  int Broadcast(double * MyVals, int Count, int Root) const;

  //! Epetra_MpiComm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param Values InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  int Broadcast(int * MyVals, int Count, int Root) const;

  //! Epetra_MpiComm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param Values InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  int Broadcast(long * MyVals, int Count, int Root) const;

  //! Epetra_MpiComm Broadcast function.
  /*! Take list of input values from the root processor and sends to all other processors.
    \param Values InOut
           On entry, the root processor contains the list of values.  On exit,
	   all processors will have the same list of values.  Note that values must be
	   allocated on all processor before the broadcast.
    \param Count In
           On entry, contains the length of the list of Values.
    \param Root In
           On entry, contains the processor from which all processors will receive a copy of Values.
  */

  int Broadcast(char * MyVals, int Count, int Root) const;
  //@}

  //! @name Gather Methods
  //@{ 
  //! Epetra_MpiComm All Gather function.
  /*! Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values, to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must by of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  int GatherAll(double * MyVals, double * AllVals, int Count) const;

  //! Epetra_MpiComm All Gather function.
  /*!Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values, to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must by of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  int GatherAll(int * MyVals, int * AllVals, int Count) const;

  //! Epetra_MpiComm All Gather function.
  /*!Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values, to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must by of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  int GatherAll(long * MyVals, long * AllVals, int Count) const;

  //! Epetra_MpiComm All Gather function.
  /*!Take list of input values from all processors in the communicator and creates an ordered contiguous list of
    those values on each processor.
    \param MyVals In
           On entry, contains the list of values, to be sent to all processors.
    \param AllVals Out
           On exit, contains the list of values from all processors. Must by of size NumProc*Count.
    \param Count In
           On entry, contains the length of the list of MyVals.
  */

  int GatherAll(long long * MyVals, long long * AllVals, int Count) const;
  //@}

  //! @name Sum Methods
  //@{ 
  //! Epetra_MpiComm Global Sum function.
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

  int SumAll(double * PartialSums, double * GlobalSums, int Count) const;

  //! Epetra_MpiComm Global Sum function.
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
  int SumAll(int * PartialSums, int * GlobalSums, int Count) const;

  //! Epetra_MpiComm Global Sum function.
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
  int SumAll(long * PartialSums, long * GlobalSums, int Count) const;

  //! Epetra_MpiComm Global Sum function.
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
  int SumAll(long long * PartialSums, long long * GlobalSums, int Count) const;
  //@}

  //! @name Max/Min Methods
  //@{ 
  //! Epetra_MpiComm Global Max function.
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
  int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const;

  //! Epetra_MpiComm Global Max function.
  /*!Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
					 using these Partial Maxs, the max across all processors will be computed.
    \param GlobalMaxs Out
           On exit, contains the list of maxs computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const;

  //! Epetra_MpiComm Global Max function.
  /*!Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
					 using these Partial Maxs, the max across all processors will be computed.
    \param GlobalMaxs Out
           On exit, contains the list of maxs computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const;

  //! Epetra_MpiComm Global Max function.
  /*!Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
					 using these Partial Maxs, the max across all processors will be computed.
    \param GlobalMaxs Out
           On exit, contains the list of maxs computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MaxAll(long long * PartialMaxs, long long * GlobalMaxs, int Count) const;

  //! Epetra_MpiComm Global Min function.
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
  int MinAll(double * PartialMins, double * GlobalMins, int Count) const;

  //! Epetra_MpiComm Global Min function.
  /*!Take list of input values from all processors in the communicator, computes the min and returns the
    min to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
					 using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MinAll(int * PartialMins, int * GlobalMins, int Count) const;

  //! Epetra_MpiComm Global Min function.
  /*!Take list of input values from all processors in the communicator, computes the min and returns the
    min to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
					 using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MinAll(long * PartialMins, long * GlobalMins, int Count) const;

  //! Epetra_MpiComm Global Min function.
  /*!Take list of input values from all processors in the communicator, computes the min and returns the
    min to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
					 using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MinAll(long long * PartialMins, long long * GlobalMins, int Count) const;
  //@}

  //! @name Parallel Prefix Methods
  //@{ 
  //! Epetra_MpiComm Scan Sum function.
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
  int ScanSum(double * MyVals, double * ScanSums, int Count) const;

  //! Epetra_MpiComm Scan Sum function.
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
  int ScanSum(int * MyVals, int * ScanSums, int Count) const;

  //! Epetra_MpiComm Scan Sum function.
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
  int ScanSum(long * MyVals, long * ScanSums, int Count) const;

  //! Epetra_MpiComm Scan Sum function.
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
  int ScanSum(long long * MyVals, long long * ScanSums, int Count) const;
  //@}

  //! @name Attribute Accessor Methods
  //@{ 
  
  //! Extract MPI Communicator from a Epetra_MpiComm object.
  MPI_Comm Comm() const {return(MpiCommData_->Comm_);};

  //! Return my process ID. 
  /*! In MPI mode returns the rank of the calling process.  In serial mode
    returns 0.
  */
  int MyPID() const {return(MpiCommData_->rank_);};
  
  //! Returns total number of processes. 
  /*! In MPI mode returns the size of the MPI communicator.  In serial mode
    returns 1.
  */
  int NumProc() const {return(MpiCommData_->size_);};
  //@}

  //! @name Gather/Scatter and Directory Constructors
  //@{ 
  //! Create a distributor object.
  Epetra_Distributor * CreateDistributor() const;
  //! Create a directory object for the given Epetra_BlockMap.
  Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const;
  //@}

  //! @name MPI-specific Methods
  //@{ 
  //! Acquire an MPI tag from the Epetra range of 24050-24099, increment tag. 
  int GetMpiTag() const {int tag = MpiCommData_->curTag_++; if (tag > MpiCommData_->maxTag_) tag = MpiCommData_->minTag_; return(tag);};

  //! Get the MPI Communicator (identical to Comm() method; used when we know we are MPI. 
  MPI_Comm GetMpiComm() const {return(MpiCommData_->Comm_);};
  //@}
  //! @name Print object to an output stream
  //@{ 
  //! Print method that implements Epetra_Object virtual Print method
  inline void Print(ostream & os) const {
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel 
{
  int numThreads = omp_get_num_threads();
  int threadNum = omp_get_thread_num();
#pragma omp critical
                os << "\n   Thread/TotalThreads " << threadNum << "/" << numThreads << " on Processor/TotalProcs "
		<<  MyPID() << "/" << NumProc();
}
#else
                os << "::Processor "<< MyPID()<<" of " << NumProc() << " total processors.";
#endif
  return;}
  //! Print method that implements Epetra_Comm virtual PrintInfo method
  void PrintInfo(ostream & os) const {Epetra_MpiComm::Print(os);return;};
  //@}

  //! @name Expert Users and Developers Only
  //@{ 

	//! Returns the reference count of MpiCommData.
	/*! (Intended for testing purposes.) */
	int ReferenceCount() const {return(MpiCommData_->ReferenceCount());};

	//! Returns a pointer to the MpiCommData instance this MpiComm uses. 
	/*! (Intended for developer use only for testing purposes.) */
	const Epetra_MpiCommData * DataPtr() const {return(MpiCommData_);};

  //@}

	//! Assignment Operator
	Epetra_MpiComm & operator=(const Epetra_MpiComm & Comm);
  
 private:
  
  int CheckInput(double * ptr, int count) const {if ((ptr==0) && (count>0)) return(-1); return(0);};
  int CheckInput(int * ptr, int count) const {if ((ptr==0) && (count>0)) return(-1); return(0);};
  int CheckInput(long * ptr, int count) const {if ((ptr==0) && (count>0)) return(-1); return(0);};
  int CheckInput(char * ptr, int count) const {if ((ptr==0) && (count>0)) return(-1); return(0);};
  int CheckInput(long long * ptr, int count) const {if ((ptr==0) && (count>0)) return(-1); return(0);};

  void CleanupData();
  Epetra_MpiCommData * MpiCommData_;

};
#endif /* EPETRA_MPICOMM_H */
