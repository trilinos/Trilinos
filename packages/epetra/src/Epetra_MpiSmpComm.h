
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_MPISMPCOMM_H
#define EPETRA_MPISMPCOMM_H
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiDistributor.h"
#include <mpi.h>
class Epetra_Distributor;
#include "Epetra_MpiSmpCommData.h"

//! Epetra_MpiSmpComm:  The Epetra MPI Shared Memory Parallel Communication Class.
/*! The Epetra_MpiSmpComm class is an implementation of Epetra_Comm that encapsulates the general
  information and services needed for other Epetra classes to run on a parallel computer using MPI and
  shared memory threads.
  \warning This is an experimental class that marginally supported nested share memory parallelism within
  MPI processes.
*/

class Epetra_MpiSmpComm: public Epetra_Object, public virtual Epetra_Comm {
    
  public:

  //@{ \name Constructor/Destructor Methods
  //! Epetra_MpiSmpComm MPI Constructor.
  /*! Creates a Epetra_MpiSmpComm instance for use with MPI.  If no specialized
    MPI communicator is needed, this constuctor can be called with the
    argument MPI_COMM_WORLD.  
  */
  Epetra_MpiSmpComm(MPI_Comm comm);


  //! Epetra_MpiSmpComm Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_MpiSmpComm instance.
  */
  Epetra_MpiSmpComm(const Epetra_MpiSmpComm& Comm);

	//! Clone method.
	Epetra_Comm * Clone() const {
		return(dynamic_cast<Epetra_Comm *>(new Epetra_MpiSmpComm(*this)));
	};

  //! Epetra_MpiSmpComm Destructor.
  /*! Completely deletes a Epetra_MpiSmpComm object.  
    \warning Note:  All objects that depend
    on a Epetra_MpiSmpComm instance should be destroyed prior to calling this
    function.
  */
  virtual ~Epetra_MpiSmpComm();
  //@}

  //@{ \name Barrier Methods
  //! Epetra_MpiSmpComm Barrier function.
  /*!Causes each processor in the communicator to wait until all processors
    have arrived.
  */
  void Barrier() const;
  //@}

  //@{ \name Broadcast Methods
  //! Epetra_MpiSmpComm Broadcast function.
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

  //! Epetra_MpiSmpComm Broadcast function.
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
  //@}

  //@{ \name Gather Methods
  //! Epetra_MpiSmpComm All Gather function.
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

  //! Epetra_MpiSmpComm All Gather function.
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
  //@}

  //@{ \name Sum Methods
  //! Epetra_MpiSmpComm Global Sum function.
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

  //! Epetra_MpiSmpComm Global Sum function.
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
  //@}

  //@{ \name Max/Min Methods
  //! Epetra_MpiSmpComm Global Max function.
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


  //! Epetra_MpiSmpComm Global Max function.
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

  //! Epetra_MpiSmpComm Global Min function.
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


  //! Epetra_MpiSmpComm Global Min function.
  /*!Take list of input values from all processors in the communicator, computes the max and returns the
    max to all processors.
    \param PartialMins In
           On entry, contains the list of values, usually partial mins computed locally;
					 using these Partial Mins, the min across all processors will be computed.
    \param GlobalMins Out
           On exit, contains the list of mins computed across all processors.
    \param Count In
           On entry, contains the length of the list of values.
  */
  int MinAll(int * PartialMins, int * GlobalMins, int Count) const;
  //@}

  //@{ \name Parallel Prefix Methods
  //! Epetra_MpiSmpComm Scan Sum function.
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


  //! Epetra_MpiSmpComm Scan Sum function.
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
  //@}

  //@{ \name Attribute Accessor Methods
  
  //! Extract MPI Communicator from a Epetra_MpiSmpComm object.
  MPI_Comm Comm() const {return(MpiSmpCommData_->Comm_);};

  //! Return my process ID. 
  /*! In MPI mode returns the rank of the calling process.  In serial mode
    returns 0.
  */
  int MyPID() const {return(MpiSmpCommData_->rank_);};
  
  //! Returns total number of processes. 
  /*! In MPI mode returns the size of the MPI communicator.  In serial mode
    returns 1.
  */
  int NumProc() const {return(MpiSmpCommData_->size_);};
  //@}

  //@{ \name Gather/Scatter and Directory Constructors
  //! Create a distributor object.
  Epetra_Distributor * CreateDistributor() const;
  //! Create a directory object for the given Epetra_BlockMap.
  Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const;
  //@}

  //@{ \name MPI-specific Methods
  //! Acquire an MPI tag from the Epetra range of 24050-24099, increment tag. 
  int GetMpiTag() const {int tag = MpiSmpCommData_->curTag_++; if (tag>MpiSmpCommData_->maxTag_) tag = MpiSmpCommData_->minTag_; return(tag);};
  //! Acquire an MPI tag from the Epetra range of 24050-24099, increment tag. 
  MPI_Comm GetMpiComm() const {return(MpiSmpCommData_->Comm_);};
  //@}
  //@{ \name Experimental SMP cluster methods (not rigorously implemented)
  //! Epetra_MpiSmpComm Node Barrier function.
  /*! A no-op for a serial communicator.  For MPI, it 
    causes each process on a given node in the communicator to wait until all processes on that node
    have arrived.

    This function can be used to select a subset of MPI processes that are associated with a group of threaded
    processes and synchronize only with this subset.
  */
  void NodeBarrier() const;

  //! Return my thread ID. 
  /*! If SetMyThreadID was called to set a thread value, this function returns the thread ID of the calling process.  
      Otherwise returns 0.
  */
  int MyThreadID() const {return(MpiSmpCommData_->ThreadID_);};

  //! Return my node ID. 
  /*! If SetMyNodeD was called to set a node value, this function returns the thread ID of the calling process.  
      Otherwise returns the same value as MyPID().
  */
  int MyNodeID() const {return(MpiSmpCommData_->NodeID_);};
  
  //! Set number of threads on this node. 
  /*! Sets the number of threads on the node that owns the calling process.  By default the number of threads is 1.
  */
  int SetNumThreads(int NumThreads) {MpiSmpCommData_->NumThreads_ = NumThreads; return(0);};
  
  //! Get number of threads on this node. 
  /*! Sets the number of threads on the node that owns the calling process.  By default the number of threads is 1.
  */
  int NumThreads() const {return(MpiSmpCommData_->NumThreads_);};
   
  //! Set my thread ID. 
  /*! Sets the thread ID for the calling process.  Can be used to facilitate threaded programming across an MPI
      application by allowing multiple MPI processes to be considered threads of a virtual shared memory process.
      Threads and nodes should be used together.  By default the thread ID is zero.
  */
  int SetMyThreadID(int ThreadID) {MpiSmpCommData_->ThreadID_ = ThreadID; return(0);};
  
  //! Set my node ID. 
  /*! Sets the node ID for the calling process.  Can be used to facilitate threaded programming across an MPI
      application by associating several MPI processes with a single node. By default, each MPI process is
      associated with a single node with the same ID.
  */
  int SetMyNodeID(int NodeID) {MpiSmpCommData_->NodeID_ = NodeID; return(0);};
  //@}

  //@{ \name Print object to an output stream
  //! Print method that implements Epetra_Object virtual Print method
  inline void Print(ostream & os) const {
  os << "::Processor "<< MyPID()<<" of " << NumProc() << " total processors" << endl; 
  os << "::Thread "<< MyThreadID()<<" of " << NumThreads() << " on node " << MyNodeID() << endl; 
  return;}

  //! Print method that implements Epetra_Comm virtual PrintInfo method
  void PrintInfo(ostream & os) const {Epetra_MpiSerialComm::Print(os);return;};

  //@}

	//! Assignment Operator
	Epetra_MpiSmpComm & operator=(const Epetra_MpiSmpComm & Comm);
  
 private:

	void CleanupData();
	Epetra_MpiSmpCommData * MpiSmpCommData_;

};
#endif /* EPETRA_MPISMPCOMM_H */
