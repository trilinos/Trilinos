#ifndef _Epetra_ESI_IndexSpace_h_
#define _Epetra_ESI_IndexSpace_h_

#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"

#include "Epetra_ESI_CHK_ERR.h"
#include "Epetra_ESI_Object.h"

namespace epetra_esi {

/** Petra's ESI IndexSpace implementation.
This class implements these ESI interfaces:
<ul>
<li>   esi::Object
<li>   esi::IndexSpace
</ul>
Note that although this class is templated, it may only be instantiated on
the type int.
*/

template<class Ordinal>
class IndexSpace : public virtual esi::IndexSpace<Ordinal>,
                      public virtual epetra_esi::Object,
                      public virtual Epetra_Map
{
 public:
  /** Constructor. */
  IndexSpace(Ordinal globalSize, Ordinal localSize, Ordinal indexBase,
                const Epetra_Comm& comm);

  /** Destructor. */
  virtual ~IndexSpace();


  //
  //esi::IndexSpace functions.
  //

  /** Get the global size of the basis that this map describes.
    @param globalSize Output. Global (across all processors) basis size.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalSize(Ordinal& globalSize);


  /** Get the local size of the basis that this map describes.
    @param localSize Output. Local (to 'this' processor) basis size.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getLocalSize(Ordinal& localSize)
    { localSize = NumMyPoints(); return(0); }


  /** Get the number of 'colors' represented in this map. epetra_esi::IndexSpace doesn't
    currently do coloring, so the output value will always be zero.
    @param colorSetSize Output. Always 0.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalColorSetSize(Ordinal& colorSetSize)
    { colorSetSize = 0; return(0); }


  /** Get the local 'colors' represented in this map. epetra_esi::IndexSpace doesn't
    currently do coloring, so this function is a no-op.
    @param localColors Output. Not referenced by this implementation.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getLocalColors(Ordinal*)
    { return(0); }


  /** Get the locally owned identifiers. This actually obtains the locally-
     owned set of global identifiers. i.e., global equation numbers.
   @param localIdentifiers Output. List of length getLocalSize(), allocated by
               the caller.
   @return error-code, 0 if successful.
  */
  virtual esi::ErrorCode getLocalIdentifiers(Ordinal* localIdentifiers)
    { return( MyGlobalElements(localIdentifiers) ); }


  /** Get the global number of partitions of the basis, represented by this
     map. This is the number of processors for this implementation.
    @param numPartitions Output. num-processors.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalPartitionSetSize(Ordinal& numPartitions)
    { numPartitions = Comm().NumProc(); return(0); }


  /** Get the rank of 'this' partition. This is the local processor rank, or
    the local mpi rank.
    @param localRank On exit, local mpi rank.
    @return error-code, 0 if successful.
  */
  virtual esi::ErrorCode getLocalPartitionRank(Ordinal& localRank)
    { localRank = Comm().MyPID(); return(0); }


  /** Offset into the global basis, of the first locally owned equation.
    @param globalOffset Output.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode getLocalPartitionOffset(Ordinal& globalOffset)
    { globalOffset = MinMyGID(); return(0); }


  /** Get a list containing the number-of-local-identifiers on all processors.
   @param partitionSizes Output. Caller-allocated list of length
       getGlobalPartitionSetSize(). On exit, the i-th entry in this list is
       the number of equations held on processor i.
   @param error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalPartitionSizes(Ordinal* partitionSizes)
    { Ordinal mySize = NumMyPoints();
      return( Comm().GatherAll(&mySize, partitionSizes, 1) ); }


  /** Get a list containing the offsets of all processors' first locally 
     owned equations.
    @param partitionOffsets Output. Caller-allocated list of length 
      getGlobalPartitionSetSize().
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalPartitionOffsets(Ordinal* partitionOffsets)
    { Ordinal myOffset = MinMyGID();
      return( Comm().GatherAll(&myOffset, partitionOffsets, 1) ); }

 private:
#ifdef EPETRA_MPI
  MPI_Comm comm_;
#endif
};

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_IndexSpace.cpp"
#endif

#endif

