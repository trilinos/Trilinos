#ifndef _Trilinos_ESI_Broker_h_
#define _Trilinos_ESI_Broker_h_

class Epetra_Comm;
class Epetra_CrsGraph;

#ifdef EPETRA_SERIAL
#define MPI_Comm int
#define MPI_SUCCESS 0
#endif
#include "ESI_Broker.h"

//forward declarations for esi interfaces.

class esi::Object;
template<class T> class esi::IndexSpace;
template<class Scalar,class Ordinal> class esi::Vector;

#include "Epetra_Array.h"

/** Petra_ESI implementation of ESI_Broker, a broker and manager for
   the pieces of a linear system, which appear as (mostly) ESI interfaces.

  This class is intended to be given to the FEI class FEI_Implementation
  as a constructor argument. The FEI implementation will then request various
  interfaces from us, and use them for passing data into as it assembles a
  linear system from the finite-element data it receives from the application.

  Important usage notes for Trilinos_ESI_Broker:
    - The parameters function is used to select the preconditioner and solver.
    Parameters should be space-separated key-value string pairs.
    - Interfaces for maps and vectors are not available until after
    'setGlobalOffsets' has been called. Interfaces for matrices are not
    available until after setMatrixStructure has been called.
    - The solver interface is not available until after a matrix has been
    constructed.
*/

class Trilinos_ESI_Broker : public virtual ESI_Broker {
 public:
  /** Default constructor. Intended for use by a component framework.
      (Components need to have argument-free constructors.) If this class
      is constructed this way, then one of the first methods called after
      construction should be Trilinos_ESI_Broker::setMPIComm.
  */
  Trilinos_ESI_Broker();

  /** Constructor which takes an MPI Communicator as an argument. This is the
      constructor that will be used by SIERRA when using this class with
      SNL's FEI implementation.
  */
  Trilinos_ESI_Broker(MPI_Comm comm);

  /** Destructor. */
  ~Trilinos_ESI_Broker();

  /** This method need only be called if the default constructor is used to
      construct this class. In that case, this should be the first method
      called after construction.
  */
  int setMPIComm(MPI_Comm comm);

   bool hasESIinterfaces() { return(true); };

   int getLibraryName(char*& libName);

   int parameters(int numParams, char** paramStrings);

   int clone(ESI_Broker*& esi_lsmgr);

   int setGlobalOffsets(int len, int* nodeOffsets,
                        int* eqnOffsets, int* blkEqnOffsets);

   int setMultCREqns(int numCRs, 
		     const int* numNodesPerCR,
		     const int* const* nodeNumbers,
		     const int* const* eqnNumbers,
		     const int* fieldIDs,
		     const int* multiplierEqnNumbers);

   int setMatrixStructure(int** ptColIndices,
			  int* ptRowLengths,
			  int** blkColIndices,
			  int* blkRowLengths,
			  int* ptRowsPerBlkRow);

   int getInterfaceInstance(const char* instanceName,
			    const char* interfaceName,
			    void*& objectPtr);

   int setInterfaceInstance(const char* instanceName,
			    const char* interfaceName,
			    void* objectPtr);

 private:
   bool stringsMatch(const char* str1, const char* str2)
     {
       int len1 = strlen(str1);   int len2 = strlen(str2);

       if (len1 != len2) return(false);

       if (strncmp(str1, str2, len1) == 0) return(true);
       else return(false);
     };

   int openDebugOutput(const char* path, const char* fileName);

   bool getParam(const char* flag, int numParams,
                 char** paramStrings, char* result);

   int getIndexSpaceInstance(const char* mapName, void*& mapPtr);
   int getVectorInstance(const char* mapName,
			 const char* vecName, void*& vecPtr);
   int getMatrixInstance(const char* mapName,
			 const char* matName, void*& matPtr);
   int getSolverInstance(const char* slvName, void*& slvPtr);
   int getPreconditionerInstance(const char* pcName, void*& pcPtr);

   char* libName_;

   MPI_Comm comm_;
   bool haveMPIComm_;
   Epetra_Comm* petra_comm_;
   Epetra_CrsGraph* petra_graph_;

   //we'll manage lists of maps, matrices, and vectors (both rhs and soln)
   Epetra_Array<esi::IndexSpace<int>*>       ispaces_;
   Epetra_Array<char*>                        ispaceNames_;
   Epetra_Array<esi::Object*>                  matrices_;
   Epetra_Array<char*>                        matrixNames_;
   Epetra_Array<esi::Vector<double,int>*>      vecs_;
   Epetra_Array<char*>                        vecNames_;

   bool setGlobalOffsetsCalled_;
   bool setMatrixStructureCalled_;

   Epetra_Array<esi::Object*> solvers_;
   Epetra_Array<char*> solverNames_;

   int localProc_, numProcs_;

   int globalSize_;
   int localSize_;
   int localOffset_;

   Epetra_Array<int> rowLengths_;

   Epetra_Array<char*> params_;

   int outputLevel_;

   bool debugOutput_;
   FILE* debugFile_;
};

#endif

