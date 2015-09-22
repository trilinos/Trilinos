#ifndef _FEData_h_
#define _FEData_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>
#include <fei_fstream.hpp>

#include <fei_mpi.h>

#define dbgOut() if (debugOutputLevel_ > 0) *dbgOStreamPtr_

#include <fei_defs.h>
#include <fei_FiniteElementData.hpp>

#include <snl_fei_Utils.hpp>

/** Simply a test harness to use for checking the FEI layer's usage of the
    FiniteElementData interface.
*/

class FEData : public virtual FiniteElementData {
 public:
  /** Constructor. */
  FEData(MPI_Comm comm)
    :
    comm_(comm), numProcs_(1), localProc_(0),
    debugOutputLevel_(0),
    dbgPath_(NULL),
    dbgOStreamPtr_(NULL),
    dbgFileOpened_(false)
    {
#ifndef FEI_SER
      if (MPI_Comm_rank(comm_, &localProc_) != MPI_SUCCESS) MPI_Abort(comm_,-1);
      if (MPI_Comm_size(comm_, &numProcs_) != MPI_SUCCESS) MPI_Abort(comm_,-1);
#endif
      setDebugLog(0, ".");
    }

  virtual ~FEData()
    {
      if (dbgFileOpened_ == true) { dbgFStreamPtr_->close(); }

      delete [] dbgPath_;
      delete dbgOStreamPtr_;
    }

  /** For setting argc/argv style parameters.
     @param numParams Number of strings in the params argument
     @param params A list of strings which will usually contain space-separated
        key-value pairs. Example: "debugOutput /usr/users/me/work_dir"
  */
  int parameters(int numParams, char** params);

  /** Supply the FiniteElementData implementation with an object (created and
      owned by the caller) that can be used to obtain various information about
      problem layout, shared finite-element nodes, etc.
      For details, see the documentation for the Lookup interface.
      @param lookup Input. Reference to an implementation of the
      Lookup interface
  */
  int setLookup(Lookup& lookup)
    {
      dbgOut() << "setLookup" << FEI_ENDL;
      return(0);
    }


  /** For describing the general structure of the finite-element problem that
      is to be assembled.
  */
  int describeStructure(int numElemBlocks,
			const int* numElemsPerBlock,
			const int* numNodesPerElem,
			const int* elemMatrixSizePerBlock,
			int totalNumNodes,
			int numSharedNodes,
			int numMultCRs)
    {
      dbgOut() << "describeStructure" << FEI_ENDL
	<< "   numElemBlocks: " << numElemBlocks << FEI_ENDL;
      for(int i=0; i<numElemBlocks; i++) {
	dbgOut() << "   elem-block " << i << ": " << FEI_ENDL
	  << "      number of elements: " << numElemsPerBlock[i] << FEI_ENDL
	  << "      nodes per element:  " << numNodesPerElem[i] << FEI_ENDL;
	dbgOut() << "         elemMatrixSizePerBlock: "
	    << elemMatrixSizePerBlock[i] << FEI_ENDL;
      }
      return(0);
    }

   /** For passing element-connectivity arrays.
   */
   int setConnectivity(int elemBlockID,
		       int elemID,
		       int numNodes,
		       const int* nodeNumbers,
		       const int* numDofPerNode,
           const int* dof_ids)
     {
       dbgOut() << "setConnectivity" << FEI_ENDL
	 << "   elemBlockID: " << elemBlockID << ", elemID: " << elemID 
		<< ", numNodes: " << numNodes << FEI_ENDL
	 << "   nodeNumbers: ";
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << nodeNumbers[i] << " ";
       }
       dbgOut() << FEI_ENDL;

       dbgOut() << "   numDOFPerNode: ";
       for(int j=0; j<numNodes; j++) {
	 dbgOut() << numDofPerNode[j] << " ";
       }
       dbgOut() << FEI_ENDL;

       return(0);
     }


   int setElemMatrix(int elemBlockID,
		     int elemID,
		     int numNodes,
		     const int* nodeNumbers,
		     const int* dofPerNode,
         const int* dof_ids,
		     const double *const * coefs)
     {
       dbgOut() << "setElemMatrix" << FEI_ENDL
	 << "   elemBlockID: " << elemBlockID << ", elemID: " << elemID << FEI_ENDL
	 << "   numNodes: " << numNodes << FEI_ENDL;
       int i;
       dbgOut() << "   nodeNumbers: ";
       for(i=0; i<numNodes; i++) {
	 dbgOut() << nodeNumbers[i] << " ";
       }
       dbgOut() << FEI_ENDL << "   dofPerNode: ";
       int numRows = 0;
       for(i=0; i<numNodes; i++) {
	 dbgOut() << dofPerNode[i] << " ";
	 numRows += dofPerNode[i];
       }
       dbgOut() << FEI_ENDL << "   coefs:" << FEI_ENDL;
       for(i=0; i<numRows; i++) {
	 dbgOut() << "      ";
	 for(int j=0; j<numRows; j++) {
	   dbgOut() << coefs[i][j] << " ";
	 }
	 dbgOut() << FEI_ENDL;
       }

       return(0);
     }


   int setElemVector(int elemBlockID,
		     int elemID,
		     int numNodes,
		     const int* nodeNumbers,
		     const int* dofPerNode,
         const int* dof_ids,
		     const double* coefs)
     {
       dbgOut() << "setElemVector" << FEI_ENDL
	 << "   elemBlockID: " << elemBlockID << ", elemID: " << elemID << FEI_ENDL
	 << "   numNodes: " << numNodes << FEI_ENDL;
       int i;
       dbgOut() << "   nodeNumbers: ";
       for(i=0; i<numNodes; i++) {
	 dbgOut() << nodeNumbers[i] << " ";
       }
       dbgOut() << FEI_ENDL << "   dofPerNode: ";
       int numRows = 0;
       for(i=0; i<numNodes; i++) {
	 dbgOut() << dofPerNode[i] << " ";
	 numRows += dofPerNode[i];
       }
       dbgOut() << FEI_ENDL << "   coefs:" << FEI_ENDL << "      ";
       for(i=0; i<numRows; i++) {
	 dbgOut() << coefs[i] << " ";
       }
       dbgOut() << FEI_ENDL;
       return(0);
     }

   int setDirichletBCs(int numBCs,
		       const int* nodeNumbers,
		       const int* dofOffsets,
		       const double* values)
     {
       dbgOut() << "setDirichletBCs" << FEI_ENDL
	 << "   numBCs: " << numBCs << FEI_ENDL;
       for(int i=0; i<numBCs; i++) {
	 dbgOut() << "     nodeNumber: " << nodeNumbers[i] << ", "
	   << "dof-offset: " << dofOffsets[i] << ", value: " << values[i]<<FEI_ENDL;
       }

       return(0);
     }

   int sumIntoMatrix(int numRowNodes,
		     const int* rowNodeNumbers,
		     const int* rowDofOffsets,
		     const int* numColNodesPerRow,
		     const int* colNodeNumbers,
		     const int* colDofOffsets,
		     const double* coefs)
     {
       dbgOut() << "sumIntoMatrix, numRowNodes: " << numRowNodes << FEI_ENDL;
       int offset = 0;
       for(int i=0; i<numRowNodes; i++) {
	 dbgOut() << "   rowNodeNumber " << rowNodeNumbers[i]
	   << ", rowDofOffset " << rowDofOffsets[i] << FEI_ENDL;
	 for(int j=0; j<numColNodesPerRow[i]; j++) {
	   dbgOut() << "      colNodeNumber " << colNodeNumbers[offset]
		    << ", colDofOffset " << colDofOffsets[offset]
		    << ", value: " << coefs[offset]<<FEI_ENDL;
	   offset++;
	 }
       }

       return(0);
     }

   int sumIntoRHSVector(int numNodes,
		     const int* nodeNumbers,
		     const int* dofOffsets,
		     const double* coefs)
     {
       dbgOut() << "sumIntoRHSVector, numNodes: " << numNodes << FEI_ENDL;
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << "   nodeNumber " << nodeNumbers[i]
	   << ", dof-offset " << dofOffsets[i] << ", value: " << coefs[i]<<FEI_ENDL;
       }

       return(0);
     }

   int putIntoRHSVector(int numNodes,
		     const int* nodeNumbers,
		     const int* dofOffsets,
		     const double* coefs)
     {
       dbgOut() << "putIntoRHSVector, numNodes: " << numNodes << FEI_ENDL;
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << "   nodeNumber " << nodeNumbers[i]
	   << ", dof-offset " << dofOffsets[i] << ", value: " << coefs[i]<<FEI_ENDL;
       }

       return(0);
     }

   int loadComplete()
     {
       dbgOut() << "loadComplete" << FEI_ENDL;

       return(0);
     }

   /** Function called to request the launching of the linear solver.
    @param solveStatus Output, should indicate the status of the solve. A
    successful solve is usually indicated by a value of 0.
    @param iterations Output, how many iterations were performed.
    @return error-code, 0 if convergence tolerance was achieved within the
     specified maximum number of iterations. If error return is non-zero, the
    calling application will be expected to check solveStatus, and consult the
    solver-library's documentation to figure out exactly what happened.
   */
   int launchSolver(int& solveStatus, int& iterations)
     {
       dbgOut() << "launchSolver" << FEI_ENDL;

       solveStatus = 0;
       iterations = 0;

       return(0);
     }

   int reset()
     {
       dbgOut() << "reset" << FEI_ENDL;
       return(0);
     }

   int resetRHSVector()
     {
       dbgOut() << "resetRHSVector" << FEI_ENDL;
       return(0);
     }

   int resetMatrix()
     {
       dbgOut() << "resetMatrix" << FEI_ENDL;
       return(0);
     }

   int deleteConstraints()
     {
       dbgOut() << "deleteConstraints" << FEI_ENDL;
       return(0);
     }

   int getSolnEntry(int nodeNumber,
		    int dofOffset,
		    double& value)
     {
       dbgOut() << "getSolnEntry, nodeNumber: " << nodeNumber 
	 << ", dofOffset: " << dofOffset << FEI_ENDL;

       value = -999.99;

       return(0);
    }

   int getMultiplierSoln(int CRID, double& lagrangeMultiplier)
     {
       lagrangeMultiplier = -999.99;
       return(0);
     }

   /** Pass nodal data that probably doesn't mean anything to the FEI
     implementation, but may mean something to the linear solver. Examples:
     geometric coordinates, nullspace data, etc.
    @param fieldID Identifier for the field that describes this data. Lists of
       field identifiers and field sizes defined for the finite-element problem
       may be obtained from the Lookup interface that is supplied to the
       ESI_Broker by the FEI implementation.
    @param nodeNumbers List of nodes for which data is being supplied.
    @param numNodes
    @param data List of length numNodes * (size of field 'fieldID')
   */
   int putNodalFieldData(int fieldID,
			 int fieldSize,
			 int numNodes,
			 const int* nodeNumbers,
			 const double* coefs)
     {
       dbgOut() << "putNodalFieldData, fieldID: " << fieldID << ", fieldSize: "
	 << fieldSize << FEI_ENDL;
       int offset = 0;
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << "   nodeNumber " << nodeNumbers[i] << ", coefs: ";
	 for(int j=0; j<fieldSize; j++) {
	   dbgOut() << coefs[offset++] << " ";
	 }
	 dbgOut() << FEI_ENDL;
       }
       
       return(0);
     }

   int setMultiplierCR(int CRID,
		       int numNodes,
		       const int* nodeNumbers,
		       const int* dofOffsets,
		       const double* coefWeights,
		       double rhsValue)
     {
       dbgOut() << "setMultiplierCR, CRID: " << CRID << ", numNodes: " << numNodes << FEI_ENDL;
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << "   nodeNumber " << nodeNumbers[i] << ", dof-offset "
		  << dofOffsets[i] << ", coefWeight: " << coefWeights[i] <<FEI_ENDL;
       }

       dbgOut() << "   rhsValue: " << rhsValue << FEI_ENDL;

       return(0);
     }

   int setPenaltyCR(int CRID,
		    int numNodes,
		    const int* nodeNumbers,
		    const int* dofOffsets,
		    const double* coefWeights,
		    double penaltyValue,
		    double rhsValue)
     {
       dbgOut() << "setPenaltyCR, CRID: " << CRID << ", numNodes: " << numNodes << FEI_ENDL;
       for(int i=0; i<numNodes; i++) {
	 dbgOut() << "   nodeNumber " << nodeNumbers[i] << ", dof-offset "
		  << dofOffsets[i] << ", coefWeight: " << coefWeights[i] <<FEI_ENDL;
       }

       dbgOut() << "   penaltyValue: " << penaltyValue << FEI_ENDL;
       dbgOut() << "   rhsValue: " << rhsValue << FEI_ENDL;

       return(0);
     }

 private:
   int setDebugLog(int debugOutputLevel, const char* path);

   MPI_Comm comm_;
   int numProcs_, localProc_;

   int debugOutputLevel_;
   char* dbgPath_;
   FEI_OSTREAM* dbgOStreamPtr_;
   bool dbgFileOpened_;
   FEI_OFSTREAM* dbgFStreamPtr_;
};

#endif // _FEData_h_
