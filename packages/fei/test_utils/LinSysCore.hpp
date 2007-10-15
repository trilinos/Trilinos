#ifndef _TEST_LinSysCore_h_
#define _TEST_LinSysCore_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

class Lookup;

#include <fei_sstream.hpp>
#include <fei_fstream.hpp>
#include <fei_mpi.h>
#include <fei_defs.h>

#include <fei_Data.hpp>
#include <fei_LinearSystemCore.hpp>

#include <feiArray.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_Lookup.hpp>
#include <fei_SSVec.hpp>
#include <fei_SSMat.hpp>
#include <snl_fei_Utils.hpp>

/** A simple test harness to check whether the FEI layer is calling the
LinearSystemCore layer correctly.
*/

class TEST_LinSysCore : public virtual LinearSystemCore {
 public:
  TEST_LinSysCore(MPI_Comm comm)
    : comm_(comm), numProcs_(1), localProc_(0),
    firstLocalEqn_(0), lastLocalEqn_(0), numGlobalEqns_(0),
    A_(NULL), b_(0, 1), rhsIDs_(0, 1), currentRHS_(0), x_(),
    lookup_(NULL),
    debugOutputPath_(NULL)
    {
#ifndef FEI_SER
      if (MPI_Comm_rank(comm_, &localProc_) != MPI_SUCCESS) abort();
      if (MPI_Comm_size(comm_, &numProcs_)  != MPI_SUCCESS) abort();
#endif
    }

  ~TEST_LinSysCore()
    {
      delete A_;
      for(int j=0; j<b_.length(); j++) delete b_[j];

      delete [] debugOutputPath_;
    }

  /** For cloning a LinearSystemCore instance. Caller recieves a pointer to a
     new instantiation of the implementing object. */
    LinearSystemCore* clone() { return(new TEST_LinSysCore(comm_)); }


  /** For setting argc/argv style parameters.
     @param numParams Number of strings in the params argument
     @param params A list of strings which will usually contain space-separated
        key-value pairs. Example: "debugOutput /usr/users/me/work_dir"
  */

    int parameters(int numParams, char** params)
      {
	const char* param = snl_fei::getParamValue("debugOutput",
							  numParams,params);

	if (param != NULL) {
	  if (debugOutputPath_ != NULL) delete [] debugOutputPath_;
	  debugOutputPath_ = new char[strlen(param)+1];

	  strcpy(debugOutputPath_, param);
	}

	return(0);
      }


  /** Supply the LinearSystemCore implementation with an object (created and
     owned by the caller) that can be used to obtain various information about
     problem layout, shared finite-element nodes, etc.
     For details, see the documentation for the Lookup interface.
     @param lookup Input. Reference to an implementation of the Lookup interface
  */
    int setLookup(Lookup& lookup)
      {
	lookup_ = &lookup;
	if (lookup_ == NULL) return(-1);

	for(int p=0; p<numProcs_; p++) {
	  MPI_Barrier(comm_);
	  if (p != localProc_) continue;

	  int numSharedNodes = lookup_->getNumSharedNodes();
	  if (numSharedNodes<0) return(-1);
/*        const int* shNodeNums = lookup_->getSharedNodeNumbers(); */

/* 	  for(int i=0; i<numSharedNodes; i++) { */
/* 	    int numSubdomains = lookup_->getNumSubdomains(shNodeNums[i]); */
/* 	    FEI_COUT << p << "   shNode " << shNodeNums[i] << ", numSubdomains: " */
/* 	      << numSubdomains << FEI_ENDL; */
/* 	  } */
	}

	return(0);
      }


  /** Supply LinearSystemCore with global offset information for the problem
      being assembled.
      @param len Length of the following list arguments. This will be numProcs+1
      @param nodeOffsets The FEI implementation assigns a global 0-based
          numbering to the finite-element nodes in the problem. Each processor
          is given ownership of a contiguous subset of these node-numbers.
          nodeOffsets[i] gives the first local node-number for processor i.
          nodeOffsets[len-1] gives the total number of nodes.
      @param eqnOffsets eqnOffsets[i] gives the first local equation number for
         processor i, eqnOffsets[len-1] gives the global number of equations.
      @param blkEqnOffsets Contains the same kind of information as eqnOffsets,
          but for 'block-equations'. A block-equation contains all of the
          point-equations present at a finite-element node. Special case: if
          this problem contains Lagrange Multiplier constraints, they will
          be equations that don't correspond to any node, and there will only be
          one of these equations mapped to a block-equation.
  */
    int setGlobalOffsets(int len, int* nodeOffsets,
			 int* eqnOffsets, int* blkEqnOffsets)
      {
	if (len != numProcs_+1) {
	  FEI_CERR << "TEST_LinSysCore::setGlobalOffsets, bad length ("<<len<<")."
	       << FEI_ENDL;
	  abort();
	}

	firstLocalEqn_ = eqnOffsets[localProc_];
	lastLocalEqn_ = eqnOffsets[localProc_+1]-1;
	numGlobalEqns_ = eqnOffsets[numProcs_];

	if (lastLocalEqn_ - firstLocalEqn_ >= 0) {
	  A_ = new SSMat;
	  for(int i=0; i<rhsIDs_.length(); i++) {
	    b_.append(new SSVec);
	  }
	}

	return(0);
      }


   /** For passing element-connectivity arrays.
    @param elemBlock Identifier for the element-block that these elements
        belong to.
    @param numElements Length of the elemIDs list.
    @param numNodesPerElem Length of each row in the connNodes table.
    @param elemIDs Identifiers for each element for which connectivities are
       being supplied.
    @param connNodes Table, with one row for each element. Each row is a list of
       the nodes that are connected to that element.
   */
    int setConnectivities(GlobalID elemBlock,
			  int numElements,
			  int numNodesPerElem,
			  const GlobalID* elemIDs,
			  const int* const* connNodes) { return(0); }


   /** For passing element-stiffness arrays.
    @param elemBlock Identifier for the element-block that these elements
        belong to.
    @param numElems Length of the elemIDs list.
    @param elemIDs Identifiers for each element for which a stiffness array is
       being supplied.
    @param stiff List of 'numElems' tables, each table is of size 
         'numEqnsPerElem' X 'numEqnsPerElem'.
    @param numEqnsPerElem
    @param eqnIndices Table, with 'numElems' rows, each row being a list of
        'numEqnsPerElem' scatter indices (0-based global matrix row/column
        indices).
   */
    int setStiffnessMatrices(GlobalID elemBlock,
			     int numElems,
			     const GlobalID* elemIDs,
			     const double *const *const *stiff,
			     int numEqnsPerElem,
			     const int *const * eqnIndices)
      { return(0); }


   /** For passing element-load vectors.
    @param elemBlock Identifier for the element-block that these elements
        belong to.
    @param numElems Length of the elemIDs list.
    @param elemIDs Identifiers for each element for which a load vector is
       being supplied.
    @param load Table with 'numElems' rows, each row is of length 
         'numEqnsPerElem'.
    @param numEqnsPerElem
    @param eqnIndices Table, with 'numElems' rows, each row being a list of
        'numEqnsPerElem' scatter indices (0-based global equation numbers).
   */
    int setLoadVectors(GlobalID elemBlock,
		       int numElems,
		       const GlobalID* elemIDs,
		       const double *const * load,
		       int numEqnsPerElem,
		       const int *const * eqnIndices)
      { return(0); }


   /** Supply LinearSystemCore with information defining the structure of the
       sparse matrix to be assembled. Implementers of LinearSystemCore may safely
       assume that this function will not be called until after the
       function 'setGlobalOffsets' has been called. Using the information 
       provided via setGlobalOffsets, the number-of-local-equations can be
       trivially calculated. After setMatrixStructure has been called, there
       should be enough information to instantiate internal linear-algebra
       entities, such as vectors, matrix, etc.
     @param ptColIndices Table, with num-local-eqns rows, and the i-th row is of
             length ptRowLengths[i].
     @param ptRowLengths
     @param blkColIndices Table, with num-local-blkEqns rows, and the i-th row
             is of length blkRowLengths[i].
     @param blkRowLengths
     @param ptRowsPerBlkRow The i-th local block-equation corresponds to
          ptRowsPerBlkRow[i] point-equations.
   */
    int setMatrixStructure(int** ptColIndices,
			   int* ptRowLengths,
			   int** blkColIndices,
			   int* blkRowLengths,
			   int* ptRowsPerBlkRow)
      {
	int numLocalEqns = lastLocalEqn_ - firstLocalEqn_ + 1;
	for(int i=0; i<numLocalEqns; ++i) {
	  int row = firstLocalEqn_ + i;
	  for(int j=0; j<ptRowLengths[i]; ++j) {
	    A_->putCoef(row, ptColIndices[i][j], 0.0);
	  }
	}
	return(0);
      }


   /** Specify which global equation numbers correspond to Lagrange Multiplier
      equations. This function won't be called if there are no Lagrange
      Multiplier constraints in the problem. If this function is called, it is
      guaranteed to be called after 'setGlobalOffsets' and before
      'setMatrixStructure'. The primary purpose of this function is to give
      LinearSystemCore implementers the opportunity to deal with constraints in a
      special way, rather than assembling everything into one matrix. If the
      problem being assembled does have Lagrange constraints, then the FEI
      implementation will request an ESI_MatrixRowWriteAccess interface for
      "C_Matrix" from LinearSystemCore. If that is not available, then the FEI
      implementation will request "A_Matrix" and assemble everything into that.
     @param numCRs number of constraint relations
     @param numNodesPerCR number of constrained node in each constraint relation
     @param nodeNumbers Table of constrained nodes. 'numCRs' rows, with the i-th
        row being of length numNodesPerCR[i].
     @param eqnNumbers Table, same dimensions as 'nodeNumbers'. These are the
        global 0-based matrix column indices of the constraint coefficients.
     @param multiplierEqnNumbers Equation numbers that the Lagrange Multipliers
        correspond to.
   */
    int setMultCREqns(int multCRSetID,
		      int numCRs, int numNodesPerCR,
		      int** nodeNumbers, int** eqnNumbers,
		      int* fieldIDs,
		      int* multiplierEqnNumbers) { return(0); }

   /** Specify which nodes and equation numbers correspond to penalty
      constraints. This function is included for completeness, but hasn't
     yet been proven to be useful or necessary, and will probably not be
     included in the successor to LinearSystemCore (ESI_LSManager).
   */
    int setPenCREqns(int penCRSetID,
		     int numCRs, int numNodesPerCR,
		     int** nodeNumbers, int** eqnNumbers,
		     int* fieldIDs) { return(0); }


   /** Provides point-entry data, as well as block-entry data. This is the
    primary assembly function, through which the FEI implementation provides
   the local equation contributions of all element contributions.
   */
    int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
			    int numPtCols, const int* ptCols,
			    int numBlkRows, const int* blkRows,
			    int numBlkCols, const int* blkCols,
			    const double* const* values)
      {
	if (A_ == NULL) {
	  FEI_CERR << "TEST_LinSysCore::sumIntoSystemMatrix A_ == NULL." << FEI_ENDL;
	  return(-1);
	}

	int err = 0;
	for(int i=0; i<numPtRows; i++) {
	  err = A_->sumInRow(ptRows[i], ptCols, values[i], numPtCols);
	  if (err != 0) {
	    FEI_CERR << "TEST_LinSysCore::sumIntoSystemMatrix A_->sumInRow failed."
	      << " row " << ptRows[i] << FEI_ENDL;
	    return(-1);
	  }
	}
	return(0);
      }

   /** Purely point-entry version for accumulating coefficient data into the 
      matrix. This will be called when a matrix contribution fills only part of
      a block-equation. e.g., when a penalty constraint is being applied to a
  single solution field on a node that has several solution fields.
   (A block-equation contains all solution field equations at a node.)
   */
    int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
			    int numPtCols, const int* ptCols,
			    const double* const* values)
      {
	if (A_ == NULL) {
	  FEI_CERR << "TEST_LinSysCore::sumIntoSystemMatrix A_ == NULL." << FEI_ENDL;
	  return(-1);
	}

	int err = 0;
	for(int i=0; i<numPtRows; i++) {
	  err = A_->sumInRow(ptRows[i], ptCols, values[i], numPtCols);
	  if (err != 0) {
	    FEI_CERR << "TEST_LinSysCore::sumIntoSystemMatrix A_->sumInRow failed."
	      << " row " << ptRows[i] << FEI_ENDL;
	    return(-1);
	  }
	}
	return(0);
      }

   /** Point-entry matrix data as for 'sumIntoSystemMatrix', but in this case
     the data should be "put" into the matrix (i.e., overwrite any coefficients
     already present) rather than being "summed" into the matrix.
   */
    int putIntoSystemMatrix(int numPtRows, const int* ptRows,
			    int numPtCols, const int* ptCols,
			    const double* const* values)
      {
	if (A_ == NULL) {
	  FEI_CERR << "TEST_LinSysCore::putIntoSystemMatrix A_ == NULL." << FEI_ENDL;
	  return(-1);
	}

	int err = 0;
	for(int i=0; i<numPtRows; i++) {
	  err = A_->putRow(ptRows[i], ptCols, values[i], numPtCols);
	  if (err != 0) {
	    FEI_CERR << "TEST_LinSysCore::putIntoSystemMatrix A_->sumInRow failed."
	      << " row " << ptRows[i] << FEI_ENDL;
	    return(-1);
	  }
	}
	return(0);
      }

   /** Get the length of a row of the matrix.
       @param row Global 0-based equation number
       @param length Output. Length of the row.
       @return error-code non-zero if any error occurs, e.g., row is not local.
   */
    int getMatrixRowLength(int row, int& length)
      {
	length = -1;
	SSVec* mrow = A_->getRow(row);
	if (mrow == NULL) return(-1);
	length = mrow->length();
	return(0);
      }

   /** Obtain the coefficients and indices for a row of the matrix.
       @param row Global 0-based equation number
       @param coefs Caller-allocated array, length 'len', to be filled with
       coefficients
       @param indices Caller-allocated array, length 'len', to be filled with
       indices. (These indices will be global 0-based equation numbers.)
       @param len Length of the caller-allocated coefs and indices arrays
       @param rowLength Output. Actual length of this row. Not referenced if
       row is not in the local portion of the matrix.
       @return error-code non-zero if any error occurs, e.g., row is not local.
   */
    int getMatrixRow(int row, double* coefs, int* indices,
                            int len, int& rowLength)
      {
	rowLength = -1;
	SSVec* mrow = A_->getRow(row);
	if (mrow == NULL) return(-1);
	rowLength = mrow->length();
	for(int i=0; i<rowLength; ++i) {
	  coefs[i] = mrow->coefs()[i];
	  indices[i] = mrow->indices()[i];
	}
	return(0);
      }

   /** For accumulating coefficients into the rhs vector */

    int sumIntoRHSVector(int num, const double* values,
			 const int* indices)
      {
	if (b_.length() <= 0) return(-1);

	SSVec& b = *(b_[currentRHS_]);

	int err = 0;
	for(int i=0; i<num; i++) {
	  err = b.addEntry(indices[i], values[i]);
	  if (err != 0) return(-1);
	}

	return(0);
      }

   /** For putting coefficients into the rhs vector */
    int putIntoRHSVector(int num, const double* values,
			 const int* indices)
      {
	if (b_.length() <= 0) return(-1);

	SSVec& b = *(b_[currentRHS_]);

	int err = 0;
	for(int i=0; i<num; i++) {
	  err = b.putEntry(indices[i], values[i]);
	  if (err != 0) return(-1);
	}

	return(0);
      }

   /** For getting coefficients out of the rhs vector */
    int getFromRHSVector(int num, double* values,
			 const int* indices)
      {
	SSVec& b = *(b_[currentRHS_]);
	feiArray<int>& inds = b.indices();
	feiArray<double>& coefs = b.coefs();
	for(int i=0; i<num; ++i) {
	  int offset = snl_fei::binarySearch(indices[i], inds);
	  if (offset > -1) {
	    values[i] = coefs[offset];
	  }
	  else {
	    values[i] = -99.9;
	  }
	}
	return(0);
      }

   /** The FEI implementation calls this function to signal the linsyscore
      object that data-loading is finished.
   */
    int matrixLoadComplete()
      {
       std::cout << "debugOutputPath: " << debugOutputPath_ << std::endl;
	dumpMatrixVector(debugOutputPath_);
	return(0);
      }

   /** Pass nodal data that probably doesn't mean anything to the FEI
     implementation, but may mean something to the linear solver. Examples:
     geometric coordinates, nullspace data, etc.
    @param fieldID Identifier for the field that describes this data. Lists of
       field identifiers and field sizes defined for the finite-element problem
       may be obtained from the Lookup interface that is supplied to the
       LinearSystemCore by the FEI implementation.
    @param nodeNumbers List of nodes for which data is being supplied.
    @param numNodes
    @param data List of length numNodes * (size of field 'fieldID')
   */
    int putNodalFieldData(int fieldID, int fieldSize,
                                  int* nodeNumbers, int numNodes,
                                  const double* data) { return(-1); }


   /** For setting the scalar 's' (usually 0.0) throughout the matrix rhs 
    vector.
   */
    int resetMatrixAndVector(double s)
      {
	if (A_ != NULL) A_->logicalClear();

	for(int i=0; i<b_.length(); i++) b_[i]->logicalClear();

	return(0);
      }

   /** For setting the scalar 's' (usually 0.0) throughout the matrix.
   */
    int resetMatrix(double s)
      {
	if (A_ != NULL) A_->logicalClear();
	return(0);
      }

   /** For setting the scalar 's' (usually 0.0) throughout the rhs vector.
   */
    int resetRHSVector(double s)
      {
	for(int i=0; i<b_.length(); i++) b_[i]->logicalClear();
	return(0);
      }

   /** The FEI implementation calls this function to inform LinearSystemCore
       of equations that need to have essential (Dirichlet) boundary conditions
      enforced on them. The intent is that the LinearSystemCore implementation
      will  perform the column-modification b[i] = gamma[i]/alpha[i] if i == 
     globalEqn[i], otherwise b[i] -= gamma[i]/alpha[i] * A(i,globalEqn[i]) if 
     i != globalEqn[i]. After this operation is performed, all of row
     globalEqn[i] and column globalEqn[i] should be set to 0.0, except for the
     diagonal position. (Naturally the implementer is free to enforce the
     boundary condition another way if they wish.)
     @param globalEqn List, of length 'len', of global 0-based equation numbers.
     @param alpha List, of length 'len', of coefficients. When the solution to
        the linear system is later requested, the solution value for
       globalEqn[i] should be gamma[i]/alpha[i].
     @param gamma
     @param len
   */
    int enforceEssentialBC(int* globalEqn,
			   double* alpha,
			   double* gamma,
			   int len)
      {
	SSVec& b = *(b_[currentRHS_]);
	SSMat& A = *A_;
	feiArray<int>& rowNumbers = A.getRowNumbers();
	feiArray<SSVec*>& rows = A.getRows();
	int err;

	for(int i=0; i<len; i++) {
	  int bceqn = globalEqn[i];
	  double value = gamma[i]/alpha[i];

	  //first, do the column operations.
	  //b[bceqn] -= A(row, bceqn)*value for all row != bceqn
	  //then A(row, bceqn) = 0.0;

	  for(int row=firstLocalEqn_; row<= lastLocalEqn_; row++) {
	    if (row == bceqn) continue;

	    double* AvalPtr = NULL;
	    err = A.coefficientPointer(row, bceqn, AvalPtr);
	    if (err != 0) continue;

	    double rhsterm = *AvalPtr*value * -1.0;
	    err = b.addEntry(row, rhsterm);

	    *AvalPtr = 0.0;
	  }

	  //Now we'll set row bceqn of A to be all zero except for 1.0 on
	  //the diagonal, if row bceqn is local.

	  if (firstLocalEqn_ > bceqn || lastLocalEqn_ < bceqn) continue;
	  int index = snl_fei::binarySearch(bceqn, rowNumbers);
	  if (index < 0) continue;

	  feiArray<int>& colIndices = A.getRows()[index]->indices();
	  int rowLength = rows[index]->length();
	  for(int col=0; col<rowLength; col++) {
	    double* Aptr = NULL;
	    err = A.coefficientPointer(bceqn, colIndices[col], Aptr);
	    if (err != 0) continue;

	    if (bceqn == colIndices[col]) *Aptr = 1.0;
	    else *Aptr = 0.0;
	  }
	}

	return(0);
      }

   /** The FEI implementation calls this function to inform LinearSystemCore
     that certain local equations need column-modifications made due to
     essential boundary-conditions being enforced on other processors. The
     column modification is roughly this: b(globalEqns[i]) -= A(globalEqns[i],
       colIndices[i][j]) * coefs[i][j], for i in [0..numEqns-1] and j in [0..
       colIndLen[i]-1]. (Note that A(globalEqns[i], colIndices[i][j]) should be
       set = 0.0 after the appropriate value has been accumulated into b also.)
     @param numEqns Length of 'globalEqns'
     @param globalEqns Equations that are local to this processor.
     @param colIndices Table, with 'numEqns' rows, and the i-th row is of length
          colIndLen[i]. The i-th row contains column indices in equation
        globalEqns[i]. These column indices correspond to equations (rows) that
        are owned by other processors, and those other processors are imposing
        essential boundary conditions on those equations.
     @param colIndLen List of length 'numEqns'.
     @param coefs This table holds the gamma/alpha coeficients that are the
        value of the boundary conditions being enforced on each of the remote
        equations in the 'colIndices' table.
   */
    int enforceRemoteEssBCs(int numEqns, int* globalEqns,
			    int** colIndices, int* colIndLen,
			    double** coefs)
      {
	SSVec& b = *(b_[currentRHS_]);
	SSMat& A = *A_;
	feiArray<int>& rowNumbers = A.getRowNumbers();
	int err;

	for(int i=0; i<numEqns; i++) {
	  int row = globalEqns[i];
	  int index = snl_fei::binarySearch(row, rowNumbers);
	  if (index < 0) continue;

	  for(int j=0; j<colIndLen[i]; j++) {
	    double* Aptr = NULL;
	    err = A.coefficientPointer(row, colIndices[i][j], Aptr);
	    if (err != 0) continue;
	    double value = *Aptr * coefs[i][j] * -1.0;
	    err = b.addEntry(row, value);
	    *Aptr = 0.0;
	  }
	}

	return(0);
      }

   /** This function is called to inform LinearSystemCore that natural (or
    Neumann) or mixed boundary conditions are to be enforce on some equations.
    Basically, these can be enforced by performing this operation: 
    A(globalEqn[i], globalEqn[i]) += alpha[i]/beta[i], and b(globalEqn[i]) += 
    gamma[i]/beta[i], for all i in [0 .. len-1]. (Note that alpha[i]==0.0
    implies a Natural or Neumann boundary condition, while gamma[i]==0.0 implies
    a mixed boundary condition.
    @param globalEqn List, length 'len', of equation on which to impose
       boundary conditions.
    @param alpha List, length 'len', of coefficients described above.
    @param beta
    @param gamma
    @param len
   */
    int enforceOtherBC(int* globalEqn, double* alpha,
		       double* beta, double* gamma, int len)
      {
	SSVec& b = *(b_[currentRHS_]);
	SSMat& A = *A_;
	feiArray<int>& rowNumbers = A.getRowNumbers();
	int err;

	for(int i=0; i<len; i++) {
	  int row = globalEqn[i];
	  if (firstLocalEqn_ > row || lastLocalEqn_ < row) continue;

	  err = b.addEntry(row, gamma[i]/beta[i]);

	  int index = snl_fei::binarySearch(row, rowNumbers);
	  if (index < 0) continue;
	  double* Aptr = NULL;
	  err = A.coefficientPointer(row, row, Aptr);
	  if (err != 0) continue;
	  *Aptr += alpha[i]/beta[i];
	}

	return(0);
      }

   /** The FEI implementation calls this function to request a pointer to the
     internal 'A-matrix' data.
     @param data See Data class documentation.
   */
    int getMatrixPtr(Data& data) { return(-1); }

   /** LinearSystemCore's internal 'A-matrix' should be replaced with a scaled
     copy of the incoming data.
     @param scalar coefficient by which to scale the incoming data.
     @param data See documentation for Data class.
   */
    int copyInMatrix(double scalar, const Data& data) { return(-1); }

   /** The FEI implementation calls this function to request a scaled copy of
     the internal 'A-matrix' data. The FEI implementation will then be
    responsible for deciding when this matrix data should be destroyed. The
    LinearSystemCore implementation should not keep a reference to the pointer
    that was handed out.
    @param scalar
    @param data See documentation for Data class.
   */
    int copyOutMatrix(double scalar, Data& data) { return(-1); }

   /** A scaled copy of the incoming data should be added to the internal
    'A-matrix'.
    @param scalar
    @param data See documentation for Data class.
   */
    int sumInMatrix(double scalar, const Data& data) { return(-1); }

   /** Same semantics as getMatrixPtr, but applied to rhs vector. */
    int getRHSVectorPtr(Data& data) { return(-1); }

   /** Same semantics as copyInMatrix, but applied to rhs vector. */
    int copyInRHSVector(double scalar, const Data& data) { return(-1); }

   /** Same semantics as copyOutMatrix, but applied to rhs vector. */
    int copyOutRHSVector(double scalar, Data& data) { return(-1); }

   /** Same semantics as sumInMatrix, but applied to rhs vector. */
    int sumInRHSVector(double scalar, const Data& data) { return(-1); }

   /** Utility function for destroying the matrix in a Data container. The 
    caller (owner of 'data') can't destroy the matrix because they don't know
    what concrete type it is and can't get to its destructor. The contents of 
    'data' is a matrix previously passed out via 'copyOutMatrix'.
    @param data See documentation for Data class.
   */
    int destroyMatrixData(Data& data) { return(-1); }

   /** Utility function for destroying the vector in a Data container. The 
 caller (owner of 'data') can't destroy the vector because they don't know what 
    concrete type it is and can't get to its destructor. The contents of 'data'
    is a vector previously passed out via 'copyOutRHSVector'.
    @param data See documentation for Data class.
   */
    int destroyVectorData(Data& data) { return(-1); }

   /** Indicate the number of rhs-vectors being assembled/solved for.
     This function will be called by the FEI implementation at or near the
     beginning of the problem assembly. If numRHSs is greater than 1, then
     calls to 'getMemberInterface' requesting an interface to an rhs vector
     will use the 'objName' argument "b_Vector_n", where n is in [0 .. 
     numRHSs-1]. If there is only one rhs-vector, then 'objName' will be
     simply "b_Vector".
    @param numRHSs Length of the rhsIDs list.
    @param rhsIDs Caller-supplied integer identifiers for the rhs vectors. This
       argument will probably be removed, as it is obsolete (a carry-over from
       LinearSystemCore).
   */
    int setNumRHSVectors(int numRHSs, const int* rhsIDs)
      {
	rhsIDs_.resize(numRHSs);

	for(int i=0; i<numRHSs; i++) rhsIDs_[i] = rhsIDs[i];
	currentRHS_ = 0;
	return(0);
      }

   /** Set the 'current context' to the rhs-vector corresponding to rhsID. 
    Subsequent data received via 'sumIntoRHSVector' should be directed into
    this rhs-vector. Any other function-calls having to do with the rhs-vector
    should also effect this rhs-vector.
    @param rhsID
   */
    int setRHSID(int rhsID)
      {
	int index = rhsIDs_.find(rhsID);
	if (index < 0) {
	  FEI_CERR << "TEST_LinSysCore::setRHSID rhsID " << rhsID << " not found."
	       << FEI_ENDL;
	  return(-1);
	}

	currentRHS_ = index;

	return(0);
      }

   /** The FEI implementation will call this function to supply initial-guess
     data that should be used as the starting 'x-vector' if an iterative
    solution is to be performed.
    @param eqnNumbers Global 0-based equation numbers for which the initial
      guess should be set.
    @param values The initial guess data.
    @param len Number of equations for which an initial guess is being supplied.
   */
    int putInitialGuess(const int* eqnNumbers, const double* values,
                                int len) { return(0); }

   /** The FEI implementation will call this function to request all local
    solution values.
    @param answers Solution coefficients.
    @param len This should equal the number of local equations. If it is less,
      the LinearSystemCore implementation should simply pass out the first 'len'
      local solution values. If it is more, then just pass out numLocalEqns
      solution values.
   */
    int getSolution(double* answers, int len)
      {
	for(int i=0; i<len; ++i) answers[i] = -99.9;
	return(0); 
      }

   /** The FEI implementation will call this function to request a single
    solution value.
    @param eqnNumber Global 0-based equation number.
    @param answer
   */
    int getSolnEntry(int eqnNumber, double& answer)
      {
	answer = -99.9;
	return(0);
      }

   /** This will be called to request that LinearSystemCore form the residual
    vector r = b - A*x, and pass the coefficients for r back out in the 'values'
    list.
    @param values
    @param len This should equal num-local-eqns.
   */
    int formResidual(double* values, int len) { return(-1); }

   /** Function called to request the launching of the linear solver.
    @param solveStatus Output, should indicate the status of the solve. A
    successful solve is usually indicated by a value of 0.
    @param iterations Output, how many iterations were performed.
    @return error-code, 0 if convergence tolerance was achieved within the
     specified maximum number of iterations. If error return is non-zero, the
    calling application will be expected to check solveStatus, and consult the
    solver-library's documentation to figure out exactly what happened.
   */
    int launchSolver(int& solveStatus, int& iterations) { return(-1); }

   /** This function's intent is to provide a file-name to be used by 
    LinearSystemCore in writing the linear system into disk files. Format is
    not specified. Implementers may choose to augment this name in the style
    of writing 3 files: A_name, x_name, b_name, or some other convention.
    This function is ill-defined, obsolete, and will probably never be called
    by the FEI implementation.
   */
    int writeSystem(const char* name) { return(-1); }

 private:
    int dumpMatrixVector(char* path)
      {
	if (A_ == NULL) return(-1);
	if (b_.length() <= 0) return(-1);

	if (path == NULL) return(-1);

	FEI_OSTRINGSTREAM filename;

	filename<< path<<"/A_TLSC.mtx."<<numProcs_<<"."<<localProc_;

	FEI_OFSTREAM ofstr(filename.str().c_str());
	if (ofstr.bad()) {
	  FEI_CERR << "TEST_LinSysCore::matrixLoadComplete: failed to open "
	    << filename << FEI_ENDL;
	  return(-1);
	}

	ofstr.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

	ofstr << numGlobalEqns_ << " " << numGlobalEqns_ << FEI_ENDL;

	feiArray<int>& rowNumbers = A_->getRowNumbers();
	feiArray<SSVec*>& rows = A_->getRows();

	int i;
	for(i=0; i<rowNumbers.length(); i++) {
	  feiArray<int>& cols = rows[i]->indices();
	  feiArray<double>& rowCoefs = rows[i]->coefs();

	  for(int j=0; j<cols.length(); j++) {
	    ofstr << rowNumbers[i] << " " << cols[j] << " "
		  << rowCoefs[j] << FEI_ENDL;
	  }
	}

	ofstr.close();

        FEI_OSTRINGSTREAM bname;
	bname<< path<<"/b_TLSC."<<numProcs_<<"."<<localProc_;
	ofstr.open(bname.str().c_str());

	if (ofstr.bad()) {
	  FEI_CERR << "TEST_LinSysCore::matrixLoadComplete: failed to open "
	    << filename << FEI_ENDL;
	  return(-1);
	}

	feiArray<int>& eqns = b_[0]->indices();
	feiArray<double>& entries = b_[0]->coefs();

	ofstr << numGlobalEqns_ << FEI_ENDL;

	for(i=0; i<eqns.length(); i++) {
	  ofstr << eqns[i] << " " << entries[i] << FEI_ENDL;
	}

	return(0);
      }

    MPI_Comm comm_;
    int numProcs_, localProc_;

    int firstLocalEqn_, lastLocalEqn_, numGlobalEqns_;

    SSMat* A_;
    feiArray<SSVec*> b_;
    feiArray<int> rhsIDs_;
    int currentRHS_;
    SSVec x_;

    Lookup* lookup_;

    char* debugOutputPath_;
};

#endif // _TEST_LinSysCore_h_
