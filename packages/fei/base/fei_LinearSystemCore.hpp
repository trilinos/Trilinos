#ifndef _fei_LinearSystemCore_hpp_
#define _fei_LinearSystemCore_hpp_

class Data;
class Lookup;

#include <fei_defs.h>
#include <cstdlib>

/**  
  This is the original internal FEI interface to solver-libraries --
  the destination for the data being assembled into a linear system by the FEI 
  implementation.

  When creating a specific FEI implementation, i.e., a version that
  supports a specific underlying linear solver library, the main
  task that must be performed is the implementation of this interface,
  LinearSystemCore (Note: the LinearSystemCore interface is being
  replaced by the ESI_Broker interface.).

  To date (as of August 2001), implementations of this interface exist for
  coupling the following solver libraries to the FEI implementation:
  <ul>
  <li>Aztec
  <li>HYPRE
  <li>ISIS++
  <li>PETSc
  <li>Prometheus
  <li>SPOOLES
  </ul>

  An implementation of LinearSystemCore holds and manipulates all
  solver-library-specific stuff, such as matrices/vectors, 
  solvers/preconditioners, etc. An instance of this class is owned and used by
  the class that implements the public FEI spec. i.e., when element
  contributions, etc., are received from the finite-element application, the
  data is ultimately passed to this
  class for assembly into the sparse matrix and associated vectors. This class
  will also be asked to launch any underlying solver, and finally to
  return the solution.

  Terminology notes:
  - I generally refer to connectivities as nodes or nodeIDs, which are
    finite-element entities. However, the LinearSystemCore interface functions
    don't receive nodeIDs -- instead, nodeNumbers are provided by the FEI
    implementation code layer. See note 3 below.

  - I refer to the rows of a matrix or entries in a vector as equations
    or equation-numbers, which are linear-algebra entities.

  - An equation can be mapped to a node/field pair. i.e., given a globally
    unique node number, and a solution field on that node, a global
    equation number can be obtained, and vice versa.
  
  - The term 'block' is used in two different ways here. A block-entry matrix
    is a matrix made up of dense sub-blocks, with each sub-block containing a
    number of equations. An element-block is a finite-element term, denoting a
    group of elements which are homogeneous in topology. i.e., all elements in
    an element-block have the same number of nodes, and the same numbers of
    fields per node.

   Key services provided by the FEI implementation layer:

   Mappings:
     Map from node/field pairs to equation numbers.
     Map from lagrange multipliers to equation numbers.
     Map from element-dofs to equation numbers.

   Partitioning and communication:
     Decide which processor shared equations belong to.
     Move element-contributions for shared equations from
     the sharing/contributing processor to the owning processor.
     i.e., the LinearSystemCore object will only be given local
     equation data.

NOTES:

1. The LinearSystemCore object is given a 'Lookup' interface, from which it can
  obtain information about the structure of the finite-element problem. e.g.,
  it can look up the number of fields, the fields' sizes, the number of
  element-blocks, etc., etc. For details on the Lookup interface, see the file
  Lookup.h.

2.  Except for making calls on the Lookup object, LinearSystemCore is reactive.
  All of the functions below are member functions of LinearSystemCore, which
  are called by FEI code. 'set' and 'put' functions are for information to be
  provided TO LinearSystemCore by the FEI implementation layer, while 'get'
  functions are for information to be requested FROM LinearSystemCore by the
  FEI implementation layer.

3.  Node-numbers that are given to LinearSystemCore (through setConnectivities,
  etc.,) are 0-based and contiguous. Each processor owns a contiguous block
  of node-numbers. However, non-local node-numbers may of course appear in
  connectivity lists. These node-numbers are not the same as the 'nodeIDs'
  that appear in the FEI public interface, and throughout the FEI implementation
  code. nodeIDs are provided by the application, and may not be assumed to be
  0-based or contiguous. They are arbitrary node identifiers which need only
  be globally unique. nodeIDs never appear in LinearSystemCore functions.

4.  Not all LinearSystemCore functions are necessary for assembling a
  linear system. For instance, many implementations will ignore the data that
  is passed in by 'setConnectivities', 'setStiffnessMatrices', etc.

5.  Some data is passed into LinearSystemCore redundantly. i.e., matrix
  coefficient data appears in the 'setStiffnessMatrices' function as well as
  the 'sumIntoSystemMatrix/sumIntoSystemBlkMatrix' function. The difference is,
  only local data is supplied through sumIntoSystemMatrix. Additionally, the
  FEI implementation layer performs the necessary communication to ensure that
  data from shared finite-element nodes is moved onto the 'owning' processor
  before being passed to LinearSystemCore::sumIntoSystemMatrix.
  The finite-element application doesn't have the concept of node ownership,
  instead processors can equally share nodes. The FEI implementation designates
  a single processor as the owner of these nodes, and moves any shared 
  stiffness data onto the appropriate owning processor.
  The data supplied through setStiffnessMatrices is the un-modified stiffness
  data supplied by the application, including any shared-non-local portions.

  Similarly for setLoadVectors and sumIntoRHSVector. sumIntoRHSVector only
  provides local data, while setLoadVectors supplies the un-modified element-
  load vectors from the application, including any shared portions.

*/  

class LinearSystemCore {
 public:
  /** Default constructor, typically overridden by the implementing class, which
   probably requires an MPI Communicator, and possibly other arguments.
  */
   LinearSystemCore(){};

  /** Destructor, should call through to the implementation's destructor and
    destroy all allocated memory, internal objects, etc. Exceptions: objects
    created in reponse to calls to the functions 'copyOutMatrix' and
    'copyOutRHSVector' are not destroyed here. The caller is assumed to have
   taken responsibility for those matrix/vector copies.
  */
   virtual ~LinearSystemCore() {};


  /** For cloning a LinearSystemCore instance. Caller recieves a pointer to a
     new instantiation of the implementing object. */
   virtual LinearSystemCore* clone() = 0;


  /** For setting argc/argv style parameters.
     @param numParams Number of strings in the params argument
     @param params A list of strings which will usually contain space-separated
        key-value pairs. Example: "debugOutput /usr/users/me/work_dir"
  */

   virtual int parameters(int numParams, const char*const * params) = 0;


  /** Supply the LinearSystemCore implementation with an object (created and
     owned by the caller) that can be used to obtain various information about
     problem layout, shared finite-element nodes, etc.
     For details, see the documentation for the Lookup interface.
     @param lookup Input. Reference to an implementation of the Lookup interface
  */
   virtual int setLookup(Lookup& lookup) = 0;


   /** Query a named property (such as timing statistics, etc.) from the solver
       library.
       @param name Input. Name of the property for which a value is being
       requested.
       @pararm value Output. Requested property's value.
       @return error-code 0 if successful. -1 probably indicates that the
       named property is not recognized.
   */
   virtual int getProperty(const char* /*name*/, double& /*value*/)
     {
       return(-1);
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
   virtual int setGlobalOffsets(int len, int* nodeOffsets,
                                 int* eqnOffsets, int* blkEqnOffsets) = 0;


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
   virtual int setConnectivities(GlobalID elemBlock,
                                  int numElements,
                                  int numNodesPerElem,
                                  const GlobalID* elemIDs,
                                  const int* const* connNodes) = 0;


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
   virtual int setStiffnessMatrices(GlobalID elemBlock,
                                     int numElems,
                                     const GlobalID* elemIDs,
                                     const double *const *const *stiff,
                                     int numEqnsPerElem,
                                     const int *const * eqnIndices) = 0;


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
   virtual int setLoadVectors(GlobalID elemBlock,
                               int numElems,
                               const GlobalID* elemIDs,
                               const double *const * load,
                               int numEqnsPerElem,
                               const int *const * eqnIndices) = 0;


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
   virtual int setMatrixStructure(int** ptColIndices,
                                   int* ptRrowLengths,
                                   int** blkColIndices,
                                   int* blkRowLengths,
                                   int* ptRowsPerBlkRow) = 0;


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
   virtual int setMultCREqns(int multCRSetID,
                              int numCRs, int numNodesPerCR,
                              int** nodeNumbers, int** eqnNumbers,
                              int* fieldIDs,
                              int* multiplierEqnNumbers) = 0;

   /** Specify which nodes and equation numbers correspond to penalty
      constraints. This function is included for completeness, but hasn't
     yet been proven to be useful or necessary, and will probably not be
     included in the successor to LinearSystemCore (ESI_LSManager).
   */
   virtual int setPenCREqns(int penCRSetID,
                              int numCRs, int numNodesPerCR,
                              int** nodeNumbers, int** eqnNumbers,
                              int* fieldIDs) = 0;


   /** Provides point-entry data, as well as block-entry data. This is the
    primary assembly function, through which the FEI implementation provides
   the local equation contributions of all element contributions.
   */
   virtual int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                                    int numPtCols, const int* ptCols,
                                    int numBlkRows, const int* blkRows,
                                    int numBlkCols, const int* blkCols,
                                    const double* const* values) = 0;

   /** Purely point-entry version for accumulating coefficient data into the 
      matrix. This will be called when a matrix contribution fills only part of
      a block-equation. e.g., when a penalty constraint is being applied to a
  single solution field on a node that has several solution fields.
   (A block-equation contains all solution field equations at a node.)
   */
   virtual int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                                    int numPtCols, const int* ptCols,
                                    const double* const* values) = 0;

   /** Point-entry matrix data as for 'sumIntoSystemMatrix', but in this case
     the data should be "put" into the matrix (i.e., overwrite any coefficients
     already present) rather than being "summed" into the matrix.
   */
   virtual int putIntoSystemMatrix(int numPtRows, const int* ptRows,
                                    int numPtCols, const int* ptCols,
                                    const double* const* values) = 0;

   /** Get the length of a row of the matrix.
       @param row Global 0-based equation number
       @param length Output. Length of the row.
       @return error-code non-zero if any error occurs, e.g., row is not local.
   */
   virtual int getMatrixRowLength(int row, int& length) = 0;

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
   virtual int getMatrixRow(int row, double* coefs, int* indices,
                            int len, int& rowLength) = 0;

   /** For accumulating coefficients into the rhs vector */

   virtual int sumIntoRHSVector(int num, const double* values,
                                 const int* indices) = 0;
   /** For putting coefficients into the rhs vector */
   virtual int putIntoRHSVector(int num, const double* values,
                                 const int* indices) = 0;
   /** For getting coefficients out of the rhs vector */
   virtual int getFromRHSVector(int num, double* values,
                                 const int* indices) = 0;

   /** The FEI implementation calls this function to signal the linsyscore
      object that data-loading is finished.
   */
   virtual int matrixLoadComplete() = 0;

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
   virtual int putNodalFieldData(int fieldID, int fieldSize,
                                  int* nodeNumbers, int numNodes,
                                  const double* data) = 0;


   /** For setting the scalar 's' (usually 0.0) throughout the matrix rhs 
    vector.
   */
   virtual int resetMatrixAndVector(double s) = 0;

   /** For setting the scalar 's' (usually 0.0) throughout the matrix.
   */
   virtual int resetMatrix(double s) = 0;

   /** For setting the scalar 's' (usually 0.0) throughout the rhs vector.
   */
   virtual int resetRHSVector(double s) = 0;

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
   virtual int enforceEssentialBC(int* globalEqn, double* alpha,
                                   double* gamma, int len) = 0;

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
   virtual int enforceRemoteEssBCs(int numEqns, int* globalEqns,
                                          int** colIndices, int* colIndLen,
                                          double** coefs) = 0;

   /** The FEI implementation calls this function to request a pointer to the
     internal 'A-matrix' data.
     @param data See Data class documentation.
   */
   virtual int getMatrixPtr(Data& data) = 0;

   /** LinearSystemCore's internal 'A-matrix' should be replaced with a scaled
     copy of the incoming data.
     @param scalar coefficient by which to scale the incoming data.
     @param data See documentation for Data class.
   */
   virtual int copyInMatrix(double scalar, const Data& data) = 0;

   /** The FEI implementation calls this function to request a scaled copy of
     the internal 'A-matrix' data. The FEI implementation will then be
    responsible for deciding when this matrix data should be destroyed. The
    LinearSystemCore implementation should not keep a reference to the pointer
    that was handed out.
    @param scalar
    @param data See documentation for Data class.
   */
   virtual int copyOutMatrix(double scalar, Data& data) = 0;

   /** A scaled copy of the incoming data should be added to the internal
    'A-matrix'.
    @param scalar
    @param data See documentation for Data class.
   */
   virtual int sumInMatrix(double scalar, const Data& data) = 0;

   /** Same semantics as getMatrixPtr, but applied to rhs vector. */
   virtual int getRHSVectorPtr(Data& data) = 0;

   /** Same semantics as copyInMatrix, but applied to rhs vector. */
   virtual int copyInRHSVector(double scalar, const Data& data) = 0;

   /** Same semantics as copyOutMatrix, but applied to rhs vector. */
   virtual int copyOutRHSVector(double scalar, Data& data) = 0;

   /** Same semantics as sumInMatrix, but applied to rhs vector. */
   virtual int sumInRHSVector(double scalar, const Data& data) = 0;

   /** Utility function for destroying the matrix in a Data container. The 
    caller (owner of 'data') can't destroy the matrix because they don't know
    what concrete type it is and can't get to its destructor. The contents of 
    'data' is a matrix previously passed out via 'copyOutMatrix'.
    @param data See documentation for Data class.
   */
   virtual int destroyMatrixData(Data& data) = 0;

   /** Utility function for destroying the vector in a Data container. The 
 caller (owner of 'data') can't destroy the vector because they don't know what 
    concrete type it is and can't get to its destructor. The contents of 'data'
    is a vector previously passed out via 'copyOutRHSVector'.
    @param data See documentation for Data class.
   */
   virtual int destroyVectorData(Data& data) = 0;

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
   virtual int setNumRHSVectors(int numRHSs, const int* rhsIDs) = 0;

   /** Set the 'current context' to the rhs-vector corresponding to rhsID. 
    Subsequent data received via 'sumIntoRHSVector' should be directed into
    this rhs-vector. Any other function-calls having to do with the rhs-vector
    should also effect this rhs-vector.
    @param rhsID
   */
   virtual int setRHSID(int rhsID) = 0;

   /** The FEI implementation will call this function to supply initial-guess
     data that should be used as the starting 'x-vector' if an iterative
    solution is to be performed.
    @param eqnNumbers Global 0-based equation numbers for which the initial
      guess should be set.
    @param values The initial guess data.
    @param len Number of equations for which an initial guess is being supplied.
   */
   virtual int putInitialGuess(const int* eqnNumbers, const double* values,
                                int len) = 0;

   /** The FEI implementation will call this function to request all local
    solution values.
    @param answers Solution coefficients.
    @param len This should equal the number of local equations. If it is less,
      the LinearSystemCore implementation should simply pass out the first 'len'
      local solution values. If it is more, then just pass out numLocalEqns
      solution values.
   */
   virtual int getSolution(double* answers, int len) = 0;

   /** The FEI implementation will call this function to request a single
    solution value.
    @param eqnNumber Global 0-based equation number.
    @param answer
   */
   virtual int getSolnEntry(int eqnNumber, double& answer) = 0;

   /** This will be called to request that LinearSystemCore form the residual
    vector r = b - A*x, and pass the coefficients for r back out in the 'values'
    list.
    @param values
    @param len This should equal num-local-eqns.
   */
   virtual int formResidual(double* values, int len) = 0;

   /** Function called to request the launching of the linear solver.
    @param solveStatus Output, should indicate the status of the solve. A
    successful solve is usually indicated by a value of 0.
    @param iterations Output, how many iterations were performed.
    @return error-code, 0 if convergence tolerance was achieved within the
     specified maximum number of iterations. If error return is non-zero, the
    calling application will be expected to check solveStatus, and consult the
    solver-library's documentation to figure out exactly what happened.
   */
   virtual int launchSolver(int& solveStatus, int& iterations) = 0;

   /** This function's intent is to provide a file-name to be used by 
    LinearSystemCore in writing the linear system into disk files. Format is
    not specified. Implementers may choose to augment this name in the style
    of writing 3 files: A_name, x_name, b_name, or some other convention.
    This function is ill-defined, obsolete, and will probably never be called
    by the FEI implementation.
   */
   virtual int writeSystem(const char* name) = 0;

   virtual double* getMatrixBeginPointer()
   { return NULL; }

   virtual int getMatrixOffset(int /*row*/, int /*col*/)
   { return -1; }

};

#endif

